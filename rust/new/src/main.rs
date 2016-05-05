#[macro_use]
extern crate ndarray;

use ndarray::{
    ArrayView,
    ArrayViewMut,
    OwnedArray,
    Ix,
};

use std::f64;
use std::f64::consts::PI;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

type Ix2 = (Ix, Ix);

fn main() {
    let nx = 400;        // columns
    let ny = 400;        // rows
    let nx2 = nx + 1;
    let ny2 = ny + 1;

    let rho: f64 = 1500.0;      // density
    let vp0: f64 = 2000.0;      // velocity

    // array for velocity initialized with ones
    /* old style
    let mut vp = OwnedArray::from_elem((ny2, nx2), 1.0);

    for (_, elt) in vp.indexed_iter_mut() {
        // multiply each cell of array (which is 1.0) by vp0
        *elt *= vp0
    }
     */

    // new style for initialization
    let vp= OwnedArray::from_elem((ny2, nx2), 1.0f64).map(|& x| x*vp0*vp0);

    println!("vp[:5, :5] is:\n{}", vp.slice(s![..5, ..5]));

    // define source and receiver positions
    let isx = 200;
    let isy = 100;
    let irx = 200;
    let iry = 300;

    // source and its parameters
    let nt = 3000;
    let fc: f64 = 1000.0;
    let mut dx: f64 = vp0/(30.0*fc);
    let mut dy = dx;
    let dt: f64 = dx/vp0*1.0/((2.0 as f64).sqrt()*2.0*PI);

    println!("nt: {}\nfc: {}\ndx: {}\ndy: {}\ndt: {}", nt, fc, dx, dy, dt);

    // define source function in time
    let tsour = 1.0/fc;
    let mut t = OwnedArray::from_elem((nt), 0.0f64);

    for (i, elt) in t.indexed_iter_mut() {
        *elt = (i as f64)*dt;
    }
    // println!("time :\n{}", t.slice(s![..10]));
    let t0 = tsour*1.5;
    let tau = t.map(|&elt| PI*(elt - t0)/t0);
    let a = 4.0;
    let fs = tau.map(|&x| (1.0 - a*x*x)*((-2.0*x*x).exp()));

    // println!("source function: {}", fs.slice(s![..10]));
    for g in 0..10 {
      println!("{}, {}", fs.get((g)).unwrap(), fs.get((nt-1-g)).unwrap());
    }

    // define PML properties
    let npml = 30;
    let pmlfac: f64 = 50.0;
    let pmlexp: i32 = 2;
    let mut qx = OwnedArray::from_elem((ny2, nx2), 0f64);
    let mut qy = OwnedArray::from_elem((ny2, nx2), 0f64);

    for i in 0..npml {
        for j in 0..nx2 {
            unsafe{
                *qy.uget_mut((i, j)) = pmlfac*(((npml-i) as f64 ).powi(pmlexp));
                *qy.uget_mut((ny2-i-1, j)) = pmlfac*(((npml-i) as f64).powi(pmlexp));
            }
        }
    };
    for i in 0..ny2 {
        for j in 0..npml {
            unsafe{
                *qx.uget_mut((i, nx2-j-1)) = pmlfac*(((npml-j) as f64).powi(pmlexp));
                *qx.uget_mut((i, j)) = pmlfac*(((npml-j) as f64).powi(pmlexp));
            }
        }
    };

    // initialize fields
    let mut px = OwnedArray::from_elem((ny2, nx2), 0f64);
    let mut py = OwnedArray::from_elem((ny2, nx2), 0f64);
    let mut ux = OwnedArray::from_elem((ny2, nx2), 0f64);
    let mut uy = OwnedArray::from_elem((ny2, nx2), 0f64);

    // receiver trace
    let mut sfd = OwnedArray::from_elem((nt), 0f64);

    let mut diffop: f64;
    let mut pmlop: f64;

    //some optimization tricks
    dx = 1.0/dx;
    dy = 1.0/dy;
    let one_over_rho = 1.0/rho;

    let mut count = 0;
    // Main loop
    'time: for k in 1..nt {
        /* some tests to understand how things work
        println!("source: {:?}", fs.get((k)));
        let tmp = 0.5*fs.get((k)).unwrap();
        println!("1/2: {}", tmp);
        if k==10 {
            break 'main;
        }
        */

        // insert source in safe mode
        // *(px.get_mut((isy, isx))).unwrap() += dt*0.5*(fs.get((k))).unwrap();
        // *(py.get_mut((isy, isx))).unwrap() += dt*0.5*(fs.get((k))).unwrap();

        // insert source in unsafe mode
        unsafe{
            *px.uget_mut((isy, isx)) = px.uget((isy, isx)) + dt*0.5*fs.uget((k));
            *py.uget_mut((isy, isx)) = px.uget((isy, isx)) + dt*0.5*fs.uget((k));
        }
        // staggered grid, loop over space
        'space: for i in 0..ny {
            for j in 0..nx {
                unsafe{
                    // update px
                    diffop = (ux.uget((i+1, j)) - ux.uget((i, j)))*dx;
                    pmlop  = qx.uget((i+1, j+1)) * px.uget((i+1, j+1));
                    *px.uget_mut((i+1, j+1)) = px.uget((i+1, j+1)) - dt*(pmlop + rho*vp.uget((i+1, j+1))*diffop);

                    // update py
                    diffop = (uy.uget((i, j+1)) - uy.uget((i, j)))*dy;
                    pmlop  = qy.uget((i+1, j+1)) * py.uget((i+1, j+1));
                    *py.uget_mut((i+1, j+1)) = *py.uget((i+1, j+1)) - dt*(pmlop + rho*vp.uget((i+1, j+1))*diffop);

                    //update ux
                    diffop = (px.uget((i+1, j+1)) - px.uget((i, j+1)) + py.uget((i+1, j+1)) -py.uget((i, j+1)))*dx;
                    pmlop  = 0.5*(qx.uget((i+1, j+1))+qx.uget((i, j+1)))*ux.uget((i,j));
                    *ux.uget_mut((i,j)) = ux.uget((i,j)) - dt*one_over_rho*(pmlop + diffop);

                    //update uy
                    diffop = (px.uget((i+1, j+1)) - px.uget((i+1, j)) + py.uget((i+1, j+1)) -py.uget((i+1, j)))*dy;
                    pmlop  = 0.5*(qy.uget((i+1, j+1))+qy.uget((i+1, j)))*uy.uget((i,j));
                    *uy.uget_mut((i,j)) = uy.uget((i,j)) - dt*one_over_rho*(pmlop + diffop);

                    count += 1;
                }
            };
        };

        unsafe{
            *sfd.uget_mut((k)) = px.uget((iry, irx)) + py.uget((iry, irx));
        }
        //println!("{}", sfd.get((k)).unwrap())
        count += 1;

    }
    println!("Done. (count ! {})", count);

    /* uncomment to output the trace at receiver postion */

    let path = Path::new("trace.txt");
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}",
                           display,
                           Error::description(&why)),
        Ok(file) => file,
    };

    for i in 0..sfd.len(){
        let s = sfd.get((i)).unwrap().to_string()+"\n";
        match file.write_all(s.as_bytes()) {
            Err(why) => {
                panic!("couldn't write to {}: {}", display,
                       Error::description(&why))
            },
            //Ok(_) => println!("successfully wrote to {}", display),
            Ok(_) => (),
        }
    }

}
