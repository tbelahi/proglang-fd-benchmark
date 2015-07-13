use std::f64;
use std::thread;
use std::iter;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

fn main() {
    //println!("Hello, world!);

    const NX:usize = 400; //columns
    const NY:usize = 400;  // rows
    const NX2:usize = NX+1; //columns
    const NY2:usize = NY+1;  // rows

    let rho:f64  = 1500.0;
    let vp0: f64 = 2000.0;

    // mutable float64 pseudo array of size (rows = NY2, columns = NX2)
    // vp end up in the heap, no allocation problem
    let row = iter::repeat(vp0).take(NX2).collect();
    let mut vp: Vec<Vec<f64>> = iter::repeat(row).take(NY2).collect();

    for i in 0..NY2 {
        for j in 0..NX2 {
           // mutability is ugly but today I am not going to care too much
           vp[i][j] = vp[i][j]*vp[i][j];
        };
    };

    //println!("{}", vp[NY2-1][NX2-1]);

    //define source and receiver position
    let isx = 200;
    let isy = 100;
    let irx = 200;
    let iry = 300;

    //source and parameters
    const NT: usize = 3000;
    let fc = 1000.0;
    let dx: f64  = vp0/(30.0*fc);
    let dy = dx;
    let dt: f64 = dx/vp0*1.0/(f64::consts::SQRT_2*2.0*f64::consts::PI);

    //println!("{} {} {} {} {}", nt, fc, dx, dy, dt);

    // define source function in time
    let tsour = 1.0/fc;
    let mut t: Vec<f64> = iter::repeat(0f64).take(NT).collect();
    for i in 0..NT {
        t[i] = (i as f64)*dt;
    };
    //let lsour = tsour/dt;
    let t0 = tsour*1.5;
    let mut tau: Vec<f64> = iter::repeat(0f64).take(NT).collect();
    for i in 0..NT {
        tau[i] = f64::consts::PI*(t[i]-t0)/t0;
    };
    let a = 4.0;
    let mut fs: Vec<f64> = iter::repeat(0f64).take(NT).collect();
    for i in 0..NT {
        fs[i] = (1.0 - a*tau[i]*tau[i])*((-2.0*tau[i]*tau[i]).exp())
    };

    // define PML properties
    let npml = 30;
    let pmlfac: f64 = 50.0;
    let pmlexp: i32 = 2;
    let qxrow = iter::repeat(0f64).take(NX2).collect();
    let qyrow = iter::repeat(0f64).take(NX2).collect();
    let mut qx: Vec<Vec<f64>> = iter::repeat(qxrow).take(NY2).collect();
    let mut qy: Vec<Vec<f64>> = iter::repeat(qyrow).take(NY2).collect();

    let mut i= 0;
    let mut j= 0;
    while i<npml {
        while j<NX2 {
            qy[i][j] = pmlfac*(((npml-i) as f64 ).powi(pmlexp));
            qy[NY2-i-1][j] = pmlfac*(((npml-i) as f64).powi(pmlexp));
            i += 1;
            j += 1;
        }
    };
    let mut i= 0;
    let mut j= 0;
    while i<NY2 {
        while j<npml {
            qx[i][NX2-j-1] = pmlfac*(((npml-j) as f64).powi(pmlexp));
            qx[i][j] = pmlfac*(((npml-j) as f64).powi(pmlexp));
            i += 1;
            j += 1;
        }
    };

    // initialize fields

    let pxrow = iter::repeat(0f64).take(NX2).collect();
    let uxrow = iter::repeat(0f64).take(NX2).collect();
    let pyrow = iter::repeat(0f64).take(NX2).collect();
    let uyrow = iter::repeat(0f64).take(NX2).collect();
    let mut px: Vec<Vec<f64>> = iter::repeat(pxrow).take(NY2).collect();
    let mut py: Vec<Vec<f64>> = iter::repeat(pyrow).take(NY2).collect();
    let mut ux: Vec<Vec<f64>> = iter::repeat(uxrow).take(NY2).collect();
    let mut uy: Vec<Vec<f64>> = iter::repeat(uyrow).take(NY2).collect();

    let mut sfd: Vec<f64> = iter::repeat(0f64).take(NT).collect();

    let mut diffop = 0f64;
    let mut pmlop  = 0f64;

    // Main Loop
    let mut k = 1;
    while k<NT {

        //if k%100 == 0 { println!("timestep : {}", k)}
        // inject source
        px[isy][isx] = px[isy][isx] + dt*0.5*fs[k];
        py[isy][isx] = py[isy][isx] + dt*0.5*fs[k];

        //staggered grid, loop over space
        let mut i = 0;
        let mut j = 0;
        while i<NY {
            while j<NX{
                //update px
                diffop = (ux[i+1][j] - ux[i][j])/dx;
                pmlop  = qx[i+1][j+1]*px[i+1][j+1];
                px[i+1][j+1] = px[i+1][j+1] - dt*(pmlop + rho*vp[i+1][j+1]*diffop);

                //update py
                diffop = (uy[i][j+1] - uy[i][j])/dy;
                pmlop  = qy[i+1][j+1]*py[i+1][j+1];
                py[i+1][j+1] = py[i+1][j+1] - dt*(pmlop + rho*vp[i+1][j+1]*diffop);

                //update ux
                diffop = (px[i+1][j+1] - px[i][j+1] + py[i+1][j+1] -py[i][j+1])/dx;
                pmlop  = 0.5*(qx[i+1][j+1]+qx[i][j+1])*ux[i][j];
                ux[i][j] = ux[i][j] - dt/rho*(pmlop + diffop);

                //update uy
                diffop = (px[i+1][j+1] - px[i+1][j] + py[i+1][j+1] -py[i+1][j])/dy;
                pmlop  = 0.5*(qy[i+1][j+1]+qy[i][j+1])*uy[i][j];
                uy[i][j] = uy[i][j] - dt/rho*(pmlop + diffop);

                i += 1;
                j += 1;

            };
        };

        sfd[k] = px[iry][irx] + py[iry][irx];
        k += 1;
    };

    /* uncomment to output the trace at receiver postion

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
        let s = sfd[i].to_string()+"\n";
        match file.write_all(s.as_bytes()) {
            Err(why) => {
                panic!("couldn't write to {}: {}", display,
                                                Error::description(&why))
            },
            Ok(_) => println!("successfully wrote to {}", display),
        }
    }
    */

}
