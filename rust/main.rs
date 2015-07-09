use std::f64;
use std::thread;

// there is a crate ndarray that would make my life easier !!

fn main() {
    println!("Hello, world! Euh, Thomas !");

    thread::sleep_ms(2000);


    const NX:usize = 300; //columns
    const NY:usize = 300;  // rows
    const NX2:usize = NX+1; //columns
    const NY2:usize = NY+1;  // rows

    let rho:f64  = 1500.0;
    let vp0: f64 = 2000.0;

    // mutable float64 array of size (rows = NY, columns = NZ)
    let mut vp = Box::new([[vp0; NX2]; NY2]);

    for i in 0..NY2 {
        for j in 0..NX2 {
           // mutability is ugly but today I am not going to care too much
           vp[i][j] = vp[i][j]*vp[i][j];
        };
    };



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
    let mut t = Box::new([0f64; NT]);
    for i in 0..NT {
        t[i] = (i as f64)*dt;
    };
    let lsour = tsour/dt;
    let t0 = tsour*1.5;
    let mut tau = Box::new([0f64; NT]);
    for i in 0..NT {
        tau[i] = f64::consts::PI*(t[i]-t0)/t0;
    };
    let a = 4.0;
    let mut fs = Box::new([0f64;NT]);
    for i in 0..NT {
        fs[i] = (1.0 - a*tau[i]*tau[i])*((-2.0*tau[i]*tau[i]).exp())
    };

    // define PML properties
    let npml = 30;
    let pmlfac = 50.0;
    let pmlexp = 2;
/*    let mut qx = Box::new([[0f64; NX2]; NY2]);
    let mut qy = Box::new([[0f64; NX2]; NY2]);
    for i in 0..npml {
        for j in 0..NX2 {
            qy[i][j] = pmlfac*((npml-i) as f64).powi(pmlexp) as f64;
            qy[NY2-i][j] = pmlfac*((npml-i) as f64).powi(pmlexp) as f64;
        }
    };
    for j in 0..npml {
        for i in 0..NY2 {
            qx[i][j] = pmlfac*((npml-i) as f64).powi(pmlexp) as f64;
            qx[i][NX2-j] = pmlfac*((npml-i) as f64).powi(pmlexp) as f64;
        }
    };

    // initialize fields
    let mut px = [[0f64; NX2]; NY2]; // pressure field
    let mut py = [[0f64; NX2]; NY2];
    let mut ux = [[0f64; NX2]; NY2]; // particle velocity field
    let mut uy = [[0f64; NX2]; NY2];

    let mut sfd = [0f64; NT];

    let mut diffop = 0f64;
    let mut pmlop  = 0f64;

    // Main Loop
    for k in 1..NT {
        // inject source
        px[isy][isx] = px[isy][isx] + dt*0.5*fs[k];
        py[isy][isx] = py[isy][isx] + dt*0.5*fs[k];

        //staggered grid, loop over space
        for i in 0..NY {
            for j in 0..NX{
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

            };
        };

        sfd[k] = px[iry][irx] + py[iry][irx];
    };
*/

}
