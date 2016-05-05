Speed test
==========

### Objective
The goal of this project is to develop a little benchmark of various programming languages when it comes to model the acoustic wave equation, I use the staggered grid approach to model the acoustic wave equation in a finite-differences time domain scheme.

### Languages
As of now, the programming languages that are compared are:

 - Fortran
 - Julia
 - Python
 - Matlab
 - Rust 

### Usage
Just run the `comparison.sh` script and look at the output.

### requirements
several tools need to be on your PATH for comparison.sh to be able to run:

 - ifort
 - gfortran
 - python + numpy + matplotlib + numba
 - julia (+ PyCall)
 - matlab
 - cargo

#### Results
See by yourself. There is no big surprise. The Julia implementation is faster (on my machine) than Fortran90 compiled with gfortran thanks to the use of SIMD instructions, I guess it also sets appart Julia and "jitted" python. This is the main taking home point I got from this little exercise. Rust implementation is working fine. Maybe it can be optimized and don't forget that it should be compiled with `cargo build --release` which provides a tremendous improvement compared to `cargo  build`. Julia and Rust are pretty close performance wise, which should not be surprising since they both target LLVM.However, writting Julia (scientific) code is much easier and straightforward than writting arrays operiation code in Rust. But Rust reached version 1.0 and is stable now, Julia should really focus on getting to the 1.0 stable API version and will then become a very compelling choice. Overall Fortran + ifort is the speed winner (on intel hardware) and the complexity of writting Fortran code is not higher than writting Julia or Python code. 

One other important message is: avoid sliced arrays with Julia. They are slow, more than 2x slower than similar matlab code and 3.5x slower than equivalent python code.

