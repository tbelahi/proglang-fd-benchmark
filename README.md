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
 - julia
 - matlab
 - cargo

#### Results
See by yourself. There is no big surprise. The Julia implementation is faster (on my machine) than Fortran90 compiled with gfortran thanks to the use of SIMD instructions. This is the main taking home point I got from this little exercise. Rust implementation is working fine. Maybe it can be optimized. Already figuring out that it should be compiled with `cargo build --release` made a tremendous improvement. Overall Fortran + ifort is the speed winner.

update : it seems that between the first version of Julia I used (v0.3.3) and the version I used recently (v0.3.9) a lot of improvement has been made on speeding up the use of slices of arrays. My vectorized code used to be slower than the code with unfolded loop combined with simd instruction. It is not the case anymore. This is a strong case for Julia. Great performance out of the box without so much effort put into it.
