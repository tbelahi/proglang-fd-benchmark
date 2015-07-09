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
 - Rust (work in progress)

### Usage
Just run the `comparison.sh` script and look at the output.

#### Results
See by yourself. There is no big surprise. The Julia implementation is faster (on my machine) than Fortran90 compiled with gfortran thanks to the use of SIMD instructions. This is the main taking home point I got from this little exercise. My poor knowledge of Rust lead to a super slow implementation. I guess the way I construct a 2D array with collection of vecs does not help at all.
