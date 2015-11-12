#!/bin/sh

#############
## FORTRAN ##
#############
cd fortran

echo "Fortran (ifort -O3): "
ifort fdtd2Ds.f90 -o fdtd.exe -O3
time ./fdtd.exe

echo "Fortran (gfortran): "
gfortran fdtd2Ds.f90 -o fdtd_gfortran.exe -O3
time ./fdtd_gfortran.exe

echo "Fortran (ifort -04 -... agresive)"
ifort fdtd2Ds.f90 -o fdtd_agressive.exe -O4 -check nobounds -xAVX -ftz -shared-intel -mcmodel=medium
time ./fdtd_agressive.exe

cd ..

###########
## JULIA ##
###########
cd julia

echo "julia :"
julia fdtd2ds.jl

#echo "julia (profiling version) :"
# profiler get confused by @simd/@inbounds in
# fdtd2ds.jl so I use fdtd2ds_prof.jl to spot the hotpots
#julia-b24213b893/bin/julia fdtd2ds_prof.jl

cd ..

############
## PYTHON ##
############
cd python

echo "python vectorised :"
python fdtd_wave_equation.py

echo "python loop+numba(jit) :"
python fdtd_wave_equation_jitted.py

cd ..

############
## MATLAB ##
############
cd matlab

echo "matlab :"
matlab -nodesktop -nosplash -r run

cd ..

##########
## RUST ##
##########
cd rust

echo "rust :"
cargo build --release
time ./target/release/fdtd2ds

cd ..

