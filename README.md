# Plane wave scattering by a sphere

Computes the electromagnetic fields scattered by a dielectric or metallic sphere using Mie theory. The sphere is illuminated by an x-polarized plane wave propagating in the positive z-direction. The relative permittivity and permeability of the sphere can both specified as complex values. The background medium is vacuum (hardcoded) but can be changed by editing the corresponding permittivity and permeability in *src/main.cpp*.

This program may be used to compute the fields scattered by sphere on: 1) a 0D, 1D, 2D and 3D rectangular grid, or 2) on the surface of a sphere of points. It may also be used to compute the absorption, extinction and scattering cross sections as well as the radiation pattern. To perform these operations, refer to the sections below.

The infinite series in Mie theory are truncated such that only `n` terms are taken into account, where `n=1` corresponds to dipolar radiation, `n=2` to dipolar + quadrupolar radiation, and so on. Therefore, the higher the value of `n`, the more accurate is the result but the more time it takes to perform the computation. It is the responsibility of the user to select the appropriate  value of `n` and to check for convergence.

The program is written in C++ and is parallelized so as to use all available cores. Once executed, it generates a file containing the relevant data.

Note that the time harmonic function `exp(jwt)` is assumed, which implies that *negative values* of the imaginary parts of the permittivity and the permeability correspond to loss.

## Obtaining and compiling the program

Download the archive and extract it or simply use git:

`git clone https://github.com/kachourim/MieScattering.git`

Then, go into the directory and compile the code:

`cd MieScattering`

`make -C src`

This program has been successfully tested on Debian 10 and Ubuntu 18.04. Note that the program uses fortran libraries and thus requires the GNU fortran compiler `gfortran`. All necessary dependencies may be installed with:

`sudo apt install make g++ gfortran`


(For older versions of `g++`, it may be required to add `-std=c++11` to the variable `FLAGS` in the makefile.)

## Mie computation on a 0D, 1D, 2D or 3D rectangular grid

Run the command: `./Mie 1 lambda n r ers_r ers_i mrs_r mrs_i Lx Ly Lz Cx Cy Cz res filename`, where the parameters are defined by

Parameter	|	Description	
--- | --- 
lambda	 	| Free space wavelength [m]
n 			| Number of Mie coefficients
r			| Sphere radius [m]
ers_r		| Sphere relative permittivity (real part)
ers_i		| Sphere relative permittivity (imag part)
mrs_r		| Sphere relative permeability (real part)
mrs_i		| Sphere relative permeability (imag part)
Lx 			| Grid length along x [m] 			
Ly 			| Grid length along y [m] 				
Lz 			| Grid length along z [m] 
Cx 			| Grid center along x 			
Cy 			| Grid center along y				
Cz 			| Grid center along z
res 		| Grid resolution [m]
filename	| Name of output data file

The figure below compares the result obtained from this program to the same one obtained in COMSOL for the parameters `r = 250e-9 m`, `lambda = 1000e-9 m`, `ers = 5 - 2i` and `mrs = -4 - 6i`.

<img src="/images/2Dgrid.png" width="550">


## Mie computation on the surface of a sphere

Run the command: `./Mie 2 lambda n r ers_r ers_i mrs_r mrs_i R res filename`, where the parameters are defined by

Parameter	|	Description	
--- | --- 
lambda	 	| Free space wavelength [m]
n 			| Number of Mie coefficients
r			| Sphere radius [m]
ers_r		| Sphere relative permittivity (real part)
ers_i		| Sphere relative permittivity (imag part)
mrs_r		| Sphere relative permeability (real part)
mrs_i		| Sphere relative permeability (imag part)
R 			| Points sphere radius [m] 			
res 		| Points sphere angular resolution [°]
filename	| Name of output data file


<img src="/images/Sphere.png" width="550">

## Compute the scattering, absorption and extinction cross sections (fixed wavelength, varying radius)

Run the command: `./Mie 3 lambda n ers_r ers_i mrs_r mrs_i rstart rstop npts filename`, where the parameters are defined by

Parameter	|	Description	
--- | --- 
lambda	 	| Free space wavelength [m]
n 			| Number of Mie coefficients
ers_r		| Sphere relative permittivity (real part)
ers_i		| Sphere relative permittivity (imag part)
mrs_r		| Sphere relative permeability (real part)
mrs_i		| Sphere relative permeability (imag part)
rstart 		| Smallest sphere radius [m] 			
rstop 		| Largest sphere radius [m] 				
npts 		| Number of steps between rstart and rstop [m] 
filename	| Name of output data file

The figure below compares the result obtained from this program to the same one obtained in COMSOL (dots) for the parameters `r = 250e-9 m`, `lambda = 1000e-9 m`, `ers = 5 - 2i` and `mrs = -4 - 6i`.

<img src="/images/CS.png" width="500">


## Compute the scattering, absorption and extinction cross sections (fixed radius, varying wavelength)

Run the command: `./Mie 4 r n ers_r ers_i mrs_r mrs_i wstart wstop npts filename`, where the parameters are defined by

Parameter	|	Description	
--- | --- 
r	 		| Sphere radius [m]
n 			| Number of Mie coefficients
ers_r		| Sphere relative permittivity (real part)
ers_i		| Sphere relative permittivity (imag part)
mrs_r		| Sphere relative permeability (real part)
mrs_i		| Sphere relative permeability (imag part)
wstart 		| Smallest free-space wavelength [m] 			
wstop 		| Largest free-space wavelength [m] 				
npts 		| Number of steps between wstart and wstop [m] 
filename	| Name of output data file

Note that this command does the same as the one of the previous section, except that the radius of the sphere is fixed, while the wavelength changes instead of the opposite.

## Radiation pattern in the xy-, xz- or yz-plane

Run the command: `./Mie 5 lambda n r ers_r ers_i mrs_r mrs_i R res plane filename`, where the parameters are defined by

Parameter	|	Description	
--- | --- 
lambda	 	| Free space wavelength [m]
n 			| Number of Mie coefficients
r			| Sphere radius [m]
ers_r		| Sphere relative permittivity (real part)
ers_i		| Sphere relative permittivity (imag part)
mrs_r		| Sphere relative permeability (real part)
mrs_i		| Sphere relative permeability (imag part)
R 			| Points sphere radius [m] 			
res 		| Points sphere angular resolution [°]
plane 		| 1=xy-plane, 2=xz-plane or 3=yz-plane
filename	| Name of output data file

The figure below compares the result obtained from this program to the same one obtained in COMSOL (dots) for the parameters `r = 250e-9 m`, `lambda = 1000e-9 m`, `ers = 5 - 2i` and `mrs = -4 - 6i`.

<img src="/images/RadiationPattern.png" width="400">


## Python script

A python scripts entitled *RunAndPlot.py* is also provided. It may be used to run the computation and plot the corresponding results. The purpose of this python script is to illustrate how to use the program and not really to serve as a front end for it.

This script requires `numpy` and `matplotlib`, if you don't have them, then run:

`sudo apt install python3-numpy python3-matplotlib`

You may edit this script and then run it using the command `python RunAndPlot.py SEL`, where SEL takes the following values:

SEL	|	Description	
--- | --- 
0	| Run the computation only
1 	| Plot the data only
2	| Run the computation and plot the data


