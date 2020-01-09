#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

#include "Mie_Functions.h"

#define pi 3.1415926535897932384626433
const complex<double> j(0.0,1.0);


int main(int argc, char **argv)
{
	int n,  i, SEL, npts, plane;
    double r = 0.0, lam0 = 0.0, R, rstart, rstop, wstart, wstop, IN, ers_r, ers_i = 0.0, mrs_r, mrs_i = 0.0;
    double Lx, Ly, Lz, Cx, Cy, Cz, res;
	char *fname = NULL;
    complex<double>  *F, ers, mrs;
	vector<double> rv, wv, FF;
	vector<complex<double>> CSsca, CSext;
    vector<vector<double>> pts;
	vector<vector<complex<double>>> vF;

	
	
	complex<double> erb(1,0); 		// permittivity of background medium 
	complex<double> mrb(1,0);		// permeabilitty of background medium 
	
	
	
	if (argc == 1)
	{
		printf("Choose what you what to do by specifying the following arguments:\n");
		printf("-----------------------------------------------------------------\n");
		printf("Mie computation on a rectangular grid: 1 lam n r ers_r ers_i mrs_r mrs_i Lx Ly Lz Cx Cy Cz res filename\n");
		printf("Mie computation on the surface of a sphere: 2 lam n r ers_r ers_i mrs_r mrs_i R angular res filename\n");
		printf("Cross sections (fixed wavelength): 3 lam n ers_r ers_i mrs_r mrs_i rstart rstop npts filename\n");
		printf("Cross sections (fixed radius): 4 r n ers_r ers_i mrs_r mrs_i wstart wstop npts filename\n");
		printf("Radiation pattern in the xy-, xz-, or yz-plane: 5 lam n r ers_r ers_i mrs_r mrs_i R angular res plane filename\n");
		printf("-----------------------------------------------------------------\n");
		return 1;
	}
	
	
    // --------------------------------------------------
	// Define parameters
	// --------------------------------------------------
	
	SEL			= int(atof(argv[1]));				// Select computation

	if (SEL == 1)
	{						
		lam0	= atof(argv[2]);   					// Free space wavelength
		n   	= atof(argv[3]);          			// Number of Mie coefficients
		r		= atof(argv[4]);         			// Sphere Radius
		ers_r	= atof(argv[5]);         			// Sphere relative permittivity (real part)
		ers_i	= atof(argv[6]);         			// Sphere relative permittivity (imag. part)
		mrs_r	= atof(argv[7]);         			// Sphere relative permeabilitty (real part)
		mrs_i	= atof(argv[8]);         			// Sphere relative permeabilitty (imag. part)
		Lx  	= atof(argv[9]);					// Grid length along x
		Ly  	= atof(argv[10]);					// Grid length along y
		Lz  	= atof(argv[11]);					// Grid length along z
		Cx  	= atof(argv[12]);					// Grid center along x
		Cy  	= atof(argv[13]);					// Grid center along y
		Cz  	= atof(argv[14]);					// Grid center along z
		res 	= atof(argv[15]);					// Grid resolution [m]
		fname 	= argv[16];							// File name
	}
	else if (SEL == 2)
	{
		lam0	= atof(argv[2]);   					// Free space wavelength
		n   	= atof(argv[3]);          			// Number of Mie coefficients
		r		= atof(argv[4]);         			// Sphere Radius
		ers_r	= atof(argv[5]);         			// Sphere relative permittivity (real part)
		ers_i	= atof(argv[6]);         			// Sphere relative permittivity (imag. part)
		mrs_r	= atof(argv[7]);         			// Sphere relative permeabilitty (real part)
		mrs_i	= atof(argv[8]);         			// Sphere relative permeabilitty (imag. part)		
		R  		= atof(argv[9]);					// Points sphere radius
		res  	= atof(argv[10]);					// Angular resolution [°]
		fname 	= argv[11];							// File name
	}
	else if (SEL == 3)
	{
		lam0	= atof(argv[2]);   					// Free space wavelength
		n   	= atof(argv[3]);          			// Number of Mie coefficients
		ers_r	= atof(argv[4]);         			// Sphere relative permittivity (real part)
		ers_i	= atof(argv[5]);         			// Sphere relative permittivity (imag. part)
		mrs_r	= atof(argv[6]);         			// Sphere relative permeabilitty (real part)
		mrs_i	= atof(argv[7]);         			// Sphere relative permeabilitty (imag. part)
		rstart  = atof(argv[8]);					// Smallest sphere radius [m]
		rstop  	= atof(argv[9]);					// Largest sphere radius [m]
		npts 	= atof(argv[10]);					// Number of steps between rstart and rstop
		fname 	= argv[11];							// File name
	}
	else if (SEL == 4)
	{
		r		= atof(argv[2]);   					// Radius [m]
		n   	= atof(argv[3]);          			// number of Mie coefficients
		ers_r	= atof(argv[4]);         			// Sphere relative permittivity (real part)
		ers_i	= atof(argv[5]);         			// Sphere relative permittivity (imag. part)
		mrs_r	= atof(argv[6]);         			// Sphere relative permeabilitty (real part)
		mrs_i	= atof(argv[7]);         			// Sphere relative permeabilitty (imag. part)
		wstart  = atof(argv[8]);					// Smallest free-space wavelength [m]
		wstop  	= atof(argv[9]);					// Largest free-space wavelength [m]
		npts 	= atof(argv[10]);					// Number of steps between wstart and wstop
		fname 	= argv[11];							// File name
	}
	else if (SEL == 5)
	{
		lam0	= atof(argv[2]);   					// Free space wavelength
		n   	= atof(argv[3]);          			// Number of Mie coefficients
		r		= atof(argv[4]);         			// Sphere Radius
		ers_r	= atof(argv[5]);         			// Sphere relative permittivity (real part)
		ers_i	= atof(argv[6]);         			// Sphere relative permittivity (imag. part)
		mrs_r	= atof(argv[7]);         			// Sphere relative permeabilitty (real part)
		mrs_i	= atof(argv[8]);         			// Sphere relative permeabilitty (imag. part)
		R  		= atof(argv[9]);					// Points sphere radius
		res  	= atof(argv[10]);					// Angular resolution [°]
		plane  	= int(atof(argv[11]));				// Plane (1:xy, 2:xz, 3:yz) 
		fname 	= argv[12];							// File name
	}
	
	ers = ers_r + j*ers_i;
	mrs = mrs_r + j*mrs_i;
	
    // --------------------------------------------------
	// Generate points
	// --------------------------------------------------
    
	// IN: incident + scattered fields (IN = 1), scattered fields only (IN = 0)
	
    if (SEL == 1) {pts = RectPtsGen(Lx, Ly, Lz, Cx, Cy, Cz, res); IN = 1.0;}
	if (SEL == 2) {pts = SpherePtsGen(R, res); IN = 1.0;}
	if (SEL == 3) rv  = linspace(rstart, rstop, npts);
	if (SEL == 4) wv  = linspace(wstart, wstop, npts);
	if (SEL == 5) {pts = CirclePtsGen(plane, R, res); IN = 0.0;}
    
      
	if (SEL == 1 || SEL == 2 || SEL == 5)
	{
		printf("Number of points to be evaluated: %d\n", int(pts.size()));
    
		for(i = 0; i < int(pts.size()); i++) vF.push_back({0,0,0,0,0,0});
		
		// Computes fields (in parallel)
		#pragma omp parallel for
		for(i = 0; i < int(pts.size()); i++)
		{
			complex<double> F0[6];
			  
			F = Mie_fields(IN, pts[i][0], pts[i][1], pts[i][2], n, r, lam0, erb, mrb, ers, mrs, F0);
			
			vF[i][0] = F[0]; vF[i][1] = F[1]; vF[i][2] = F[2];
			vF[i][3] = F[3]; vF[i][4] = F[4]; vF[i][5] = F[5];      
		}
		
		// Save to file
		FieldsToFileXYZ(fname, pts, vF, r, lam0, n, erb, mrb, ers, mrs);
	}
	
	if (SEL == 3)
	{
		for(i = 0; i < int(rv.size()); i++) {CSsca.push_back(0); CSext.push_back(0);}
		
		// Computes scattering cross section (in parallel)
		#pragma omp parallel for
		for(i = 0; i < int(rv.size()); i++)	
		{
			CSsca[i] = Csca(n, rv[i], lam0, erb, mrb, ers, mrs);
			CSext[i] = Cext(n, rv[i], lam0, erb, mrb, ers, mrs);
		}
		
		// Save to file		
		CsToFileR(fname, rv, CSsca, CSext, lam0, n, erb, mrb, ers, mrs);
	}
	
	
	if (SEL == 4)
	{	
		for(i = 0; i < int(wv.size()); i++) {CSsca.push_back(0); CSext.push_back(0);}
		
		// Computes scattering cross section (in parallel)
		#pragma omp parallel for
		for(i = 0; i < int(wv.size()); i++)	
		{
			CSsca[i] = Csca(n, r, wv[i], erb, mrb, ers, mrs);
			CSext[i] = Cext(n, r, wv[i], erb, mrb, ers, mrs);
		}
		// Save to file
		CsToFileW(fname, wv, CSsca, CSext, r, n, erb, mrb, ers, mrs);
	}
	
	return 0;
}


