#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

#define pi 3.1415926535897932384626433
#define eta0 376.7303134617707

const complex<double> j(0.0,1.0);


extern "C" 
{
    void zbesj_( double *zr, double *zi, double *fnu, int *kode, int *n,
                 double *Jr, double *Ji, int *nz, int *ierr );
	
	void zbesy_(double *zr, double *zi, double *fnu, int *kode, int *n,
                 double *Jr, double *Ji, int *nz, double *Wr, double *Wi, int *ierr);
	
	void dxlegf_(double *dnu1, int *nudiff, int *mu1, int *mu2, double *theta, int *id, double *pqa, double *ipqa, int* ierr);
}


int ConfigReader(const char *filename, double D[])
{
	int i = 0;
    char buf[100]; 

    FILE* ptr = fopen(filename,"r"); 
    if (ptr==NULL) 
    { 
        printf("Config file not found!\n"); 
        return 1; 
    } 
  

	while(fgets(buf, 100, (FILE*) ptr)) {
		if (buf[0] != '#' && buf[0] != '\n')
		{
	    	D[i] = atof(buf);
			i++;
		}
		if(i==7) break;
	}

	fclose(ptr);
  
    return 0; 
}


vector<vector<double>> RectPtsGen(double Lx, double Ly, double Lz, double Cx, double Cy, double Cz, double res)
{
	double x, y, z;
	vector<vector<double>> pts;
    
	for(z = -Lz/2.0 + Cz; z <= Lz/2.0 + Cz; z += res)
	{
		for(y = -Ly/2.0 + Cy; y <= Ly/2.0 + Cy; y += res)
		{
			for(x = -Lx/2.0 + Cx; x <= Lx/2.0 + Cx; x += res)
			{
				pts.push_back({x,y,z});
			}
		}
	}

	return pts;
}



vector<vector<double>> SpherePtsGen(double R, double res)
{
	double  phi, theta;
	vector<vector<double>> pts;
	
	res *= pi/180.0;

	pts.push_back({0.0,0.0,R});
	pts.push_back({0.0,0.0,-R});

	for(theta = res; theta < pi; theta += res)
	{
		for(phi = 0; phi < 2*pi; phi += res)
		{
			pts.push_back({R*sin(theta)*cos(phi),R*sin(theta)*sin(phi),R*cos(theta)});
		}
	}

	return pts;
}


vector<vector<double>> CirclePtsGen(int plane, double R, double res)
{
	double  phi, theta;
	vector<vector<double>> pts;
	
	res *= pi/180.0;

	if (plane == 1)	// xy-plane
	{
		theta = pi/2.0;
		
		for(phi = 0; phi < 2.0*pi; phi += res)	pts.push_back({R*sin(theta)*cos(phi),R*sin(theta)*sin(phi),R*cos(theta)});
	}
	else if (plane == 2)	// xz-plane
	{
		phi = 0.0;
		
		for(theta = 0; theta < 2.0*pi; theta += res)	pts.push_back({R*sin(theta)*cos(phi),R*sin(theta)*sin(phi),R*cos(theta)});
	}
	else if (plane == 3)	// yz-plane
	{
		phi = pi/2.0;
		
		for(theta = 0; theta < 2.0*pi; theta += res)	pts.push_back({R*sin(theta)*cos(phi),R*sin(theta)*sin(phi),R*cos(theta)});
	}

	return pts;
}



complex<double> J1(double fnu, complex<double> z)
{
	int   n, nz, ierr, kode;
    double  zr, zi, Jr, Ji;
	
	n = 1;
	kode = 1.0; 		// leave as is
    fnu += 0.5; 		// order
    zr = real(z); zi = imag(z); // real and imag of input

    zbesj_(&zr, &zi, &fnu, &kode, &n, &Jr, &Ji, &nz, &ierr);
	
	return sqrt(pi*z/2.0)*(Jr + j*Ji);
}

double Y1(double fnu, double x, double k)
{
	int   n, nz, ierr, kode;
    double  zr, zi, Yr, Yi, Wr, Wi;
	
	n = 1;
	kode = 1.0; 		// leave as is
    fnu += 0.5; 		// order
    zr = k*x; zi = 0.0; // real and imag of input

	zbesy_(&zr, &zi, &fnu, &kode, &n, &Yr, &Yi, &nz, &Wr, &Wi, &ierr);
	
	return sqrt(pi*k*x/2.0)*Yr;
}



complex<double> J1d(double fnu, complex<double> z)
{
	int   n, nz, ierr, kode;
    double  fnu1, fnu2, zr, zi, Jr1, Ji1, Ji2, Jr2;
	
	n = 1;
	kode = 1.0; 		// leave as is
    fnu1 = fnu + 0.5; 	// order
	fnu2 = fnu + 1.5; 	// order
    zr = real(z); zi = imag(z); // real and imag of input

    zbesj_(&zr, &zi, &fnu1, &kode, &n, &Jr1, &Ji1, &nz, &ierr);
	zbesj_(&zr, &zi, &fnu2, &kode, &n, &Jr2, &Ji2, &nz, &ierr);
		
	return sqrt(pi/(z*2.0))*((Jr1+j*Ji1) + fnu*(Jr1+j*Ji1) - z*(Jr2+j*Ji2));
}

double Y1d(double fnu, double x, double k)
{
	int   n, nz, ierr, kode;
    double  fnu1, fnu2, zr, zi, Yr1, Yi, Yr2, Wr, Wi;
	
	n = 1;
	kode = 1.0; 		// leave as is
    fnu1 = fnu + 0.5; 	// order
	fnu2 = fnu + 1.5; 	// order
    zr = k*x; zi = 0.0; // real and imag of input

	zbesy_(&zr, &zi, &fnu1, &kode, &n, &Yr1, &Yi, &nz, &Wr, &Wi, &ierr);
	zbesy_(&zr, &zi, &fnu2, &kode, &n, &Yr2, &Yi, &nz, &Wr, &Wi, &ierr);
		
	return sqrt(pi/(k*x*2.0))*(Yr1 + fnu*Yr1 - k*x*Yr2);
}

complex<double> H1(double fnu, complex<double> z)
{
	int   n, nz, ierr, kode;
    double  zr, zi, Jr, Ji, Yr, Yi, Wr, Wi;


	n = 1;
	kode = 1.0; 		// leave as is
    fnu += 0.5; 		// order
    zr = real(z); zi = imag(z); // real and imag of input

    zbesj_(&zr, &zi, &fnu, &kode, &n, &Jr, &Ji, &nz, &ierr);	
	zbesy_(&zr, &zi, &fnu, &kode, &n, &Yr, &Yi, &nz, &Wr, &Wi, &ierr);
	
	return sqrt(pi*z/2.0)*((Jr+j*Ji) + j*(Yr+j*Yi));
}

complex<double> H2(double fnu, complex<double> z)
{
	int   n, nz, ierr, kode;
    double  zr, zi, Jr, Ji, Yr, Yi, Wr, Wi;


	n = 1;
	kode = 1.0; 		// leave as is
    fnu += 0.5; 		// order
    zr = real(z); zi = imag(z); // real and imag of input

    zbesj_(&zr, &zi, &fnu, &kode, &n, &Jr, &Ji, &nz, &ierr);	
	zbesy_(&zr, &zi, &fnu, &kode, &n, &Yr, &Yi, &nz, &Wr, &Wi, &ierr);
	
	return sqrt(pi*z/2.0)*((Jr+j*Ji) - j*(Yr+j*Yi));
}

complex<double> H1d(double fnu, complex<double> z)
{
	int   n, nz, ierr, kode;
    double  zr, zi, Jr1, Jr2, Ji1, Ji2, Yr1, Yr2, Yi1, Yi2, Wr, Wi, fnu1, fnu2;

	n = 1;
	kode = 1.0; 		// leave as is
    fnu1 = fnu + 0.5; 	// order
	fnu2 = fnu + 1.5; 	// order
    zr = real(z); zi = imag(z); // real and imag of input

	zbesj_(&zr, &zi, &fnu1, &kode, &n, &Jr1, &Ji1, &nz, &ierr);
    zbesj_(&zr, &zi, &fnu2, &kode, &n, &Jr2, &Ji2, &nz, &ierr);	
	zbesy_(&zr, &zi, &fnu1, &kode, &n, &Yr1, &Yi1, &nz, &Wr, &Wi, &ierr);
	zbesy_(&zr, &zi, &fnu2, &kode, &n, &Yr2, &Yi2, &nz, &Wr, &Wi, &ierr);
		
	return sqrt(pi/(z*2.0))*((fnu+1.0)*(Jr1+j*Ji1) - z*(Jr2+j*Ji2) + j*(fnu+1.0)*(Yr1+j*Yi1) - j*z*(Yr2+j*Yi2));
}

complex<double> H2d(double fnu, complex<double> z)
{
	int   n, nz, ierr, kode;
    double  zr, zi, Jr1, Jr2, Ji1, Ji2, Yr1, Yr2, Yi1, Yi2, Wr, Wi, fnu1, fnu2;

	n = 1;
	kode = 1.0; 		// leave as is
    fnu1 = fnu + 0.5; 	// order
	fnu2 = fnu + 1.5; 	// order
    zr = real(z); zi = imag(z); // real and imag of input

	zbesj_(&zr, &zi, &fnu1, &kode, &n, &Jr1, &Ji1, &nz, &ierr);
    zbesj_(&zr, &zi, &fnu2, &kode, &n, &Jr2, &Ji2, &nz, &ierr);	
	zbesy_(&zr, &zi, &fnu1, &kode, &n, &Yr1, &Yi1, &nz, &Wr, &Wi, &ierr);
	zbesy_(&zr, &zi, &fnu2, &kode, &n, &Yr2, &Yi2, &nz, &Wr, &Wi, &ierr);
		
	return sqrt(pi/(z*2.0))*((fnu+1.0)*(Jr1+j*Ji1) - z*(Jr2+j*Ji2) - j*(fnu+1.0)*(Yr1+j*Yi1) + j*z*(Yr2+j*Yi2));
}




double P0(double dnu1, double theta)
{
	double pqa, ipqa;
	int  nudiff, mu1, mu2, id, ierr;
	
	mu1   = 0;  	// value of m
	
	// to compute the positive legendre function
	id     = 3;
	mu2    = mu1;
	nudiff = 0;
	
	
	if (theta > pi/2 && fmod(dnu1,2.0) == 0)
	{
		theta = pi - theta;
		dxlegf_(&dnu1, &nudiff, &mu1, &mu2, &theta, &id, &pqa, &ipqa, &ierr);
		return pqa;
	}
	else if (theta > pi/2 && fmod(dnu1,2.0) == 1)
	{
		theta = pi - theta;
		dxlegf_(&dnu1, &nudiff, &mu1, &mu2, &theta, &id, &pqa, &ipqa, &ierr);
		return -pqa;
	}
	else
	{
		dxlegf_(&dnu1, &nudiff, &mu1, &mu2, &theta, &id, &pqa, &ipqa, &ierr);
		return pqa;	
	}
}



double P1(double dnu1, double theta)
{
	double pqa, ipqa;
	int  nudiff, mu1, mu2, id, ierr;
	
	mu1    = 1;  	// value of m
	
	// to compute the positive legendre function
	id     = 3;
	mu2    = mu1;
	nudiff = 0;
	
	
	if (theta > pi/2 && fmod(dnu1,2.0) == 0)
	{
		theta = pi - theta;
		dxlegf_(&dnu1, &nudiff, &mu1, &mu2, &theta, &id, &pqa, &ipqa, &ierr);
		return -pqa;
	}
	else if (theta > pi/2 && fmod(dnu1,2.0) == 1)
	{
		theta = pi - theta;
		dxlegf_(&dnu1, &nudiff, &mu1, &mu2, &theta, &id, &pqa, &ipqa, &ierr);
		return pqa;
	}
	else
	{
		dxlegf_(&dnu1, &nudiff, &mu1, &mu2, &theta, &id, &pqa, &ipqa, &ierr);
		return pqa;	
	}
}



double P1d(double n, double theta)
{
	return -P1(n,theta)/tan(theta) - n*(n+1.0)*P0(n,theta);
}




complex<double> Csca(double N, double r, double lam0, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs)
{
	// From Bohren and Huffman
	double n;
	complex<double> CS(0,0), a_n, b_n, m, k0;
	
	erb = conj(erb);
	mrb = conj(mrb);
	ers = conj(ers);
	mrs = conj(mrs);
	
	m = sqrt(ers*mrs)/sqrt(erb*mrb);
	k0 = 2.0*pi/lam0*sqrt(erb*mrb);

	for(n = 1.0; n < N+1.0; n++)
	{
		a_n = (mrb*m*J1(n,m*k0*r)*J1d(n,k0*r) - mrs*J1(n,k0*r)*J1d(n,m*k0*r))/(mrb*m*J1(n,m*k0*r)*H1d(n,k0*r) - mrs*H1(n,k0*r)*J1d(n,m*k0*r));
		b_n = (mrs*J1(n,m*k0*r)*J1d(n,k0*r) - mrb*m*J1(n,k0*r)*J1d(n,m*k0*r))/(mrs*J1(n,m*k0*r)*H1d(n,k0*r) - mrb*m*H1(n,k0*r)*J1d(n,m*k0*r));

		CS += (2*n + 1)*(pow(abs(a_n),2) + pow(abs(b_n),2));
	}
	
	return CS*2.0*pi/pow(k0,2);
}



complex<double> Cext(double N, double r, double lam0, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs)
{
	// From Bohren and Huffman
	double n;
	complex<double> CS(0,0), a_n, b_n, m, k0;
	
	erb = conj(erb);
	mrb = conj(mrb);
	ers = conj(ers);
	mrs = conj(mrs);
	
	m = sqrt(ers*mrs)/sqrt(erb*mrb);
	k0 = 2.0*pi/lam0*sqrt(erb*mrb);

	for(n = 1.0; n < N+1.0; n++)
	{
		a_n = (mrb*m*J1(n,m*k0*r)*J1d(n,k0*r) - mrs*J1(n,k0*r)*J1d(n,m*k0*r))/(mrb*m*J1(n,m*k0*r)*H1d(n,k0*r) - mrs*H1(n,k0*r)*J1d(n,m*k0*r));
		b_n = (mrs*J1(n,m*k0*r)*J1d(n,k0*r) - mrb*m*J1(n,k0*r)*J1d(n,m*k0*r))/(mrs*J1(n,m*k0*r)*H1d(n,k0*r) - mrb*m*H1(n,k0*r)*J1d(n,m*k0*r));

		CS += (2*n + 1)*real(a_n + b_n);
	}
	
	return CS*2.0*pi/pow(k0,2);
}



complex<double> CscaPEC(double N, double r, double lam0, complex<double> erb, complex<double> mrb)
{
	double n;
	complex<double> CS(0,0), k0;
	
	k0 = 2.0*pi/lam0*sqrt(erb*mrb);

	for(n = 1.0; n < N+1.0; n++)	CS += pow(-1,n)*(2*n + 1)/(H2(n,k0*r)*H2d(n,k0*r));
	
	return pow(abs(CS),2)*pi/pow(k0,2);
}



complex<double> *Mie_fields(double IN, double x, double y, double z, double N, double r, double lam0, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs, complex<double> Fields[])
{
	double n, p1, p1d, theta, phi, R;
	complex<double>  A1, A2, A3, A4, B1, B2, a_n, b_n, c_n, d_n, e_n, kd, k0, alpha, eta, etas;
	complex<double> Eri, Ethetai, Ephii, Hri, Hthetai, Hphii;
    complex<double> Ers, Ethetas, Ephis, Hrs, Hthetas, Hphis;
    complex<double> Er, Etheta, Ephi, Hr, Htheta, Hphi;


	complex<double> Er0i(0,0), Etheta1i(0,0),  Etheta2i(0,0), Ephi1i(0,0),  Ephi2i(0,0);
	complex<double> Hr0i(0,0), Htheta1i(0,0),  Htheta2i(0,0), Hphi1i(0,0),  Hphi2i(0,0);
	complex<double> Er0s(0,0), Etheta1s(0,0), Etheta2s(0,0), Ephi1s(0,0), Ephi2s(0,0);
	complex<double> Hr0s(0,0), Htheta1s(0,0), Htheta2s(0,0), Hphi1s(0,0), Hphi2s(0,0);

    complex<double> E0(1,0);
	
	alpha 	= sqrt(ers/mrs);
	eta     = eta0*sqrt(mrb/erb);
	etas    = eta0*sqrt(mrs/ers);
	k0 		= 2.0*pi/lam0*sqrt(erb*mrb);
	kd		= 2.0*pi/lam0*sqrt(ers*mrs);
	R 		= sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	phi   	= -atan2(y,-x) + pi;
	theta 	= acos(z/R);


    if (theta == 0)		theta = 1e-20;
    if (theta == pi)	theta = pi - pi*1e-3;

	if (R < r)
		{				
			for(n = 1.0; n < N+1.0; n++)
			{
				// Define Mie coefficients 
				a_n = pow(j,-n)*(2*n + 1)/(n*(n + 1));
				
				p1  = P1(n,theta);
				p1d = P1d(n,theta);
						
				d_n	= -j*alpha*mrs/(alpha*H2d(n,k0*r)*J1(n,kd*r) - H2(n,k0*r)*J1d(n,kd*r))*a_n;
				e_n	= j*(kd/k0)*a_n/(alpha*H2(n,k0*r)*J1d(n,kd*r) - H2d(n,k0*r)*J1(n,kd*r));
				
				// Scattered fields 
				Er0s     	+=  n*(n+1)*d_n*J1(n,kd*R)*p1;
				Etheta1s 	+=  e_n*J1(n,kd*R)		*p1;
				Etheta2s 	+=  d_n*J1d(n,kd*R)		*p1d;
				Ephi1s   	+=  e_n*J1(n,kd*R) 		*p1d;
				Ephi2s   	+=  d_n*J1d(n,kd*R)		*p1;

				Hr0s     	+=  n*(n+1)*e_n*J1(n,kd*R)*p1;
				Htheta1s 	+=  d_n*J1(n,kd*R)		*p1;
				Htheta2s 	+=  e_n*J1d(n,kd*R)		*p1d;
				Hphi1s   	+=  d_n*J1(n,kd*R) 		*p1d;
				Hphi2s   	+=  e_n*J1d(n,kd*R)		*p1;		
			}
		
			// total fields 
			Er     	= + E0*cos(phi)/(j*pow(kd*R,2))*Er0s;
			Etheta 	= - E0/(kd*R*sin(theta))*cos(phi)*Etheta1s + E0/(j*kd*R)*cos(phi)*Etheta2s;
			Ephi   	= + E0/(kd*R)*sin(phi)*Ephi1s - E0/(j*kd*R*sin(theta))*sin(phi)*Ephi2s;

			Hr     	= + E0*sin(phi)/(j*etas*pow(kd*R,2))*Hr0s;
			Htheta	= - E0/(etas*kd*R*sin(theta))*sin(phi)*Htheta1s + E0/(j*etas*kd*R)*sin(phi)*Htheta2s;
			Hphi   	= - E0/(etas*kd*R)*cos(phi)*Hphi1s + E0/(j*etas*kd*R*sin(theta))*cos(phi)*Hphi2s;
			
			// converting to cartesian 
			Fields[0]  = Er*sin(theta)*cos(phi) +	Etheta*cos(theta)*cos(phi)	-	Ephi*sin(phi);
			Fields[1]  = Er*sin(theta)*sin(phi) +	Etheta*cos(theta)*sin(phi)	+	Ephi*cos(phi);
			Fields[2]  = Er*cos(theta)			- 	Etheta*sin(theta);

			Fields[3]  = Hr*sin(theta)*cos(phi) + 	Htheta*cos(theta)*cos(phi) 	-	Hphi*sin(phi);
			Fields[4]  = Hr*sin(theta)*sin(phi) + 	Htheta*cos(theta)*sin(phi) 	+	Hphi*cos(phi);
			Fields[5]  = Hr*cos(theta)          - 	Htheta*sin(theta);
		}
		else
		{
			for(n = 1.0; n < N+1.0; n++)
			{
				// Define Mie coefficients 
				a_n = pow(j,-n)*(2*n + 1)/(n*(n + 1));
				
				p1  = P1(n,theta);
				p1d = P1d(n,theta);
				
				b_n	= (-alpha*J1d(n,k0*r)*J1(n,kd*r) + J1(n,k0*r)*J1d(n,kd*r))/(alpha*H2d(n,k0*r)*J1(n,kd*r) - H2(n,k0*r)*J1d(n,kd*r))*a_n;
				c_n	= (-alpha*J1(n,k0*r)*J1d(n,kd*r) + J1d(n,k0*r)*J1(n,kd*r))/(alpha*H2(n,k0*r)*J1d(n,kd*r) - H2d(n,k0*r)*J1(n,kd*r))*a_n;
				
				A1 	= a_n*J1(n,k0*R) *p1;
				A2 	= a_n*J1d(n,k0*R)*p1d;
				A3 	= a_n*J1(n,k0*R) *p1d;
				A4 	= a_n*J1d(n,k0*R)*p1;

				B1 	= b_n*H2(n,k0*R)*p1;
				B2 	= c_n*H2(n,k0*R)*p1;

				// Incident fields 
				Er0i     	+= A1*n*(n+1);
				Etheta1i  	+= A1;
				Etheta2i 	+= A2;
				Ephi1i    	+= A3;
				Ephi2i   	+= A4;

				Hr0i     	+= A1*n*(n+1);
				Htheta1i  	+= A1;
				Htheta2i 	+= A2;
				Hphi1i    	+= A3;
				Hphi2i   	+= A4;

				// Scattered fields 
				Er0s     	+=  n*(n+1)*B1;
				Etheta1s 	+=  B2;
				Etheta2s 	+=  b_n*H2d(n,k0*R)	*p1d;
				Ephi1s   	+=  c_n*H2(n,k0*R) 	*p1d;
				Ephi2s   	+=  b_n*H2d(n,k0*R)	*p1;

				Hr0s     	+=  n*(n+1)*B2;
				Htheta1s 	+=  B1;
				Htheta2s 	+=  c_n*H2d(n,k0*R)	*p1d;
				Hphi1s   	+=  b_n*H2(n,k0*R) 	*p1d;
				Hphi2s   	+=  c_n*H2d(n,k0*R)	*p1;
				
			}
			
			// Incident fields 
			Eri     	= + E0*cos(phi)/(j*pow(k0*R,2))*Er0i;
			Ethetai		= - E0/(k0*R*sin(theta))*cos(phi)*Etheta1i + E0/(j*k0*R)*cos(phi)*Etheta2i;
			Ephii   	= + E0/(k0*R)*sin(phi)*Ephi1i - E0/(j*k0*R*sin(theta))*sin(phi)*Ephi2i;

			Hri    		= + E0*sin(phi)/(j*eta*pow(k0*R,2))*Hr0i;
			Hthetai		= - E0/(eta*k0*R*sin(theta))*sin(phi)*Htheta1i + E0/(j*eta*k0*R)*sin(phi)*Htheta2i;
			Hphii   	= - E0/(eta*k0*R)*cos(phi)*Hphi1i + E0/(j*eta*k0*R*sin(theta))*cos(phi)*Hphi2i;

			// Scattered fields 
			Ers     	= + E0*cos(phi)/(j*pow(k0*R,2))*Er0s;
			Ethetas		= - E0/(k0*R*sin(theta))*cos(phi)*Etheta1s + E0/(j*k0*R)*cos(phi)*Etheta2s;
			Ephis   	= + E0/(k0*R)*sin(phi)*Ephi1s - E0/(j*k0*R*sin(theta))*sin(phi)*Ephi2s;

			Hrs     	= + E0*sin(phi)/(j*eta*pow(k0*R,2))*Hr0s;
			Hthetas		= - E0/(eta*k0*R*sin(theta))*sin(phi)*Htheta1s + E0/(j*eta*k0*R)*sin(phi)*Htheta2s;
			Hphis   	= - E0/(eta*k0*R)*cos(phi)*Hphi1s + E0/(j*eta*k0*R*sin(theta))*cos(phi)*Hphi2s;

			// total fields 
			Er     	= IN*Eri     	+ Ers;
			Etheta 	= IN*Ethetai	+ Ethetas;
			Ephi   	= IN*Ephii   	+ Ephis;

			Hr     	= IN*Hri     	+ Hrs;
			Htheta 	= IN*Hthetai	+ Hthetas;
			Hphi   	= IN*Hphii   	+ Hphis;
			
			// converting to cartesian 
			Fields[0]  = Er*sin(theta)*cos(phi) +	Etheta*cos(theta)*cos(phi)	-	Ephi*sin(phi);
			Fields[1]  = Er*sin(theta)*sin(phi) +	Etheta*cos(theta)*sin(phi)	+	Ephi*cos(phi);
			Fields[2]  = Er*cos(theta)			- 	Etheta*sin(theta);

			Fields[3]  = Hr*sin(theta)*cos(phi) + 	Htheta*cos(theta)*cos(phi) 	-	Hphi*sin(phi);
			Fields[4]  = Hr*sin(theta)*sin(phi) + 	Htheta*cos(theta)*sin(phi) 	+	Hphi*cos(phi);
			Fields[5]  = Hr*cos(theta)          - 	Htheta*sin(theta);
		}

	return Fields;
}




complex<double> *Mie_fieldsPEC(double IN, double x, double y, double z, double N, double r, double lam0, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs, complex<double> Fields[])
{
	double n, p1, p1d, theta, phi, R;
	complex<double> A1, A2, A3, A4, B1, B2, a_n, b_n, c_n, k0;
	complex<double> Eri, Ethetai, Ephii, Hri, Hthetai, Hphii;
    complex<double> Ers, Ethetas, Ephis, Hrs, Hthetas, Hphis;
    complex<double> Er, Etheta, Ephi, Hr, Htheta, Hphi;


	complex<double> Er0i(0,0), Etheta1i(0,0),  Etheta2i(0,0), Ephi1i(0,0),  Ephi2i(0,0);
	complex<double> Hr0i(0,0), Htheta1i(0,0),  Htheta2i(0,0), Hphi1i(0,0),  Hphi2i(0,0);
	complex<double> Er0s(0,0), Etheta1s(0,0), Etheta2s(0,0), Ephi1s(0,0), Ephi2s(0,0);
	complex<double> Hr0s(0,0), Htheta1s(0,0), Htheta2s(0,0), Hphi1s(0,0), Hphi2s(0,0);


    complex<double> E0(1,0);

	k0 		= 2.0*pi/lam0*sqrt(erb*mrb);
	R 		= sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	phi   	= -atan2(y,-x) + pi;
	theta 	= acos(z/R);

    if (theta == 0)		theta = 1e-20;
    if (theta == pi)	theta = pi - pi*1e-3;

    if (R < r)
    {
        Fields[0] = 0.0;
        Fields[1] = 0.0;
        Fields[2] = 0.0;
		
        Fields[3] = 0.0;
        Fields[4] = 0.0;
        Fields[5] = 0.0;
    }
    else
    {
        for(n = 1.0; n < N+1.0; n++)
		{
            // Define Mie coefficients 
            a_n = pow(j,-n)*(2*n + 1)/(n*(n + 1));
            b_n	= -a_n*J1d(n,k0*r)/H2d(n,k0*r);
            c_n = -a_n*J1(n,k0*r)/H2(n,k0*r);

            p1  = P1(n,theta);
            p1d = P1d(n,theta);

            A1 	= a_n*J1(n,k0*r) *p1;
            A2 	= a_n*J1d(n,k0*r)*p1d;
            A3 	= a_n*J1(n,k0*r) *p1d;
            A4 	= a_n*J1d(n,k0*r)*p1;

            B1 	= b_n*H2(n,k0*r)*p1;
            B2 	= c_n*H2(n,k0*r)*p1;

            // Incident fields 
            Er0i     	+= A1*n*(n+1);
            Etheta1i  	+= A1;
            Etheta2i 	+= A2;
            Ephi1i    	+= A3;
            Ephi2i   	+= A4;

            Hr0i     	+= A1*n*(n+1);
            Htheta1i  	+= A1;
            Htheta2i 	+= A2;
            Hphi1i    	+= A3;
            Hphi2i   	+= A4;

            // Scattered fields 
            Er0s     	+=  n*(n+1)*B1;
            Etheta1s 	+=  B2;
            Etheta2s 	+=  b_n*H2d(n,k0*r)	*p1d;
            Ephi1s   	+=  c_n*H2(n,k0*r) 	*p1d;
            Ephi2s   	+=  b_n*H2d(n,k0*r)	*p1;

            Hr0s     	+=  n*(n+1)*B2;
            Htheta1s 	+=  B1;
            Htheta2s 	+=  c_n*H2d(n,k0*r)	*p1d;
            Hphi1s   	+=  b_n*H2(n,k0*r) 	*p1d;
            Hphi2s   	+=  c_n*H2d(n,k0*r)	*p1;
        }
        
		// Incident fields 
        Eri     	= + E0*cos(phi)/(j*pow(k0*R,2))*Er0i;
        Ethetai		= - E0/(k0*R*sin(theta))*cos(phi)*Etheta1i + E0/(j*k0*R)*cos(phi)*Etheta2i;
        Ephii   	= + E0/(k0*R)*sin(phi)*Ephi1i - E0/(j*k0*R*sin(theta))*sin(phi)*Ephi2i;

        Hri    		= + E0*sin(phi)/(j*eta0*pow(k0*R,2))*Hr0i;
        Hthetai		= - E0/(eta0*k0*R*sin(theta))*sin(phi)*Htheta1i + E0/(j*eta0*k0*R)*sin(phi)*Htheta2i;
        Hphii   	= - E0/(eta0*k0*R)*cos(phi)*Hphi1i + E0/(j*eta0*k0*R*sin(theta))*cos(phi)*Hphi2i;

        // Scattered fields 
        Ers     	= + E0*cos(phi)/(j*pow(k0*R,2))*Er0s;
        Ethetas		= - E0/(k0*R*sin(theta))*cos(phi)*Etheta1s + E0/(j*k0*R)*cos(phi)*Etheta2s;
        Ephis   	= + E0/(k0*R)*sin(phi)*Ephi1s - E0/(j*k0*R*sin(theta))*sin(phi)*Ephi2s;

        Hrs     	= + E0*sin(phi)/(j*eta0*pow(k0*R,2))*Hr0s;
        Hthetas		= - E0/(eta0*k0*R*sin(theta))*sin(phi)*Htheta1s + E0/(j*eta0*k0*R)*sin(phi)*Htheta2s;
        Hphis   	= - E0/(eta0*k0*R)*cos(phi)*Hphi1s + E0/(j*eta0*k0*R*sin(theta))*cos(phi)*Hphi2s;
    
        // total fields 
        Er     	= IN*Eri     	+ Ers;
        Etheta 	= IN*Ethetai	+ Ethetas;
        Ephi   	= IN*Ephii   	+ Ephis;

        Hr     	= IN*Hri     	+ Hrs;
        Htheta 	= IN*Hthetai	+ Hthetas;
        Hphi   	= IN*Hphii   	+ Hphis;
        
        // converting to cartesian 
        Fields[0]  = Er*sin(theta)*cos(phi) +	Etheta*cos(theta)*cos(phi)	-	Ephi*sin(phi);
        Fields[1]  = Er*sin(theta)*sin(phi) +	Etheta*cos(theta)*sin(phi)	+	Ephi*cos(phi);
        Fields[2]  = Er*cos(theta)			- 	Etheta*sin(theta);

        Fields[3]  = Hr*sin(theta)*cos(phi) + 	Htheta*cos(theta)*cos(phi) 	-	Hphi*sin(phi);
        Fields[4]  = Hr*sin(theta)*sin(phi) + 	Htheta*cos(theta)*sin(phi) 	+	Hphi*cos(phi);
        Fields[5]  = Hr*cos(theta)          - 	Htheta*sin(theta);
    }

	return Fields;
}


double * force(vector<vector<complex<double>>> F, vector<vector<double>> pts, double res){
    
    complex<double> Ex, Ey, Ez, Hx, Hy, Hz, Ecx, Ecy, Ecz, Hcx, Hcy, Hcz;
	static double FF[3];
    double m0, e0, c0, D, fx = 0.0, fy = 0.0, fz = 0.0, theta, phi, R;
    int i;
   
    c0 = 299792458.0;
    m0 = pi*4e-7;
    e0 = 1.0/(m0*pow(c0,2));
	R 	= sqrt(pow(pts[0][0],2) + pow(pts[0][1],2) + pow(pts[0][2],2));
    D = pow(R*res*pi/180.0,2);
    
    for(i = 0; i < int(pts.size()); i++)
    {
		phi   	= -atan2(pts[i][1],-pts[i][0]) + pi;
		theta 	= acos(pts[i][2]/R);
		
		Ex = F[i][0]; Ecx = conj(Ex);
		Ey = F[i][1]; Ecy = conj(Ey);
		Ez = F[i][2]; Ecz = conj(Ez);
		Hx = F[i][3]; Hcx = conj(Hx);
		Hy = F[i][4]; Hcy = conj(Hy);
		Hz = F[i][5]; Hcz = conj(Hz);
   
		fx += 0.5*real((0.5*((Ecx*Ex - Ecy*Ey - Ecz*Ez)*e0 + (Hcx*Hx - Hcy*Hy - Hcz*Hz)*m0))*sin(theta)*cos(phi)+(Ecy*Ex*e0 + Hcy*Hx*m0)*sin(theta)*sin(phi)+(Ecz*Ex*e0+ Hcz*Hx*m0)*cos(theta))*sin(theta);
		fy += 0.5*real((Ecx*Ey*e0 + Hcx*Hy*m0)*sin(theta)*cos(phi)+(0.5*((-Ecx*Ex + Ecy*Ey - Ecz*Ez)*e0 - (Hcx*Hx- Hcy*Hy + Hcz*Hz)*m0))*sin(theta)*sin(phi)+(Ecz*Ey*e0 + Hcz*Hy*m0)*cos(theta))*sin(theta);    
		fz += 0.5*real((Ecx*Ez*e0 + Hcx*Hz*m0)*sin(theta)*cos(phi)+(Ecy*Ez*e0 + Hcy*Hz*m0)*sin(theta)*sin(phi)+(0.5*((-Ecx*Ex- Ecy*Ey + Ecz*Ez)*e0 - (Hcx*Hx + Hcy*Hy - Hcz*Hz)*m0))*cos(theta))*sin(theta);
    }
    
    FF[0] = fx*D;
    FF[1] = fy*D;
    FF[2] = fz*D;
    
    return FF;
}



double * forceTETM(vector<vector<complex<double>>> TM, vector<vector<complex<double>>> TE, vector<vector<double>> A, double dA, double R, double F[]){
    
    complex<double> Ex, Ey, Ez, Hx, Hy, Hz, Ecx, Ecy, Ecz, Hcx, Hcy, Hcz;
    
    double m0, e0, c0, D, fx = 0.0, fy = 0.0, fz = 0.0, theta, phi;
    int i;
   
    c0 = 299792458.0;
    m0 = pi*4e-7;
    e0 = 1.0/(m0*pow(c0,2));
    D = pow(R*dA*pi/180.0,2);
    
    for(i = 0; i < int(A.size()); i++)
    {
		Ex = TM[i][0] + TE[i][0]; Ecx = conj(Ex);
		Ey = TM[i][1] + TE[i][1]; Ecy = conj(Ey);
		Ez = TM[i][2] + TE[i][2]; Ecz = conj(Ez);
		Hx = TM[i][3] + TE[i][3]; Hcx = conj(Hx);
		Hy = TM[i][4] + TE[i][4]; Hcy = conj(Hy);
		Hz = TM[i][5] + TE[i][5]; Hcz = conj(Hz);
		
		phi   = A[i][0]*pi/180.0;
		theta = A[i][1]*pi/180.0;
   
		fx += 0.5*real((0.5*((Ecx*Ex - Ecy*Ey - Ecz*Ez)*e0 + (Hcx*Hx - Hcy*Hy - Hcz*Hz)*m0))*sin(theta)*cos(phi)+(Ecy*Ex*e0 + Hcy*Hx*m0)*sin(theta)*sin(phi)+(Ecz*Ex*e0+ Hcz*Hx*m0)*cos(theta))*sin(theta);
		fy += 0.5*real((Ecx*Ey*e0 + Hcx*Hy*m0)*sin(theta)*cos(phi)+(0.5*((-Ecx*Ex + Ecy*Ey - Ecz*Ez)*e0 - (Hcx*Hx- Hcy*Hy + Hcz*Hz)*m0))*sin(theta)*sin(phi)+(Ecz*Ey*e0 + Hcz*Hy*m0)*cos(theta))*sin(theta);    
		fz += 0.5*real((Ecx*Ez*e0 + Hcx*Hz*m0)*sin(theta)*cos(phi)+(Ecy*Ez*e0 + Hcy*Hz*m0)*sin(theta)*sin(phi)+(0.5*((-Ecx*Ex- Ecy*Ey + Ecz*Ez)*e0 - (Hcx*Hx + Hcy*Hy - Hcz*Hz)*m0))*cos(theta))*sin(theta);
    }
    
    F[0] = fx*D;
    F[1] = fy*D;
    F[2] = fz*D;
    
    return F;
}


void CsToFileR(const char *Name, vector<double> rv, vector<complex<double>> CSsca, vector<complex<double>> CSext, double lam, double n, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs)
{
	int i;
	FILE *f = fopen(Name, "w");

	fprintf(f, "# wavelength = %.e, number of Mie coefficients = %.2e, erb = %.2e+%.2ej, mrb = %.2e+%.2ej, ers = %.2e+%.2ej, mrs = %.2e+%.2ej\n", lam, n, real(erb), imag(erb), real(mrb), imag(mrb), real(ers), imag(ers), real(mrs), imag(mrs));
    fprintf(f, "# r Csca Cext\n");
  
	for(i = 0; i < int(rv.size()); i++)	fprintf(f, "%.10e\t%.10e\t%.10e\n", rv[i], abs(CSsca[i]), abs(CSext[i]));
      
    fclose(f);
}

void CsToFileW(const char *Name, vector<double> wv, vector<complex<double>> CSsca, vector<complex<double>> CSext, double r, double n, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs)
{
	int i;
	FILE *f = fopen(Name, "w");

	fprintf(f, "# radius = %.e, number of Mie coefficients = %.2e, erb = %.2e+%.2ej, mrb = %.2e+%.2ej, ers = %.2e+%.2ej, mrs = %.2e+%.2ej\n", r, n, real(erb), imag(erb), real(mrb), imag(mrb), real(ers), imag(ers), real(mrs), imag(mrs));
    fprintf(f, "# lam Csca Cext\n");
  
	for(i = 0; i < int(wv.size()); i++)	fprintf(f, "%.10e\t%.10e\t%.10e\n", wv[i], abs(CSsca[i]), abs(CSext[i]));
      
    fclose(f);
}

void FieldsToFileXYZ(const char *Name, vector<vector<double>> A, vector<vector<complex<double>>> vT, double r, double lam0, double n, complex<double> erb, complex<double> mrb, complex<double> ers, complex<double> mrs)
{
	int i;
	FILE *f = fopen(Name, "w");

	fprintf(f, "# Sphere radius = %.2e, wavelength = %.e, number of Mie coefficients = %.2e, erb = %.2e+%.2ej, mrb = %.2e+%.2ej, ers = %.2e+%.2ej, mrs = %.2e+%.2ej\n", r, lam0, n, real(erb), imag(erb), real(mrb), imag(mrb), real(ers), imag(ers), real(mrs), imag(mrs));
    fprintf(f, "# x y z Exr Exi Eyr Eyi Ezr Ezi Hxr Hxi Hyr Hyi Hzr Hzi\n");
  
	 
	for(i = 0; i < int(A.size()); i++)
    {
		fprintf(f, "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", A[i][0], A[i][1], A[i][2], real(vT[i][0]), imag(vT[i][0]), real(vT[i][1]), imag(vT[i][1]), real(vT[i][2]), imag(vT[i][2]), real(vT[i][3]), imag(vT[i][3]), real(vT[i][4]), imag(vT[i][4]), real(vT[i][5]), imag(vT[i][5]));
	}
   
      
    fclose(f);
}

void FieldsToFileXYZ_PEC(const char *Name, vector<vector<double>> A, vector<vector<complex<double>>> vT, double r, double lam, double n)
{
	int i;
	FILE *f = fopen(Name, "w");

	fprintf(f, "# Sphere radius = %.2e, wavelength = %.e, number of Mie coefficients = %.2e\n", r, lam, n);
    fprintf(f, "# x y z Exr Exi Eyr Eyi Ezr Ezi Hxr Hxi Hyr Hyi Hzr Hzi\n");
  
	 
	for(i = 0; i < int(A.size()); i++)
    {
		fprintf(f, "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", A[i][0], A[i][1], A[i][2], real(vT[i][0]), imag(vT[i][0]), real(vT[i][1]), imag(vT[i][1]), real(vT[i][2]), imag(vT[i][2]), real(vT[i][3]), imag(vT[i][3]), real(vT[i][4]), imag(vT[i][4]), real(vT[i][5]), imag(vT[i][5]));
	}
   
      
    fclose(f);
}


void FieldsToFile(const char *Name, vector<vector<double>> A, vector<vector<complex<double>>> vT)
{
	int i;
	FILE *f = fopen(Name, "w");

    fprintf(f, "# phi theta Exr Exi Eyr Eyi Ezr Ezi Hxr Hxi Hyr Hyi Hzr Hzi\n");
  
	 
	for(i = 0; i < int(A.size()); i++)
    {
		fprintf(f, "%f\t%f\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", A[i][0], A[i][1], real(vT[i][0]), imag(vT[i][0]), real(vT[i][1]), imag(vT[i][1]), real(vT[i][2]), imag(vT[i][2]), real(vT[i][3]), imag(vT[i][3]), real(vT[i][4]), imag(vT[i][4]), real(vT[i][5]), imag(vT[i][5]));
	}
   
      
    fclose(f);
}



void ForcesToFile(const char *Name, vector<vector<double>> F, vector<double> r, double D[])
{
	int i;
	FILE *f = fopen(Name, "w");

	fprintf(f, "# Wavelength = %.3e, r_start = %.3e, r_stop = %.3e, r_step = %.3e, n = %.3e, angular resolution = %.3e, R = %.3e\n", D[0], D[1], D[2], D[3], D[4], D[5], D[6]);
    fprintf(f, "# r Fx Fy Fz\n");
  
	 
	for(i = 0; i < int(F.size()); i++)
    {
		fprintf(f, "%.10e\t%.10e\t%.10e\t%.10e\n", r[i], F[i][0], F[i][1], F[i][2]);
	}
   
      
    fclose(f);
}

vector<double> range(double start, double stop, double step)
{
    static vector<double> A;
    double i;

    for(i = start; i <= stop; i = i + step)
    {
        A.push_back(i);
    }
    return A;
}


vector<double> linspace(double start, double stop, int N)
{
    static vector<double> A;
    int i;


    for(i = 0; i < N; i++)
    {
        A.push_back((stop-start)/(N-1)*i + start);
    }
    return A;
}












