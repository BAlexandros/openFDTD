#include <iostream>
#include <fstream>
#include <cmath>

double cc;
double mu0;
double eps0;

enum class Source_t { Gaussian, Sinusoidal, Ricker };

double gaussian(double t, double x, double f);
double sinusoidal(double t, double x, double f);
double ricker(double t, double x, double f);

int main(void)
{
  cc   = 2.99792458e+8;
  mu0  = M_PI*4.0*1e-7;
  eps0 = 1.0 / (cc*cc*mu0);

  double E0     = 1.0;
  double freq   = 9.00e6;
  double lambda = cc/freq;
  double T      = 1.0 / freq;
  Source_t st   = Source_t::Gaussian;

  size_t Nx = 60;
  size_t Nt = 200;
  
  double Sc = 1.0;
  double dx = lambda/20.0;
  double dt = Sc*dx/cc;

  size_t tfsfBL = Nx/4;

  double ce = dt / eps0 / dx; 
  double ch = dt / mu0 / dx;
 
  double* Ce = new double[Nx]();
  double* Ch = new double[Nx]();

  for(size_t i = 0; i < Nx; i++){
    Ce[i] = ce;
    Ch[i] = ch;
  }

  double* Ez = new double[Nx]();
  double* Hy = new double[Nx]();

  std::ofstream out;
  out.open("./results.dat");

  // Time loop
  for (size_t n = 0; n < Nt; n++){

    out << "\n\n";

    // Magnetic field loop
    for (size_t i = 0; i < Nx - 1; i++){
      Hy[i] = Hy[i] + Ch[i]*(Ez[i+1] - Ez[i]);
    }

    // TFSF correction for Hy
    if (st == Source_t::Gaussian){
      Hy[tfsfBL-1] -= Ch[tfsfBL-1]*gaussian(n*dt,0,freq);
    } 
    else if (st == Source_t::Sinusoidal){
      Hy[tfsfBL-1] -= Ch[tfsfBL-1]*sinusoidal(n*dt,0,freq);
    } 
    else if (st == Source_t::Ricker){
      Hy[tfsfBL-1] -= Ch[tfsfBL-1]*ricker(n*dt,0,freq);
    }

    // 1st order Mur temporary variables
    double EzBL = Ez[1];
    double EzBR = Ez[Nx-2];

    // Electric field loop
    for (size_t i = 1; i < Nx - 1; i++){
      Ez[i] = Ez[i] + Ce[i]*(Hy[i] - Hy[i-1]);
    }

    // Update Mur Boundaries 
    Ez[0]    = EzBL + (cc*dt-dx)*(Ez[1]   - Ez[0]  )/(cc*dt+dx);
    // Ez[Nx-1] = EzBR + (cc*dt-dx)*(Ez[Nx-2]-Ez[Nx-1])/(cc*dt+dx); 

    // TFSF correction for Ex
    if (st == Source_t::Gaussian){
      Ez[tfsfBL] += Ce[tfsfBL]*gaussian((n+0.5)*dt,-dx/2.0,freq)/sqrt(mu0/eps0);;
    } 
    else if (st == Source_t::Sinusoidal){
      Ez[tfsfBL] += Ce[tfsfBL]*sinusoidal((n+0.5)*dt,-dx/2.0,freq)/sqrt(mu0/eps0);
    } 
    else if (st == Source_t::Ricker){
      Ez[tfsfBL] += Ce[tfsfBL]*ricker((n+0.5)*dt,-dx/2.0,freq)/sqrt(mu0/eps0);
    }

    // Print to file
    for (size_t i = 0; i < Nx; i++){
      out << i*dx << " " << Ez[i] << " " << Hy[i] << "\n";
    }
   
  }

  out.close();

  delete[] Ce;
  delete[] Ch;
  delete[] Ez;
  delete[] Hy;

  return 0;
}

double gaussian(double t, double x, double f){
  double tau = 0.5/f;
  double t0  = 5.0*tau;
  double c   = (t - x/cc - t0)/tau;
  return exp(-c*c);
}

double sinusoidal(double t, double x, double f){
  return sin(2.0*M_PI*f*(t - x/cc));
}

double ricker(double t, double x, double f){
  double d = 1.5/f;
  double c = M_PI*f*(t - x/cc - d);
  return (1.0 - 2.0*c*c)*exp(-c*c);
}
