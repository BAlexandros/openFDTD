#include <iostream>
#include <fstream>
#include <cmath>

double cc;
double mu0;
double eps0;

int main(void)
{
  cc   = 2.99792458e+8;
  mu0  = M_PI*4.0*1e-7;
  eps0 = 1.0 / (cc*cc*mu0);

  double E0     = 1.0;
  double freq   = 9.00e6;
  double lambda = cc/freq;
  double T      = 1.0 / freq;

  size_t Nx = 100;
  size_t Nt = 250;
  
  double Sc = 1.0;
  double dx = lambda/15.0;
  double dt = Sc*dx/cc;

  double ce = dt / eps0 / dx; 
  double ch = dt / mu0 / dx;
 
  double* Ce = new double[Nx]();
  double* Ch = new double[Nx]();

  for(size_t i = 0; i < Nx; i++){
    Ce[i] = ce;
    Ch[i] = ch;
  }

  double* Ex = new double[Nx]();
  double* Hy = new double[Nx]();

  std::ofstream out;
  out.open("./results.dat");

  // Time loop
  for (size_t n = 0; n < Nt; n++){

    out << "\n\n";

    // 1st order Mur
    double ExBL = Ex[1];
    double ExBR = Ex[Nx-2];

    // Electric field loop
    for (size_t i = 1; i < Nx - 1; i++){
      Ex[i] = Ex[i] - Ce[i]*(Hy[i] - Hy[i-1]);
    }
   
    Ex[0] = ExBL + (cc*dt-dx)*(Ex[1] - Ex[0])/(cc*dt+dx);
    Ex[Nx-1] = ExBR + (cc*dt-dx)*(Ex[Nx-2]-Ex[Nx-1])/(cc*dt+dx);

        
    // Update source
    Ex[Nx/2] = sin(2.0*M_PI*freq*n*dt);
  
    // Magnetic field loop
    for (size_t i = 0; i < Nx - 1; i++){
      Hy[i] = Hy[i] - Ch[i]*(Ex[i+1] - Ex[i]);
    }

    for (auto i = 0; i < Nx; i++){
      out << i*dx << " " << Ex[i] << " " << Hy[i] << "\n";
    }
    
   
  }

  out.close();

  delete[] Ce;
  delete[] Ch;
  delete[] Ex;
  delete[] Hy;

  return 0;
}
