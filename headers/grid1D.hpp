#ifndef __GRID1D__HEADER__
#define __GRID1D__HEADER__

#include <iostream>
#include <fstream>
#include <cmath>
#include "materials.hpp"
#include <vector>

enum class Boundary_t { PEC, Mur1};
enum class Source_t { Gaussian, Sinusoidal, Ricker };

struct Grid1Dsettings {
    // Grid information
    size_t Nx;      // Total spatial steps
    size_t Nt;      // Total time steps
    double dx;      // Spatial step
    double Sc = 1;      // Courant 
    size_t tfsfL;  // Left  TFSF boundary
    Boundary_t leftBound = Boundary_t::Mur1;
    Boundary_t rightBound = Boundary_t::Mur1;
    
    // Source information
    double E0 = 1.0;      // Signal amplitude
    double f  = 1.0e6;       // Peak frequency
    Source_t sourcetype = Source_t::Gaussian;    // Source type
    std::string fname = "results.dat";
};

class Grid1D
{
  private:
    double cc   = 2.99792458e+8;
    double mu0  = M_PI*4.0*1e-7;
    double eps0 = 1.0 / (cc*cc*mu0);

    // Grid information
    size_t Nx;      // Total spatial steps
    size_t Nt;      // Total time steps
    double dx;      // Spatial step
    double dt;      // Time step
    double Sc;      // Courant 
    size_t tfsfL;  // Left  TFSF boundary
    Boundary_t lbt;
    Boundary_t rbt;
    double EzBL;
    double EzBR;
    std::vector<Material> material_list = { materials["Air"] };
    
    // Source information
    double E0;      // Signal amplitude
    double f;       // Peak frequency
    double T;       // Period
    double lambda;  // Wavelength
    Source_t st;    // Source type

    // FDTD arrays
    double* epsr;   // Relative eps
    double* mur;    // Relative mu
    double* Ce;     // Coefficient for E
    double* Ch;     // Coefficient for H
    double* Ez;     // E field values
    double* Hy;     // M field values
    
    // File to save
    std::string fname = "out.dat";
    std::ofstream handle;

  public:
    
    Grid1D(Grid1Dsettings gs){
      E0 = gs.E0;
      f  = gs.f;
      lambda = cc/f;
      T  = 1.0/f;
      st = gs.sourcetype;

      Nx = gs.Nx;
      Nt = gs.Nt;
      Sc = gs.Sc;
      dx = lambda/20.0;
      dt = Sc*dx/cc;
      tfsfL = gs.tfsfL;

      lbt = gs.leftBound;
      rbt = gs.rightBound;

      fname = gs.fname;

      epsr = new double[Nx]();
      mur  = new double[Nx]();

      init_vacuum();

      Ce = new double[Nx]();
      Ch = new double[Nx]();

      Ez = new double[Nx]();
      Hy = new double[Nx]();
      
    }

    void init_vacuum(){
      for (size_t i = 0; i < Nx; i++){
        epsr[i] = 1.0;
        mur[i]  = 1.0;
      }
      return;
    }

    void add_material(double iL, double iR,
                      std::string mname){
      for (size_t i = iL; i < iR; i++){
        epsr[i] = materials[mname].epsr;
        mur[i]  = materials[mname].mur;
      }
      material_list.push_back(materials[mname]);
      return;
    }

    void update_coefs(){
      for (size_t i = 0; i < Nx; i++){
        Ce[i] = dt / (eps0*epsr[i]) / dx;
        Ch[i] = dt / (mu0*mur[i])   / dx;
      }
      return;
    }

    void update_magnetic(){
      for (size_t i = 0; i < Nx - 1; i++){
        Hy[i] = Hy[i] + Ch[i]*(Ez[i+1] - Ez[i]);
      }
      return;
    }

    void update_tfsf_magnetic(size_t n){
      // TFSF correction for Hy
      if (st == Source_t::Gaussian){
        Hy[tfsfL-1] -= Ch[tfsfL-1]*gaussian(n*dt,0,f);
      } 
      else if (st == Source_t::Sinusoidal){
        Hy[tfsfL-1] -= Ch[tfsfL-1]*sinusoidal(n*dt,0,f);
      } 
      else if (st == Source_t::Ricker){
        Hy[tfsfL-1] -= Ch[tfsfL-1]*ricker(n*dt,0,f);
      }
      return;
    }

    void update_electric(){
      for (size_t i = 1; i < Nx - 1; i++){
        Ez[i] = Ez[i] + Ce[i]*(Hy[i] - Hy[i-1]);
      }
      return;
    }

    void update_tfsf_electric(size_t n){
      // TFSF correction for Ez
      if (st == Source_t::Gaussian){
        Ez[tfsfL] += Ce[tfsfL]*gaussian((n+0.5)*dt,-dx/2.0,f)/sqrt(mu0/eps0);;
      } 
      else if (st == Source_t::Sinusoidal){
        Ez[tfsfL] += Ce[tfsfL]*sinusoidal((n+0.5)*dt,-dx/2.0,f)/sqrt(mu0/eps0);
      } 
      else if (st == Source_t::Ricker){
        Ez[tfsfL] += Ce[tfsfL]*ricker((n+0.5)*dt,-dx/2.0,f)/sqrt(mu0/eps0);
      }
      return;
    }
    
    void setup_ABC(){
      EzBL = Ez[1];
      EzBR = Ez[Nx-2];
      return;
    }

    void update_ABC(){
      if (lbt == Boundary_t::Mur1){
        Ez[0]    = EzBL + (cc*dt-dx)*(Ez[1]   - Ez[0]  )/(cc*dt+dx);
      }
      if (rbt == Boundary_t::Mur1){
        Ez[Nx-1] = EzBR + (cc*dt-dx)*(Ez[Nx-2]-Ez[Nx-1])/(cc*dt+dx);
      }
    }
    
    void open_result_file(){
      handle.open(fname);
      return;
    }

    void close_result_file(){
      handle.close();
      return;
    }

    void save_fields(){
      // Print to file
      handle << "\n\n";
      for (size_t i = 0; i < Nx; i++){
        handle << i*dx << " " << Ez[i] << " " << Hy[i] << "\n";
      }
      return;
    }

    size_t tsteps(){
      return Nt;
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

    ~Grid1D(){
      if (!epsr) { delete[] epsr; }
      if (!mur)  { delete[] mur;  }
      if (!Ce)   { delete[] Ce;   }
      if (!Ch)   { delete[] Ch;   }
      if (!Ez)   { delete[] Ez;   }
      if (!Hy)   { delete[] Hy;   }
    }
};

#endif
