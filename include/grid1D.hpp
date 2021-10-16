#ifndef __GRID1D__HEADER__
#define __GRID1D__HEADER__

#include <iostream>
#include <fstream>
#include <cmath> 
#include <vector>
#include "materials.hpp"

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
    std::vector<Material> material_list = { materialdb["Air"] };
    
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
    
    Grid1D(Grid1Dsettings gs);

    void init_vacuum();

    void add_material(double iL, double iR,
                      std::string mname);

    void update_coefs();

    void update_magnetic();

    void update_tfsf_magnetic(size_t n);

    void update_electric();

    void update_tfsf_electric(size_t n);
    
    void setup_ABC();

    void update_ABC();
    
    void open_result_file();

    void close_result_file();

    void save_fields();

    size_t tsteps();

    double gaussian(double t, double x, double f);

    double sinusoidal(double t, double x, double f);

    double ricker(double t, double x, double f);

    ~Grid1D();
};

#endif
