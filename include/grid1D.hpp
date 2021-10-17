#ifndef __GRID1D__HEADER__
#define __GRID1D__HEADER__

#include <iostream>
#include <fstream>
#include <cmath> 
#include <vector>
#include <complex>
#include "materials.hpp"

enum class Boundary_t { PEC, Mur1};
enum class Source_t { Gaussian, Sinusoidal, Ricker };

struct Grid1Dsettings {

    // Source information
    double E0 = 1.0;      // Signal amplitude
    double f  = 1.0e6;       // Peak frequency
    Source_t sourcetype = Source_t::Gaussian;    // Source type
    std::string field_fname = "data/fields.dat";
    std::string spectral_fname = "data/spectrum.dat";

    // Grid information
    size_t Nx;                      // Total spatial steps
    double Nl = 20;                 // Cell resolution
    double dmin = 10000;            // Finest detail in grid
    double Nd   = 1;                // Sampled points in dmin
    double Sc = 1;                  // Courant 
    size_t tfsfL;                   // Left  TFSF boundary
    Boundary_t leftBound = Boundary_t::Mur1;
    Boundary_t rightBound = Boundary_t::Mur1;

    size_t Nf = 50;
    double fm = 1.0e7;
    double df = 1.0e7;
    
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
    double Nl;      // Cell resolution
    double dmin;    // Finest detail in grid
    double Nd;      // Sampled points in dmin
    double dx;      // Spatial step
    double dt;      // Time step
    double Sc;      // Courant 
    size_t tfsfL;   // Left TFSF boundary
    Boundary_t lbt; // Left  Boundary condition
    Boundary_t rbt; // Right Boundary condition
    double EzBL;    // Temporary variable to facilitate ABC
    double EzBR;    // Temporary variable to facilitate ABC
    std::vector<Material> material_list = { materialdb["Air"] };
                    // List of materials inside the grid
    
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
    std::complex<double>* Rf;
    std::complex<double>* Tf;
    std::complex<double>* Sf;
    
    // Fourier transforms
    std::complex<double> * Ks;     // The kernels
    size_t   Nf;     // Number of frequencies
    double   fm;     // The max frequency
    double   df;     // The frequency resolution

    // File to save
    std::string field_fname;
    std::string spectral_fname;
    std::ofstream fhandle;
    std::ofstream shandle;

  public:
    
    Grid1D(Grid1Dsettings gs);

    // Initialize the grid to vacuum
    void init_vacuum();

    // Add a material spanning from iL to iR
    void add_material(double iL, double iR,
                      std::string mname);

    // Update the Ce and Ch coefficients of the grid
    void update_coefs();

    // Update the fields
    void update_magnetic();
    void update_electric();

    // Update the TFSF boundary
    void update_tfsf_magnetic(size_t n);
    void update_tfsf_electric(size_t n);
    
    // Store and use values for ABC
    void setup_ABC();
    void update_ABC();
    
    // Fourier Transforms
    void init_kernels();
    void update_kernels(double n);
    void finalize_kernels();

    // Field data output
    void open_result_file();
    void close_result_file();
    void save_results();
    void save_spectrum();

    // Return the number of total time steps Nt
    size_t tsteps();

    // The available sources for the 1D grid
    double gaussian(double t, double x, double f);
    double sinusoidal(double t, double x, double f);
    double ricker(double t, double x, double f);

    ~Grid1D();
};

#endif
