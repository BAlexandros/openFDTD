#ifndef __GRID1D__HEADER__
#define __GRID1D__HEADER__

#include <iostream>
#include <fstream>
#include <cmath> 
#include <vector>
#include <complex>
#include "materials.hpp"

enum class Boundary_t { PEC, Mur1};
enum class Source_t { Gaussian, Sinusoidal, Ricker, Custom };

struct GridMat {
  int x1;
  int x2;
  std::string matname;
};

class Grid1D
{
  private:
    double cc   = 2.99792458e+8;
    double mu0  = M_PI*4.0*1e-7;
    double eps0 = 1.0 / (cc*cc*mu0);

    // Grid information
    int Nx;      // Total spatial steps
    int Nt;      // Total time steps
    int Nl;      // Cell resolution
    double dm;      // Finest detail in grid
    int Nd;      // Sampled points in dmin
    double dx;      // Spatial step
    double dt;      // Time step
    double Sc;      // Courant 
    int tfsfL;   // Left TFSF boundary
    Boundary_t lbt; // Left  Boundary condition
    Boundary_t rbt; // Right Boundary condition
    double EzBL;    // Temporary variable to facilitate ABC
    double EzBR;    // Temporary variable to facilitate ABC
    
    // Source information
    double E0;      // Signal amplitude
    double f;       // Peak frequency
    double T;       // Period
    double lambda;  // Wavelength
    Source_t st;    // Source type
    double (Grid1D::*sourcefunc)(double,double); // Source function
    double (*othersource)(double,double);        // Custom source

    // FDTD arrays
    double* epsr;   // Relative eps
    double* mur;    // Relative mu
    double* Ce;     // Coefficient for E
    double* Ch;     // Coefficient for H
    double* Ez;     // E field values
    double* Hy;     // M field values

    // Fourier transforms
    std::complex<double> *Ks;   // The kernels
    std::complex<double> *Rf;   // Reflection Fourier array
    std::complex<double> *Tf;   // Transmission Fourier array
    std::complex<double> *Sf;   // Source Fourier array
    int Nf;                  // Number of frequencies
    double fm;                  // The max frequency
    double df;                  // The frequency resolution

    // File to save
    std::string   field_fname;    // filename to save field
    std::string   spectrum_fname; // filename to save spectrum
    std::ofstream fhandle;        // field file handle
    std::ofstream shandle;        // specrrum file handle

    // Simulation settings
    bool parallelism_enabled  = false;
    bool save_field_progress  = false;
    bool save_refl_trans_spec = false;

  public:
    
    Grid1D();

    std::vector<GridMat> material_list;
                    // List of materials inside the grid

    // Initializers
    void init_coefs();   // Initialize update coefficients
    void init_fields();  // Initialize E and H fields
    void init_kernels(); // Initialize complex kernels
    void init_fourier(); // Initialize Fourier transformations
    void init_ABC();     // Initialize Absorbing Boundary vars

    // Add a material spanning from iL to iR
    void add_material(int iL, int iR,
                      std::string mname);

    // Update the Ce and Ch coefficients of the grid
    void update_magnetic();
    void update_electric();
    void update_ABC();
    void update_tfsf_magnetic(int n);
    void update_tfsf_electric(int n);
    void update_fourier(double n);
    
    // Fourier Transforms
    void finalize_kernels();

    // Field data output
    void open_result_file();
    void close_result_file();
    void save_cur_field();
    void save_spectrum();

    // Execute the simulation
    void run_simulation();

    // The available sources for the 1D grid
    double gaussian(double q, double m);
    double sinusoidal(double q, double m);
    double ricker(double q, double m);

    // Parallelism flag
    void enable_parallelism();
    void disable_parallelism();

    // Parameter setters
    void setE0(double);
    void setFreq(double);
    void setSourceType(Source_t);
    void setSc(double);
    void setNx(int);
    void setNt(int);
    void setNl(int);
    void setNd(int);
    void setdm(double);
    void setFieldProgFname(std::string);
    void setSpectrumFname(std::string);
    void recordFieldProgTrue();
    void recordFieldProgFalse();
    void recordRefTranSpectTrue();
    void recordRefTranSpectFalse();

    void makeFieldAnimation();

    ~Grid1D();
};

#endif
