#include "../include/grid1D.hpp"

Grid1D::Grid1D(Grid1Dsettings gs){
  E0 = gs.E0;
  f  = gs.f;
  lambda = cc/f;
  T  = 1.0/f;
  st = gs.sourcetype;

  Nl = gs.Nl;
  dmin = gs.dmin;
  Nd   = gs.Nd;
  Sc = gs.Sc;
  dx = (lambda/Nl < dmin/Nd) ? lambda/Nl : dmin/Nd;
  dt = Sc*dx/cc;
  tfsfL = gs.tfsfL;
  Nx = gs.Nx;
  Nt = 6.0*Nx*dx/cc/dt;
  std::cout << "Nt: " << Nt << "\n";

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

void Grid1D::init_vacuum(){
  for (size_t i = 0; i < Nx; i++){
    epsr[i] = 1.0;
    mur[i]  = 1.0;
  }
  std::cout << "dx = " << dx << "\n";
  return;
}

void Grid1D::add_material(double iL, double iR,
    std::string mname){
  double erm = materialdb[mname].epsr;
  double mrm = materialdb[mname].mur;
  for (size_t i = iL; i < iR; i++){
    epsr[i] = erm;
    mur[i]  = mrm;
  }
  material_list.push_back(materialdb[mname]);
  double lambda_mat = cc / sqrt(erm*mrm) / f;
  double dxnew = lambda_mat / Nl;
  dx = (dxnew < dx) ? dxnew : dx;
  dt = Sc*dx/cc;
  std::cout << "dx after " << mname << " added: " << dx << "\n";
  return;
}

void Grid1D::update_coefs(){
  for (size_t i = 0; i < Nx; i++){
    Ce[i] = dt / (eps0*epsr[i]) / dx;
    Ch[i] = dt / (mu0*mur[i])   / dx;
  }
  return;
}

void Grid1D::update_magnetic(){
  for (size_t i = 0; i < Nx - 1; i++){
    Hy[i] = Hy[i] + Ch[i]*(Ez[i+1] - Ez[i]);
  }
  return;
}

void Grid1D::update_tfsf_magnetic(size_t n){
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

void Grid1D::update_electric(){
  for (size_t i = 1; i < Nx - 1; i++){
    Ez[i] = Ez[i] + Ce[i]*(Hy[i] - Hy[i-1]);
  }
  return;
}

void Grid1D::update_tfsf_electric(size_t n){
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

void Grid1D::setup_ABC(){
  EzBL = Ez[1];
  EzBR = Ez[Nx-2];
  return;
}

void Grid1D::update_ABC(){
  if (lbt == Boundary_t::Mur1){
    Ez[0]    = EzBL + (cc*dt-dx)*(Ez[1]   - Ez[0]  )/(cc*dt+dx);
  }
  if (rbt == Boundary_t::Mur1){
    Ez[Nx-1] = EzBR + (cc*dt-dx)*(Ez[Nx-2]-Ez[Nx-1])/(cc*dt+dx);
  }
}

void Grid1D::open_result_file(){
  handle.open(fname);
  return;
}

void Grid1D::close_result_file(){
  handle.close();
  return;
}

void Grid1D::save_fields(){
  // Print to file
  handle << "\n\n";
  for (size_t i = 0; i < Nx; i++){
    handle << i*dx << " " << Ez[i] << " " << Hy[i] << "\n";
  }
  return;
}

size_t Grid1D::tsteps(){
  return Nt;
}

double Grid1D::gaussian(double t, double x, double f){
  double tau = 0.1/f;
  double t0  = 3.0*tau;
  double c   = (t - x/cc - t0)/tau;
  return exp(-c*c);
}

double Grid1D::sinusoidal(double t, double x, double f){
  return sin(2.0*M_PI*f*(t - x/cc));
}

double Grid1D::ricker(double t, double x, double f){
  double d = 2.5/f;
  double c = M_PI*f*(t - x/cc - d);
  return (1.0 - 2.0*c*c)*exp(-c*c);
}

Grid1D::~Grid1D(){
  if (!epsr) { delete[] epsr; }
  if (!mur)  { delete[] mur;  }
  if (!Ce)   { delete[] Ce;   }
  if (!Ch)   { delete[] Ch;   }
  if (!Ez)   { delete[] Ez;   }
  if (!Hy)   { delete[] Hy;   }
}

