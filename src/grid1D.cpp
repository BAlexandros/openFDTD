#include "../include/grid1D.hpp"

Grid1D::Grid1D(){
  st = Source_t::Gaussian;
  sourcefunc = &Grid1D::gaussian;
  Nl = 20;
  Nd = 20;
  dm = 1e10;
  dx = 1e10;
  tfsfL = 30;
  Sc = 1.0;
  lbt = Boundary_t::Mur1;
  rbt = Boundary_t::Mur1;
  field_fname = "";
  spectrum_fname = "";

}

void Grid1D::init_fields(){
  Ez = new double[Nx]();
  Hy = new double[Nx]();
  return;
}

void Grid1D::init_coefs(){
  epsr = new double[Nx]();
  mur  = new double[Nx]();
  Ce = new double[Nx]();
  Ch = new double[Nx]();
  #pragma omp parallel for if(parallelism_enabled) \
                           firstprivate(dx,dt)\
                           shared(epsr,mur,Ce,Ch)
  for (int i = 0; i < Nx; i++){
    epsr[i] = 1.0;
    mur[i]  = 1.0;
    Ce[i] = dt / (eps0*epsr[i]) / dx;
    Ch[i] = dt / (mu0*mur[i])   / dx;
  }
  #pragma omp parallel for if(parallelism_enabled) \
                           firstprivate(dx,dt)\
                           shared(material_list,epsr,mur,Ce,Ch)
  for (int i = 0; i < material_list.size(); i++){
    GridMat cur_mat = material_list[i];
    for (int x = cur_mat.x1; x < cur_mat.x2; x++){
      epsr[x] = materialdb[cur_mat.matname].epsr;
      mur[x]  = materialdb[cur_mat.matname].mur;
    }
    for (int x = cur_mat.x1; x < cur_mat.x2; x++){
      Ce[x] = dt / (eps0*epsr[x]) / dx;
      Ch[x] = dt / (mu0*mur[x])   / dx;
    }
  }
  return;
}

void Grid1D::init_kernels(){
  fm = 0.5/dt;      std::cout<< fm << "\n";
  df = 1.0/Nt/dt;   std::cout<< df << "\n";
  Nf = Nt;       std::cout<< Nt << "\n";
  Ks = new std::complex<double>[Nf]();

  std::complex<double> j(0,1);
  for (int fi = 0; fi < Nf; fi++){
    Ks[fi] = exp(-j*2.0*M_PI*dt*(fi*df));
  }
  return;
}

void Grid1D::init_fourier(){
  Rf = new std::complex<double>[Nf]();
  Tf = new std::complex<double>[Nf]();
  Sf = new std::complex<double>[Nf]();
  return;
}

void Grid1D::update_fourier(double n){

  #pragma omp parallel for if(parallelism_enabled) \
                           shared(Rf,Tf,Sf,Ks,Ez)
  for (int fi = 0; fi < Nf; fi++){
    Rf[fi] += std::pow(Ks[fi],n)*Ez[1];
    Tf[fi] += std::pow(Ks[fi],n)*Ez[Nx-2];
    Sf[fi] += std::pow(Ks[fi],n)*(this->*sourcefunc)(n,0);
  }
  return;
}

void Grid1D::finalize_kernels(){
  for (int fi = 0; fi < Nf; fi++){
    Rf[fi] *= dt;
    Tf[fi] *= dt;
    Sf[fi] *= dt;
  }
}

void Grid1D::add_material(int iL, int iR,
    std::string mname){
  material_list.push_back(GridMat{iL,iR,mname});
  return;
}

void Grid1D::update_magnetic(){
  #pragma omp simd if(parallelism_enabled) 
  for (int i = 0; i < Nx - 1; i++){
    Hy[i] = Hy[i] + Ch[i]*(Ez[i+1] - Ez[i]);
  }
  return;
}

void Grid1D::update_tfsf_magnetic(int n){
  // TFSF correction for Hy
  Hy[tfsfL-1] -= Ch[tfsfL-1]*(this->*sourcefunc)(n,0);
  return;
}

void Grid1D::update_electric(){
  #pragma omp simd if(parallelism_enabled)
  for (int i = 1; i < Nx - 1; i++){
    Ez[i] = Ez[i] + Ce[i]*(Hy[i] - Hy[i-1]);
  }
  return;
}

void Grid1D::update_tfsf_electric(int n){
  Ez[tfsfL] += Ce[tfsfL]*(this->*sourcefunc)(n+0.5,-0.5)/sqrt(mu0*mur[tfsfL]/(eps0*epsr[tfsfL]));
  return;
}

void Grid1D::init_ABC(){
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
  fhandle.open("./data/"+field_fname);
  shandle.open("./data/"+spectrum_fname);
  return;
}

void Grid1D::close_result_file(){
  fhandle.close();
  shandle.close();
  return;
}

void Grid1D::save_cur_field(){
  fhandle << "\n\n";
  for (int i = 0; i < Nx; i++){
    fhandle << i*dx << " " << Ez[i] << " " << Hy[i] << "\n";
  }
  return;
}

void Grid1D::save_spectrum(){
  for (int i = 0; i < Nf; i++){
    double Ri = pow(std::abs(Rf[i]/Sf[i]),2);
    double Ti = pow(std::abs(Tf[i]/Sf[i]),2);
    shandle << i*df << " " <<  Ri     \
                          << " " <<  Ti     \
                          << " " <<  Ri+Ti  \
                          << "\n";
  }
}

void Grid1D::run_simulation(){

  init_coefs();
  std::cout << "Coefficients initialized\n";
  init_fields();
  std::cout << "Fields initialized\n";
  init_kernels();
  std::cout << "Kernels initialized\n";
  init_fourier();
  std::cout << "Fourier initialized\n";
  open_result_file();

  for (int n = 0; n < Nt; n++){

    update_magnetic();
    update_tfsf_magnetic(n);

    init_ABC();
    update_electric();
    update_ABC();
    update_tfsf_electric(n);

    update_fourier(n);

    save_cur_field();
  }

  // finalize_kernels();
  save_spectrum();

  close_result_file();

  makeFieldAnimation();

  std::cout << "Simulation Complete\n";
  return;
}

double Grid1D::gaussian(double q, double m){
  double tau = 0.5/f;
  double t0  = 6.0*tau;
  double c   = (q*dt - m*dx/cc - t0)/tau;
  return exp(-c*c);
}

double Grid1D::sinusoidal(double q, double m){
  return sin(2.0*M_PI*f*(q*dt - m*dx/cc));
}

double Grid1D::ricker(double q, double m){
  double d = 2.5/f;
  double c = M_PI*f*(q*dt - m*dx/cc - d);
  return (1.0 - 2.0*c*c)*exp(-c*c);
}

void Grid1D::enable_parallelism(){
  parallelism_enabled = true;
  return;
}

void Grid1D::disable_parallelism(){
  parallelism_enabled = false;
  return;
}

void Grid1D::setE0(double E0_){
  E0 = E0_;
  return;
}

void Grid1D::setFreq(double f_){
  f = f_;
  T = 1.0/f;
  lambda = cc/f;
  dx = (lambda/Nl < dx) ? lambda/Nl : dx;
  return;
}

void Grid1D::setSourceType(Source_t st_){
  st = st_;
  if (st == Source_t::Gaussian){
    sourcefunc = &Grid1D::gaussian;
  } else if
     (st == Source_t::Ricker){
    sourcefunc = &Grid1D::ricker;
  } else if
     (st == Source_t::Sinusoidal){
    sourcefunc = &Grid1D::sinusoidal;
  } 
  return;
}

void Grid1D::setSc(double Sc_){
  Sc = Sc_;
  dt = Sc*dx/cc;
  return;
}

void Grid1D::setNx(int Nx_){
  Nx = Nx_;
  return;
}

void Grid1D::setNl(int Nl_){
  Nl = Nl_;
  dx = (lambda/Nl < dx) ? lambda/Nl : dx;
  dt = Sc*dx/cc;
  return;
}

void Grid1D::setdm(double dm_){
  dm = dm_;
  dx = (dm/Nd < dx) ? dm/Nd : dx;
  dt = Sc*dx/cc;
  return;
}

void Grid1D::setNd(int Nd_){
  Nd = Nd_;
  dx = (dm/Nd < dx) ? dm/Nd : dx;
  dt = Sc*dx/cc;
  return;
}

void Grid1D::setNt(int Nt_){
  Nt = Nt_;
  return;
}

void Grid1D::setFieldProgFname(std::string fn){
  field_fname = fn;
  return;
}

void Grid1D::setSpectrumFname(std::string fn){
  spectrum_fname = fn;
  return;
}

void Grid1D::makeFieldAnimation(){
  std::ofstream gnuplotscript;
  gnuplotscript.open("tools/plot.g");
  gnuplotscript << "set terminal gif animate delay 4 optimize\n";
  gnuplotscript << "set output './gallery/field.gif'\n";
  gnuplotscript << "stats './data/field.dat' nooutput\n";
  gnuplotscript << "set yr [" << -2*E0 << ":" << 2*E0 << "]\n";
  gnuplotscript << "set xlabel \"Spatial step\"\n";
  gnuplotscript << "set ylabel \"V/m\"\n";
  gnuplotscript << "unset key\n";
  gnuplotscript << "do for [i=1:int(STATS_blocks)] {\n";
  gnuplotscript << "    set title sprintf(\"n = %d\",i+1)\n";
  gnuplotscript << "    p './data/field.dat' index (i-1) u 0:2 w l t 'Ez'\n";
  gnuplotscript << "}";
  gnuplotscript.close();
  system("gnuplot tools/plot.g");
  
  return;
}

Grid1D::~Grid1D(){
  if (!epsr) { delete[] epsr; }
  if (!mur)  { delete[] mur;  }
  if (!Ce)   { delete[] Ce;   }
  if (!Ch)   { delete[] Ch;   }
  if (!Ez)   { delete[] Ez;   }
  if (!Hy)   { delete[] Hy;   }
  if (!Ks)   { delete[] Ks;   }
  if (!Rf)   { delete[] Rf;   }
  if (!Tf)   { delete[] Tf;   }
  if (!Sf)   { delete[] Tf;   }
}
