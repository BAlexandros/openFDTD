#include "../include/grid1D.hpp"

int main(void)
{
  Grid1Dsettings gs;
  gs.Nx     = 140;
  gs.Sc     = 1.0;
  gs.tfsfL  = 30;
  gs.f      = 1.0e9;
  gs.leftBound  = Boundary_t::Mur1;
  gs.rightBound = Boundary_t::Mur1;
  gs.sourcetype = Source_t::Ricker;
  
  Grid1D g(gs);
  g.init_vacuum();
  g.add_material(80,120,"Glass");
  g.add_material(85,115,"Glycerin");

  g.update_coefs();
  g.open_result_file();

  g.init_kernels();

  for (size_t n = 0; n < g.tsteps(); n++){

    g.update_magnetic();
    g.update_tfsf_magnetic(n);

    g.setup_ABC();
    g.update_electric();
    g.update_ABC();
    g.update_tfsf_electric(n);

    g.update_kernels(n);

    g.save_results();
  }

  g.finalize_kernels();
  g.save_spectrum();

  g.close_result_file();
  
  return 0;
}

