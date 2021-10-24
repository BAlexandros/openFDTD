#include "../include/grid1D.hpp"

int main(void)
{
  Grid1D g;
  g.enable_parallelism();
  g.setE0(1.0);
  g.setFreq(1.0e9);
  g.setNx(140);
  g.setNl(100);
  g.setSc(1);
  g.setNt(15000);
  g.setSourceType(Source_t::Ricker);
  g.add_material(80,120,"Glass");
  g.add_material(85,115,"Glycerin");
  g.setFieldProgFname("field.dat");
  g.setSpectrumFname("spectrum.dat");

  g.run_simulation();
  
  return 0;
}

