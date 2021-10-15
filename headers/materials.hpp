#ifndef __MATERIALS_HEADER__
#define __MATERIALS_HEADER__

#include <map>

struct Material {
  std::string name;
  double      epsr;
  double      mur;

  Material(const std::string _name = "Air", double _er = 1, double _mr = 1)
           :name(_name), epsr(_er), mur(_mr){}
};

extern std::map<std::string, Material> materials;

std::map<std::string, Material> materials = {
  {"Air",   Material("Air"  , 1.0,  1.0)},
  {"PTFE",  Material("PTFE" , 2.1,  1.0)},
  {"Glass", Material("Glass", 10.0, 1.0)},
};

#endif
