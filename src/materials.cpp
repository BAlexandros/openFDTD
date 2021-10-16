#include "../include/materials.hpp"

std::map<std::string, Material> materialdb = {
  {"Air",   Material("Air"  , 1.0,  1.0) },
  {"PTFE",  Material("PTFE" , 2.1,  1.0) },
  {"Glass", Material("Glass", 10.0, 1.0) },
};

