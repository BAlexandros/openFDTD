#include "../include/materials.hpp"

std::map<std::string, Material> materialdb = {
  {"Air",       Material("Air",       1.0,  1.0) },
  {"Glass",     Material("Glass",     10.0, 1.0) },
  {"Teflon",    Material("Teflon",    2.1,  1.0) },
  {"Glycerine", Material("Glycerine", 42.5, 1.0) },
  {"Water",     Material("Water",     80.4, 1.0) },
};

