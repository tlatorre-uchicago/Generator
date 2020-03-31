//____________________________________________________________________________
/*!

\program gtestAlgorithms

\brief   Program used for testing / debugging GENIE's Algorithms & AlgFactory

\author  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

\created October 26, 2006

\cpright Copyright (c) 2003-2020, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________
//
#include <stdio>
#include "Tools/Flux/GFLUKAAtmoFlux.h"

typedef int testFunction(char *err);

int testGetTotalFlux(void)
{
  GAtmoFlux *atmo_flux_driver;
  double emin, emax;

  GFLUKAAtmoFlux *fluka_flux = new GFLUKAAtmoFlux;
  atmo_flux_driver = dynamic_cast<GAtmoFlux *>(fluka_flux);

  emin = -1;
  emax = 1e9;

  // Configure GAtmoFlux options (common to all concrete atmospheric flux drivers)
  // set min/max energy:
  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);
  // set flux files:
  atmo_flux_driver->AddFluxFile(12, "fmax20_i0403z.sno_nue");
  atmo_flux_driver->LoadFluxData();

  if (atmo_flux_driver->GetTotalFlux() == atmo_flux_driver->GetTotalFlux(emin,emax))
    return 0;

  return 1;
}

struct tests {
    testFunction *test;
    char *name;
} tests[] = {
    {testGetTotalFlux, "testGetTotalFlux"},
};

int main(int argc, char **argv)
{
  int i;
  char err[256];
  int retval = 0;
  struct tests test;

  for (i = 0; i < LEN(tests); i++) {
    test = tests[i];

    if (!test.test(err)) {
      printf("[\033[92mok\033[0m] %s\n", test.name);
    } else {
      printf("[\033[91mfail\033[0m] %s: %s\n", test.name, err);
      retval = 1;
    }
  }

  return retval;
}
