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
#include <cstdio>
#include "Tools/Flux/GFLUKAAtmoFlux.h"
#include "Tools/Flux/GAtmoFlux.h"
#include "Framework/Conventions/Units.h"
#include <stdlib.h> /* For getenv(). */

/* Macro to compute the size of a static C array.
 *
 * See https://stackoverflow.com/questions/1598773. */
#define LEN(x) ((sizeof(x)/sizeof(0[x]))/((size_t)(!(sizeof(x) % sizeof(0[x])))))

typedef int testFunction(char *err);

using namespace genie;
using namespace genie::flux;

int testGetTotalFlux(char *err)
{
  GAtmoFlux *atmo_flux_driver;
  double emin, emax;
  char filename[256];

  fprintf(stderr, "blah\n");

  sprintf(filename, "%s/src/contrib/test/fmax20_i0403z.sno_nue", getenv("GENIE"));

  fprintf(stderr, "creating atmoflux\n");

  GFLUKAAtmoFlux *fluka_flux = new GFLUKAAtmoFlux;
  atmo_flux_driver = dynamic_cast<GAtmoFlux *>(fluka_flux);

  emin = -1;
  emax = 1e9;

  fprintf(stderr, "forcing min energy\n");
  // Configure GAtmoFlux options (common to all concrete atmospheric flux drivers)
  // set min/max energy:
  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);
  // set flux files:
  fprintf(stderr, "adding flux file\n");
  atmo_flux_driver->AddFluxFile(12, filename);
  fprintf(stderr, "loading data\n");
  atmo_flux_driver->LoadFluxData();
  fprintf(stderr, "done loading data\n");

  if (atmo_flux_driver->GetTotalFlux() != atmo_flux_driver->GetTotalFlux(emin,emax)) {
    sprintf(err, "GetTotalFlux(%.2f,%.2f) = %f which is not equalt to the expected total flux = %f", emin, emax, atmo_flux_driver->GetTotalFlux(emin,emax), atmo_flux_driver->GetTotalFlux());
    return 1;
  }

  return 0;
}

struct tests {
    testFunction *test;
    const char *name;
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
