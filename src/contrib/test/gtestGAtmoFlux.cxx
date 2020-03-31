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
#include "Tools/Flux/GBGLRSAtmoFlux.h"
#include "Tools/Flux/GAtmoFlux.h"
#include "Framework/Conventions/Units.h"
#include <stdlib.h> /* For getenv(). */
#include "TH3D.h"

/* Macro to compute the size of a static C array.
 *
 * See https://stackoverflow.com/questions/1598773. */
#define LEN(x) ((sizeof(x)/sizeof(0[x]))/((size_t)(!(sizeof(x) % sizeof(0[x])))))

typedef int testFunction(char *err);

int isclose(double a, double b, double rel_tol, double abs_tol)
{
    /* Returns 1 if a and b are "close". This algorithm is taken from Python's
     * math.isclose() function.
     *
     * See https://www.python.org/dev/peps/pep-0485/. */
    return fabs(a-b) <= fmax(rel_tol*fmax(fabs(a),fabs(b)),abs_tol);
}

using namespace genie;
using namespace genie::flux;

/* Tests the GetTotalFlux() function. */
int testGetTotalFlux(char *err)
{
  GAtmoFlux *atmo_flux_driver;
  double emin, emax, value, expected;
  char filename[256];

  sprintf(filename, "%s/src/contrib/test/fmax20_i0403z.sno_nue", getenv("GENIE"));

  GBGLRSAtmoFlux *bartol_flux = new GBGLRSAtmoFlux;
  atmo_flux_driver = dynamic_cast<GAtmoFlux *>(bartol_flux);

  emin = -1;
  emax = 1e9;

  // Configure GAtmoFlux options (common to all concrete atmospheric flux drivers)
  // set min/max energy:
  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);
  // set flux files:
  atmo_flux_driver->AddFluxFile(12, filename);
  atmo_flux_driver->LoadFluxData();

  /* Test that GetTotalFlux() is the same as GetFlux(12) since we are only
   * including a single neutrino flavour. */
  if (atmo_flux_driver->GetTotalFlux() != atmo_flux_driver->GetFlux(12)) {
    sprintf(err, "GetTotalFlux() = %f which is not equal to GetFlux(12) = %f", atmo_flux_driver->GetTotalFlux(), atmo_flux_driver->GetFlux(12));
    goto err;
  }

  delete atmo_flux_driver;

  return 0;

err:
  delete atmo_flux_driver;

  return 1;
}

/* Tests the GetTotalFluxInEnergyRange() function. */
int testGetTotalFluxInEnergyRange(char *err)
{
  GAtmoFlux *atmo_flux_driver;
  double emin, emax, value, expected;
  char filename[256];

  sprintf(filename, "%s/src/contrib/test/fmax20_i0403z.sno_nue", getenv("GENIE"));

  GBGLRSAtmoFlux *bartol_flux = new GBGLRSAtmoFlux;
  atmo_flux_driver = dynamic_cast<GAtmoFlux *>(bartol_flux);

  // set flux files:
  atmo_flux_driver->AddFluxFile(12, filename);
  atmo_flux_driver->LoadFluxData();

  emin = -1;
  emax = 1e9;

  // Configure GAtmoFlux options (common to all concrete atmospheric flux drivers)
  // set min/max energy:
  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);

  value = atmo_flux_driver->GetTotalFlux();
  expected = atmo_flux_driver->GetTotalFluxInEnergyRange();

  if (value != expected) {
    sprintf(err, "GetTotalFlux(%.2f,%.2f) = %f which is not equal to the expected total flux = %f", emin, emax, value, expected);
    return 1;
  }

  /* Now set emin and emax both above the bounds and make sure we get 0. */
  emin = 1e9;
  emax = 1e10;

  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);

  value = atmo_flux_driver->GetTotalFluxInEnergyRange();
  expected = 0;

  if (value != expected) {
    sprintf(err, "GetTotalFlux(%.1e,%.1e) = %f, but expected %f!", emin, emax, value, expected);
    return 1;
  }

  /* Now set emin and emax both below the bounds and make sure we get 0. */
  emin = 0;
  emax = 0.01;

  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);

  value = atmo_flux_driver->GetTotalFluxInEnergyRange();
  expected = 0;

  if (value != expected) {
    sprintf(err, "GetTotalFlux(%.1e,%.1e) = %f, but expected %f!", emin, emax, value, expected);
    return 1;
  }

  /* Now we test when both emin and emax are in the same bin. */
  emin = 0.106;
  emax = 0.11;

  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);

  value = atmo_flux_driver->GetTotalFluxInEnergyRange();
  expected = atmo_flux_driver->GetFlux(12,emin)*(emax-emin);

  if (!isclose(value,expected,1e-5,0)) {
  sprintf(err, "GetTotalFlux(%.3f,%.3f) = %f, but expected %f!", emin, emax, value, expected);
    return 1;
  }

  /* Now we test when emin and emax are just past the low and high bin edges. */
  emin = 0.10 + 1e-10;
  emax = 10.0 - 1e-10;

  atmo_flux_driver->ForceMinEnergy(emin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy(emax * units::GeV);

  value = atmo_flux_driver->GetTotalFluxInEnergyRange();
  expected = atmo_flux_driver->GetTotalFlux();

  if (!isclose(value,expected,1e-5,0)) {
  sprintf(err, "GetTotalFlux(%.3f,%.3f) = %f, but expected %f!", emin, emax, value, expected);
    return 1;
  }

  delete atmo_flux_driver;

  return 0;

err:
  delete atmo_flux_driver;

  return 1;
}

struct tests {
    testFunction *test;
    const char *name;
} tests[] = {
    {testGetTotalFlux, "testGetTotalFlux"},
    {testGetTotalFluxInEnergyRange, "testGetTotalFluxInEnergyRange"},
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
