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
  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
    int neutrino_code = file_iter->first;
    string filename   = file_iter->second;
    atmo_flux_driver->AddFluxFile(neutrino_code, filename);
  }
  atmo_flux_driver->LoadFluxData();
  // configure flux generation surface:
  atmo_flux_driver->SetRadii(gOptRL, gOptRT);
  // set rotation for coordinate tranformation from the topocentric horizontal
  // system to a user-defined coordinate system:
  if(!gOptRot.IsIdentity()) {
     atmo_flux_driver->SetUserCoordSystem(gOptRot);
  }


int main(int argc, char **argv)
{
  testReconfigInCommonPool();
  testReconfigInOwnedModules();

  return 0;
}
