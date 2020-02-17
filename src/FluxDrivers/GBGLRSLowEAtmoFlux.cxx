//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Anthony LaTorre <tlatorre at uchicago dot edu>
         University of Chicago

 based on GBGLRSAtmoFlux.cxx by 

 Author: Christopher Backhouse <c.backhouse1@physics.ox.ac.uk>
         Oxford University
*/
//____________________________________________________________________________

#include <fstream>
#include <cassert>

#include <TH3D.h>
#include <TMath.h>

#include "FluxDrivers/GBGLRSLowEAtmoFlux.h"
#include "Messenger/Messenger.h"
#include "Conventions/Constants.h"

#include "FluxDrivers/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GBGLRSLowEAtmoFlux,genie::flux::GBGLRSLowEAtmoFlux)

using std::ifstream;
using std::ios;

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

GBGLRSLowEAtmoFlux::GBGLRSLowEAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
     << "Instantiating the GENIE BGLRS + low energy atmospheric neutrino flux driver";

  this->Initialize();
  this->SetBinSizes();
}

GBGLRSLowEAtmoFlux::~GBGLRSLowEAtmoFlux()
{

}

/* Generate the correct cos(theta) and energy bin sizes.
 *
 * The cos(theta) bins are equivalent to np.linspace(-1,1,21) and the energy
 * bins are equivalent to np.logspace(-2,1,61). */
void GBGLRSLowEAtmoFlux::SetBinSizes(void)
{
  unsigned int i;
  double costheta, logE, dcostheta, logEmin, dlogE;

  costheta = kBGLRSLowE3DCosThetaMin;

  fPhiBins       = new double [2];
  fCosThetaBins  = new double [kBGLRSLowE3DNumCosThetaBins + 1];
  fEnergyBins    = new double [kBGLRSLowE3DNumLogEvBins + 1];

  fPhiBins[0] = 0;
  fPhiBins[1] = 2.*kPi;

  dcostheta = (kBGLRSLowE3DCosThetaMax - kBGLRSLowE3DCosThetaMin)/(double) kBGLRSLowE3DNumCosThetaBins;
     
  logEmin = TMath::Log10(kBGLRSLowE3DEvMin);
  dlogE = 1.0/(double) kBGLRSLowE3DNumLogEvBinsPerDecade;

  costheta = kBGLRSLowE3DCosThetaMin;
  for (i = 0; i <= kBGLRSLowE3DNumCosThetaBins; i++) {
    if (i > 0)
      costheta += dcostheta;
    fCosThetaBins[i] = costheta;
    if (i != kBGLRSLowE3DNumCosThetaBins) {
      LOG("Flux", pDEBUG)
        << "FLUKA 3d flux: CosTheta bin " << i+1
        << ": lower edge = " << fCosThetaBins[i];
    } else {
      LOG("Flux", pDEBUG)
        << "FLUKA 3d flux: CosTheta bin " << kBGLRSLowE3DNumCosThetaBins
        << ": upper edge = " << fCosThetaBins[kBGLRSLowE3DNumCosThetaBins];
    }
  }
     
  logE = logEmin;
  for (i = 0; i <= kBGLRSLowE3DNumLogEvBins; i++) {
    if (i > 0)
      logE += dlogE;
    fEnergyBins[i] = TMath::Power(10.0, logE);
    if (i != kBGLRSLowE3DNumLogEvBins) {
      LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: Energy bin " << i+1
         << ": lower edge = " << fEnergyBins[i];
    } else {
      LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: Energy bin " << kBGLRSLowE3DNumLogEvBins
         << ": upper edge = " << fEnergyBins[kBGLRSLowE3DNumLogEvBins];
    }
  }

  fNumPhiBins      = 1;
  fNumCosThetaBins = kBGLRSLowE3DNumCosThetaBins;
  fNumEnergyBins   = kBGLRSLowE3DNumLogEvBins; 
  fMaxEv = fEnergyBins[fNumEnergyBins];
}

bool GBGLRSLowEAtmoFlux::FillFluxHisto(int nu_pdg, string filename)
{
  unsigned int i;
  int ibin;
  double energy, costheta, flux;
  double junkd; // throw away 4th column

  LOG("Flux", pNOTICE) 
    << "Loading BGLRS low energy flux for neutrino: " << nu_pdg 
    << " from file: " << filename;

  TH3D* histo = 0;
  std::map<int,TH3D*>::iterator myMapEntry = fRawFluxHistoMap.find(nu_pdg);
  if (myMapEntry != fRawFluxHistoMap.end()){
      histo = myMapEntry->second;
  }

  if (!histo) {
     LOG("Flux", pERROR) << "Null flux histogram!";
     return false;
  }

  ifstream flux_stream(filename.c_str(), ios::in);

  if (!flux_stream) {
     LOG("Flux", pERROR) << "Could not open file: " << filename;
     return false;
  }

  // throw away comment line
  flux_stream.ignore(99999, '\n');

  i = 0;
  while (!flux_stream.eof()) {
    i += 1;

    flux = 0.0;
    flux_stream >> energy >> costheta >> flux >> junkd;
    if (flux > 0) {
      LOG("Flux", pINFO)
        << "Flux[Ev = " << energy 
        << ", cos8 = " << costheta << "] = " << flux;
      ibin = histo->FindBin((Axis_t) energy, (Axis_t) costheta, (Axis_t) kPi);
      histo->SetBinContent(ibin, (Stat_t)(flux));
    } else {
      LOG("Flux", pERROR) << "Flux on line " << i << " is " << flux << " which is negative!";
    }
  }
  return true;
}
