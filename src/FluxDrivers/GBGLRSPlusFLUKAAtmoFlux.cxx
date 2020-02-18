//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Anthony LaTorre <tlatorre at uchicago dot edu>
         University of Chicago

 based on GBGLRSAtmoFlux.cxx by Christopher Backhouse.
*/
//____________________________________________________________________________

#include <fstream>
#include <cassert>

#include <TH3D.h>
#include <TMath.h>

#include "FluxDrivers/GBGLRSPlusFLUKAAtmoFlux.h"
#include "Messenger/Messenger.h"
#include "Conventions/Constants.h"

#include "FluxDrivers/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GBGLRSPlusFLUKAAtmoFlux,genie::flux::GBGLRSPlusFLUKAAtmoFlux)

using std::ifstream;
using std::ios;
using std::string;
using std::istringstream;

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

GBGLRSPlusFLUKAAtmoFlux::GBGLRSPlusFLUKAAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
     << "Instantiating the GENIE BGLRS + FLUKA atmospheric neutrino flux driver";

  this->Initialize();
  this->SetBinSizes();
}

GBGLRSPlusFLUKAAtmoFlux::~GBGLRSPlusFLUKAAtmoFlux()
{
  if (fPhiBins) delete fPhiBins;
  if (fCosThetaBins) delete fCosThetaBins;
  if (fEnergyBins) delete fEnergyBins;
}

/* Generate the correct cos(theta) and energy bin sizes.
 *
 * The cos(theta) bins are equivalent to np.linspace(-1,1,21) and the energy
 * bins are equivalent to np.logspace(-2,1,61). */
void GBGLRSPlusFLUKAAtmoFlux::SetBinSizes(void)
{
  unsigned int i;
  double logE, dcostheta, logEmin, dlogE;

  if (fPhiBins) delete fPhiBins;
  if (fCosThetaBins) delete fCosThetaBins;
  if (fEnergyBins) delete fEnergyBins;

  fPhiBins       = new double [2];
  fCosThetaBins  = new double [kBGLRSPlusFLUKA3DNumCosThetaBins + 1];
  fEnergyBins    = new double [kBGLRSPlusFLUKA3DNumLogEvBins + 1];

  fPhiBins[0] = 0;
  fPhiBins[1] = 2.*kPi;

  dcostheta = (kBGLRSPlusFLUKA3DCosThetaMax - kBGLRSPlusFLUKA3DCosThetaMin)/(double) kBGLRSPlusFLUKA3DNumCosThetaBins;

  logEmin = TMath::Log10(kBGLRSPlusFLUKA3DEvMin);
  dlogE = 1.0/(double) kBGLRSPlusFLUKA3DNumLogEvBinsPerDecade;

  fCosThetaBins[0] = kBGLRSPlusFLUKA3DCosThetaMin;
  for (i = 0; i <= kBGLRSPlusFLUKA3DNumCosThetaBins; i++) {
    if (i > 0)
      fCosThetaBins[i] = fCosThetaBins[i-1] + dcostheta;
    if (i != kBGLRSPlusFLUKA3DNumCosThetaBins) {
      LOG("Flux", pDEBUG)
        << "BGLRSPlusFLUKA 3d flux: CosTheta bin " << i+1
        << ": lower edge = " << fCosThetaBins[i];
    } else {
      LOG("Flux", pDEBUG)
        << "BGLRSPlusFLUKA 3d flux: CosTheta bin " << kBGLRSPlusFLUKA3DNumCosThetaBins
        << ": upper edge = " << fCosThetaBins[kBGLRSPlusFLUKA3DNumCosThetaBins];
    }
  }

  logE = logEmin;
  for (i = 0; i <= kBGLRSPlusFLUKA3DNumLogEvBins; i++) {
    if (i > 0)
      logE += dlogE;
    fEnergyBins[i] = TMath::Power(10.0, logE);
    if (i != kBGLRSPlusFLUKA3DNumLogEvBins) {
      LOG("Flux", pDEBUG)
         << "BGLRSPlusFLUKA 3d flux: Energy bin " << i+1
         << ": lower edge = " << fEnergyBins[i];
    } else {
      LOG("Flux", pDEBUG)
         << "BGLRSPlusFLUKA 3d flux: Energy bin " << kBGLRSPlusFLUKA3DNumLogEvBins
         << ": upper edge = " << fEnergyBins[kBGLRSPlusFLUKA3DNumLogEvBins];
    }
  }

  fNumPhiBins      = 1;
  fNumCosThetaBins = kBGLRSPlusFLUKA3DNumCosThetaBins;
  fNumEnergyBins   = kBGLRSPlusFLUKA3DNumLogEvBins;
  fMaxEv = fEnergyBins[fNumEnergyBins];
}

bool GBGLRSPlusFLUKAAtmoFlux::FillFluxHisto(int nu_pdg, string filename)
{
  unsigned int i;
  int ibin;
  double energy, costheta, flux;
  string str;
  TH3D* histo;

  LOG("Flux", pNOTICE)
    << "Loading BGLRS low energy flux for neutrino: " << nu_pdg
    << " from file: " << filename;

  std::map<int,TH3D*>::iterator myMapEntry = fRawFluxHistoMap.find(nu_pdg);
  if (myMapEntry != fRawFluxHistoMap.end()){
      histo = myMapEntry->second;
  }

  if (!histo) {
     LOG("Flux", pERROR) << "Null flux histogram!";
     return false;
  }

  ifstream flux_stream(filename.c_str(), ios::in);

  if (!flux_stream || !flux_stream.good()) {
     LOG("Flux", pERROR) << "Error opening file: " << filename;
     return false;
  }

  i = 0;
  while (getline(flux_stream,str)) {
    i += 1;

    /* Skip comments. */
    if (str.length() > 0 && str.at(0) == '#') continue;

    /* Loop over each line and parse the energy, cos(theta), and flux.
     *
     * See https://stackoverflow.com/questions/21334173/read-from-file-in-c-till-end-of-line. */
    istringstream ss(str);

    flux = 0.0;
    ss >> energy >> costheta >> flux;
    if (flux > 0) {
      LOG("Flux", pDEBUG)
        << "Flux[Ev = " << energy
        << ", cos = " << costheta << "] = " << flux;
      ibin = histo->FindBin((Axis_t) energy, (Axis_t) costheta, (Axis_t) kPi);
      histo->SetBinContent(ibin, (Stat_t)(flux));
    } else {
      LOG("Flux", pERROR) << "Flux on line " << i << " is " << flux << " which is negative!";
    }
  }

  return true;
}
