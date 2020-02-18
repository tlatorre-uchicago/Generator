//____________________________________________________________________________
/*!

\class   genie::flux::GBGLRSPlusFLUKAAtmoFlux

\brief   A flux driver for a combined version of the Bartol Atmospheric
         Neutrino Flux and the Battistoni flux at low energies.

\ref     Please note that this class expects to read flux files formatted as
         described in the above BGLRS flux page.
         Each file contains 5 columns:
         - neutrino energy (GeV) at bin centre
         - neutrino cos(zenith angle) at bin centre
         - neutrino flux dN/dE (#neutrinos /m^2 /sec /sr)

\author  Anthony LaTorre <tlatorre at uchicago dot edu>
         University of Chicago

based heavily on GBGLRSAtmoFlux by Christopher Backhouse.

\created February 17, 2020

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_BGLRS_PLUS_FLUKA_ATMO_FLUX_H_
#define _GENIE_BGLRS_PLUS_FLUKA_ATMO_FLUX_H_

#include "FluxDrivers/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos(zenith) and energy bins in flux simulation
const unsigned int kBGLRSPlusFLUKA3DNumCosThetaBins           = 20;
const double       kBGLRSPlusFLUKA3DCosThetaMin               = -1.0;
const double       kBGLRSPlusFLUKA3DCosThetaMax               =  1.0;
/* Energy bins are equivalent to np.logspace(-2,1,61). */
const unsigned int kBGLRSPlusFLUKA3DNumLogEvBins              = 60;
const unsigned int kBGLRSPlusFLUKA3DNumLogEvBinsPerDecade     = 20;
const double       kBGLRSPlusFLUKA3DEvMin                     = 0.01; // GeV

class GBGLRSPlusFLUKAAtmoFlux: public GAtmoFlux {

public :
  GBGLRSPlusFLUKAAtmoFlux();
 ~GBGLRSPlusFLUKAAtmoFlux();

  // Most implementation is derived from the base GAtmoFlux
  // The concrete driver is only required to implement a function for
  // loading the input data files

private:

  void SetBinSizes   (void);
  bool FillFluxHisto (int nu_pdg, string filename);
};

} // flux namespace
} // genie namespace

#endif // _GBARTOL_ATMO_FLUX_I_H_

