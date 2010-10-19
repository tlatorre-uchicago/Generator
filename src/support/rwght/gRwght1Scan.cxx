//____________________________________________________________________________
/*!

\program grwght1scan

\brief   Generates weights given an input GHEP event file and for a given
         systematic parameter (supported by the ReWeight package).
         It outputs a ROOT file containing a tree with an entry for every 
         input event. Each such tree entry contains a TArrayF of all computed 
         weights and a TArrayF of all used tweak dial values. 
         Is a RAL/T2K analysis program.

\syntax  grwght1scan \
           -f filename [-n n1[,n2]] -s systematic -t n_twk_diall_values
           [-p neutrino_codes]

         where 
         [] is an optional argument

         -f Specifies a GHEP input file.
         -n Specifies an event range.
            Examples:
            - Type `-n 50,2350' to process all 2301 events from 50 up to 2350.
              Note: Both 50 and 2350 are included.
            - Type `-n 1000' to process the first 1000 events;
              from event number 0 up to event number 999.
            This is an optional argument. By default GENIE will process all events.
         -t Specified the number of tweak dial values between -1 and 1 
            (must be odd so as to include al -1, 0 and 1 / if it is an even
             number it will be incremented by 1)
         -s Specifies the systematic param to tweak.
            See $GENIE/src/ReWeight/GSyst.h for a list of parameters and
            their corresponding label, which is what should be input here.
         -p If set, grwght1scan reweights *only* the specified neutrino 
            species. The input is a comma separated list of PDG codes.
            This is an optional argument. By default GENIE will reweight
            interactions of all neutrino species.

\author  Jim Dobson
         Imperial College London

         Costas Andreopoulos, Jelena Ilic, Nick Grant
         STFC, Rutherford Appleton Laboratory

\created June 10, 2010

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <cassert>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TArrayF.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSyst.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"  
#include "ReWeight/GReWeightNuXSecDIS.h"    
#include "Utils/CmdLnArgParser.h"

using std::string;

using namespace genie;
using namespace genie::rew;

void GetCommandLineArgs (int argc, char ** argv);
void GetEventRange      (Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast);

string      gOptInpFilename; ///< filename for input event tree
Long64_t    gOptNEvt1;       ///< range of events to process (1st input, if any)
Long64_t    gOptNEvt2;       ///< range of events to process (2nd input, if any)
GSyst_t     gOptSyst;        ///< input systematic param
int         gOptInpNTwk;     ///< # of tweaking dial values between [-1,1]
PDGCodeList gOptNu(false);   ///< neutrinos to consider

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // Get the input event sample

  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");
  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("RewScan1", pNOTICE) 
    << "Input tree header: " << *thdr;

  if(!tree){
    LOG("RewScan1", pFATAL) 
      << "Can't find a GHEP tree in input file: "<< file.GetName();
    gAbortingInErr = true;
    exit(1);
  }

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  Long64_t nev_in_file = tree->GetEntries();

  // The tweaking dial takes N values between [-1,1]

  const int   n_points      = gOptInpNTwk;
  const float twk_dial_min  = -1.0;
  const float twk_dial_max  =  1.0;
  const float twk_dial_step = (twk_dial_max - twk_dial_min) / (n_points-1);

  // Work-out the range of events to process
  Long64_t nfirst = 0;
  Long64_t nlast  = 0;
  GetEventRange(nev_in_file, nfirst, nlast);

  Long64_t nev = (nlast - nfirst + 1);

  //
  // Summarize
  //

  LOG("RewScan1", pNOTICE) 
    << "\n"
    << "\n** grwght1scan: Will start processing events promptly."
    << "\nHere is a summary of inputs: "
    << "\n - Input event file: " << gOptInpFilename 
    << "\n - Processing: " << nev << " events in the range [" << nfirst << ", " << nlast << "]"
    << "\n - Systematic parameter to tweak: " << GSyst::AsString(gOptSyst)
    << "\n - Number of tweak dial values in [-1,1] : " << gOptInpNTwk
    << "\n - Neutrino species to reweight : " << gOptNu
    << "\n\n";


  // Declare the weights and twkdial arrays 
  const int n_events = (const int) nev;
  float weights  [n_events][n_points]; 
  float twkdials [n_events][n_points]; 

  // Create a GReWeight object and add to it a set of weight calculators

  GReWeight rw;
  rw.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  rw.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  rw.AdoptWghtCalc( "xsec_coh",        new GReWeightNuXSecCOH       );
  rw.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
  rw.AdoptWghtCalc( "nuclear_qe",      new GReWeightFGM             );
  rw.AdoptWghtCalc( "nuclear_dis",     new GReWeightDISNuclMod      );
  rw.AdoptWghtCalc( "hadro_res_decay", new GReWeightResonanceDecay  );
  rw.AdoptWghtCalc( "hadro_fzone",     new GReWeightFZone           );
  rw.AdoptWghtCalc( "hadro_intranuke", new GReWeightINuke           );
  rw.AdoptWghtCalc( "hadro_agky",      new GReWeightAGKY            );
  rw.AdoptWghtCalc( "xsec_ccqe_vec",   new GReWeightNuXSecCCQEvec   );
  rw.AdoptWghtCalc( "xsec_ncres",      new GReWeightNuXSecNCRES     );
  rw.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );

  // Get GSystSet and include the (single) input systematic parameter

  GSystSet & syst = rw.Systematics();
  syst.Init(gOptSyst);

  string syst_name = GSyst::AsString(gOptSyst);

  // Twk dial loop
  for(int ith_dial = 0; ith_dial < n_points; ith_dial++){  

     // Set non-default values and re-configure.    
     double twk_dial = twk_dial_min + ith_dial * twk_dial_step;  
     LOG("RewScan1", pNOTICE) 
       << "Reconfiguring systematic: " << syst_name 
       << " - Setting tweaking dial to: " << twk_dial;
     syst.Set(gOptSyst, twk_dial);
     rw.Reconfigure();

     // Event loop
     for(int i = nfirst; i <= nlast; i++) {

          // Get next event
          tree->GetEntry(i);
          EventRecord & event = *(mcrec->event);
          LOG("RewScan1", pDEBUG) << event;

          // Reset arrays
          weights  [i][ith_dial] = -99999.0;
          twkdials [i][ith_dial] = twk_dial;

          // Reweight this event?
          int nupdg = event.Probe()->Pdg();
          bool do_reweight = gOptNu.ExistsInPDGCodeList(nupdg);

          // Calculate weight
          double wght=1.;
          if(do_reweight) {
	     wght = rw.CalcWeight(event);
          }

          // Print/store
          LOG("RewScan1", pDEBUG) 
              << "Overall weight = " << wght;
          weights[i][ith_dial] = wght;

          if(i%100 == 0) {
              LOG("RewScan1", pNOTICE) 
                 << "***** Processed "<< i+1 << " events";
           }

           // Clean-up
           mcrec->Clear();

      } // evt loop
  } // twk_dial loop

  // Close event file
  file.Close();

  //
  // Save weights 
  //

  // Make an output tree for saving the weights. As only considering 
  // varying a single systematic use this for name of tree.
  TString wght_filename; 
  wght_filename.Form("weights_%s.root", syst_name.c_str());
  TFile * wght_file = new TFile(wght_filename.Data(), "RECREATE");
  TTree * wght_tree = new TTree(syst_name.c_str(),   "weights tree");
  int branch_eventnum = 0;
  TArrayF * branch_weight_array   = new TArrayF(n_points);
  TArrayF * branch_twkdials_array = new TArrayF(n_points);  
  wght_tree->Branch("eventnum", &branch_eventnum);
  wght_tree->Branch("weights",  &branch_weight_array);
  wght_tree->Branch("twkdials", &branch_twkdials_array);

  for(int i = nfirst; i <= nlast; i++) {
    branch_eventnum = i; 
    for(int ith_dial = 0; ith_dial < n_points; ith_dial++){  
        LOG("RewScan1", pDEBUG)
          << "Filling tree with wght = " << weights[i][ith_dial] 
          << ", twk dial = "<< twkdials[i][ith_dial];
       branch_weight_array   -> AddAt (weights [i][ith_dial], ith_dial);
       branch_twkdials_array -> AddAt (twkdials[i][ith_dial], ith_dial);
    } // twk_dial loop
    wght_tree->Fill();
  } 

  wght_file->cd();
  wght_tree->Write();
  delete wght_tree; 
  wght_tree = 0; 
  wght_file->Close();

  LOG("RewScan1", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("RewScan1", pINFO) 
     << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if(parser.OptionExists('f')) {
    LOG("RewScan1", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("RewScan1", pFATAL) 
        << "Unspecified input filename - Exiting";
    gAbortingInErr = true;
    exit(1);
  }

  // range of event numbers to process
  if ( parser.OptionExists('n') ) {
    //
    LOG("gevdump", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("gevdump", pFATAL) << "Invalid syntax";
         gAbortingInErr = true;
         exit(1);
      }
      // User specified a comma-separated set of values n1,n2.
      // Use [n1,n2] as the event range to process.
      gOptNEvt1 = vecn[0];
      gOptNEvt2 = vecn[1];
    } else {
      // User specified a single number n.
      // Use [0,n] as the event range to process.
      gOptNEvt1 = -1;
      gOptNEvt2 = parser.ArgAsLong('n');
    }  
  } else {
    LOG("gevdump", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt1 = -1;
    gOptNEvt2 = -1;
  }
  LOG("RewScan1", pDEBUG) 
    << "Input event range: " << gOptNEvt1 << ", " << gOptNEvt2;

  // get the number of tweak dials to scan
  if(parser.OptionExists('t')) {
    LOG("RewScan1", pINFO) 
       << "Reading number of tweak dial values";
    gOptInpNTwk = parser.ArgAsInt('t');
    if(gOptInpNTwk % 2 == 0)
    { 
      gOptInpNTwk+=1;
    }
    if(gOptInpNTwk < 3)
    { 
      LOG("RewScan1", pFATAL)
	 << "Specified number of tweak dial is too low, min value is 3 - Exiting";
      gAbortingInErr = true;
      exit(1);
    } 
  } else {
     LOG("RewScan1", pFATAL) 
       << "Unspecified number of tweak dials - Exiting";
     gAbortingInErr = true;
     exit(1);
  }

  // get the systematics
  if(parser.OptionExists('s')) {
   LOG("RewScan1", pINFO) 
      << "Reading input systematic parameter";
   string systematic = parser.ArgAsString('s');
   gOptSyst = GSyst::FromString(systematic);
   if(gOptSyst == kNullSystematic) {
      LOG("RewScan1", pFATAL) << "Unknown systematic: " << systematic;
      gAbortingInErr = true;
      exit(1);
   }
  } else {
    LOG("RewScan1", pFATAL) 
       << "You need to specify a systematic param using -s";
    gAbortingInErr = true;
    exit(1);
  }

  // which species to reweight?
  if(parser.OptionExists('p')) {
   LOG("RewScan1", pINFO) 
      << "Reading input list of neutrino codes";
   vector<int> vecpdg = parser.ArgAsIntTokens('p',",");
   if(vecpdg.size()==0) {
      LOG("RewScan1", pFATAL) 
         << "Empty list of neutrino codes!?";
      gAbortingInErr = true;
      exit(1);
   }
   vector<int>::const_iterator it = vecpdg.begin();
   for( ; it!=vecpdg.end(); ++it) {
     gOptNu.push_back(*it);
   }
  } else {
    LOG("RewScan1", pINFO) 
       << "Considering all neutrino species";
    gOptNu.push_back (kPdgNuE      );
    gOptNu.push_back (kPdgAntiNuE  );
    gOptNu.push_back (kPdgNuMu     );
    gOptNu.push_back (kPdgAntiNuMu );
    gOptNu.push_back (kPdgNuTau    );
    gOptNu.push_back (kPdgAntiNuTau);
  }
}
//_________________________________________________________________________________
void GetEventRange(Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast)
{
  nfirst = 0;
  nlast  = 0;

  if(gOptNEvt1>=0 && gOptNEvt2>=0) {
    // Input was `-n N1,N2'.
    // Process events [N1,N2].
    // Note: Incuding N1 and N2.
    nfirst = gOptNEvt1;
    nlast  = TMath::Min(nev_in_file-1, gOptNEvt2);
  }
  else 
  if(gOptNEvt1<0 && gOptNEvt2>=0) {
    // Input was `-n N'.
    // Process first N events [0,N). 
    // Note: Event N is not included.
    nfirst = 0;
    nlast  = TMath::Min(nev_in_file-1, gOptNEvt2-1);
  }
  else
  if(gOptNEvt1<0 && gOptNEvt2<0) {
    // No input. Process all events.
    nfirst = 0;
    nlast  = nev_in_file-1;    
  }

  assert(nfirst < nlast && nfirst >= 0 && nlast <= nev_in_file-1);
}
//_________________________________________________________________________________
