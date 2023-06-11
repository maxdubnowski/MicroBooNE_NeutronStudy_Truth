#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TVector3.h>

#include "Constants.h"

using namespace std;
using namespace Constants;

void NeutrinoSelectionFilter::Loop() {

  //--------------------------------------------------//

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TH1D::SetDefaultSumw2();
   TH2D::SetDefaultSumw2();

   //--------------------------------------------------//

   // Output file

   TString FileName = "MicroBooNE_Truth.root";
   TFile* OutputFile = new TFile(FileName,"recreate");
   std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

   //--------------------------------------------------//

   // Plot declaration

   TH1D* TrueVertexXPlot[NInte];
   TH1D* TrueVertexYPlot[NInte];
   TH1D* TrueVertexZPlot[NInte];

   TH1D* TrueMuonCosThetaPlot[NInte];

   //--------------------------------------------------//

   // Loop over the interaction processes

   for (int inte = 0; inte < NInte; inte++) {

     //--------------------------------------------------//

     // 1D analysis

     TrueVertexXPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
     TrueVertexYPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
     TrueVertexZPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

     TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

     //--------------------------------------------------//

   } // End of the loop over the interaction processes

   //--------------------------------------------------//

   // Loop over the events

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     //--------------------------------------------------//

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;

      //--------------------------------------------------//

      // MC weight to scale events to data pot

      double event_weight = (data_pot / mc_pot) * weightSplineTimesTune;

      //--------------------------------------------------//

      int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0, TruePi0Counter = 0;
      int TrueNeutronCounter = 0, TrueHeavierMesonCounter = 0;
      
      std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
      std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();
      std::vector<int> VectorTrueNeutronIndex; VectorTrueNeutronIndex.clear();

      int NMCParticles = mc_pdg->size();

      for (int WhichMCParticle = 0; WhichMCParticle < NMCParticles; WhichMCParticle++) {

	// MC truth information for the final-state primary particles

	// CC numu events

	if ( ccnc == 0 && nu_pdg == 14) {

	  TVector3 MCParticle(mc_px->at(WhichMCParticle),mc_py->at(WhichMCParticle),mc_pz->at(WhichMCParticle));
	  double MCParticleMomentum = MCParticle.Mag();
	  int MCParticlePdg = mc_pdg->at(WhichMCParticle);

	  if ( MCParticlePdg == MuonPdg && MCParticleMomentum >= ArrayNBinsMuonMomentum[0] ) 
	    { TrueMuonCounter++;  VectorTrueMuonIndex.push_back(WhichMCParticle); }

	  if ( MCParticlePdg == ProtonPdg && MCParticleMomentum >= ArrayNBinsProtonMomentum[0] ) 
	    { TrueProtonCounter++; VectorTrueProtonIndex.push_back(WhichMCParticle); }

	  if ( fabs(MCParticlePdg) == AbsChargedPionPdg && MCParticleMomentum >= ChargedPionMomentumThres ) 
	    { TrueChargedPionCounter++; }

	  if (MCParticlePdg == NeutralPionPdg) { TruePi0Counter++; }


	} // End of the demand stable final state particles and primary interactions

      } // end of the loop over the MCParticles

      //--------------------------------------------------//      

      // CC1p0pi signal events

      bool CC1p0pi = false;

      if (TrueMuonCounter == 1 && TrueProtonCounter == 1 && TrueChargedPionCounter == 0 && TruePi0Counter == 0 && TrueHeavierMesonCounter == 0) {

	CC1p0pi = true;

      }

      if (!CC1p0pi) { continue; }

      //--------------------------------------------------//

      // Define muon / proton vectors & deltapt

      //--------------------------------------------------//

      // Use only true CC1p0pi events inside fiducial volume of interest

      //--------------------------------------------------//

      // Classify events based on interaction


      int genie_mode = -1;

      if (interaction == 0) { genie_mode = 1; } // QE
      else if (interaction == 10) { genie_mode = 2; } // MEC
      else if (interaction == 1) { genie_mode = 3; } // RES
      else if (interaction == 2) { genie_mode = 4; } // DIS
      else { genie_mode = 5; } // COH / other

      //--------------------------------------------------//

      // Fill in plots

      // All events 


      // For a specific interaction
      

      //--------------------------------------------------//

   } // End of the loop over the events

   //--------------------------------------------------//

   // Divide by bin width with Reweight function

   //--------------------------------------------------//

} // End of the program
