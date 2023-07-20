#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TVector3.h>

#include "../Constants.h"

using namespace std;
using namespace Constants;

void Reweight(TH1D* h); //reweights the histogram to produce the absolute cross section
bool inFV(TVector3 vector); //Checks if the vertex is within the fiducial volume
void xSec(TH1D* h);
void xSec(TH2D* h);


void NeutrinoSelectionFilter::Loop() {



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TH1D::SetDefaultSumw2();
   TH2D::SetDefaultSumw2();


   TFile* OutputFile;
   TString FileName;
   // Output file

   if (XSecMode != true){
     FileName = "MicroBooNE_Truth.root";
     OutputFile = new TFile(FileName,"recreate");
     std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;
   }

   else {
     FileName = "XSection_MicroBooNE_Truth.root";
     OutputFile = new TFile(FileName,"recreate");
     std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;
   }


   // Plot declaration
   
   TH1D* TrueVertexXPlot[NInte];
   TH1D* TrueVertexYPlot[NInte];
   TH1D* TrueVertexZPlot[NInte];   
   TH1D* TrueMuonCosThetaPlot[NInte];
   TH1D* TrueDeltaPtPlot[NInte];
   TH1D* TrueNeutronMultiplicity[NInte];
   TH1D* TrueDeltaAlphaTPlot[NInte];
   TH1D* TruePMissingMagnitudePlot[NInte][NNeut+1];
   TH1D* TruePMissingDirectionPlot[NInte][NNeut+1];

   TH2D* TruePMissingMagvsDirPlot[NInte][NNeut+1];
   TH2D* TrueNuEvsRecoNuEPlot[NInte][NNeut+1];



   // Loop over the interaction processes
   for (int inte = 0; inte < NInte; inte++) {
     // 1D analysis    
     TrueVertexXPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX-20,MaxVertexX+20);
     TrueVertexYPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY-20,MaxVertexY+20);
     TrueVertexZPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ-80,MaxVertexZ+80);
     TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
     TrueDeltaPtPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
     TrueNeutronMultiplicity[inte] = new TH1D(InteractionLabels[inte]+"TrueNeutronMultiplicityPlot",";Neutron Multiplicity;#frac{d#sigma}{dN}  [cm^{2} N^{-1} Ar^{-1}]",7,-0.5, 6.5);
     TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
    
     TruePMissingDirectionPlot[inte][NNeut] = new TH1D(InteractionLabels[inte]+"TruePMissingDirectionPlot_AllNeutron",";cos(#theta_{miss});#frac{d#sigma}{dcos(#theta_{miss})}  [cm^{2} Ar^{-1}]",10,-1,1);
     TruePMissingMagnitudePlot[inte][NNeut] = new TH1D(InteractionLabels[inte]+"TruePMissingMagnitudePlot_AllNeutron", ";P_{miss}  [GeV/c]; #frac{d#sigma}{dP_{miss}}  [cm^{2} GeV^{-1}c Ar^{-1}]", 10, 0,1);

     


     //2D analysis
      TruePMissingMagvsDirPlot[inte][NNeut] = new TH2D(InteractionLabels[inte]+"TruePMissingMagvsDirPlot_AllNeutron",";cos(#theta_{miss}) ;Magnitude of p_{miss}  [GeV/c]" , 20,-1,1,20,0,1.5);
      TrueNuEvsRecoNuEPlot[inte][NNeut] = new TH2D(InteractionLabels[inte]+"TrueNuEvsRecoNuEPlot_AllNeutron",";Reco #nu_{#mu} Energy [GeV] ;True #nu_{#mu} Energy  [GeV]" , 20,0,2,20,0,2);

     //Splitting into the number of neutrons
     for(int neut=0; neut <NNeut; neut++){
       TruePMissingDirectionPlot[inte][neut] = new TH1D(InteractionLabels[inte]+Form("TruePMissingDirectionPlot_Neutron%d",neut),";cos(#theta_{miss});#frac{d#sigma}{dcos(#theta_{miss})}  [cm^{2} Ar^{-1}]",10,-1,1);
       TruePMissingMagnitudePlot[inte][neut] = new TH1D(InteractionLabels[inte]+Form("TruePMissingMagnitudePlot_Neutron%d",neut),";P_{miss}  [GeV/c];#frac{d#sigma}{dP_{miss}}  [cm^{2} GeV^{-1}c Ar^{-1}]", 10, 0,1);
       TruePMissingMagvsDirPlot[inte][neut] = new TH2D(InteractionLabels[inte]+Form("TruePMissingMagvsDirPlot_Neutron%d",neut),";cos(#theta_{miss}) ;Magnitude of p_{miss}  [GeV/c]" , 20,-1,1,20,0,1.5);

       TrueNuEvsRecoNuEPlot[inte][neut] = new TH2D(InteractionLabels[inte]+Form("TrueNuEvsRecoNuEPlot_Neutron%d",neut),";Reco #nu_{#mu} Energy [GeV] ;True #nu_{#mu} Energy  [GeV]" , 20,0,2,20,0,2);
     }


     
   } // End of the loop over the interaction processes



   // Loop over the events
   int numberBeforeCut =0; int numberAfterCut =0;
   for (Long64_t jentry=0; jentry< nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;



      // MC weight to scale events to data pot
      //if (fabs(weightSplineTimesTune) != weightSplineTimesTune){ continue;}
      if (weightSplineTimesTune <= 0 || weightSplineTimesTune > 30) { 
	cout << weightSplineTimesTune << endl;
	continue;
      } // bug fix weight
      
      double event_weight = (data_pot / mc_pot_highStat) * weightSplineTimesTune;
      
      //std::cout << weightSpline << "   " << weightTune << "   " << weightSplineTimesTune << std::endl; //the weights were not set


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

	  if ( MCParticlePdg == MuonPdg && (MCParticleMomentum >= ArrayNBinsMuonMomentum[0] && MCParticleMomentum <= ArrayNBinsMuonMomentum[NBinsMuonMomentum] ) ) 
	    { TrueMuonCounter++;  VectorTrueMuonIndex.push_back(WhichMCParticle); }

	  if ( MCParticlePdg == ProtonPdg && ( MCParticleMomentum >= ArrayNBinsProtonMomentum[0] && MCParticleMomentum <= ArrayNBinsProtonMomentum[NBinsProtonMomentum] ) ) 
	    { TrueProtonCounter++; VectorTrueProtonIndex.push_back(WhichMCParticle); }

	  if ( fabs(MCParticlePdg) == AbsChargedPionPdg && MCParticleMomentum >= ChargedPionMomentumThres ) 
	    { TrueChargedPionCounter++; }

	  if (MCParticlePdg == NeutronPdg)
	    { TrueNeutronCounter++; VectorTrueNeutronIndex.push_back(WhichMCParticle);  }
	  
	  if (MCParticlePdg == NeutralPionPdg) { TruePi0Counter++; }
	} // End of the demand stable final state particles and primary interactions
      } // end of the loop over the MCParticles







      // CC1p0pi signal events

      bool CC1p0pi = false;
      if (TrueMuonCounter == 1 && TrueProtonCounter == 1 && TrueChargedPionCounter == 0 && TruePi0Counter == 0 && TrueHeavierMesonCounter == 0){ CC1p0pi = true; }
      if (!CC1p0pi) { continue; }

     





      // Define vectors and quantities which will be used throughout the plots
      TVector3 muonMomentumVector(mc_px->at(VectorTrueMuonIndex[0]), mc_py->at(VectorTrueMuonIndex[0]), mc_pz->at(VectorTrueMuonIndex[0]));
      TVector3 muonTransverseMom(muonMomentumVector.X(), muonMomentumVector.Y(), 0 );
      TVector3 protonMomentumVector(mc_px->at(VectorTrueProtonIndex[0]), mc_py->at(VectorTrueProtonIndex[0]), mc_pz->at(VectorTrueProtonIndex[0]));
      TVector3 deltaPtVector(muonMomentumVector.X() + protonMomentumVector.X(), muonMomentumVector.Y() + protonMomentumVector.Y() , 0);
      double transverseMissingMom = deltaPtVector.Mag();


      TVector3 VertexLocation(mc_vx->at(VectorTrueMuonIndex[0]), mc_vy->at(VectorTrueMuonIndex[0]), mc_vz->at(VectorTrueMuonIndex[0]) );
      TVector3 VertexLocationProton(mc_vx->at(VectorTrueProtonIndex[0]), mc_vy->at(VectorTrueProtonIndex[0]), mc_vz->at(VectorTrueProtonIndex[0]) );
      //std::cout << (VertexLocation-VertexLocationProton).X() << "   " << (VertexLocation-VertexLocationProton).Y()<< "   " << (VertexLocation-VertexLocationProton).Z() << std::endl;

      double MuonCosTheta = muonMomentumVector.CosTheta();
      double ProtonCosTheta = protonMomentumVector.CosTheta();
      
      double protonKE = mc_E->at(VectorTrueProtonIndex[0]) - ProtonMass_GeV;
      double CalEne = mc_E->at(VectorTrueMuonIndex[0]) + protonKE +0.04;
      TVector3 nuMomentumVector(0,0, CalEne);
      TVector3 pMissing = nuMomentumVector - protonMomentumVector - muonMomentumVector;
      double pMissingMagnitude = pMissing.Mag();
      double pMissingDirection = pMissing.CosTheta();

      double cosAlpha = (-1.0*(muonTransverseMom.Dot(deltaPtVector)) ) / (muonTransverseMom.Mag() * deltaPtVector.Mag());
      double deltaAlphaT = TMath::ACos(cosAlpha)* (180.0 / 3.1415926);
      
     
      double trueNuEne = nu_e;
      




      numberBeforeCut++;
      // Use only true CC1p0pi events inside fiducial volume of interest
      if (!inFV(VertexLocation)){ continue; }
      //if ((pMissingMagnitude > 0.3) && (pMissingDirection >-0.2)){ continue; }
      numberAfterCut++;
   


      // Classify events based on interaction
      int genie_mode = -1;

      if (interaction == 0) { genie_mode = 1; } // QE
      else if (interaction == 10) { genie_mode = 2; } // MEC
      else if (interaction == 1) { genie_mode = 3; } // RES
      else if (interaction == 2) { genie_mode = 4; } // DIS
      else { genie_mode = 5; } // COH / other









      // Fill in plots
      // All events 
      TrueVertexXPlot[0]->Fill(VertexLocation.X(), event_weight);
      TrueVertexYPlot[0]->Fill(VertexLocation.Y(), event_weight);
      TrueVertexZPlot[0]->Fill(VertexLocation.Z(), event_weight);
      TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta, event_weight);
      
      TrueDeltaPtPlot[0]->Fill(transverseMissingMom, event_weight);
      TrueNeutronMultiplicity[0]->Fill(TrueNeutronCounter , event_weight);
      TrueDeltaAlphaTPlot[0]->Fill(deltaAlphaT , event_weight);
      
      
      // // For a specific interaction
      
      TrueVertexXPlot[genie_mode]->Fill(VertexLocation.X(), event_weight);	
      TrueVertexYPlot[genie_mode]->Fill(VertexLocation.Y(), event_weight);	
      TrueVertexZPlot[genie_mode]->Fill(VertexLocation.Z(), event_weight);	
      TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta, event_weight);


      TrueDeltaPtPlot[genie_mode]->Fill(transverseMissingMom, event_weight);
      TrueNeutronMultiplicity[genie_mode]->Fill(TrueNeutronCounter , event_weight);
      TrueDeltaAlphaTPlot[genie_mode]->Fill(deltaAlphaT , event_weight);


      TruePMissingMagnitudePlot[0][NNeut]->Fill(pMissingMagnitude , event_weight);
      TruePMissingDirectionPlot[0][NNeut]->Fill(pMissingDirection , event_weight);
      TruePMissingMagvsDirPlot[0][NNeut]->Fill(pMissingDirection, pMissingMagnitude, event_weight);      

      TruePMissingMagnitudePlot[genie_mode][NNeut]->Fill(pMissingMagnitude , event_weight);
      TruePMissingDirectionPlot[genie_mode][NNeut]->Fill(pMissingDirection , event_weight);
      TruePMissingMagvsDirPlot[genie_mode][NNeut]->Fill(pMissingDirection, pMissingMagnitude, event_weight);      
     

      TrueNuEvsRecoNuEPlot[0][NNeut]->Fill(CalEne, trueNuEne, event_weight);
      TrueNuEvsRecoNuEPlot[genie_mode][NNeut]->Fill(CalEne, trueNuEne, event_weight);
      
      
      //for number of neutrons
      if (TrueNeutronCounter == 0){
      	TruePMissingMagnitudePlot[0][0]->Fill(pMissingMagnitude , event_weight);
      	TruePMissingDirectionPlot[0][0]->Fill(pMissingDirection , event_weight);
      	TruePMissingMagnitudePlot[genie_mode][0]->Fill(pMissingMagnitude , event_weight);
      	TruePMissingDirectionPlot[genie_mode][0]->Fill(pMissingDirection , event_weight);
      	TruePMissingMagvsDirPlot[genie_mode][0]->Fill(pMissingDirection, pMissingMagnitude, event_weight);      
	TrueNuEvsRecoNuEPlot[0][0]->Fill(CalEne, trueNuEne, event_weight);
	TrueNuEvsRecoNuEPlot[genie_mode][0]->Fill(CalEne, trueNuEne, event_weight);

      }

      if (TrueNeutronCounter ==1){
      	TruePMissingMagnitudePlot[0][1]->Fill(pMissingMagnitude , event_weight);
      	TruePMissingDirectionPlot[0][1]->Fill(pMissingDirection , event_weight);
      	TruePMissingMagnitudePlot[genie_mode][1]->Fill(pMissingMagnitude , event_weight);
      	TruePMissingDirectionPlot[genie_mode][1]->Fill(pMissingDirection , event_weight);
      	TruePMissingMagvsDirPlot[genie_mode][1]->Fill(pMissingDirection, pMissingMagnitude, event_weight);      
	TrueNuEvsRecoNuEPlot[0][1]->Fill(CalEne, trueNuEne, event_weight);
	TrueNuEvsRecoNuEPlot[genie_mode][1]->Fill(CalEne, trueNuEne, event_weight);

      }


      if (TrueNeutronCounter >=2){
      	TruePMissingMagnitudePlot[0][2]->Fill(pMissingMagnitude , event_weight);
      	TruePMissingDirectionPlot[0][2]->Fill(pMissingDirection , event_weight);
      	TruePMissingMagnitudePlot[genie_mode][2]->Fill(pMissingMagnitude , event_weight);
      	TruePMissingDirectionPlot[genie_mode][2]->Fill(pMissingDirection , event_weight);
      	TruePMissingMagvsDirPlot[genie_mode][2]->Fill(pMissingDirection, pMissingMagnitude, event_weight);      
	TrueNuEvsRecoNuEPlot[0][2]->Fill(CalEne, trueNuEne, event_weight);
	TrueNuEvsRecoNuEPlot[genie_mode][2]->Fill(CalEne, trueNuEne, event_weight);

      }





   } // End of the loop over the events
   

   std::cout << "Before: " << numberBeforeCut << " ||  After: " << numberAfterCut << "  || Ratio: " << ((double)numberAfterCut) / ((double)numberBeforeCut) << std::endl;

 



   //Save the files without turning to cross sections
   if (XSecMode != true){
     OutputFile->cd();
     OutputFile->Write();
     OutputFile->Close();
     std::cout << std::endl;
     std::cout << "File " << FileName << " has been created" << std::endl;
     std::cout << std::endl;
     std::cout << std::endl << "--------------------------------------------------------------" << std::endl <<std::endl;
   }

   else {
     // Convert Histograms to Cross sections
     for (int inte = 0; inte < NInte; inte++){
       xSec(TrueDeltaPtPlot[inte]);
       // xSec(TrueVertexXPlot[inte]);
       // xSec(TrueVertexYPlot[inte]);
       // xSec(TrueVertexZPlot[inte]);
       xSec(TrueMuonCosThetaPlot[inte]);
       xSec(TrueNeutronMultiplicity[inte]);
       xSec(TrueDeltaAlphaTPlot[inte]);
       
       for (int neut; neut <=NNeut; neut++){
	 xSec(TruePMissingMagnitudePlot[inte][neut]);
	 xSec(TruePMissingDirectionPlot[inte][neut]);
	 xSec(TruePMissingMagvsDirPlot[inte][neut]);
       }       
     }
     
     //Save the files with turning to cross sections
     OutputFile->cd();
     OutputFile->Write();
     OutputFile->Close();
     std::cout << std::endl;
     std::cout << "File " << FileName  << " has been created" << std::endl;
     std::cout << std::endl;
     std::cout << std::endl << "--------------------------------------------------------------" << std::endl <<std::endl;
     
   }


} // End of the program



void Reweight(TH1D* h){
  int NBins = h->GetXaxis()->GetNbins();
 
  for(int i=0; i<NBins; i++){
    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = CurrentEntry/ h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = CurrentError/ h->GetBinWidth(i+1);
    
    h->SetBinContent(i+1, NewEntry);
    h->SetBinError(i+1,NewError);
  }
}


bool inFV(TVector3 vector){
  if((vector.X() < MaxVertexX) && (vector.X() > MinVertexX) && (vector.Y() < MaxVertexY) && (vector.Y() > MinVertexY) && (vector.Z() < MaxVertexZ) && (vector.Z() > MinVertexZ) ) return true;
  
  else return false;
}


void xSec(TH1D* h){
  int NBins = h->GetXaxis()->GetNbins();
  double toXsec = integrated_flux * number_targets;
  for(int i=0; i<NBins; i++){
    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = CurrentEntry/ (h->GetBinWidth(i+1) * toXsec);

    double CurrentError = h->GetBinError(i+1);
    double NewError = CurrentError/ (h->GetBinWidth(i+1) * toXsec);
    
    h->SetBinContent(i+1, NewEntry);
    h->SetBinError(i+1,NewError);
  }  
}


void xSec(TH2D* h){
  int NXBins = h->GetXaxis()->GetNbins();
  int NYBins = h->GetYaxis()->GetNbins();

  double toXsec = integrated_flux * number_targets;
  for(int i=0; i<NXBins; i++){
    for(int j=0; j<NYBins; j++){
      double CurrentEntry = h->GetBinContent(i+1,j+1);
      double NewEntry = CurrentEntry/ ((h->GetXaxis()->GetBinWidth(i+1))*(h->GetYaxis()->GetBinWidth(j+1)) * toXsec);
      
      double CurrentError = h->GetBinError(i+1,j+1);
      double NewError = CurrentError/ ((h->GetXaxis()->GetBinWidth(i+1))*(h->GetYaxis()->GetBinWidth(j+1)) * toXsec);
      
      h->SetBinContent(i+1, j+1, NewEntry);
      h->SetBinError(i+1, j+1, NewError);
    }
  }  
}


