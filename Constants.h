#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

namespace Constants {

        //--------------------------------------------------//

        // Run 1 constants

        double data_pot = 1.62E20;
	double mc_pot = 5.1303E+19;  // 1.31E+21;
	double integrated_flux = 1.19E11; // cm^-2                                                                                                                                                      
	double number_targets = 1.05E30; // Argon nuclei, not nucleons    

	static const double ChargedPionMomentumThres = 0.07; // GeV/c

	static const int NNeut = 4;
	bool XSecMode = true;
        //--------------------------------------------------//

        // Pdg codes & masses

        static const int NuMuPdg = 14, MuonPdg = 13, ProtonPdg = 2212, AbsChargedPionPdg = 211, NeutralPionPdg = 111;
        static const int ElectronPdg = 11, PhotonPdg = 22, NeutronPdg = 2112, KaonPdg = 321;
        static const int DeuteriumPdg = 1000010020, HeliumPdg = 1000020040, ArgonPdg = 1000180400;

        static const double MuonMass = 106, ProtonMass = 938.272, NeutronMass = 939.565; // MeV
        static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
        static const double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.); // GeV^2

        //--------------------------------------------------//

	// Binning

	static const int NBinsDeltaPT = 13; static const double ArrayNBinsDeltaPT[NBinsDeltaPT+1] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9};
	static const int NBinsDeltaPtx = 11; static const double ArrayNBinsDeltaPtx[NBinsDeltaPtx+1] = {-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55};
	static const int NBinsDeltaPty = 12; static const double ArrayNBinsDeltaPty[NBinsDeltaPty+1] = {-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45};
	static const int NBinsDeltaAlphaT = 7; static const double ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT+1] = { 0.,22.,44.,66.,88.,110.,145.,180. };
	static const int NBinsDeltaPhiT = 12; static const double ArrayNBinsDeltaPhiT[NBinsDeltaPhiT+1] = {0.,12.5,25.,37.5,50.,60.,75.,90.,106.,126.,145.,162.,180.};
	static const int NBinsECal = 9; static const double ArrayNBinsECal[NBinsECal+1] = {0.2,0.35,0.5,0.65,0.8,0.95,1.1,1.25,1.4,1.6}; 
	static const int NBinsQ2 = 8; static const double ArrayNBinsQ2[NBinsQ2+1] = { 0.,0.08,0.18,0.28,0.39,0.5,0.65,0.8,1.}; 	
	static const int NBinsMuonMomentum = 10; static const double ArrayNBinsMuonMomentum[NBinsMuonMomentum+1] = { 0.1,0.2,0.3,0.4,0.5,0.64,0.77,0.9,1.,1.1,1.2};
	static const int NBinsMuonPhi = 10; static const double ArrayNBinsMuonPhi[NBinsMuonPhi+1] = { -180.,-144.,-108.,-72.,-36.,0.,36.,72.,108.,144.,180.};
	static const int NBinsMuonCosTheta = 18;
	static const double ArrayNBinsMuonCosTheta[NBinsMuonCosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};
	static const int NBinsProtonMomentum = 10; static const double ArrayNBinsProtonMomentum[NBinsProtonMomentum+1] = {0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.};
	static const int NBinsProtonPhi = 10; static const double ArrayNBinsProtonPhi[NBinsProtonPhi+1] = { -180.,-144.,-108.,-72.,-36.,0.,36.,72.,108.,144.,180.}; 
	static const int NBinsProtonCosTheta = 11;
	static const double ArrayNBinsProtonCosTheta[NBinsProtonCosTheta+1] = { -1.,-0.73,-0.43,-0.18,0.05,0.2,0.37,0.54,0.7,0.8,0.9,1.};

	static const int NBinsVertexX = 15; static const double MinVertexX = 10., MaxVertexX = 246.;
	static const int NBinsVertexY = 15; static const double MinVertexY = -105., MaxVertexY = 105.;
	static const int NBinsVertexZ = 25; static const double MinVertexZ = 10., MaxVertexZ = 1026.;

	TString RecoLabelXAxisVertexX = ";Vertex x [cm] ;Weighted Events / 5.1303E+19 POT";
	TString RecoLabelXAxisVertexY = ";Vertex y [cm] ;Weighted Events / 5.1303E+19 POT";
	TString RecoLabelXAxisVertexZ = ";Vertex z [cm] ;Weighted Events / 5.1303E+19 POT";
	
        //--------------------------------------------------//

	// Labels for 1D plots
	
	static TString LabelXAxisECal = ";E^{Cal} [GeV] ;#frac{d#sigma}{dE^{cal}}  [cm^{2} GeV^{-1} Ar^{-1}] "; static TString LabelXAxisTrueECal = ";True E^{Cal} [GeV]";
	static TString LabelXAxisMuonMomentum = ";p_{#mu} [GeV/c]"; static TString LabelXAxisTrueMuonMomentum = ";True p_{#mu} [GeV/c]";
	static TString LabelXAxisProtonMomentum = ";p_{p} [GeV/c]"; static TString LabelXAxisTrueProtonMomentum = ";True p_{p} [GeV/c]";
	static TString LabelXAxisMuonPhi = ";#phi_{#mu} [deg]"; static TString LabelXAxisTrueMuonPhi = ";True #phi_{#mu} [deg]";
	static TString LabelXAxisProtonPhi = ";#phi_{p} [deg]"; static TString LabelXAxisTrueProtonPhi = ";True #phi_{p} [deg]";
	static TString LabelXAxisMuonCosTheta = ";cos#theta_{#mu};#frac{d#sigma}{dcos#theta_{#mu}}  [cm^{2} Ar^{-1}]"; static TString LabelXAxisTrueMuonCosTheta = ";True cos#theta_{#mu}";
	static TString LabelXAxisProtonCosTheta = ";cos#theta_{p}"; static TString LabelXAxisTrueProtonCosTheta = ";True cos#theta_{p}";
	static TString LabelXAxisDeltaPtx = ";#deltap_{T,x} [GeV/c]"; static TString LabelXAxisTrueDeltaPtx = ";True #deltap_{T,x} [GeV/c]";
	static TString LabelXAxisDeltaPty = ";#deltap_{T,y} [GeV/c]"; static TString LabelXAxisTrueDeltaPty = ";True #deltap_{T,y} [GeV/c]";	
	static TString LabelXAxisDeltaPT = ";#deltap_{T} [GeV/c];#frac{d#sigma}{d#deltap_{T}}  [cm^{2} GeV^{-1}c Ar^{-1}]"; static TString LabelXAxisTrueDeltaPT = ";True #deltap_{T} [GeV/c]";
	static TString LabelXAxisDeltaAlphaT = ";#delta#alpha_{T} [deg];#frac{d#sigma}{d#delta#alpha_{T}}  [cm^{2} deg^{-1} Ar^{-1}]"; static TString LabelXAxisTrueDeltaAlphaT = ";True #delta#alpha_{T} [deg]";
	static TString LabelXAxisDeltaPhiT = ";#delta#phi_{T} [deg]"; static TString LabelXAxisTrueDeltaPhiT = ";True #delta#phi_{T} [deg]";

        //--------------------------------------------------//

	// Constants, Cuts & Thresholds

	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};
	const int NInte = InteractionLabels.size();

        //--------------------------------------------------//
	
}
#endif
