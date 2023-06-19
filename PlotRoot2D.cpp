#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TPad.h>
#include "Constants.h"
using namespace std;
using namespace Constants;

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

void PlotRoot2D() {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			
	double LegendTextSize = 0.03;

	TString OutFilePath = "/uboone/app/users/maxd/MicroBooNE_NeutronStudy_Truth/";


	
	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; std::vector<int> Colors;
	if (XSecMode == false){
	  Names.push_back(OutFilePath+"MicroBooNE_Truth.root"); 
	  Labels.push_back(" ");
	  Colors.push_back(kBlack);	
	}

	else {
	  Names.push_back(OutFilePath+"XSection_MicroBooNE_Truth.root"); 
	  Labels.push_back(" ");
	  Colors.push_back(kBlack);
	}

	const int NSamples = Names.size();
	std::vector<TFile*> Files; Files.resize(NSamples);



	// Plots to overlay

	std::vector<TString> PlotNames;
	for (int neut =0; neut <3; neut++){
	  // PlotNames.push_back(Form("TruePMissingMagvsDirPlot_Neutron%d",neut));
	  // PlotNames.push_back(Form("QETruePMissingMagvsDirPlot_Neutron%d",neut));
	  // PlotNames.push_back(Form("MECTruePMissingMagvsDirPlot_Neutron%d",neut));
	  // PlotNames.push_back(Form("RESTruePMissingMagvsDirPlot_Neutron%d",neut));
	  // PlotNames.push_back(Form("DISTruePMissingMagvsDirPlot_Neutron%d",neut));
	}

	PlotNames.push_back("TrueNuEvsRecoNuEPlot_AllNeutron");



	const int NPlots = PlotNames.size();




	// Loop over the samples to open the files and the TTree

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples





	// Loop over the plots to be compared


	
	for (int iPlot = 0; iPlot < NPlots; iPlot++) {
	  // Loop over the samples to open the files and to get the corresponding plot
	  for (int iSample = 0; iSample < NSamples; iSample++) {	
	
	    std::vector<TH2D*> Histos; Histos.resize(NSamples);
	    TString CanvasName = "Canvas_" + PlotNames[iPlot];
	    TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	    PlotCanvas->cd();
	    PlotCanvas->SetTopMargin(0.12);
	    PlotCanvas->SetLeftMargin(0.18);
	    PlotCanvas->SetBottomMargin(0.15);		
	    PlotCanvas->Draw();	
	    
	    
	    Histos[iSample] = (TH2D*)(Files[iSample]->Get(PlotNames[iPlot]));
	    
	    // Histos[iSample]->SetLineWidth(4);
	    // Histos[iSample]->SetLineColor( Colors.at(iSample) );	
	    Histos[iSample]->SetTitle(Labels[iSample]);
	    Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
	    //Histos[iSample]->GetXaxis()->SetTitle("cos(#theta_{miss})");
	    Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
	    
	    // Histos[iSample]->GetXaxis()->SetNdivisions(8);
	    Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
	    Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
	    Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
	    Histos[iSample]->GetXaxis()->CenterTitle();						
	    
	    Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
	    //Histos[iSample]->GetYaxis()->SetTitle("Magnitude p_{missing}");
	    Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
	    //Histos[iSample]->GetYaxis()->SetNdivisions(6);
	    
	    Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
	    Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
	    Histos[iSample]->GetYaxis()->SetTitleOffset(1.3);
	    Histos[iSample]->GetYaxis()->SetTickSize(0);
	    Histos[iSample]->GetYaxis()->CenterTitle();	
	    
	    PlotCanvas->cd()->SetLogz();
	    
	    
	    PlotCanvas->cd();
	    Histos[iSample]->Draw("colz");
	    Histos[0]->Draw("colz");	
	    
	    // leg->AddEntry(Histos[iSample],Labels[iSample],"l");
	    
	    
	    PlotCanvas->cd();
	    if (XSecMode==false){ PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"NoCut_MicroBooNE_Truth.pdf"); }
	    else {PlotCanvas->SaveAs("myPlotsXsec/"+PlotNames[iPlot]+"XSectionNoCut_MicroBooNE_Truth.pdf"); }
	 
	  } // End of the loop over the samples grabing the plots	
	} // End of the loop over the plots
	
	
	
	

	

	
	

	

	
	
	
} // End of the program
