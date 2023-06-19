#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include "Constants.h"

using namespace std;
using namespace Constants;

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

void PlotRoot() {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			
	double LegendTextSize = 0.04;

	TString OutFilePath = "/uboone/app/users/maxd/MicroBooNE_NeutronStudy_Truth/";

	//------------------------------//

	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; std::vector<int> Colors;
	if (XSecMode == false){
	  Names.push_back(OutFilePath+"MicroBooNE_Truth.root"); 
	  Labels.push_back("All");
	  Colors.push_back(kBlack);	
	  
	  Names.push_back(OutFilePath+"MicroBooNE_Truth.root"); 
	  Labels.push_back("QE");
	  Colors.push_back(kCyan+2);	
	  
	  Names.push_back(OutFilePath+"MicroBooNE_Truth.root"); 
	  Labels.push_back("MEC");
	  Colors.push_back(kGreen+2);	
	  
	  Names.push_back(OutFilePath+"MicroBooNE_Truth.root"); 
	  Labels.push_back("RES");
	  Colors.push_back(kOrange+2);	
	  
	  Names.push_back(OutFilePath+"MicroBooNE_Truth.root"); 
	  Labels.push_back("DIS");
	  Colors.push_back(kBlue+2);	
	}

	else {
	  Names.push_back(OutFilePath+"XSection_MicroBooNE_Truth.root"); 
	  Labels.push_back("All");
	  Colors.push_back(kBlack);	
	  
	  Names.push_back(OutFilePath+"XSection_MicroBooNE_Truth.root"); 
	  Labels.push_back("QE");
	  Colors.push_back(kCyan+2);	
	  
	  Names.push_back(OutFilePath+"XSection_MicroBooNE_Truth.root"); 
	  Labels.push_back("MEC");
	  Colors.push_back(kGreen+2);	
	  
	  Names.push_back(OutFilePath+"XSection_MicroBooNE_Truth.root"); 
	  Labels.push_back("RES");
	  Colors.push_back(kOrange+2);	
	  
	  Names.push_back(OutFilePath+"XSection_MicroBooNE_Truth.root"); 
	  Labels.push_back("DIS");
	  Colors.push_back(kBlue+2);

	}

	const int NSamples = Names.size();
	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;
	///*
	//Choose which plots to display in the files

	 // PlotNames.push_back("TrueVertexXPlot");
	 // PlotNames.push_back("TrueVertexYPlot");
	 // PlotNames.push_back("TrueVertexZPlot");

	//PlotNames.push_back("TrueMuonCosThetaPlot");
	//PlotNames.push_back("TrueDeltaPtPlot");
	
	// PlotNames.push_back("TrueNeutronMultiplicityPlot");
	 //PlotNames.push_back("TrueDeltaAlphaTPlot");
	
	// for (int neut=0; neut < 2; neut++){
	//   PlotNames.push_back(Form("TruePMissingDirectionPlot_Neutron%d",neut));
	//   PlotNames.push_back(Form("TruePMissingMagnitudePlot_Neutron%d",neut));
	// }

	

	const int NPlots = PlotNames.size();







	// Loop over the samples to open the files and the TTree
	for (int iSample = 0; iSample < NSamples; iSample++) {
		Files[iSample] = new TFile(Names[iSample],"readonly");
	} // End of the loop over the samples



	

	// Loop over the plots to be compared
	for (int iPlot = 0; iPlot < NPlots; iPlot++) {

		TString CanvasName = "Canvas_" + PlotNames[iPlot];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetTopMargin(0.12);
		PlotCanvas->SetLeftMargin(0.2);
		PlotCanvas->SetBottomMargin(0.15);		
		PlotCanvas->Draw();	

		TLegend* leg  = new TLegend(0.35,0.91, 0.85, 0.99 ); //On the right side of plot and above border
		//TLegend* leg  = new TLegend(0.2,0.7,0.55,0.83); //Original from Afro
		leg->SetBorderSize(0);
		leg->SetNColumns(5);
		leg->SetTextSize(LegendTextSize); 
		leg->SetTextFont(FontStyle);						

		// Loop over the samples to open the files and to get the corresponding plot

		std::vector<TH1D*> Histos; Histos.resize(NSamples);

		for (int iSample = 0; iSample < NSamples; iSample++) {	
		  
		 
		  Histos[iSample] = (TH1D*)(Files[iSample]->Get(InteractionLabels[iSample] +PlotNames[iPlot]));
		  
		  Histos[iSample]->SetLineWidth(4);
		  Histos[iSample]->SetLineColor( Colors.at(iSample) );	
		  
		  Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
		  Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
		  Histos[iSample]->GetXaxis()->SetNdivisions(8);
		  Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
		  Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
		  Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
		  Histos[iSample]->GetXaxis()->CenterTitle();						

		  Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
		  Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
		  Histos[iSample]->GetYaxis()->SetNdivisions(6);
		  Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
		  //Histos[iSample]->GetYaxis()->SetTitle("Weighted Events / 1.62E20 POT ");
		  Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
		  Histos[iSample]->GetYaxis()->SetTitleOffset(1.2);
		  Histos[iSample]->GetYaxis()->SetTickSize(0);
		  Histos[iSample]->GetYaxis()->CenterTitle();	

		  double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());			
		  Histos[iSample]->GetYaxis()->SetRangeUser(0.,1.1*imax);
		  Histos[0]->GetYaxis()->SetRangeUser(0.,1.1*imax);			

		  PlotCanvas->cd();
		  Histos[iSample]->Draw("hist same"); // e");
		  Histos[0]->Draw("hist same"); // e");	

		  leg->AddEntry(Histos[iSample],Labels[iSample],"l");
			
		
		} // End of the loop over the samples grabing the plots	
		
		PlotCanvas->cd();
		leg->Draw();
		
		if (XSecMode==false){ PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"NoCut_MicroBooNE_Truth.pdf"); }
		else {PlotCanvas->SaveAs("myPlotsXsec/"+PlotNames[iPlot]+"XSection_NoCut_MicroBooNE_Truth.pdf"); }
		
	} // End of the loop over the plots
	
	
} // End of the program
