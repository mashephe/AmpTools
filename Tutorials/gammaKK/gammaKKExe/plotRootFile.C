
/************************************************************
 *
 * 2012/08/04 Kei Moriya (Indiana University)
 *
 * This program takes in a ROOT file created by
 * plotResults.cc and plots the histograms.
 *
 * The histograms are for J/psi -> gamma + K + K.
 * For data, acc MC, gen MC, we have
 * - m12, m23, m13
 * - photon angles
 * - Kaon angles
 *
 ************************************************************/


#include <iostream>

int plotRootFile(string filename=""){

  gROOT->ForceStyle();
  gStyle->SetTitleSize(0.050,"XYZ"); 
  gStyle->SetTitleOffset(0.850,"XYZ");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetOptStat(0);

  if(filename==""){
    cout << endl;
    cout << "----------------------------------------------" << endl;
    cout << endl;
    cout << "plotRootFile:" << endl;
    cout << "Usage: plotRootFile (\"input ROOT file name\")" << endl;
    cout << endl;
    cout << "----------------------------------------------" << endl;
    cout << endl;
    abort();
  }

  // Read in input ROOT file
  TFile *infile = new TFile(filename.c_str());

  // Plot histograms for
  // data, acc MC, gen MC.
  // Histograms are
  // - m12, m23, m13
  // - photon angles
  // - Kaon angles
  // for each amplitude, and total.

  const Int_t NUMTYPES = 3;
  const Int_t NUMAMPS  = 2;

  const Int_t mycolors[] = {kRed, kBlue, kGreen, kMagenta, kCyan};

  char datatype[20];
  char ampnum[20];

  char canvasname[200];
  char canvastitle[200];

  char hname[200];
  char hnewname[200];

  char text[200];

  // Get PhotonCosTheta.
  // All histogram names are determined in plotResults.cc
  TH1F *hdataPhotonCosTheta[NUMTYPES];
  TH1F *hMCPhotonCosTheta[NUMTYPES][NUMAMPS];
  TH1F *htotMCPhotonCosTheta[NUMTYPES];
  TH1F *hsumMCPhotonCosTheta[NUMTYPES];

  for(Int_t itype=0;itype<3;itype++){ // for data, acc MC, gen MC
    if(itype==0)      sprintf(datatype,"dat");
    else if(itype==1) sprintf(datatype,"acc");
    else if(itype==2) sprintf(datatype,"gen");
    
    // Data
    if(itype==0){
      sprintf(hname,"hPhotonCosTheta%s",datatype);
      sprintf(hnewname,"_hPhotonCosTheta%s",datatype);
      hdataPhotonCosTheta[itype] = (TH1F*)infile->Get(hname)->Clone(hnewname);
    }
    // MC
    else{

      sprintf(hname,"hPhotonCosTheta%s",datatype,iamp+1);
      sprintf(hnewname,"_hPhotonCosTheta%s",datatype);
      htotMCPhotonCosTheta[itype] = (TH1F*)infile->Get(hname)->Clone(hnewname);

      for(Int_t iamp=0;iamp<NUMAMPS;iamp++){ // for each amplitude
	sprintf(hname,"hPhotonCosTheta%s%d",datatype,iamp+1);
	sprintf(hnewname,"_hPhotonCosTheta%s%d",datatype,iamp+1);
	hMCPhotonCosTheta[itype][iamp] = (TH1F*)infile->Get(hname)->Clone(hnewname);
	
	// Create total intensity histogram
	if(iamp==0){
	  sprintf(hname,"hPhotonCosTheta%s",datatype);
	  sprintf(hnewname,"_hPhotonCosTheta%s",datatype);
	  hsumMCPhotonCosTheta[itype] = (TH1F*)hMCPhotonCosTheta[itype][iamp]->Clone(hnewname);
	}
	else{
	  // Add to total
	  hsumMCPhotonCosTheta[itype]->Add(hMCPhotonCosTheta[itype][iamp]);
	}
      } // end of loop over amps
    } // end of if MC
  } // end of data, acc MC, gen MC

  TCanvas *cPhotonCosTheta[NUMTYPES];
  TLegend *legendPhotonCosTheta[NUMTYPES];
  for(Int_t itype=0;itype<3;itype++){ // for data, acc MC, gen MC
    if(itype==0)      sprintf(datatype,"dat");
    else if(itype==1) sprintf(datatype,"acc");
    else if(itype==2) sprintf(datatype,"gen");

    // Create canvas
    sprintf(canvasname,"c_PhotonCosTheta_%s",datatype);
    sprintf(canvastitle,"PhotonCosTheta for %s",datatype);
    cPhotonCosTheta[itype] = new TCanvas(canvasname, canvastitle,0,0,1200,750);

    // Draw
    legendPhotonCosTheta[itype] = new TLegend(0.36,0.12,0.85,0.32);
    legendPhotonCosTheta[itype]->SetBorderSize(0);
    legendPhotonCosTheta[itype]->SetFillStyle(0);

    // Data
    if(itype==0){
      hdataPhotonCosTheta[itype]->SetLineColor(kBlack);
      hdataPhotonCosTheta[itype]->SetLineWidth(2);
      hdataPhotonCosTheta[itype]->SetLineStyle(1);
 
      hdataPhotonCosTheta[itype]->SetMinimum(0);
      hdataPhotonCosTheta[itype]->GetXaxis()->CenterTitle();
      hdataPhotonCosTheta[itype]->GetYaxis()->CenterTitle();

      hdataPhotonCosTheta[itype]->Draw();

      sprintf(text,"total for %s",datatype);
      legendPhotonCosTheta[itype]->AddEntry(hdataPhotonCosTheta[itype],text,"LP");

      // Draw acc MC on top
      hsumMCPhotonCosTheta[1]->SetLineColor(kViolet);
      hsumMCPhotonCosTheta[1]->SetLineWidth(2);
      hsumMCPhotonCosTheta[1]->SetLineStyle(2);
      hsumMCPhotonCosTheta[1]->Draw("same");

      sprintf(text,"sum for acc MC");
      legendPhotonCosTheta[itype]->AddEntry(hsumMCPhotonCosTheta[1],text,"LP");

      // Draw each amplitude contribution for acc MC
      for(Int_t iamp=0;iamp<NUMAMPS;iamp++){ // for each amplitude
	hMCPhotonCosTheta[1][iamp]->SetLineColor(mycolors[iamp]);
	hMCPhotonCosTheta[1][iamp]->SetLineWidth(1);
	hMCPhotonCosTheta[1][iamp]->SetLineStyle(2);
	hMCPhotonCosTheta[1][iamp]->SetMarkerColor(mycolors[iamp]);
	hMCPhotonCosTheta[1][iamp]->SetMarkerSize(0.50);
	hMCPhotonCosTheta[1][iamp]->SetMarkerStyle(20);
	hMCPhotonCosTheta[1][iamp]->Draw("same");
	
	sprintf(text,"for amp %d",iamp+1);
	legendPhotonCosTheta[itype]->AddEntry(hMCPhotonCosTheta[1][iamp],text,"LP");
      }

    }
    // MC
    else{
      htotMCPhotonCosTheta[itype]->SetLineColor(kBlack);
      htotMCPhotonCosTheta[itype]->SetLineWidth(2);
      htotMCPhotonCosTheta[itype]->SetLineStyle(1);
      htotMCPhotonCosTheta[itype]->SetMinimum(0);
      htotMCPhotonCosTheta[itype]->GetXaxis()->CenterTitle();
      htotMCPhotonCosTheta[itype]->GetYaxis()->CenterTitle();
      htotMCPhotonCosTheta[itype]->Draw();

      hsumMCPhotonCosTheta[itype]->SetLineColor(kRed);
      hsumMCPhotonCosTheta[itype]->SetLineWidth(2);
      hsumMCPhotonCosTheta[itype]->SetLineStyle(1);
      hsumMCPhotonCosTheta[itype]->SetMinimum(0);
      hsumMCPhotonCosTheta[itype]->Draw("same");

      sprintf(text,"total for %s MC",datatype);
      legendPhotonCosTheta[itype]->AddEntry(hsumMCPhotonCosTheta[itype],text,"LP");

      // If this is the generated MC, then do a fit to 1 + alpha cos^2 theta
      if(itype==2){
	TF1 *fitFunc;
	fitFunc = new TF1("fitFunc","[0] * (1. + [1] * pow(x,2.))",-1,1);
	fitFunc->SetParameter(0,htotMCPhotonCosTheta[itype]->GetMaximum());
	htotMCPhotonCosTheta[itype]->Fit(fitFunc,"EQRN");
	cout << fitFunc->GetParameter(1) << "\t+/-\t" << fitFunc->GetParError(1) << endl;

	// Draw fit function
	fitFunc->SetLineColor(kSpring);
	fitFunc->SetLineStyle(1);
	fitFunc->SetLineWidth(3);
	fitFunc->Draw("same");

	sprintf(text,"fit with 1+#alpha cos^{2}#theta_{#gamma},#alpha=%7.5f#pm%7.5f",
		fitFunc->GetParameter(1),fitFunc->GetParError(1));
	legendPhotonCosTheta[itype]->AddEntry(fitFunc,text,"LP");
      } // end of gen MC
      
      for(Int_t iamp=0;iamp<NUMAMPS;iamp++){ // for each amplitude
	hMCPhotonCosTheta[itype][iamp]->SetLineColor(mycolors[iamp]);
	hMCPhotonCosTheta[itype][iamp]->SetLineWidth(1);
	hMCPhotonCosTheta[itype][iamp]->SetLineStyle(2);
	hMCPhotonCosTheta[itype][iamp]->Draw("same");
	
	sprintf(text,"intensity for amp %d",iamp+1);
	legendPhotonCosTheta[itype]->AddEntry(hMCPhotonCosTheta[itype][iamp],text,"LP");

      }
    } // end of MC
    legendPhotonCosTheta[itype]->Draw("same");

    // Save
    sprintf(text,"figures/%s/PhotonCosTheta_%s.pdf",filename.c_str(),datatype);
    cPhotonCosTheta[itype]->SaveAs(text);
  }

}
