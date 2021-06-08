
TCanvas* c1;
void plotProjections(int canvas, TString projection);


void plotResult(){


  c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);


  _file1->cd();

  plotProjections(2,"12");
  plotProjections(3,"23");
  plotProjections(4,"13");


  _file0->cd();

  TH2F* hs12s23  = (TH2F*) gDirectory->FindObjectAny("hs12s23");

  c1->cd(1);
  hs12s23->Draw("colz");

  TLegend* tl = new TLegend(0.6,0.77,0.85,0.85);
  tl->AddEntry(hs12s23,"DATA","");
  tl->Draw("same");


  c1->cd();

}



void plotProjections(int canvas, TString projection){

  c1->cd(canvas);

  TH1F* hmdat  = (TH1F*) gDirectory->FindObjectAny("hm"+projection+"dat");
  TH1F* hmacc  = (TH1F*) gDirectory->FindObjectAny("hm"+projection+"acc");
  TH1F* hmacc1 = (TH1F*) gDirectory->FindObjectAny("hm"+projection+"acc1");
  TH1F* hmacc2 = (TH1F*) gDirectory->FindObjectAny("hm"+projection+"acc2");

  hmdat->SetTitleOffset(2.0,"Y");
  hmdat->SetTitleOffset(1.6,"X");
  hmdat->SetMarkerStyle(20);
  hmdat->SetMarkerSize(0.5);
  hmdat->SetLineColor(kBlack);

  hmacc->SetLineColor(kBlack);
  hmacc->SetLineWidth(2.0);

  hmacc1->SetLineColor(kRed);
  hmacc1->SetLineWidth(2.0);

  hmacc2->SetLineColor(kBlue);
  hmacc2->SetLineWidth(2.0);

  hmdat->Draw();
  hmacc->Draw("hist,same");
  hmacc1->Draw("hist,same");
  hmacc2->Draw("hist,same");
  hmdat->Draw("same");

  TLegend* tl = new TLegend(0.6,0.65,0.85,0.85);
  if (projection == "23") tl = new TLegend(0.2,0.65,0.45,0.85);
  tl->AddEntry(hm12dat,"DATA","p");
  tl->AddEntry(hm12acc,"FIT","l");
  tl->AddEntry(hm12acc1,"AMP1","l");
  tl->AddEntry(hm12acc2,"AMP2","l");
  tl->Draw("same");

}


