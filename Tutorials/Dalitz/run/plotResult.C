{

  _file0->cd();

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);

  c1->cd(1);

  TH1F* hm12dat  = (TH1F*) gDirectory->FindObjectAny("hm12dat");
  TH1F* hm12acc  = (TH1F*) gDirectory->FindObjectAny("hm12acc");
  TH1F* hm12acc1 = (TH1F*) gDirectory->FindObjectAny("hm12acc1");
  TH1F* hm12acc2 = (TH1F*) gDirectory->FindObjectAny("hm12acc2");

  hm12dat->SetTitleOffset(2.0,"Y");
  hm12dat->SetTitleOffset(1.6,"X");
  hm12dat->SetMarkerStyle(20);
  hm12dat->SetMarkerSize(0.5);
  hm12dat->SetLineColor(kBlack);

  hm12acc->SetLineColor(kBlack);
  hm12acc->SetLineWidth(2.0);

  hm12acc1->SetLineColor(kRed);
  hm12acc1->SetLineWidth(2.0);

  hm12acc2->SetLineColor(kBlue);
  hm12acc2->SetLineWidth(2.0);

  hm12dat->Draw();
  hm12acc->Draw("hist,same");
  hm12acc1->Draw("hist,same");
  hm12acc2->Draw("hist,same");
  hm12dat->Draw("same");

  TLegend* tl = new TLegend(0.6,0.6,0.85,0.85);
  tl->AddEntry(hm12dat,"DATA","p");
  tl->AddEntry(hm12acc,"FIT","l");
  tl->AddEntry(hm12acc1,"AMP1","l");
  tl->AddEntry(hm12acc2,"AMP2","l");
  tl->Draw("same");


  c1->cd(2);

  TH1F* hm13dat  = (TH1F*) gDirectory->FindObjectAny("hm13dat");
  TH1F* hm13acc  = (TH1F*) gDirectory->FindObjectAny("hm13acc");
  TH1F* hm13acc1 = (TH1F*) gDirectory->FindObjectAny("hm13acc1");
  TH1F* hm13acc2 = (TH1F*) gDirectory->FindObjectAny("hm13acc2");

  hm13dat->SetTitleOffset(2.0,"Y");
  hm13dat->SetTitleOffset(1.6,"X");
  hm13dat->SetMarkerStyle(20);
  hm13dat->SetMarkerSize(0.5);
  hm13dat->SetLineColor(kBlack);

  hm13acc->SetLineColor(kBlack);
  hm13acc->SetLineWidth(2.0);

  hm13acc1->SetLineColor(kRed);
  hm13acc1->SetLineWidth(2.0);

  hm13acc2->SetLineColor(kBlue);
  hm13acc2->SetLineWidth(2.0);

  hm13dat->Draw();
  hm13acc->Draw("hist,same");
  hm13acc1->Draw("hist,same");
  hm13acc2->Draw("hist,same");
  hm13dat->Draw("same");

  tl->Draw("same");



  c1->cd(3);

  TH1F* hm23dat  = (TH1F*) gDirectory->FindObjectAny("hm23dat");
  TH1F* hm23acc  = (TH1F*) gDirectory->FindObjectAny("hm23acc");
  TH1F* hm23acc1 = (TH1F*) gDirectory->FindObjectAny("hm23acc1");
  TH1F* hm23acc2 = (TH1F*) gDirectory->FindObjectAny("hm23acc2");

  hm23dat->SetTitleOffset(2.0,"Y");
  hm23dat->SetTitleOffset(1.6,"X");
  hm23dat->SetMarkerStyle(20);
  hm23dat->SetMarkerSize(0.5);
  hm23dat->SetLineColor(kBlack);

  hm23acc->SetLineColor(kBlack);
  hm23acc->SetLineWidth(2.0);

  hm23acc1->SetLineColor(kRed);
  hm23acc1->SetLineWidth(2.0);

  hm23acc2->SetLineColor(kBlue);
  hm23acc2->SetLineWidth(2.0);

  hm23dat->Draw();
  hm23acc->Draw("hist,same");
  hm23acc1->Draw("hist,same");
  hm23acc2->Draw("hist,same");
  hm23dat->Draw("same");

  c1->cd();

}
