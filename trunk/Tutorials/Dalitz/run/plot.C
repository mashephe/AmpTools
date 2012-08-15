{

  _file0->cd();

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

}
