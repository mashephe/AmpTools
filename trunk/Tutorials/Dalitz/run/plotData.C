{

  _file0->cd();

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);

  TH2F* hs12s23  = (TH2F*) gDirectory->FindObjectAny("hs12s23");
  TH1F* hm12  = (TH1F*) gDirectory->FindObjectAny("hm12");
  TH1F* hm23  = (TH1F*) gDirectory->FindObjectAny("hm23");
  TH1F* hm13  = (TH1F*) gDirectory->FindObjectAny("hm13");

  c1->cd(1);
  hs12s23->Draw("colz");

  c1->cd(2);
  hm12->Draw("");

  c1->cd(3);
  hm23->Draw("");

  c1->cd(4);
  hm13->Draw("");

  c1->cd();

}
