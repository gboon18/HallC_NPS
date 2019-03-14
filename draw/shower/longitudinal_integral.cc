{
  TFile *f = new TFile("output/longitudinal.root", "READ");
  TCanvas *c1 = new TCanvas("c1", "", 200, 10, 700, 500);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  h_10GeV->GetXaxis()->SetTitle("X_{0}");
  h_10GeV->Draw();

  TLegend *leg = new TLegend(.6, .6, .9, .9 );
  leg->AddEntry(h_10GeV, "10GeV", "l");
  leg->Draw();

  c1->Print("output/output/longitudinal.pdf");
  gPad->SetLogy(1);
  c1->Print("output/output/longitudinal_logy.pdf");
  ////////////////////////////////////////////////////////////////////

  TH1* hc_10GeV = h_10GeV->GetCumulative();

  hc_10GeV->Scale(100/h_10GeV->Integral());

  TCanvas *c2 = new TCanvas("c2", "", 200, 10, 700, 500);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(1);
  hc_10GeV->SetTitle("Cumulated energy deposition in Calorimeter");
  hc_10GeV->GetXaxis()->SetTitle("X_{0}");
  hc_10GeV->GetYaxis()->SetTitle("Energy deposition [%]");
  hc_10GeV->Draw();

  TLegend *leg_c = new TLegend(0.6, 0.1, 0.9, 0.4 );
  leg_c->AddEntry(hc_10GeV, "10GeV", "l");
  leg_c->Draw();

  c2->Print("output/output/longitudinal_integral_logy.pdf");
  c2->cd();
  gPad->SetLogy(0);
  c2->Update();
  c2->Print("output/output/longitudinal_integral.pdf");

}
