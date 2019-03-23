{
  TFile *f = new TFile("output/longitudinal.root", "READ");
  TCanvas *c1 = new TCanvas("c1", "Longitudinal energy profile", 200, 10, 700, 500);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(1);
  h->Draw();

  TLegend *leg = new TLegend(.6, .6, .9, .9 );
  leg->AddEntry(h, "Longitudinal profile", "l");
  leg->Draw();

  c1->Print("output/output/longitudinal_logy.pdf");
  ////////////////////////////////////////////////////////////////////

  TH1* hc = h->GetCumulative();

  hc->Scale(100/h->Integral());

  TCanvas *c2 = new TCanvas("c2", "Total longitudinal energy deposition", 200, 10, 700, 500);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(0);
  hc->SetTitle("Total longitudinal energy deposition in Calorimeter");
  hc->GetXaxis()->SetTitle("X_{0}");
  hc->GetYaxis()->SetTitle("Energy deposition [%]");
  hc->Draw();

  TLegend *leg_c = new TLegend(0.6, 0.1, 0.9, 0.4 );
  leg_c->AddEntry(hc, "Total longitudinal energy deposition", "l");
  leg_c->Draw();

  c2->Print("output/output/longitudinal_integral.pdf");

}
