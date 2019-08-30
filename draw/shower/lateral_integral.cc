{
  TFile *f = new TFile("output/lateral.root", "READ");
  TCanvas *c1 = new TCanvas("c1", "Lateral energy profile", 200, 10, 700, 500);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(1);
  //  h->GetYaxis()->SetRangeUser(1, 11);
  h->Draw();

  c1->Update();

  TLine *l2_1 = new TLine(21.9, 10*(gPad->GetUymin()), 21.9, 10*(gPad->GetUymax()));
  l2_1->SetLineStyle(2);
  l2_1->SetLineColor(6);
  l2_1->Draw();

  TLine *l2_2 = new TLine(21.9*2, 10*(gPad->GetUymin()), 21.9*2, 10*(gPad->GetUymax()));
  l2_2->SetLineStyle(2);
  l2_2->SetLineColor(6);
  l2_2->Draw();

  TLine *l2_3 = new TLine(21.9*3, 10*(gPad->GetUymin()), 21.9*3, 10*(gPad->GetUymax()));
  l2_3->SetLineStyle(2);
  l2_3->SetLineColor(6);
  l2_3->Draw();

  TLine *l2_4 = new TLine(21.9*4, 10*(gPad->GetUymin()), 21.9*4, 10*(gPad->GetUymax()));
  l2_4->SetLineStyle(2);
  l2_4->SetLineColor(6);
  l2_4->Draw();

  TLegend *leg = new TLegend(.6, .6, .9, .9 );
  leg->AddEntry(h, "Lateral profile", "l");
  leg->Draw();

  c1->Print("output/output/lateral_logy.pdf");
  ////////////////////////////////////////////////////////////////////

  TH1* hc = h->GetCumulative();

  hc->Scale(100/h->Integral());

  TCanvas *c2 = new TCanvas("c2", "Total lateral energy deposition", 200, 10, 700, 500);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(0);
  hc->SetTitle("Total lateral energy deposition in Calorimeter");
  hc->GetXaxis()->SetTitle("R [mm]");
  hc->GetYaxis()->SetTitle("Energy deposition [%]");
  //  hc->GetYaxis()->SetRangeUser(1, 110);
  hc->Draw();

  c2->Update();
  TLine *lc_1 = new TLine(21.9, gPad->GetUymin(), 21.9, gPad->GetUymax());
  lc_1->SetLineStyle(2);
  lc_1->SetLineColor(6);
  lc_1->Draw();

  TLine *lc_2 = new TLine(21.9*2, gPad->GetUymin(), 21.9*2, gPad->GetUymax());
  lc_2->SetLineStyle(2);
  lc_2->SetLineColor(6);
  lc_2->Draw();

  TLine *lc_3 = new TLine(21.9*3, gPad->GetUymin(), 21.9*3, gPad->GetUymax());
  lc_3->SetLineStyle(2);
  lc_3->SetLineColor(6);
  lc_3->Draw();

  TLine *lc_4 = new TLine(21.9*4, gPad->GetUymin(), 21.9*4, gPad->GetUymax());
  lc_4->SetLineStyle(2);
  lc_4->SetLineColor(6);
  lc_4->Draw();

  TLegend *leg_c = new TLegend(0.6, 0.1, 0.9, 0.4 );
  leg_c->AddEntry(hc, "Total lateral energy deposition", "l");
  leg_c->Draw();

  c2->Print("output/output/lateral_integral.pdf");
}
