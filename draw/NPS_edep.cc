{
  gStyle->SetOptStat(0);
  TFile *f = new TFile("../jlab_HallC_NPS_PbWO4/build/test.root", "READ");
  TTree *t = (TTree*)f->Get("t");

  const int array = 1116;

  double edep[array] = {0.};
  int op_pc[array] = {0};

  t->SetBranchAddress("edep", edep);
  t->SetBranchAddress("op_pc", op_pc);

  int x, y;

  TH2D* h_edep = new TH2D("h_edep", "energy deposition in 2D", 30, 0, 30, 36, 0, 36);
  TH2D* h_op = new TH2D("h_op", "number of scintillated photon arrived at the PMT", 30, 0, 30, 36, 0, 36);
  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i);

    for(int j = 0 ; j < array ; j++){
      x = int(j/36);
      y = j%36;
      h_edep->Fill(x, y, edep[j]);
      h_op->Fill(x, y, op_pc[j]);
    }
  }

  TCanvas* c1 = new TCanvas("c1", "", 200, 10, 1000, 500);
  c1->Divide(2,1);
  c1->cd(1);
  h_edep->SetTitle("Total energy depostition in each crystal");
  h_edep->Draw("colz");
  c1->cd(2);
  h_op->SetTitle("Total scintillated photon collected in each PMT");
  h_op->Draw("colz");
}
