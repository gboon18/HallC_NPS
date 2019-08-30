{
  gStyle->SetOptStat(0);
  double X0 = 8.9;//PbWO4 radiation length is 0.89 cm
  double Eb = 1e4;//Beam energy
  double percent = 100;//change the unit to %

  double totedep = 0.0;// MeV

  TFile *f = new TFile("test.root");
  TTree *t = (TTree *)f->Get("t_Flux");

  //call trees "and then" open write file
  TFile *f_write = new TFile("output/longitudinal.root", "RECREATE");

  const int array = 1080;

  int x_int, y_int;

  double fCrystal_Z = 200.5;//mm
  double Z = 23; //Crystal length in X0 unit. 22.528088989
  int    bin = 23;

  int    evtNb;
  double edep[array];
  double x[array];
  double y[array];
  double z[array];

  double edep_store[array] = {0.};

  int    evtNb_store = -1;
  vector<int>    evtNb_goodcount;

  double highestEdep = 0.;
  int    highestEdep_x = -1;
  int    highestEdep_y = -1;

  t->SetBranchAddress("evtNb", &evtNb);
  t->SetBranchAddress("edep", &edep);
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);

  TH1D *h = new TH1D("h", "", bin, 0, Z);

  TH2D *h2D_1 = new TH2D("h2D_1", "", 30, 0, 30, 36, 0, 36);
  TH2D *h2D_2 = new TH2D("h2D_2", "", 30, 0, 30, 36, 0, 36);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i); if(i%20000 == 0)cout<<""<<(double)(1.0*i/ t->GetEntries())*100<<"%"<<endl;
    if(evtNb!=evtNb_store) {memset(edep_store, 0., sizeof(edep_store)); evtNb_store = evtNb;}
    for(int j = 0 ; j < array ; j++){
      x_int = int(j/36);
      y_int = j%36;
      edep_store[j] += edep[j];
      h2D_1->Fill(x_int, y_int, edep[j]); 
    }

    t->GetEntry(i+1);
    if(evtNb!=evtNb_store){
      highestEdep = 0.;
      highestEdep_x = -1;
      highestEdep_y = -1;
      for(int k = 0 ; k < array ; k++){
	x_int = int(k/36);
	y_int = k%36;
	if(edep_store[k]>highestEdep){
	  highestEdep = edep_store[k];
	  highestEdep_x = x_int;
	  highestEdep_y = y_int;
	}
      }
      if((highestEdep_x!=0 && highestEdep_x!=1 && highestEdep_x!=28 && highestEdep_x!=29) 
	 && (highestEdep_y!=0 && highestEdep_y!=1 && highestEdep_y!=34 && highestEdep_y!=35)
	 && (highestEdep_x!=-1 && highestEdep_y!=-1)){
	evtNb_goodcount.push_back(evtNb_store);
      }   
    }
  }
  totedep = evtNb_goodcount.size()*Eb;

  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i);  if(i%20000 == 0)cout<<"Filling: "<<(double)(1.0*i/ t->GetEntries())*100<<"%"<<endl;
    //    if(find(evtNb_goodcount.begin(), evtNb_goodcount.end(), evtNb)!=evtNb_goodcount.end()){}// does not work in ROOT v5
    for(int k = 0 ; k < evtNb_goodcount.size() ; k++){
      if(evtNb == evtNb_goodcount[k]){
	for(int j = 0 ; j < array ; j++){
	  x_int = int(j/36);
	  y_int = j%36;
	  h2D_2->Fill(x_int, y_int, edep[j]);
	  h->Fill((z[j] + 0.5*fCrystal_Z)/X0, edep[j]/totedep*percent);
	}
      }
    }
  }

  cout<<"totedep     : "<<totedep<<" MeV"<<endl;
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_pattern = new TCanvas("c_pattern", "", 200, 10, 300*2, 360*2);
  TCanvas *c_pattern_2 = new TCanvas("c_pattern_2", "", 200, 10, 300*2, 360*2);
  TCanvas *c_edep    = new TCanvas("c_edep",    "", 200, 10, 300*2, 360*2);

  c_pattern->cd();
  h2D_1->Draw("colz");
  c_pattern_2->cd();
  h2D_2->Draw("colz");

  c_edep->cd();
  gPad->SetLogy(1);
  h->SetTitle("Longitudinal energy deposition in calorimeter");
  h->GetXaxis()->SetTitle("X_{0}");
  h->GetYaxis()->SetTitle("#frac{1}{E}#frac{dE}{dX_{0}} [%]");
  h->GetXaxis()->SetRangeUser(0, 100);
  h->SetLineColor(1);
  h->SetMarkerColor(1);
  h->Draw();

  TLegend *leg = new TLegend(.3, .1, .6, .4);
  leg->AddEntry(h, "Longitudinal profile", "l");
  leg->Draw();

  f_write->Write();
  // f_write->Close();
  evtNb_goodcount.clear();
}
