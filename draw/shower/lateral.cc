{
  gStyle->SetOptStat(0);
  double Eb = 1e4;//Beam energy
  double percent = 100;//change the unit to %

  double totedep = 0.0;

  TFile *f = new TFile("test.root");
  TTree *t = f->Get("t2");

  //call trees "and then" open write file
  TFile *f_write = new TFile("output/lateral.root", "RECREATE");

  const int array = 1080;

  int x, y;

  double fCrystal_Y = 20.5;//mm
  double fCrystal_X = 20.5;//mm
  double fWrapThickness = 65*1e-3;//mm
  double gap = 0.;
  double Y = fCrystal_Y + 2*fWrapThickness + gap;
  double X = fCrystal_X + 2*fWrapThickness + gap;
  double R = 309; //make a circle with radius of NPS_X
  int    bin = 309;//2R = 618.9mm, R=309.45mm

  int    evtNb;
  double edep[array];
  double x[array];
  double y[array];
  double z[array];

  int    evtNb_store = -1;
  int    evtNb_notgood = -1;
  int    evtNb_goodcount = 0;
  int    evtNb_totcount = 0;
  int    j_store = -1;
  int    x_store = -1;
  int    y_store = -1;

  double x_pos_store = -1.;
  double y_pos_store = -1.;

  t->SetBranchAddress("evtNb", &evtNb);
  t->SetBranchAddress("edep", &edep);
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);

  TH1D *h = new TH1D("h", "Lateral energy profile", bin, 0, R);

  TH2D *h2D_1 = new TH2D("h2D_1", "Pattern on the calorimeter", 30, 0, 30, 36, 0, 36);//in x axis
  TH2D *h2D_2 = new TH2D("h2D_2", "Pattern on the calorimeter except hits on the edge", 30, 0, 30, 36, 0, 36);//in x axis

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i);
    for(int j = 0 ; j < array ; j++){
      x = int(j/36);
      y = j%36;
      if(edep[j]){if(i%20000 == 0) cout<<"totEdep calculating : "<<(double)(1.0*i/t->GetEntries())*100<<"%"<<endl;
	if(evtNb!=evtNb_store){ 
	  evtNb_totcount++;
	  if ((x!=0 && x!=1 && x!=28 && x!=29) && (y!=0 && y!=1 && y!=34 && y!=35)){
	    evtNb_goodcount++;
	  }
	  else{
	    evtNb_notgood = evtNb;
	  }
	}
      }
    }
    evtNb_store = evtNb;
  }

  totedep = evtNb_goodcount*Eb;
  cout<<"Total counts  : "<<evtNb_totcount<<endl;
  cout<<"Good  counts  : "<<evtNb_goodcount<<endl;
  cout<<"totedep : "<<totedep<<" MeV"<<endl;

  /////////////////////
  evtNb_store   = -1;//
  evtNb_notgood = -1;//
  /////////////////////

  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i);
    for(int j = 0 ; j < array ; j++){
      x = int(j/36);
      y = j%36;
      if(edep[j]){if(i%20000 == 0) cout<<""<<(double)(1.0*i/t->GetEntries())*100<<"%"<<endl;
	if(evtNb!=evtNb_store){
	  if((x!=0 && x!=1 && x!=28 && x!=29) && (y!=0 && y!=1 && y!=34 && y!=35)){
	    j_store = j;
	    x_store = x;
	    y_store = y;

	    x_pos_store = x[j];
	    y_pos_store = y[j];

	    h->Fill(0, edep[j]/totedep*percent);
	    h2D_2->Fill(x, y, edep[j]); 
	  }
	  else{
	    evtNb_notgood = evtNb;
	  }
	}
	else if(evtNb != evtNb_notgood){ 
	  h->Fill(sqrt( (x_pos_store + (-x[j] + (x - x_store)*X))*(x_pos_store + (-x[j] + (x - x_store)*X)) + (y_pos_store - (y[j] + (y - y_store)*Y))*(y_pos_store - (y[j] + (y - y_store)*Y)) ), edep[j]/totedep*percent );

	  h2D_2->Fill(x, y, edep[j]); 
	}

	h2D_1->Fill(x, y, edep[j]); 
      }
    }
    evtNb_store = evtNb;
  }

  cout<<"Total counts  : "<<evtNb_totcount<<endl;
  cout<<"Good  counts  : "<<evtNb_goodcount<<endl;
  cout<<"totedep : "<<totedep<<" MeV"<<endl;

  cout<<"totedep     : "<<totedep<<" MeV"<<endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_pattern = new TCanvas("c_pattern", "Pattern on the calorimeter", 200, 10, 300*2, 360*2);
  TCanvas *c_pattern_2 = new TCanvas("c_pattern_2", "Pattern on the calorimeter except hits on the edge", 200, 10, 300*2, 360*2);
  TCanvas *c_edep    = new TCanvas("c_edep",    "Lateral energy profile", 200, 10, 300*2, 360*2);

  c_pattern->cd();
  h2D_1->Draw("colz");
  c_pattern_2->cd();
  h2D_2->Draw("colz");

  c_edep->cd();
  gPad->SetLogy(1);
  h->SetTitle("Lateral energy deposition in calorimeter");
  h->GetXaxis()->SetTitle("R [mm]");
  h->GetYaxis()->SetTitle("#frac{1}{E}#frac{dE}{dR} [%]");
  h->GetXaxis()->SetRangeUser(0, 100);
  //  h->GetYaxis()->SetRangeUser(1, 100);
  h->SetLineColor(1);
  h->SetMarkerColor(1);
  h->Draw();

  TLegend *leg = new TLegend(.4, .6, .9, .9);
  leg->AddEntry(h, "Lateral profile", "l");
  leg->Draw();

  f_write->Write();
  // f_write->Close();
}
