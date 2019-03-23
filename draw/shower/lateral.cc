{
  gStyle->SetOptStat(0);
  double Eb = 1e4;//Beam energy
  double percent = 100;//change the unit to %

  double totedep = 0.0;// MeV

  TFile *f = new TFile("test.root");
  TTree *t = f->Get("t_Flux");

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

  double edep_store[array] = {0.};

  int    evtNb_store = -1;
  vector<int>    evtNb_goodcount;
  int    x_store = -1;
  int    y_store = -1;

  double x_pos_store = -1.;
  double y_pos_store = -1.;

  double highestEdep = 0.;
  int    highestEdep_x = -1;
  int    highestEdep_y = -1;

  t->SetBranchAddress("evtNb", &evtNb);
  t->SetBranchAddress("edep", &edep);
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);

  TH1D *h = new TH1D("h", "", bin, 0, R);

  TH2D *h2D_1 = new TH2D("h2D_1", "", 30, 0, 30, 36, 0, 36);
  TH2D *h2D_2 = new TH2D("h2D_2", "", 30, 0, 30, 36, 0, 36);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i); if(i%20000 == 0)cout<<""<<(double)(1.0*i /t->GetEntries())*100<<"%"<<endl;
    if(evtNb!=evtNb_store) {memset(edep_store, 0., sizeof(edep_store)); evtNb_store = evtNb;}
    for(int j = 0 ; j < array ; j++){
      x = int(j/36);
      y = j%36;
      edep_store[j] += edep[j];
      h2D_1->Fill(x, y, edep[j]); 
    }

    t->GetEntry(i+1);
    if(evtNb!=evtNb_store){
      highestEdep = 0.;
      highestEdep_x = -1;
      highestEdep_y = -1;
      for(int k = 0 ; k < array ; k++){
  	x = int(k/36);
  	y = k%36;
  	if(edep_store[k]>highestEdep){
  	  highestEdep = edep_store[k];
  	  highestEdep_x = x;
  	  highestEdep_y = y;
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
  /////////////////////
  evtNb_store   = -1;//
  /////////////////////
  for(int i = 0 ; i < t->GetEntries() ; i++){
    t->GetEntry(i); if(i%20000 == 0)cout<<"Filling: "<<(double)(1.0*i/t->GetEntries())*100<<"%"<<endl;
    //    if(find(evtNb_goodcount.begin(), evtNb_goodcount.end(), evtNb)!=evtNb_goodcount.end()){}// does not work in ROOT v5
    for(int k = 0 ; k < evtNb_goodcount.size() ; k++){
      if(evtNb == evtNb_goodcount[k]){
  	for(int j = 0 ; j < array ; j++){
  	  x = int(j/36);
  	  y = j%36;
  	  h2D_2->Fill(x, y, edep[j]);
  	  if(edep[j] && evtNb!=evtNb_store){

  	    evtNb_store = evtNb;

  	    x_store = x;
  	    y_store = y;

  	    x_pos_store = x[j];
  	    y_pos_store = y[j];

  	    h->Fill(0, edep[j]/totedep*percent);
  	  }
  	  else{
  	    h->Fill(sqrt( (x_pos_store + (-x[j] + (x - x_store)*X))*(x_pos_store + (-x[j] + (x - x_store)*X)) + (y_pos_store - (y[j] + (y - y_store)*Y))*(y_pos_store - (y[j] + (y - y_store)*Y)) ), edep[j]/totedep*percent );
  	  }
  	}
      }
    }
  }

  cout<<"totedep : "<<totedep<<" MeV"<<endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_pattern = new TCanvas("c_pattern", "1", 200, 10, 300*2, 360*2);
  TCanvas *c_pattern_2 = new TCanvas("c_pattern_2", "2", 200, 10, 300*2, 360*2);
  TCanvas *c_edep    = new TCanvas("c_edep",    "", 200, 10, 300*2, 360*2);

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
  h->SetLineColor(1);
  h->SetMarkerColor(1);
  h->Draw();

  TLegend *leg = new TLegend(.4, .6, .9, .9);
  leg->AddEntry(h, "Lateral profile", "l");
  leg->Draw();

  f_write->Write();
  //  f_write->Close();
  evtNb_goodcount.clear();
}
