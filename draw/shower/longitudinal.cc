{
  gStyle->SetOptStat(0);
  double X0 = 8.9;//PbWO4 radiation length is 0.89 cm
  double Eb = 1e4;//Beam energy
  double percent = 100;//change the unit to %

  double totedep = 0.0;

  TFile *f_10GeV = new TFile("output.root");
  TTree *t_10GeV = f_10GeV->Get("t_Flux");

  //call trees "and then" open write file
  TFile *f_write = new TFile("output/longitudinal.root", "RECREATE");

  const int array = 1080;

  int x, y;

  double fCrystal_Z = 200.5;//mm
  double Z = 23; //Crystal length in X0 unit. 22.528088989
  int    bin = 23;//

  int    evtNb_10GeV;
  double edep_10GeV[array];
  double x_10GeV[array];
  double y_10GeV[array];
  double z_10GeV[array];

  int    evtNb_store = -1;
  int    evtNb_notgood = -1;
  int    evtNb_goodcount = 0;
  int    evtNb_totcount = 0;

  t_10GeV->SetBranchAddress("evtNb", &evtNb_10GeV);
  t_10GeV->SetBranchAddress("edep", &edep_10GeV);
  t_10GeV->SetBranchAddress("x", &x_10GeV);
  t_10GeV->SetBranchAddress("y", &y_10GeV);
  t_10GeV->SetBranchAddress("z", &z_10GeV);

  TH1D *h_10GeV = new TH1D("h_10GeV", "", bin, 0, Z);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i = 0 ; i < t_10GeV->GetEntries() ; i++){
    t_10GeV->GetEntry(i);
    for(int j = 0 ; j < array ; j++){
      x = int(j/36);
      y = j%36;
      if(edep_10GeV[j]){if(i%20000 == 0) cout<<"totEdep calculating : "<<(double)(1.0*i/t_10GeV->GetEntries())*100<<"%"<<endl;
	if(evtNb_10GeV!=evtNb_store){ 
	  evtNb_totcount++;
	  if ((x!=0 && x!=1 && x!=28 && x!=29) && (y!=0 && y!=1 && y!=34 && y!=35)){
	    evtNb_goodcount++;
	  }
	  else {
	    evtNb_notgood = evtNb_10GeV;
	  }
	}
      }
    }
    evtNb_store = evtNb_10GeV;
  }

  totedep = evtNb_goodcount*Eb;
  cout<<"Total counts  : "<<evtNb_totcount<<endl;
  cout<<"Good  counts  : "<<evtNb_goodcount<<endl;
  cout<<"totedep       : "<<totedep<<" MeV"<<endl;

  /////////////////////
  evtNb_store   = -1;//
  evtNb_notgood = -1;//
  /////////////////////

  for(int i = 0 ; i < t_10GeV->GetEntries() ; i++){
    t_10GeV->GetEntry(i);
    for(int j = 0 ; j < array ; j++){
      x = int(j/36);
      y = j%36;
      if(edep_10GeV[j]){if(i%20000 == 0) cout<<"Almost done : "<<(double)(1.0*i/t_10GeV->GetEntries())*100<<"%"<<endl;
	if(evtNb_10GeV!=evtNb_store){
	  if((x!=0 && x!=1 && x!=28 && x!=29) && (y!=0 && y!=1 && y!=34 && y!=35)){
	    h_10GeV->Fill((z_10GeV[j] + 0.5*fCrystal_Z)/X0, edep_10GeV[j]/totedep*percent);
	  }
	  else{
	    evtNb_notgood = evtNb_10GeV;
	  }
	}
	else if(evtNb_10GeV != evtNb_notgood){ 
	  h_10GeV->Fill((z_10GeV[j] + 0.5*fCrystal_Z)/X0, edep_10GeV[j]/totedep*percent );
	}
      }
    }
    evtNb_store = evtNb_10GeV;
  }

  cout<<"Total counts  : "<<evtNb_totcount<<endl;
  cout<<"Good  counts  : "<<evtNb_goodcount<<endl;
  cout<<"totedep       : "<<totedep<<" MeV"<<endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_edep    = new TCanvas("c_edep",    "", 200, 10, 300*2, 360*2);

  c_edep->cd();
  gPad->SetLogy(1);
  h_10GeV->SetTitle("Longitudinal energy deposition in calorimeter");
  h_10GeV->GetXaxis()->SetTitle("X_{0} [mm]");
  h_10GeV->GetYaxis()->SetTitle("#frac{1}{E}#frac{dE}{dX_{0}} [%]");
  h_10GeV->GetXaxis()->SetRangeUser(0, 100);
  h_10GeV->SetLineColor(1);
  h_10GeV->SetMarkerColor(1);
  h_10GeV->Draw();

  TLegend *leg = new TLegend(.3, .1, .6, .4);
  leg->AddEntry(h_10GeV, "10GeV electron", "l");
  leg->Draw();

  f_write->Write();
  f_write->Close();
}
