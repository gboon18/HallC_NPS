//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/shared/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
// $Id: EventAction.cc 92322 2015-08-27 14:54:05Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"

#include "HistoManager.hh"

#include "G4Event.hh"

#include "B5HadCalorimeterHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "CrystalCoverHit.hh"
#include "CrystalFrontCoverHit.hh"
#include "PMTcoverHit.hh"

#include "FluxHit.hh"

#include "dvcsGlobals.hh"

#include "TProcessID.h"

#include "RFunction.C"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run, HistoManager* histo)
  :G4UserEventAction(),
   fRunAct(run),
   fHistoManager(histo),
   fEvtNb(-1),
   fHCHCID(-1),fHadCalEdep(),
   fPID(),
   fOP_sc(), fOP_ce(),
   fCrystCoverHCID(-1), fCrystCoverOP(), fCrystFrontCoverHCID(-1), fCrystFrontCoverOP(), fPMTcoverHCID(-1), fPMTcoverOP(),
   fFluxEne(),
   //reuse for Pseudo Sphere
   fFlxHCID(-1),
   fFluxEnergy(-999), fFluxPID(-999),
   fPos(-999),
   hit_HMS(0),

  M_targ(-999.),
  fSmeared_vertex_z(-999.),
  fVertex_z(-999.),

  fEnergyAbs(0.), fEnergyGap(0.),
  fTrackLAbs(0.), fTrackLGap(0.),
  fPrintModulo(0)                             
{
  fPrintModulo = 100; 
  fHadCalEdep.resize(1080,0.);
  fPID.resize(1080,0),
    fOP_sc.resize(1080,0), fOP_ce.resize(1080,0),
    fCrystCoverOP.resize(1080,0), fCrystFrontCoverOP.resize(1080,0), fPMTcoverOP.resize(1080,0);
  fFluxEne = 0.;

  G4int run_number = dvcsGlobals::run_number;
  dvcs_evt = new TDVCSEvent(run_number);
  calo_evt = new TCaloEvent(run_number);
  L_calo_phot = new TLorentzVector();//If I do not do it, it crashes.....

  dvcs_evt->SetCaloEvent(calo_evt);
  G4cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111111111111111Calorimter distance (cm): "<<dvcs_evt->GetGeometry()->GetCaloDist()<<G4endl;

  G4cout<<"Caloriemeter angle (deg): "<<dvcs_evt->GetGeometry()->GetCaloTheta()*TMath::RadToDeg()<<G4endl;
  G4int target_gen_proc_type = dvcsGlobals::target_gen_proc_type;
  if( target_gen_proc_type == 0 || target_gen_proc_type == 1 || target_gen_proc_type ==2)
    {
      if(target_gen_proc_type == 0) {M_targ = 0.9383;}
      else if( target_gen_proc_type == 1 ) {M_targ = 0.939565378;}
      else if( target_gen_proc_type == 2 ) {M_targ = 1.875612859;}
    }
  else
    {
      cout<<"Event Action!! Error: wrong target mass "<<endl;
      cout<<"The program is exiting "<<endl;
      exit(1);
    }

  fPx[5] = {-999.};
  fPy[5] = {-999.};
  fPz[5] = {-999.};
  fE[5]  = {-999.};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  fEvtNb = evt->GetEventID();
  if (fEvtNb%fPrintModulo == 0) 
 
    // initialisation per event
    fEnergyAbs = fEnergyGap = 0.;
  fTrackLAbs = fTrackLGap = 0.;
  
  fFluxEne = 0.;

  if (fHCHCID==-1) {

    //cc-in2p3 does not understand(or geant4.10.02.0 does not understand) "auto"
    //    auto sdManager = G4SDManager::GetSDMpointer();
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fHCHCID = sdManager->GetCollectionID("HadCalorimeter/HadCalorimeterColl");//B2/TrackerChamberSD is the SDname, TrackerHitsCollection is the collection ID.
  }

  if (fCrystCoverHCID==-1) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fCrystCoverHCID = sdManager->GetCollectionID("CrystalCover/CrystalCoverColl");
  }

  if (fCrystFrontCoverHCID==-1) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fCrystFrontCoverHCID = sdManager->GetCollectionID("CrystalFrontCover/CrystalFrontCoverColl");
  }

  if (fPMTcoverHCID==-1) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fPMTcoverHCID = sdManager->GetCollectionID("PMTcover/PMTcoverColl");
  }

  //reuse for Pseudo Sphere

  if (fFlxHCID==-1) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fFlxHCID = sdManager->GetCollectionID("Flux/FluxColl");//B2/TrackerChamberSD is the SDname, TrackerHitsCollection is the collection ID.
  }

  if( dvcsGlobals::hit_HMS_CALO_flag )
    {
      fPx[4] = 0.;
      fPy[4] = 0.;
      fPz[4] = 0.;
      fE[4] = 0.;
    }
  fHistoManager->SetEvtNb(fEvtNb);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  //
  fRunAct->fillPerEvent(fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap);
  
  //cc-in2p3 does not understand(or geant4.10.02.0 does not understand) "auto"
  //  auto hce = evt->GetHCofThisEvent();
  G4HCofThisEvent* hce = evt->GetHCofThisEvent();
  if (!hce) {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
		"Code001", JustWarning, msg);
    return;
  }
  
  // Get hits collections
  //cc-in2p3 does not understand(or geant4.10.02.0 does not understand) "auto" 
  //  auto hcHC 
  B5HadCalorimeterHitsCollection* hcHC 
    = static_cast<B5HadCalorimeterHitsCollection*>(hce->GetHC(fHCHCID));
  if (!hcHC) {
    G4ExceptionDescription msg;
    msg << "Some of hits collections of this event not found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
		"Code001", JustWarning, msg);
    return;
  }   

  CrystalCoverHitsCollection* CrystCoverHC 
    = static_cast<CrystalCoverHitsCollection*>(hce->GetHC(fCrystCoverHCID));
  if (!CrystCoverHC) {
    G4ExceptionDescription msg;
    msg << "Some of hits collections of this event not found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
		"Code001", JustWarning, msg);
    return;
  }   
  CrystalFrontCoverHitsCollection* CrystFrontCoverHC 
    = static_cast<CrystalFrontCoverHitsCollection*>(hce->GetHC(fCrystFrontCoverHCID));
  if (!CrystFrontCoverHC) {
    G4ExceptionDescription msg;
    msg << "Some of hits collections of this event not found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
		"Code001", JustWarning, msg);
    return;
  }   
  PMTcoverHitsCollection* PMTcoverHC 
    = static_cast<PMTcoverHitsCollection*>(hce->GetHC(fPMTcoverHCID));
  if (!PMTcoverHC) {
    G4ExceptionDescription msg;
    msg << "Some of hits collections of this event not found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
		"Code001", JustWarning, msg);
    return;
  }   

  // HCEnergy
  for (G4int i=0;i<1080;i++)
    {
      B5HadCalorimeterHit* hit = (*hcHC)[i];
      G4double eDep = hit->GetEdep();//total energy deposition of each crystals
      G4int PID = hit->GetPID();
      // SetEnrgyandPID must be done for all 1080 crystals.
      // If not, the edep from last events with edep>0 will be stored again. 
      fHadCalEdep[i] = eDep*1.04;//Compensate lost 4% of the total energy of the shower at the end of the crystal.
      fPID[i] = PID;

      CrystalCoverHit* CChit = (*CrystCoverHC)[i];
      CrystalFrontCoverHit* CFChit = (*CrystFrontCoverHC)[i];
      PMTcoverHit* PMTChit = (*PMTcoverHC)[i];
      G4double sc = hit->GetOPInt_sc();
      G4double ce = hit->GetOPInt_ce();
      fOP_sc[i] = sc;
      fOP_ce[i] = ce;

      //No. of OP reflected at the side of the crystal wrapper.
      G4int CrystalCoverOP = CChit->GetOPInt();  
      fCrystCoverOP[i] = CrystalCoverOP;
      //No. of OP reflected at the front of the crystal wrapper
      G4int CrystalFrontCoverOP = CFChit->GetOPInt();  
      fCrystFrontCoverOP[i] = CrystalFrontCoverOP;
      //No. of OP arrived at the PMT cover
      G4int PMTcoverOP = PMTChit->GetOPInt();  
      fPMTcoverOP[i] = PMTcoverOP;

      fHistoManager->SetEnergyandPID( i, fPID[i], fHadCalEdep[i], fOP_sc[i], fOP_ce[i], fCrystCoverOP[i], fCrystFrontCoverOP[i], fPMTcoverOP[i]);
    }
  fHistoManager->FillNtuple();

  if( dvcsGlobals::hit_HMS_CALO_flag )
    {
      G4double Ebeam = dvcsGlobals::Ebeam; 
      G4int ObjectNumber=TProcessID::GetObjectCount(); 

      for (G4int i = 0 ; i < 1080 ; i++)
	{
	  B5HadCalorimeterHit* hit = (*hcHC)[i];
	  G4double eDep = hit->GetEdep();//total energy deposition of each crystals       
	  TCaloBlock* block=calo_evt->AddBlock(i);
	  //energy and time
	  //Currently, no time info in Geant4 simulation. default as 0 for now.
	  block->AddPulse(eDep*(1e-3)*(1.04), 0.);//energy (GeV) and time (ns) //Compensate lost 4% of the total energy of the shower at the end of the crystal.
	}
      dvcs_evt->SetVertex(0., 0., fSmeared_vertex_z);//[cm] unit

      calo_evt->TriggerSim(0.5);//0.5 GeV // if too many blocks in the event
      calo_evt->DoClustering(-3., 3.);//-3ns to +3ns windows
      G4int Nb_clust = calo_evt->GetNbClusters();
      for(G4int i = 0 ; i < Nb_clust ; i++ ){

	TCaloCluster *clus = calo_evt->GetCluster(i);
	clus->Analyze(1., -3., 3.);//"true||false", "time_min(ns)", "time_max(ns)", "weight"
	*L_calo_phot = dvcs_evt->GetPhoton(i);
	if(Nb_clust) fHistoManager->SetClusterPositionandMomentum(evt->GetEventID(), clus->GetX(), clus->GetY(), L_calo_phot->E(), L_calo_phot->Px(), L_calo_phot->Py(), L_calo_phot->Pz());
      }
      if(Nb_clust >= 1 && hit_HMS){
	G4double clust_x   = calo_evt->GetCluster(0)->GetX();
	G4double clust_y   = calo_evt->GetCluster(0)->GetY();
	G4double clust_ene = calo_evt->GetCluster(0)->GetE();
	G4int clust_size   = calo_evt->GetCluster(0)->GetClusSize();
	fHistoManager->SetCaloData(clust_ene, clust_x, clust_y, clust_size);

	G4double phot_px = L_calo_phot->Px();//GeV
	G4double phot_py = L_calo_phot->Py();
	G4double phot_pz = L_calo_phot->Pz();

	fHistoManager->SetPhotRec(phot_px, phot_py, phot_pz); //Set Reconstructed photon kinematics: Px, Py, Pz

	// fHistoManager->SetClusterPositionandMomentum(evt->GetEventID(), clust_x, clust_y, clust_ene, phot_px, phot_py, phot_pz);
	fHistoManager->FillNtupleClust();

	TVector3 k(0, 0, Ebeam);
	TVector3 kp(fPx[4], fPy[4], fPz[4]);
	TVector3 qp(phot_px, phot_py, phot_pz);

	G4double Q2 = 2*k.Mag()*kp.Mag()*(1 - cos(k.Angle(kp)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))

	TVector3 q = k - kp;

	TVector3 v1=q.Cross(kp);
	TVector3 v2=q.Cross(qp);
	G4double fphi=v1.Angle(v2);
	if(q.Dot(v1.Cross(v2))<0) fphi=2.*TMath::Pi()-fphi;

	G4double cos=(q.Dot(qp))/(q.Mag()*qp.Mag());
	G4double nu=Ebeam-kp.Mag();
	G4double tM=(Q2*M_targ+2.*nu*M_targ*(nu-sqrt(nu*nu+Q2)*cos))/(sqrt(nu*nu+Q2)*cos-nu-M_targ);

	G4double xB = Q2/(2*M_targ*nu);

	fHistoManager->SetKineRec(tM, xB, Q2, fphi);

	TLorentzVector Lb, Lp, Lgp, Lkp, L_mis;
	Lb.SetPxPyPzE(0, 0, Ebeam, Ebeam);
	Lp.SetPxPyPzE(0, 0, 0, M_targ);
	Lkp.SetPxPyPzE(fPx[4], fPy[4], fPz[4], fE[4]);
	Lgp = *L_calo_phot;

	G4double mm2 = (Lb + Lp - Lgp - Lkp).M2();

	fHistoManager->SetClustW2(mm2);

	fHistoManager->FillNtupleDVCS();
      }

      hit_HMS = false;
      calo_evt->Reset();
      TProcessID::SetObjectCount(ObjectNumber);
    }

  // Get hits collections
  //cc-in2p3 does not understand(or geant4.10.02.0 does not understand) "auto" 
  //  auto hcHC 
  FluxHitsCollection* fluxHC 
    = static_cast<FluxHitsCollection*>(hce->GetHC(fFlxHCID));
  if (!fluxHC) {
    G4ExceptionDescription msg;
    msg << "Some of hits collections of this event not found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
		"Code001", JustWarning, msg);
    return;
  }   
  FluxHit* fluxhit = (*fluxHC)[0];
  G4double energy = fluxhit->GetEdep();//energy of the particle passing through the volume
  G4int fluxPID = fluxhit->GetPID();
  G4ThreeVector position = fluxhit->GetPos();
  fFluxEnergy = energy;
  fFluxPID = fluxPID;
  fPos = position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::DefinePrimaries(TLorentzVector* L_elec_init, TLorentzVector* L_elec_scat, TLorentzVector* L_phot_final, TLorentzVector* L_prot_final, G4double vertex_z)
{
  fPx[0] = L_elec_init->Px(); fPx[1] = L_elec_scat->Px(); fPx[2] = L_phot_final->Px(); fPx[3] = L_prot_final->Px();
  fPy[0] = L_elec_init->Py(); fPy[1] = L_elec_scat->Py(); fPy[2] = L_phot_final->Py(); fPy[3] = L_prot_final->Py();
  fPz[0] = L_elec_init->Pz(); fPz[1] = L_elec_scat->Pz(); fPz[2] = L_phot_final->Pz(); fPz[3] = L_prot_final->Pz();
  fE[0]  = L_elec_init->E();  fE[1]  = L_elec_scat->E();  fE[2]  = L_phot_final->E();  fE[3]  = L_prot_final->E();

  G4double HMS_angle = dvcsGlobals::HMS_angle;
  G4double Ebeam = dvcsGlobals::Ebeam;
  fVertex_z = vertex_z;//[cm]
  //Is it okay to set seed here?
  TRandom2 rand2;
  rand2.SetSeed(0);// "0" automatically generate seed values which are different every time by using TRandom3
  fSmeared_vertex_z = fVertex_z + rand2.Gaus(0, 1.0/sin(HMS_angle))*0.1; // 0.1 is for mm->cm //Hall C HMS Target Vertex Recoonstruction Accuracy : 1.0 mm

  fHistoManager->SetVertexz(fSmeared_vertex_z, fVertex_z);

  G4double beam_px = 0.;
  G4double beam_py = 0.;
  G4double beam_pz = Ebeam;

  fHistoManager->SetInitElecRec(beam_px, beam_py, beam_pz);
  fHistoManager->SetInitElecGen(fPx[0], fPy[0], fPz[0]);
  fHistoManager->SetScatElecGen(fPx[1], fPy[1], fPz[1]);
  fHistoManager->SetPhotGen(fPx[2], fPy[2], fPz[2]);

  //===========Calculate Q2, t, phi from DVCS event generator=========

  TVector3 mom_init_elec_gen(0, 0, fE[0]);//? why not (fPx[0], fPy[0], fPz[0])?

  TVector3 mom_scat_elec_gen(fPx[1], fPy[1], fPz[1]);
  TVector3 mom_final_phot_gen(fPx[2], fPy[2], fPz[2]);

  G4double Q2_gen = 2*mom_init_elec_gen.Mag()*mom_scat_elec_gen.Mag()*(1 - cos(mom_init_elec_gen.Angle(mom_scat_elec_gen)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))
  TVector3 mom_virt_phot_gen = mom_init_elec_gen - mom_scat_elec_gen;

  TVector3 v1_gen = mom_virt_phot_gen.Cross(mom_scat_elec_gen);
  TVector3 v2_gen = mom_virt_phot_gen.Cross(mom_final_phot_gen);
  G4double phi_gen = v1_gen.Angle(v2_gen);
  if(mom_virt_phot_gen.Dot(v1_gen.Cross(v2_gen))<0) phi_gen = 2.*TMath::Pi()-phi_gen;

  G4double cos_gen = (mom_virt_phot_gen.Dot(mom_final_phot_gen))/(mom_virt_phot_gen.Mag()*mom_final_phot_gen.Mag());
  G4double nu_gen = fE[0]-mom_scat_elec_gen.Mag();
  G4double xB_gen = Q2_gen/(2*M_targ*nu_gen);
  //G4double tM_gen=(Q2_gen*M_targ+2.*nu_gen*M_targ*(nu_gen-sqrt(nu_gen*nu_gen+Q2_gen)*cos_gen))/(sqrt(nu_gen*nu_gen+Q2_gen)*cos_gen-nu_gen-M_targ);
  G4double tM_gen = 2*M_targ*(M_targ - fE[3]);

  fHistoManager->SetKineGen(tM_gen, xB_gen, Q2_gen, phi_gen);


  //============Calculate R_value from DVCS event generator===========

  //  G4double HMS_angle= dvcsGlobals::HMS_angle;  //23.91*TMath::Pi()/180.; // KINEMATICS 3

  G4double mom_central= dvcsGlobals::HMS_momentum;

  //HMS_angle
  G4int run_number = dvcsGlobals::run_number;
  G4double ry = -0.01*fVertex_z*TMath::Sin(TMath::ATan2(fPx[1],fPz[1])); // 0.01 to convert cm to meters
  G4double rdp= (sqrt(fPx[1]*fPx[1]+fPy[1]*fPy[1]+fPz[1]*fPz[1])-mom_central)/mom_central;
  G4double rtheta = TMath::ATan2(-fPy[1],TMath::Sqrt(fPx[1]*fPx[1]+fPz[1]*fPz[1]));
  G4double rphi   = TMath::ATan2(fPx[1],fPz[1]) - HMS_angle;

  //G4double r_val_gen = dvcs_event->rfunction(ry,rdp,rtheta,rphi);
  G4double r_val_gen = RFunction(run_number, rtheta, rdp, rphi, ry);//Alexa's RFunction
  fHistoManager->SetRvalGen(r_val_gen);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::DefineWeights(G4double PSF, G4double X_sum, G4double X_diff, G4double X_BH)
{

  fPSF = PSF;
  fX_sum = X_sum;
  fX_diff = X_diff;
  fX_BH = X_BH;

  fHistoManager->SetWeights(fPSF, fX_sum, fX_diff, fX_BH);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::DefineHMS_elec(G4ThreeVector mom, G4double ene)
{
  fPx[4] = mom.getX();
  fPy[4] = mom.getY();
  fPz[4] = mom.getZ();
  fE[4]  = ene;

  fHistoManager->SetScatElecRec(fPx[4], fPy[4], fPz[4]);
  hit_HMS = true;

  //============Calculate R_value from Detected electron on HMS===========

  G4double HMS_angle = dvcsGlobals::HMS_angle;
  G4double pcentral= dvcsGlobals::HMS_momentum; // KINEMATICS 3

  //HMS_angle
  G4int run_number = dvcsGlobals::run_number;
  G4double ry = -0.01*fSmeared_vertex_z*TMath::Sin(TMath::ATan2(fPx[4],fPz[4]));  // 0.01 to convert cm to meters
  G4double rdp= (sqrt(fPx[4]*fPx[4]+fPy[4]*fPy[4]+fPz[4]*fPz[4])-pcentral)/pcentral;
  G4double rtheta = TMath::ATan2(-fPy[4],TMath::Sqrt(fPx[4]*fPx[4]+fPz[4]*fPz[4]));
  G4double rphi   = TMath::ATan2(fPx[4],fPz[4]) - HMS_angle;

  //G4double r_val_rec = dvcs_event->rfunction(ry,rdp,rtheta,rphi);
  G4double r_val_rec = RFunction(run_number,rtheta,rdp,rphi,ry);//Alexa's RFunction
  fHistoManager->SetRvalRec(r_val_rec);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int EventAction::GetEventNb()
{
  return fEvtNb;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
