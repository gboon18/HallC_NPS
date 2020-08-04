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
/// \file electromagnetic/TestEm4/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//20190225(start)
#include <math.h>//arccosine(acos), pi(M_PI)
//20190225(finish)

#include "PrimaryGeneratorAction.hh"

//20190417(start)
#include "EventAction.hh"
//20190417(finish)

//20190410(start)
#include "dvcsGlobals.hh"
//20190410(finish)

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//20190225(start)
//20171006(start)
// #include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"//20190225. Back to using this.
//20171006(finish)
//20190225(finish)

//20190225(start)
#include "G4ParticleTable.hh"
//20190225(finish)

//20190626(start)
//test(start)
#include <fstream>
#include <iostream>
//test(finish)
//20190626(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//20190410(start)
// PrimaryGeneratorAction::PrimaryGeneratorAction()
//20190418(start)
// PrimaryGeneratorAction::PrimaryGeneratorAction(long seed1, long seed2)
PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* evt, long seed1, long seed2)
//20190418(finish)
//20190410(finish)
  : G4VUserPrimaryGeneratorAction(),
    //20190418(start)
    fEventAction(evt),
    //20190418(finish)
    //20190419(start)
    M_targ(-999.),
    //20190419(finish)
    fParticleGun(0),
    //20190417(start)
    L_elec_init(0), L_elec_scat(0), L_phot_final(0), L_prot_final(0)
    //20190417(finish)
{
  //20190410(start)
  //let's use dvcsGlobals with [mm], [rad] and [GeV] units
  //So....delete this one
  // //20190212(start)
  // gEv=new TGenDVCS(11,0,40000,50000);//(beam-energy [GeV], target type, seed1, seed2)
  // gEv->SetTargetParam(15.,0.0723);//(target-length [cm?], target density)
  // // gEv->SetGeometry(acos(1-(3*3/2*11*3)), 3.0, 6.3*(M_PI/180.), 6000.);//(spectrometer angle [rad](cos^-1(1-Q^2/2*k*k')), spec-central momentum [GeV], calo-angle [rad], calo-face distance from targ [mm?])
  // gEv->SetGeometry(-0.5283572643, 3.0, -0.2617993878, 4000.);//(spectrometer angle [rad](cos^-1(1-Q^2/2*k*k')), spec-central momentum [GeV], calo-angle [rad], calo-face distance from targ [mm?])
  // //20190212(finish)
  //20190410(finish)

  //20190410(start)
  //let's use dvcsGlobals with [mm], [rad] and [GeV] units
  G4double BeamEnergy = dvcsGlobals::Ebeam;
  //20190412(start)
  //test. temp
  // G4double HMSangle = -1*dvcsGlobals::HMS_angle;//the reason -1 is there, is to move the HMS to the other side of the beam-line
  G4double HMSangle = dvcsGlobals::HMS_angle;
  //20190412(finish)
  G4double HMSmomentum = dvcsGlobals::HMS_momentum;
  G4double Calo_distance = 10.*dvcsGlobals::Calo_distance;//[cm] to [mm]
  //20190412(start)
  //test. temp
  // G4double Calo_angle = -1*dvcsGlobals::Calo_angle;//the reason -1 is there, is to move the calorimeter to the other side of the beam-line
  G4double Calo_angle = dvcsGlobals::Calo_angle;
  //20190412(finish)
  target_gen_proc_type = dvcsGlobals::target_gen_proc_type;
  G4double target_density = dvcsGlobals::target_density;
  G4double target_length = dvcsGlobals::target_length;

  G4cout<<"PrimaryGeneratorAction Calo_angle : "<<Calo_angle*TMath::RadToDeg()<<G4endl;

  if( target_gen_proc_type == 0 || target_gen_proc_type == 1 || target_gen_proc_type ==2)
    {
      if(target_gen_proc_type == 0) {M_targ = 0.9383;}
      else if( target_gen_proc_type == 1 ) {M_targ = 0.939565378;}
      else if( target_gen_proc_type == 2 ) {M_targ = 1.875612859;}
    }
  else
    {
      cout<<"Primary Generator Action!! Error: wrong target mass "<<endl;
    }

  gEv=new TGenDVCS(BeamEnergy, target_gen_proc_type, seed1, seed2);//(beam-energy [GeV], target type, seed1, seed2)//But I do not think the seed1 and seed2 is actually used in TGenDVCS.
  gEv->SetTargetParam(target_length, target_density);//(target-length [cm][DVCSGen takes cm to calculate bremstraulung], target density[g/cm3][geant4 default])

  gEv->SetGeometry(HMSangle, HMSmomentum, Calo_angle, Calo_distance);//(spectrometer angle [rad](cos^-1(1-Q^2/2*k*k')), spec-central momentum [GeV], calo-angle [rad], calo-face distance from targ [mm in TDVCSgen?], [mm in simulation and NPS_SOFT])
  //20190410(finish)


  //20190225(start)
  G4int n_particle = 1;
  // fParticleGun  = new G4GeneralParticleSource(/*n_particle*/);//20171006 changed from G4ParticleGun
  fParticleGun  = new G4ParticleGun(n_particle);//20190225. Back to using this
  //20190225(finish)
  /*
    fParticleGun->SetParticleEnergy(0*eV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //20190417(start)
  dvcsGlobals::hit_HMS_CALO_flag = false;
  G4bool hit_HMS = false;
  //20190417(finish)

  //20190412(start)
  //test. temp
  // G4cout<<"Event Starts : "<<anEvent->GetEventID()<<G4endl;
  //20190225(start)
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  //20190225(finish)

  //20190212(start)
  gEv->GenerateVertex();//Generates the vertex //[cm]
  gEv->ExtBrem();//Make external pre-vertex radiative corrections
  gEv->GenKin();//Computes scattered electron kinematics
  gEv->IntRCBef();//Internal radiative corrections (before vertex)
  //;//Compute the scattered electron 4-vector

  if(gEv->ComputeElectron()){//Check hor spectro accep.

    gEv->IntRCAft();//Internal radiative corrections (after vertex)

    TLorentzVector L_b(0, 0, dvcsGlobals::Ebeam, dvcsGlobals::Ebeam);
    TLorentzVector L_scat_el = *gEv->GetScatteredElectron();
    TLorentzVector L_init_phot = L_b - L_scat_el;
    G4double nu = L_init_phot.E();

    TVector3 k_v = TVector3(L_b.Vect());
    TVector3 kp_v = TVector3(L_scat_el.Vect());

    G4double Q2 = 2*k_v.Mag()*kp_v.Mag()*(1 - cos(k_v.Angle(kp_v)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))                                                                              

    //cout<<"Q2 diff = "<<(L_init_phot.M2() + Q2)<<endl;                                                                                                                           
    G4double xB = Q2/(2*M_targ*nu);
    G4double q3=TMath::Sqrt(Q2+TMath::Power(nu,2.));
    G4double q0primemax=0.5*Q2*(1.-xB)/(xB*(M_targ+nu-q3));
    G4double q0primemin=0.5*Q2*(1.-xB)/(xB*(M_targ+nu+q3));

    G4double tmax=-Q2-2.*q0primemax*(nu-q3);
    G4double tmin=-Q2-2.*q0primemin*(nu+q3);

    // G4cout<<"tmax : "<<tmax<<G4endl;
    // G4cout<<"tmin : "<<tmin<<G4endl;

    gEv->Settmin(tmax - 2.); //Use this method to change tmin (default -2 GeV)                                                                                                     
    gEv->Settmax(0.);  //Use this method to change tmax (default 0 GeV) 

    gEv->ComputeDVCS();//Compute the gamma*p->gamma p collision

    // //20190412(start)
    // gEv->Print("all");
    // //20190412(finish)

    gEv->ApplySpecVerAcc();//Rotates all final vectors around the beam axis

    //20190417(start)
    hit_HMS = false;
    //20190417(finish)

    if(gEv->HitsSpectro(gEv->GetScatteredElectron()) && gEv->HitsCalo(gEv->GetFinalPhoton())){//Checks calo acceptance

      //20190417(start)
      hit_HMS = true;
      //20190417(finish)

      //20190417(start)
      L_elec_init  = gEv->GetInitialElectron();
      L_elec_scat  = gEv->GetScatteredElectron();
      L_phot_final = gEv->GetFinalPhoton(); 
      L_prot_final = gEv->GetFinalProton(); 

      //test(start)
      // G4cout<<"scat elec angle : "<<TMath::ATan2(TMath::Sqrt(L_elec_scat->X()*L_elec_scat->X() + L_elec_scat->Y()*L_elec_scat->Y()), L_elec_scat->Z())*TMath::RadToDeg()<<G4endl;
      // ofstream outfile("hms.txt", ios::app);
      // outfile<<TMath::ATan2(TMath::Sqrt(L_elec_scat->X()*L_elec_scat->X() + L_elec_scat->Y()*L_elec_scat->Y()), L_elec_scat->Z())*TMath::RadToDeg()<<G4endl;
      // outfile.close();
      // ofstream outfile("nps.txt", ios::app);
      // outfile<<TMath::ATan2(TMath::Sqrt(L_phot_final->X()*L_phot_final->X() + L_phot_final->Y()*L_phot_final->Y()), L_phot_final->Z())*TMath::RadToDeg()<<G4endl;
      // outfile.close();
      //test(finish)
      if(hit_HMS){
	dvcsGlobals::hit_HMS_CALO_flag = true;
	// cout<<"Sum = "<<gEv->XSecSum()<<endl;                                                                                                                               
	// cout<<"Diff = "<<gEv->XSecDif()<<endl;                                                                                                                              
	// cout<<"PSF = "<<gEv->GetPSF()<<endl;                                                                                                                                

	if(target_gen_proc_type == 0)//proton
	  {
	    G4double PSF    = gEv->GetPSF();
	    G4double X_sum  = gEv->XSecSum(0);
	    G4double X_diff = gEv->XSecDif();
	    G4double X_BH   = gEv->XSecSum(1);
	    // G4cout << "Pri_Gen_Act, " << "PSF = " << gEv->GetPSF() << " ; crs_sum = " << gEv->XSecSum(0) << " ; crs_diff = " << gEv->XSecDif() << " ; crs_BH = " << gEv->XSecSum(1) << G4endl; //debug 
	    fEventAction->DefineWeights(PSF, X_sum, X_diff, X_BH);
	    // fEventAction->DefineWeights(gEv->GetPSF(), gEv->XSecSum(0), gEv->XSecDif(), gEv->XSecSum(1));
	  }
	else if( target_gen_proc_type == 1 || target_gen_proc_type == 2 )//neutron or deutron
	  {
	    fEventAction->DefineWeights(gEv->GetPSF(), 0., 0., 0.);
	    //cout << "crs_sum = plop ; crs_BH = plop" << endl; //debug                                                                                                        
	  }
	else
	  {
	    G4cout<<"Target types is not set correctly, "<<target_gen_proc_type<<G4endl;
	    G4cout<<"The program is exiting"<<G4endl;
	    exit(1);
	  }

	fEventAction->DefinePrimaries(L_elec_init, L_elec_scat, L_phot_final, L_prot_final, gEv->GetVertex()->Z());
	//20190417(finish)

	//20190408(start)
	// G4cout<<"Photon generation"<<G4endl;
	//20190408(finish)

	//20190410(start)
	//unit change [cm] to [mm]
	// 20190225(start)
	fParticleGun->SetParticleDefinition(particleTable->FindParticle(11));//11 : electron
	fParticleGun->SetParticleEnergy(gEv->GetScatteredElectron()->E()*CLHEP::GeV - 0.510999*CLHEP::MeV);
	fParticleGun->SetParticlePosition(G4ThreeVector(gEv->GetVertex()->X()*CLHEP::cm, gEv->GetVertex()->Y()*CLHEP::cm, gEv->GetVertex()->Z()*CLHEP::cm));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector((gEv->GetScatteredElectron()->Vect()).Px(), gEv->GetScatteredElectron()->Vect().Py(), gEv->GetScatteredElectron()->Vect().Pz()));

	fParticleGun->GeneratePrimaryVertex(anEvent);

	fParticleGun->SetParticleDefinition(particleTable->FindParticle(22));//22 : photon
	fParticleGun->SetParticleEnergy(gEv->GetFinalPhoton()->E()*CLHEP::GeV);
	fParticleGun->SetParticlePosition(G4ThreeVector(gEv->GetVertex()->X()*CLHEP::cm, gEv->GetVertex()->Y()*CLHEP::cm, gEv->GetVertex()->Z()*CLHEP::cm));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector((gEv->GetFinalPhoton()->Vect()).Px(), gEv->GetFinalPhoton()->Vect().Py(), gEv->GetFinalPhoton()->Vect().Pz()));

	fParticleGun->GeneratePrimaryVertex(anEvent);

	// fParticleGun->SetParticleDefinition(particleTable->FindParticle(11));//11 : electron
	// fParticleGun->SetParticleEnergy(gEv->GetScatteredElectron()->E()*CLHEP::GeV - 0.510999*CLHEP::MeV);
	// fParticleGun->SetParticlePosition(G4ThreeVector(gEv->GetVertex()->X(), gEv->GetVertex()->Y(), gEv->GetVertex()->Z()));
	// fParticleGun->SetParticleMomentumDirection(G4ThreeVector((gEv->GetScatteredElectron()->Vect()).Px(), gEv->GetScatteredElectron()->Vect().Py(), gEv->GetScatteredElectron()->Vect().Pz()));

	// fParticleGun->GeneratePrimaryVertex(anEvent);

	// fParticleGun->SetParticleDefinition(particleTable->FindParticle(22));//22 : photon
	// fParticleGun->SetParticleEnergy(gEv->GetFinalPhoton()->E()*CLHEP::GeV);
	// fParticleGun->SetParticlePosition(G4ThreeVector(gEv->GetVertex()->X(), gEv->GetVertex()->Y(), gEv->GetVertex()->Z()));
	// fParticleGun->SetParticleMomentumDirection(G4ThreeVector((gEv->GetFinalPhoton()->Vect()).Px(), gEv->GetFinalPhoton()->Vect().Py(), gEv->GetFinalPhoton()->Vect().Pz()));

	// fParticleGun->GeneratePrimaryVertex(anEvent);
	//20190225(finish)
	//20190410(finish)
	// }
      }
    }
  }
  gEv->Clear();
  //20190212(finish)

  //20190412(start)
  //test. temp
  // G4cout<<"Event Ends : "<<anEvent->GetEventID()<<G4endl;
  //20190225(start)

  // fParticleGun->SetParticleDefinition(particleTable->FindParticle(22));
  // fParticleGun->SetParticleEnergy(7*CLHEP::GeV);
  // fParticleGun->SetParticlePosition(G4ThreeVector((3*sin(15*M_PI/180.))*CLHEP::m  + 10*CLHEP::mm, 10.*CLHEP::mm, (3*cos(15*M_PI/180.))*CLHEP::m));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector((sin(15*M_PI/180.)), 0., (cos(15*M_PI/180.))));
  // fParticleGun->GeneratePrimaryVertex(anEvent);

  // fParticleGun->SetParticleDefinition(particleTable->FindParticle(22));
  // fParticleGun->SetParticleEnergy(7*CLHEP::GeV);
  // fParticleGun->SetParticlePosition(G4ThreeVector((3*sin(15*M_PI/180.))*CLHEP::m  + 150*CLHEP::mm, -150.*CLHEP::mm, (3*cos(15*M_PI/180.))*CLHEP::m));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector((sin(15*M_PI/180.)), 0., (cos(15*M_PI/180.))));
  // fParticleGun->GeneratePrimaryVertex(anEvent);

  /*
    if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
    G4int Z = 10, A = 24;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion
    = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    }    
    //create vertex
    //
    */   
  //20190225(start)
  // fParticleGun->GeneratePrimaryVertex(anEvent);//moved to up
  //20190225(finish)

}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

