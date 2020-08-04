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

#include <math.h>//arccosine(acos), pi(M_PI)

#include "PrimaryGeneratorAction.hh"

#include "EventAction.hh"

#include "dvcsGlobals.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4ParticleGun.hh"

#include "G4ParticleTable.hh"

#include <fstream>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* evt, long seed1, long seed2)
  : G4VUserPrimaryGeneratorAction(),
    fEventAction(evt),
    M_targ(-999.),
    fParticleGun(0),
    L_elec_init(0), L_elec_scat(0), L_phot_final(0), L_prot_final(0)
{
  //let's use dvcsGlobals with [mm], [rad] and [GeV] units
  G4double BeamEnergy = dvcsGlobals::Ebeam;
  G4double HMSangle = dvcsGlobals::HMS_angle;
  G4double HMSmomentum = dvcsGlobals::HMS_momentum;
  G4double Calo_distance = 10.*dvcsGlobals::Calo_distance;//[cm] to [mm]
  G4double Calo_angle = dvcsGlobals::Calo_angle;
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

  gEv->SetGeometry(HMSangle, HMSmomentum, Calo_angle, Calo_distance);//(spectrometer angle [rad](cos^-1(1-Q^2/2*k*k')), spec-central momentum [GeV], calo-angle [rad], calo-face distance from targ [mm in TDVCSgen], [mm in simulation and NPS_SOFT])

  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  dvcsGlobals::hit_HMS_CALO_flag = false;
  G4bool hit_HMS = false;
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

  gEv->GenerateVertex();//Generates the vertex //[cm]
  gEv->ExtBrem();//Make external pre-vertex radiative corrections
  gEv->GenKin();//Computes scattered electron kinematics
  gEv->IntRCBef();//Internal radiative corrections (before vertex)

  if(gEv->ComputeElectron()){//Check hor spectro accep.

    gEv->IntRCAft();//Internal radiative corrections (after vertex)

    TLorentzVector L_b(0, 0, dvcsGlobals::Ebeam, dvcsGlobals::Ebeam);
    TLorentzVector L_scat_el = *gEv->GetScatteredElectron();
    TLorentzVector L_init_phot = L_b - L_scat_el;
    G4double nu = L_init_phot.E();

    TVector3 k_v = TVector3(L_b.Vect());
    TVector3 kp_v = TVector3(L_scat_el.Vect());

    G4double Q2 = 2*k_v.Mag()*kp_v.Mag()*(1 - cos(k_v.Angle(kp_v)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))                                                                              

    G4double xB = Q2/(2*M_targ*nu);
    G4double q3=TMath::Sqrt(Q2+TMath::Power(nu,2.));
    G4double q0primemax=0.5*Q2*(1.-xB)/(xB*(M_targ+nu-q3));
    G4double q0primemin=0.5*Q2*(1.-xB)/(xB*(M_targ+nu+q3));

    G4double tmax=-Q2-2.*q0primemax*(nu-q3);
    G4double tmin=-Q2-2.*q0primemin*(nu+q3);

    gEv->Settmin(tmax - 2.); //Use this method to change tmin (default -2 GeV)                                                                                                     
    gEv->Settmax(0.);  //Use this method to change tmax (default 0 GeV) 

    gEv->ComputeDVCS();//Compute the gamma*p->gamma p collision

    gEv->ApplySpecVerAcc();//Rotates all final vectors around the beam axis

    hit_HMS = false;

    if(gEv->HitsSpectro(gEv->GetScatteredElectron()) && gEv->HitsCalo(gEv->GetFinalPhoton())){//Checks calo acceptance

      hit_HMS = true;

      L_elec_init  = gEv->GetInitialElectron();
      L_elec_scat  = gEv->GetScatteredElectron();
      L_phot_final = gEv->GetFinalPhoton(); 
      L_prot_final = gEv->GetFinalProton(); 

      if(hit_HMS){
	dvcsGlobals::hit_HMS_CALO_flag = true;

	if(target_gen_proc_type == 0)//proton
	  {
	    G4double PSF    = gEv->GetPSF();
	    G4double X_sum  = gEv->XSecSum(0);
	    G4double X_diff = gEv->XSecDif();
	    G4double X_BH   = gEv->XSecSum(1);
	    fEventAction->DefineWeights(PSF, X_sum, X_diff, X_BH);
	  }
	else if( target_gen_proc_type == 1 || target_gen_proc_type == 2 )//neutron or deutron
	  {
	    fEventAction->DefineWeights(gEv->GetPSF(), 0., 0., 0.);
	  }
	else
	  {
	    G4cout<<"Target types is not set correctly, "<<target_gen_proc_type<<G4endl;
	    G4cout<<"The program is exiting"<<G4endl;
	    exit(1);
	  }

	fEventAction->DefinePrimaries(L_elec_init, L_elec_scat, L_phot_final, L_prot_final, gEv->GetVertex()->Z());

	//unit change [cm] to [mm]
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

      }
    }
  }
  gEv->Clear();
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

