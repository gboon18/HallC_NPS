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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 67226 2013-02-08 12:07:18Z ihrivnac $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"

//20180621(start)
#include "G4TouchableHistory.hh"
#include "HistoManager.hh"
//20180621(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//20180622(start)
// SteppingAction::SteppingAction(DetectorConstruction* det,
//                                          EventAction* evt)
SteppingAction::SteppingAction(DetectorConstruction* det,
			       EventAction* evt,
			       HistoManager* histo)
  //20180622(finish)
  : G4UserSteppingAction(), 
    fDetector(det), fEventAction(evt),
    //20180622(start)
    fHistoManager(histo), fEne(), fPID(),
    //20180622(finish)
    //20190626(start)
    fTime_opening(-999.), fID_opening(-999)
    //20190626(finish)
{
  //20180622(start)
  fEne.resize(1080,0.), fPID.resize(1080,0);
  //20180622(finish)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  G4VPhysicalVolume* volume_pre 
    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4VPhysicalVolume* volume_post
    = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

  G4Track* track = aStep->GetTrack();

  //20190417(start)//HRS
  //20190626(start)//HMS
  // if(
  //    //needs improvements. Check Frederic's.
  //    //Now, particles from any directions get recorded.
  //    volume_pre->GetLogicalVolume()->GetName() == "HMS_window_log"
  //    && aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary
  //    && track->GetParticleDefinition()->GetPDGCharge() < 0
  //    && track->GetTrackID() == 1 //!!!!!!!!!!!!! Very strong colimator.
  //    ){
  //   G4double energy = track->GetTotalEnergy()/CLHEP::GeV;//MeV to GeV
  //   G4ThreeVector mom = track->GetMomentum()/CLHEP::GeV;//MeV to GeV
  //   fEventAction->DefineHMS_elec(mom, energy);
  // }

  if(
     volume_pre->GetLogicalVolume()->GetName() == "HMS_opening_log"
     && aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary
     && track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding() == 11 //electron
     ){
    fTime_opening = track->GetProperTime();
    fID_opening = track->GetTrackID();
    // G4cout<<"HMS opening : "<<fTime_opening<<", "<<fID_opening<<G4endl;
  }
  if(
     volume_pre->GetLogicalVolume()->GetName() == "HMS_exit_log"
     && aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary
     && track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding() == 11 //electron
     &&track->GetProperTime() > fTime_opening
     &&track->GetTrackID() == fID_opening
     ){
    G4double energy = track->GetTotalEnergy()/CLHEP::GeV;//MeV to GeV
    G4ThreeVector mom = track->GetMomentum()/CLHEP::GeV;//MeV to GeV
    fEventAction->DefineHMS_elec(mom, energy);
    // G4cout<<"HMS exit : "<<track->GetProperTime()<<", "<<track->GetTrackID()<<G4endl;
  }
  //20190626(finish)
  //20190417(finish)

  //20190121(start)
  //All of these are temporary
  if(
     volume_pre->GetLogicalVolume()->GetName() == "FluxOptFiber_log"
     //20190408(start)
     //Now going to use this to compare the # of particles hitting the NPS and the # of particles reconstructed
     && volume_post->GetLogicalVolume()->GetName() == "Mother_log"
     //20190408(finish)
     ){
    // G4cout<<"Photon detection"<<G4endl;
    G4int PID = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
    G4double ene = track->GetTotalEnergy();
    G4ThreeVector position = track->GetPosition();

    //20190424(start)
    G4ThreeVector mom = track->GetMomentum();
    //20190424(finish)

    //20190124(start)//getting a position of particle in local coordinate.
    const G4AffineTransform transformation = aStep->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform();
    G4ThreeVector localPosition = transformation.TransformPoint(position);
    // G4cout<<localPosition.x()<<", "<<localPosition.y()<<", "<<localPosition.z()<<G4endl;
    // G4cout<<position.x()<<", "<<position.y()<<", "<<position.z()<<G4endl;
    //20190124(finish)
    //20190408(start)
    G4int evtNb = fEventAction->GetEventNb();
    //20190424(start)
    // fHistoManager->SetFluxOptFiberEnergyandPID(evtNb, PID, ene, position, localPosition);
    fHistoManager->SetFluxOptFiberEnergyandPID(evtNb, PID, ene, mom, position, localPosition);
    //20190424(finish)
    //20190408(finish)
    fHistoManager->FillNtupleFluxOptFiber(PID);
  }
  //20190121(finish)


  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();

  //20180621(start)
  if(//20190121 : This does not work that well.
     volume_pre->GetLogicalVolume()->GetName() == "WrapFront_log"
     //     && volume_post->GetLogicalVolume()->GetName() == "Crystal_log"
     ){

    G4TouchableHistory* touchable
      = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

    //from B5HadCalorimeterSD.cc
    //crystal numbering. column by column!
    G4VPhysicalVolume* cellPhysical = touchable->GetVolume(2);
    G4int rowNo = cellPhysical->GetCopyNo();//0~29(30) in total
    //    G4cout<<"rowNo :"<<rowNo<<G4endl;
    G4VPhysicalVolume* columnPhysical = touchable->GetVolume(1);
    G4int columnNo = columnPhysical->GetCopyNo();//0~35(36) in total
    //    G4cout<<"columnNo :"<<columnNo<<G4endl;
    G4int hitID = columnNo+36*rowNo;//0~1079(1080)
    // G4cout<<"hitID :"<<hitID<<G4endl;
    
    G4Track* track = aStep->GetTrack();
    G4int PID = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
    G4double ene = track->GetTotalEnergy();
    // G4cout<<ene<<G4endl;
    // G4cout<<"you are here"<<G4endl;

    //20180623(start)
    //    fEventAction->GetFluxEnergy(ene);
    //20180623(finish)

    fHistoManager->SetFluxEnergyandPID(hitID, PID, ene);
    fHistoManager->FillNtuple2();

  }
  //20180621(finish)

  //20180913(start)
  if(
     volume_pre->GetLogicalVolume()->GetName() == "Flux_log"
     // && volume_post->GetLogicalVolume()->GetName() != "Flux_log"
     ){
    G4Track* track = aStep->GetTrack();
    G4int PID = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
    G4double ene = track->GetTotalEnergy();
    G4ThreeVector position = track->GetPosition();
    //20180913(start)//moved from EventAction.cc
    fHistoManager->SetFluxSphereEnergyandPID(PID, ene, position);
    //20190124(start)
    // fHistoManager->FillNtuple3(PID);
    fHistoManager->FillNtupleFluxSphere(PID);//Changed its name.
    //20190124(finish)
    //20180913(finish)
    // G4cout<<PID<<G4endl;
  }
  //20180913(finish)

  //20180927(start)
  //pseudo beam-dump
  //This must be removed after corrector magnet & beam-line shield after sweeping-magnet are introduced.
  //This is to remove the effect of magnetic field of the end of the sweeping-magnet 
  if(
     volume_pre->GetLogicalVolume()->GetName() == "BeamDump_log"
     // || volume_post->GetLogicalVolume()->GetName() == "BeamDump_log"//does not work if this is included.
     ){
    G4Track* track = aStep->GetTrack();
    track->SetTrackStatus(fStopAndKill);
  }
  //20180927(finish)

  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
      
  //20171013(start)
  /*
    if (volume_pre == fDetector->GetAbsorber()) fEventAction->AddAbs(edep,stepl);
    if (volume_pre == fDetector->GetGap())      fEventAction->AddGap(edep,stepl);
  */
  if (volume_pre->GetName() == "Target_log") fEventAction->AddAbs(edep,stepl);
  if (volume_pre->GetName() == "Crystal_log")      fEventAction->AddGap(edep,stepl);
  //20171013(finish)
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
