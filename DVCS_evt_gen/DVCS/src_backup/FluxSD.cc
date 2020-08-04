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
// $Id: B5HadCalorimeterSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file B5HadCalorimeterSD.cc
/// \brief Implementation of the B5HadCalorimeterSD class

#include "FluxSD.hh"
#include "FluxHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//20180814(start)
#include "G4VProcess.hh"
//20180814(finish)
//#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FluxSD::FluxSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    collectionName.insert("FluxColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FluxSD::~FluxSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FluxSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new FluxHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
  hce->AddHitsCollection(fHCID,fHitsCollection);
    
  // fill calorimeter hits with zero energy deposition
  FluxHit* hit = new FluxHit();
  fHitsCollection->insert(hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FluxSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //  G4double edep = step->GetTotalEnergyDeposit();
  //  if (edep==0.) return true;
    
  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* volumePhysical = touchable->GetVolume(0);
  G4int hitID = volumePhysical->GetCopyNo();
  FluxHit* hit = (*fHitsCollection)[hitID];

  //20171124(start)
  G4Track* track = step->GetTrack();
  G4int PID = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
  //20171124(finish)
  //20171127(start)
  //20180713(start)
  //  G4double energy = track->GetKineticEnergy();
  G4double energy = track->GetTotalEnergy();
  //20180713(finish)
  G4ThreeVector position = track->GetPosition();

  //  G4cout<<"(x, y, z) : "<<"("<<position.x()<<", "<<position.y()<<", "<<position.z()<<")"<<G4endl;
  //  G4cout<<"(x, y, z) : "<<position<<G4endl;
  //20171127(finish)

  //20180814(start)
  //record the particle only once(record when the process is "Transportation")
  //Record if the process is "Transportation" and if the volume remains the same in the next step. See the code lines below:
  //Turned out it is not usefull when the material is pseudo(and vacuum).
  //No interaction happens

  // G4StepPoint* preStepPoint = step->GetPreStepPoint();
  // const G4VProcess* CurrentProcess=preStepPoint->GetProcessDefinedStep();
  // if (CurrentProcess != 0) {
  //   const G4String & StepProcessName = CurrentProcess->GetProcessName();
  //   G4String volumePos = track->GetNextVolume()->GetName();
  //   //    if(volumePos!="World_pos")G4cout<<volumePos<<G4endl;
  //   if(StepProcessName== "Transportation" && volumePos == "Flux_pos") {
  //     // processing hit when entering the volume
  //     // set energy of the particle
  //     //      G4cout<<volumePos<<G4endl;
  //     hit->SetEdep(energy);
  //     //  G4double testing = hit->GetEdep();
  //     //  G4cout<<"testing : "<<testing<<G4endl;
  //     //20171124(start)
  //     hit->SetPID(PID);
  //     //20171124(finish)    
  //     //20171127(start)
  //     hit->SetPos(position);
  //     //20171127(finish)
  //   }
  // }
  // processing hit when entering the volume
  // set energy of the particle
  //      G4cout<<volumePos<<G4endl;
  hit->SetEdep(energy);
  //  G4double testing = hit->GetEdep();
  //  G4cout<<"testing : "<<testing<<G4endl;
  //20171124(start)
  hit->SetPID(PID);
  //20171124(finish)    
  //20171127(start)
  hit->SetPos(position);
  //20171127(finish)

  //And kill the particle after detection
  //Does not slow down the simulation significantly
  //  track->SetTrackStatus(fStopAndKill);
  //20180814(finish)

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
