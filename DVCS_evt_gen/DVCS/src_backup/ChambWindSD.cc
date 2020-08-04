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

#include "ChambWindSD.hh"
#include "ChambWindHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChambWindSD::ChambWindSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    collectionName.insert("ChambWindColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChambWindSD::~ChambWindSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChambWindSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new ChambWindHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
  hce->AddHitsCollection(fHCID,fHitsCollection);
    
  // fill calorimeter hits with zero energy deposition
  ChambWindHit* hit = new ChambWindHit();
  fHitsCollection->insert(hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ChambWindSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //  G4double edep = step->GetTotalEnergyDeposit();
  //  if (edep==0.) return true;
    
  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* volumePhysical = touchable->GetVolume(0);
  G4int hitID = volumePhysical->GetCopyNo();
  ChambWindHit* hit = (*fHitsCollection)[hitID];

  //20171124(start)
  G4Track* track = step->GetTrack();
  G4int PID = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
  //20171124(finish)
  //20171127(start)
  G4double energy = track->GetKineticEnergy();
  G4ThreeVector position = track->GetPosition();
  //  G4cout<<"(x, y, z) : "<<"("<<position.x()<<", "<<position.y()<<", "<<position.z()<<")"<<G4endl;
  //  G4cout<<"(x, y, z) : "<<position<<G4endl;
  //20171127(finish)
  // set energy of the particle
  hit->SetEdep(energy);
  //  G4double testing = hit->GetEdep();
  //  G4cout<<"testing : "<<testing<<G4endl;
  //20171124(start)
  hit->SetPID(PID);
  //20171124(finish)    
  //20171127(start)
  hit->SetPos(position);
  //20171127(finish)

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
