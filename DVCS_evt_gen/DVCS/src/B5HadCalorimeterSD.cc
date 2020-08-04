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

#include "B5HadCalorimeterSD.hh"
#include "B5HadCalorimeterHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5HadCalorimeterSD::B5HadCalorimeterSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    collectionName.insert("HadCalorimeterColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5HadCalorimeterSD::~B5HadCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5HadCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection 
      = new B5HadCalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
    if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    hce->AddHitsCollection(fHCID,fHitsCollection);
    
    // fill calorimeter hits with zero energy deposition
    for (G4int iColumn=0 ; iColumn<36 ; iColumn++)
        for (G4int iRow=0 ; iRow<30 ; iRow++)
        {
            B5HadCalorimeterHit* hit = new B5HadCalorimeterHit();
            fHitsCollection->insert(hit);
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5HadCalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
    
  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  //changed crystal numbering. column by column!
  G4VPhysicalVolume* cellPhysical = touchable->GetVolume(2);
  G4int rowNo = cellPhysical->GetCopyNo();//0~29(30) in total
  G4VPhysicalVolume* columnPhysical = touchable->GetVolume(1);
  G4int columnNo = columnPhysical->GetCopyNo();//0~35(36) in total
  G4int hitID = columnNo+36*rowNo;//0~1079(1080)

  B5HadCalorimeterHit* hit = (*fHitsCollection)[hitID];

  G4int sc = 0;
  G4int ce = 0;
  const std::vector<const G4Track*>* secondaries =
    step->GetSecondaryInCurrentStep();
  if (secondaries->size()>0) {
    for(unsigned int i=0; i<secondaries->size(); ++i) {
      if (secondaries->at(i)->GetParentID()>0) {
	if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
	   == G4OpticalPhoton::OpticalPhotonDefinition()){
	  if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
	      == "Scintillation")sc++;
	  else if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
		   == "Cerenkov")ce++;// step->GetTrack()->SetTrackStatus(fStopAndKill);}
	}
      }
    }
  }
  G4Track* track = step->GetTrack();
  G4int PID = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();

  // check if it is first touch (for storing transformation and rotation matrix of the cell for the sake of drawing the hit)
  if (hit->GetColumnID()<0)
    {
      hit->SetColumnID(columnNo);
      hit->SetRowID(rowNo);
      G4int depth = touchable->GetHistory()->GetDepth();
      G4AffineTransform transform 
	= touchable->GetHistory()->GetTransform(depth-2);
      transform.Invert();
      hit->SetRot(transform.NetRotation());
      hit->SetPos(transform.NetTranslation());
    }
  // add energy deposition
  hit->AddEdep(edep);
  hit->SetPID(PID);
  hit->AddOPInt_sc(sc);
  hit->AddOPInt_ce(ce);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
