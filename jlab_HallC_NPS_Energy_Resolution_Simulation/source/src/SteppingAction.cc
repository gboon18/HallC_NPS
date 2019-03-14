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
#include "HistoManager.hh"

#include "G4Step.hh"
#include "G4TouchableHistory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
				 EventAction* evt,
				 HistoManager* histo)
: G4UserSteppingAction(), 
  fDetector(det), fEventAction(evt),
  fHistoManager(histo)                                         
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // // get volume of the current step
  // G4VPhysicalVolume* volume 
  // = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // // collect energy and track length step by step
  // G4double edep = aStep->GetTotalEnergyDeposit();
  
  // G4double stepl = 0.;
  // if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
  //   stepl = aStep->GetStepLength();

  G4VPhysicalVolume* volume_pre 
    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4double edep = aStep->GetTotalEnergyDeposit();

  if(
     volume_pre->GetLogicalVolume()->GetName() == "Crystal_log"
     ){

    G4TouchableHistory* touchable
      = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

    G4VPhysicalVolume* cellPhysical = touchable->GetVolume(2);
    G4int rowNo = cellPhysical->GetCopyNo();//0~29(30) in total
    G4VPhysicalVolume* columnPhysical = touchable->GetVolume(1);
    G4int columnNo = columnPhysical->GetCopyNo();//0~35(36) in total
    G4int hitID = columnNo+36*rowNo;//0~1079(1080)
    G4Track* track = aStep->GetTrack();
    G4ThreeVector position = track->GetPosition();
    const G4AffineTransform transformation = aStep->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform();
    G4ThreeVector localPosition = transformation.TransformPoint(position);

    G4double eDep = aStep->GetTotalEnergyDeposit();

    G4int evtNb = fEventAction->GetEventNb();
    fHistoManager->SetFluxEnergy(evtNb, hitID, eDep, localPosition);
    fHistoManager->FillNtuple_Flux();

  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
