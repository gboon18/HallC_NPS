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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(HistoManager* histo)
  :G4UserEventAction(),
   fHistoManager(histo),

   fEvtNb(0),

   fHCHCID(-1),fHadCalEdep(),

   fPID(),

   fOP_sc(), fOP_ce(),
   fCrystCoverHCID(-1), fCrystCoverOP(), fCrystFrontCoverHCID(-1), fCrystFrontCoverOP(), fPMTcoverHCID(-1), fPMTcoverOP(),

  fPrintModulo(0)                             
{
  fPrintModulo = 100; 

  fHadCalEdep.resize(9,0.);
  fPID.resize(9,0),

    fOP_sc.resize(9,0), fOP_ce.resize(9,0),
    fCrystCoverOP.resize(9,0), fCrystFrontCoverOP.resize(9,0), fPMTcoverOP.resize(9,0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  fEvtNb = evt->GetEventID();
  if (fEvtNb%fPrintModulo == 0) 
    //    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
 
    if (fHCHCID==-1) {

      //cc-in2p3 does not understand(or geant4.10.02.0 does not understand) "auto"
      //    auto sdManager = G4SDManager::GetSDMpointer();
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      fHCHCID = sdManager->GetCollectionID("HadCalorimeter/HadCalorimeterColl");
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* evt)
{

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
  for (G4int i=0;i<9;i++)
    {
      B5HadCalorimeterHit* hit = (*hcHC)[i];
      G4double eDep = hit->GetEdep();//total energy deposition of each crystals

      fHadCalEdep[i] = eDep;

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

      fHistoManager->SetEnergy( i, fHadCalEdep[i], fOP_sc[i], fOP_ce[i], fCrystCoverOP[i], fCrystFrontCoverOP[i], fPMTcoverOP[i]);
    }
  fHistoManager->FillNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int EventAction::GetEventNb()
{
  return fEvtNb;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
