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
/// \file analysis/shared/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
// $Id: RunAction.cc 92322 2015-08-27 14:54:05Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(HistoManager* histo, G4String fileName, G4int n, long i, long j)
  : G4UserRunAction(),
    fHistoManager(histo),
    fFileName(fileName),
    fIndex(n), fIndex2(n),
    fSeed1(i), fSeed2(j), fSeed3(i)
{
  fEdep[1080] = {0.};
  fOP_sc[1080] = {0};
  fOP_ce[1080] = {0};
  fOP_cover[1080] = {0};
  fOP_frontcover[1080] = {0};
  fOP_pmtcover[1080] = {0};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4int index = fIndex2;
  long seeds[2];
  seeds[0] = fSeed3;
  seeds[1] = fSeed2;
  G4Random::setTheSeeds(seeds, index);
  G4Random::showEngineStatus();

G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    
  //histograms
  fHistoManager->Book(fFileName); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdepPerEvent(G4int id, G4double edep, G4int sc, G4int ce, G4int opc, G4int opfc, G4int oppc)
{
  fEdep[id] += edep;
  fOP_sc[id] += sc;
  fOP_ce[id] += ce;
  fOP_cover[id] += opc;
  fOP_frontcover[id] += opfc;
  fOP_pmtcover[id] += oppc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

  for(G4int i  = 0 ; i < 1080 ; i++){
    fHistoManager->SetEnergy(i, fEdep[i], fOP_sc[i], fOP_ce[i], fOP_cover[i], fOP_frontcover[i], fOP_pmtcover[i]);
  }
  fHistoManager->FillNtuple();
  
  //save histograms
  //
  fHistoManager->Save();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
