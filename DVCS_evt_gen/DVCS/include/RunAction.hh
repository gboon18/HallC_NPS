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
/// \file analysis/shared/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//
// $Id: RunAction.hh 67226 2013-02-08 12:07:18Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class HistoManager;

class RunAction : public G4UserRunAction
{
public:

  RunAction(HistoManager*, G4String, G4int, long, long);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);
    
  void fillPerEvent(G4double, G4double, G4double, G4double); 

  void AddEdepPerEvent(G4int, G4int, G4double, G4int, G4int, G4int, G4int, G4int);

private:
  HistoManager* fHistoManager;

  G4String fFileName;

  G4double fSumEAbs, fSum2EAbs;
  G4double fSumEGap, fSum2EGap;
    
  G4double fSumLAbs, fSum2LAbs;
  G4double fSumLGap, fSum2LGap;    

  G4double fEdep[1080];
  G4int    fPID[1080];//this is dummy. Too lazy to change HistoManager::SetEnergyandPID
  G4int    fOP_sc[1080];
  G4int    fOP_ce[1080];
  G4int    fOP_cover[1080];
  G4int    fOP_frontcover[1080];
  G4int    fOP_pmtcover[1080];

  G4int    fIndex, fIndex2;

  long fSeed1, fSeed2, fSeed3;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

