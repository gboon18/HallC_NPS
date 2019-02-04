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
/// \file analysis/AnaEx02/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 92322 2015-08-27 14:54:05Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TFile;
class TTree;
class TH1D;

const G4int MaxHisto = 4;
const G4int MaxNtuple = 1080;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
public:
  HistoManager();
  ~HistoManager();
   
  void Book(G4String);
  void FillNtuple();
  void Save();
  void SetEnergyandPID(G4int , G4int, G4double, G4int, G4int, G4int, G4int, G4int);//20171017/20180117
  void SetOP(G4int , G4int );//20171017
  void PrintStatistic();
        
private:
  TFile*   fRootFile;
  TTree*   fNtuple;
  TTree*   fNtuple2;
  TTree*   fNtuple_gamma;
  TTree*   fNtuple_pos;
  TTree*   fNtuple_elec;
  TTree*   fNtuple_unknown;
  G4double fEdep[MaxNtuple];
  G4int    fPID[MaxNtuple];
  G4int    fOP_sc[MaxNtuple];
  G4int    fOP_ce[MaxNtuple];
  G4int    fOP_cover[MaxNtuple];
  G4int    fOP_frontcover[MaxNtuple];
  G4int    fOP_pmtcover[MaxNtuple];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
