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
/// \file analysis/AnaEx02/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 92374 2015-08-31 08:52:09Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fRootFile(0),
   fNtuple(0)
{
  fEdep[MaxNtuple] = {0.};
  fPID[MaxNtuple] = {0};
  fOP_sc[MaxNtuple] = {0};
  fOP_ce[MaxNtuple] = {0};
  fOP_cover[MaxNtuple] = {0};
  fOP_frontcover[MaxNtuple] = {0};
  fOP_pmtcover[MaxNtuple] = {0};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  if (fRootFile) delete fRootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book(G4String fileName)
{ 
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  // G4String fileName = "AnaEx02_20171130_mapmagon_2.root";
  fRootFile = new TFile(fileName,"RECREATE");
  if (! fRootFile) {
    G4cout << " HistoManager::Book :" 
           << " problem creating the ROOT TFile "
           << G4endl;
    return;
  }
  
  fNtuple = new TTree("t","Energy deposition and OP in crystas");
  fNtuple->Branch("edep", fEdep, "energy_deposition[1080]/D");//20171017(changed from"energy_deposition[1080]/F"<--does not work!!!)
  fNtuple->Branch("sc", fOP_sc, "scintillated OP[1080]/I");
  fNtuple->Branch("ce", fOP_ce, "cerenkov OP[1080]/I");
  fNtuple->Branch("op_cover", fOP_cover, "OP on the side of the crystal wrapper[1080]/I");
  fNtuple->Branch("op_frontcover", fOP_frontcover, "OP on the front side of the crystal wrapper[1080]/I"); 
  fNtuple->Branch("op_pc", fOP_pmtcover, "OP arrived at the pmt cover[1080]/I");
  G4cout << "\n----> Output file is open in " << fileName << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::FillNtuple()
{
  fNtuple->Fill();//20171127(temp)
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{ 
  if (! fRootFile) return;
  fRootFile->Write();       // Writing the histograms to the file
  fRootFile->Close();       // and closing the tree (and the file)
  
  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetEnergyandPID(G4int id, G4int PID, G4double edep, G4int sc, G4int ce, G4int opc, G4int opfc, G4int oppc)
{
  if(edep) fEdep[id] = edep;
  fPID[id] = PID;
  fOP_sc[id] = sc;
  fOP_ce[id] = ce;
  fOP_cover[id] = opc;
  fOP_frontcover[id] = opfc;
  fOP_pmtcover[id] = oppc;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
