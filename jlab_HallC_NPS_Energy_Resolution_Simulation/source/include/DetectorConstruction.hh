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
/// \file radioactivedecay/rdecay02/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 66586 2012-12-21 10:48:39Z ihrivnac $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"//to use cm

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
class SensitiveDetector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

public:
  
  virtual G4VPhysicalVolume* Construct();
    
  virtual void ConstructSDandField();

  void SetDetectorGap(G4double value);  
  void PrintParameters();
 
private:

  const G4double     inch = 2.54*cm;
  G4Material*        fVacuumMater;
  G4Material*        fDetectorMater;
  G4LogicalVolume*   fLogicDetector;

  G4LogicalVolume*   fLogicCrystal;//in order to call from ConstructSDandField

  G4LogicalVolume*   fLogicPMT;
  G4LogicalVolume*   fLogicWrap;
  G4LogicalVolume*   fLogicPMTcover;
  G4LogicalVolume*   fLogicWrapFront;
  G4Material*        fPMTmater;
  G4Material*        fWrapMater;
  G4Material*        fPMTcoverMater;

  G4Material* fFrameMater;

  G4double           fWorld_X;
  G4double           fWorld_Y;
  G4double           fWorld_Z;

  G4double             fMom_X;//mother volume(contains temperature control box)
  G4double             fMom_Y;
  G4double             fMom_Z;
  G4double             fMom_pos_X;//mother volume position(in World volume)
  G4double             fMom_pos_Y;
  G4double             fMom_pos_Z;

  G4double             fTemp_X;//temp control box(contains Single)
  G4double             fTemp_Y;
  G4double             fTemp_Z;
  G4double             fTemp_pos_X;//temp control box position inside mother volume
  G4double             fTemp_pos_Y;
  G4double             fTemp_pos_Z;

  G4double             fSingle_X;
  G4double             fSingle_Y;
  G4double             fSingle_Z;

  G4double             gap;

  G4double             fFrame_length;

  G4double             fCrystal_X;//PbWO4
  G4double             fCrystal_Y;//PbWO4
  G4double             fCrystal_Z;//PbWO4
  G4double             fCrystal_pos_X;//In case there will be PMT someday, 
  G4double             fCrystal_pos_Y;//In case there will be PMT someday, 
  G4double             fCrystal_pos_Z;//In case there will be PMT someday, 

  G4double             fWrapThickness;
  G4double             fPMTcoverThickness;

  G4double             fPMT_radius;
  G4double             fPMT_length;
  G4double             fPMT_pos_X;
  G4double             fPMT_pos_Y;
  G4double             fPMT_pos_Z;

  G4Material*          fWorldMater;     
  G4VPhysicalVolume*   fPhysiWorld;
                
  DetectorMessenger*   fDetectorMessenger;


private:
    
  void               DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();     
  G4bool fCheckOverlaps;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

