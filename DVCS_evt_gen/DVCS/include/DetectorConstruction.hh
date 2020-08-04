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
#include "G4FieldManager.hh"

#include "SimpleField.hh"

#include "G4SystemOfUnits.hh"//to use cm

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
class B5MagneticField;
class SensitiveDetector;

class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4ChordFinder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(G4double, G4double, G4bool, G4double);
  ~DetectorConstruction();

public:
  
  virtual G4VPhysicalVolume* Construct();
    
  virtual void ConstructSDandField();

  void SetTargetLength (G4double value);
  void SetTargetRadius (G4double value);
  void SetTargetMaterial (G4String);
    
  void SetDetectorLength(G4double value);           
  void SetDetectorThickness(G4double value);  
  void SetDetectorMaterial(G4String);               

  void SetDetectorGap(G4double value);  

  void SetMagneticField(G4bool value);  
 
  void PrintParameters();
 
public:
      
  G4double GetTargetLength();
  G4double GetTargetRadius();
  G4Material* GetTargetMaterial();       
  G4LogicalVolume* GetLogicTarget();
    
  G4double GetDetectorLength();
  G4double GetDetectorThickness();
  G4Material* GetDetectorMaterial();                 
  G4LogicalVolume* GetLogicDetector();      
                       
private:
  const G4double     inch = 2.54*cm;
  G4double           fChamberOuterRadius;
  G4double           fChamberInnerRadius;
  G4double           fChamberHight;
  G4Material*        fVacuumMater;
  G4Material*        fChamberMater;
  G4Material*        fWindowMater;
  G4Material*        fKaptonMater;
  G4LogicalVolume*   lWorld;
  G4LogicalVolume*   fLogicChamber;
  G4LogicalVolume*   fLogicInnerChamber;
  G4LogicalVolume*   fLogicInnerChamber2;
  G4Material*        fChamberWindowMater;
  G4LogicalVolume*   fLogicChamberWindow;
  G4double           fWindowHight;
  G4double           fWindowSubHight;
  G4double           fWindowFrameThickness;
  G4double           fWindowFrameHight;
  G4double           fWindowClampThickness;
  G4double           fWindowClampHight;
  G4LogicalVolume*   fLogicWindowFrame;
  G4LogicalVolume*   fLogicWindowClamp;

  G4double           fWindowThickness;
  G4double           fTargetLength; 
  G4double           fTargetRadius;
  G4double           fTargetCoverLength; 
  G4double           fTargetCoverRadius;
  G4double           fTargetWindowThickness;
  G4Material*        fTargetMater;
  G4LogicalVolume*   fLogicTarget;
  G4Material*        fTargetCoverMater;
  G4LogicalVolume*   fLogicTargetCover;

  G4double           fWindowInnerJoint1InnerRadius;
  G4double           fWindowInnerJoint1OuterRadius;
  G4double           fWindowInnerJoint1Thickness;
  G4double           fWindowInnerJoint2InnerRadius;
  G4double           fWindowInnerJoint2OuterRadius;
  G4double           fWindowInnerJoint2Thickness;
  G4double           fWindowOuterJoint1InnerRadius;
  G4double           fWindowOuterJoint1OuterRadius;
  G4double           fWindowOuterJoint1Thickness;
  G4double           fWindowOuterJoint2InnerRadius;
  G4double           fWindowOuterJoint2OuterRadius;
  G4double           fWindowOuterJoint2Thickness;
  G4double           fWindowOuterJoint2_1InnerRadius;
  G4double           fWindowOuterJoint2_1OuterRadius;
  G4double           fWindowOuterJoint2_1Thickness;
  G4double           fWindowOuterJoint2_1Position;
  G4double           fWindowOuterJoint2_2InnerRadius;
  G4double           fWindowOuterJoint2_2dx;
  G4double           fWindowOuterJoint2_2dy;
  G4double           fWindowOuterJoint2_2dz;
  G4double           fWindowOuterJoint2_3InnerRadius;
  G4double           fWindowOuterJoint2_3OuterRadius;
  G4double           fWindowOuterJoint2_3Thickness;

  G4Material*        fWindowInnerJoint1Mater;
  G4Material*        fWindowInnerJoint2Mater;
  G4Material*        fWindowOuterJoint1Mater;
  G4Material*        fWindowOuterJoint2Mater;
  G4Material*        fWindowOuterJoint2_1Mater;
  G4Material*        fWindowOuterJoint2_2Mater;
  G4Material*        fWindowOuterJoint2_3Mater;

  G4LogicalVolume*   fLogicWindowInnerJoint1;
  G4LogicalVolume*   fLogicWindowInnerJoint2;
  G4LogicalVolume*   fLogicWindowOuterJoint1;
  G4LogicalVolume*   fLogicWindowOuterJoint2;
  G4LogicalVolume*   fLogicWindowOuterJoint2_1;
  G4LogicalVolume*   fLogicWindowOuterJoint2_2;
  G4LogicalVolume*   fLogicWindowOuterJoint2_3;

  G4double           fBeampipe1Innerdx1;
  G4double           fBeampipe1Innerdx2;
  G4double           fBeampipe1Innerdy1;
  G4double           fBeampipe1Innerdy2;
  G4double           fBeampipe1Outerdx1;
  G4double           fBeampipe1Outerdx2;
  G4double           fBeampipe1Outerdy1;
  G4double           fBeampipe1Outerdy2;
  G4double           fBeampipe1Length;
  G4double           fBeampipe2OuterRadius;
  G4double           fBeampipe2InnerRadius;
  G4double           fBeampipe2Length;
  G4double           fBeampipe2FrontCoverThickness;
  G4double           fBeampipe3OuterRadius;
  G4double           fBeampipe3InnerRadius;
  G4double           fBeampipe3Length;
  G4double           fBeampipe3FrontCoverThickness;

  G4Material*        fBeampipe1Mater;
  G4Material*        fBeampipe2Mater;
  G4Material*        fBeampipe3Mater;

  G4LogicalVolume*   fLogicBeampipe1;
  G4LogicalVolume*   fLogicBeampipe2;
  G4LogicalVolume*   fLogicBeampipe3;
  G4LogicalVolume*   fLogicBeampipe2FrontCover;
  G4LogicalVolume*   fLogicBeampipe3FrontCover;

  G4double           fDetectorLength;
  G4double           fDetectorThickness;
  G4Material*        fDetectorMater;
  G4LogicalVolume*   fLogicDetector;

  G4LogicalVolume*   fLogicCrystal;//in order to call from ConstructSDandField

  G4LogicalVolume*   fLogicFluxOptFiber;

  G4LogicalVolume*   fLogicMagnetic;//in order to call from ConstructSDandField//for the magnetic field in certain volume. can be used someday

  G4Material*        fSweepingMagnet_1Mater;
  G4Material*        fSweepingMagnet_2Mater;
  G4Material*        fSweepingMagnet_3Mater;
  G4Material*        fSweepingMagnet_4Mater;
  G4Material*        fSweepingMagnet_5Mater;
  G4Material*        fSweepingMagnet_6_1Mater;
  G4Material*        fSweepingMagnet_6_234Mater;
  G4LogicalVolume*   fLogicSweepingMagnet_1_1;
  G4LogicalVolume*   fLogicSweepingMagnet_1_2;
  G4LogicalVolume*   fLogicSweepingMagnet_1_3;
  G4LogicalVolume*   fLogicSweepingMagnet_1_4_1;
  G4LogicalVolume*   fLogicSweepingMagnet_1_4_2;
  G4LogicalVolume*   fLogicSweepingMagnet_2_1_1;
  G4LogicalVolume*   fLogicSweepingMagnet_2_1_2;
  G4LogicalVolume*   fLogicSweepingMagnet_2_2;
  G4LogicalVolume*   fLogicSweepingMagnet_3_1;
  G4LogicalVolume*   fLogicSweepingMagnet_3_2;
  G4LogicalVolume*   fLogicSweepingMagnet_3_3;
  G4LogicalVolume*   fLogicSweepingMagnet_4_1;
  G4LogicalVolume*   fLogicSweepingMagnet_4_2;
  G4LogicalVolume*   fLogicSweepingMagnet_4_3;
  G4LogicalVolume*   fLogicSweepingMagnet_5_1_1;
  G4LogicalVolume*   fLogicSweepingMagnet_5_1_2;
  G4LogicalVolume*   fLogicSweepingMagnet_5_2;
  G4LogicalVolume*   fLogicSweepingMagnet_6_1;
  G4LogicalVolume*   fLogicSweepingMagnet_6_2;
  G4LogicalVolume*   fLogicSweepingMagnet_6_3;
  G4LogicalVolume*   fLogicSweepingMagnet_6_4;

  G4LogicalVolume*   fLogicPMT;
  G4LogicalVolume*   fLogicWrap;
  G4LogicalVolume*   fLogicPMTcover;
  G4LogicalVolume*   fLogicWrapFront;
  G4Material*        fPMTmater;
  G4Material*        fWrapMater;
  G4Material*        fPMTcoverMater;

  G4Material* fFrameMater;

  G4LogicalVolume*   fLogicFlux;// using it for pseudo sphere

  SimpleField *           fField;
  G4Mag_UsualEqRhs*       fEquation;
  G4MagIntegratorStepper* fStepper;
  G4ChordFinder*          fChordFinder;

  G4bool             field;

  G4double      fFieldStr;

  G4double           fWorld_X;
  G4double           fWorld_Y;
  G4double           fWorld_Z;

  G4double           fNPS_distance;

  G4double             fMom_X;//mother volume#1(contains temp control box, detector, support frames)
  G4double             fMom_Y;
  G4double             fMom_Z;
  G4double             fMom_pos_X;//mother volume#1 position(to move it freely)
  G4double             fMom_pos_Y;
  G4double             fMom_pos_Z;

  G4double             fMom_theta;

  G4double             fTemp_X;//temp control box(mother volume#2)(contains detectors)
  G4double             fTemp_Y;//temp control box(mother volume#2)(contains detectors)
  G4double             fTemp_Z;//temp control box(mother volume#2)(contains detectors)
  G4double             fTemp_pos_X;//temp control box position inside mother volume#1
  G4double             fTemp_pos_Y;//temp control box position inside mother volume#1
  G4double             fTemp_pos_Z;//temp control box position inside mother volume#1

  G4double             fNPS_X;
  G4double             fNPS_Y;
  G4double             fNPS_Z;
  G4double             fNPS_pos_X;
  G4double             fNPS_pos_Y;
  G4double             fNPS_pos_Z;

  G4double             fBulk_X;
  G4double             fBulk_Y;
  G4double             fBulk_Z;

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

  //the actual size of the physical sweeping magnet
  G4double             fSweepingMagnet_X;
  G4double             fSweepingMagnet_Y;
  G4double             fSweepingMagnet_Z;

  G4double             fSweepingMagnet_angle;//magnet_center's angle with respect to the NPS: 4deg(default) or 5.5deg
  G4double             fSweepingMagnet_center_dist;//magnet_center's distance from the target : 1.6m (for WACS, some of the settings are 1.4m)

  G4double             fSweepingMagnetSTPshift_X;
  G4double             fSweepingMagnetSTPshift_Y;
  G4double             fSweepingMagnetSTPshift_Z;
  G4double             fSweepingMagnetSTPcenter_X;
  G4double             fSweepingMagnetSTPcenter_Y;
  G4double             fSweepingMagnetSTPcenter_Z;
  G4double             fSweepingMagnet_pos;
  G4double             fSweepingMagnet_arm_theta;
  G4double             fSweepingMagnet_theta;
  G4double             fSweepingMagnetField_theta;
  G4double             fSweepingMagnetShift_X;
  G4double             fSweepingMagnetShift_Y;
  G4double             fSweepingMagnetShift_Z;

  //where the field is installed
  G4double             fSweepingMagnet_Mom_X;
  G4double             fSweepingMagnet_Mom_Y;
  G4double             fSweepingMagnet_Mom_Z;

  G4Material*          fWorldMater;     
  G4VPhysicalVolume*   fPhysiWorld;
                
  DetectorMessenger*   fDetectorMessenger;

  static G4ThreadLocal B5MagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

private:
    
  void               DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();     
  G4bool fCheckOverlaps;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

