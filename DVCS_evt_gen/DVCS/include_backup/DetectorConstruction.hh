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
//20171009(start)
#include "G4FieldManager.hh"
//20171009(finish)

//20171114(start)
#include "SimpleField.hh"
//20171114(finish)

//20171123(start)
#include "G4SystemOfUnits.hh"//to use cm
//20171123(finish)

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
//20171009(start)
class B5MagneticField;
class SensitiveDetector;
//20171009(finish)

//20171114(start)
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4ChordFinder;
//20171114(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //20180118(start)
  //  DetectorConstruction();
  //20181026(start)
  //  DetectorConstruction(G4double, G4double, G4double, G4bool);
  //20190409(start)
  // DetectorConstruction(G4double, G4double, G4double, G4bool, G4double);
  DetectorConstruction(G4double, G4double, G4bool, G4double);
  //20190409(finish)
  //20181026(finish)
  //20180117(finish)
  ~DetectorConstruction();

public:
  
  virtual G4VPhysicalVolume* Construct();
    
  //20171009
  virtual void ConstructSDandField();
  //20171009

  void SetTargetLength (G4double value);
  void SetTargetRadius (G4double value);
  void SetTargetMaterial (G4String);
    
  void SetDetectorLength(G4double value);           
  void SetDetectorThickness(G4double value);  
  void SetDetectorMaterial(G4String);               

  //20180117(start)
  void SetDetectorGap(G4double value);  
  //20180117(finish)

  //20180426(start)
  void SetMagneticField(G4bool value);  
  //20180426(finish)

  //20171114(start)
  //  void ConstructField();    
  //20171114(finish)
 
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
  //20171122(start)
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
  //20180314(start)
  G4LogicalVolume*   fLogicInnerChamber2;
  //20180314(finish)
  G4Material*        fChamberWindowMater;
  G4LogicalVolume*   fLogicChamberWindow;
  //20171122(finish)
  //20171124(start)
  G4double           fWindowHight;
  //20180313(start)
  G4double           fWindowSubHight;
  G4double           fWindowFrameThickness;
  G4double           fWindowFrameHight;
  G4double           fWindowClampThickness;
  G4double           fWindowClampHight;
  G4LogicalVolume*   fLogicWindowFrame;
  G4LogicalVolume*   fLogicWindowClamp;
  //20180313(finish)

  //20171128(start)
  G4double           fWindowThickness;
  //20171128(finish)
  //20171124(finish)
  G4double           fTargetLength; 
  G4double           fTargetRadius;
  G4double           fTargetCoverLength; 
  G4double           fTargetCoverRadius;
  G4double           fTargetWindowThickness;
  G4Material*        fTargetMater;
  G4LogicalVolume*   fLogicTarget;
  //20171006(start)
  G4Material*        fTargetCoverMater;
  G4LogicalVolume*   fLogicTargetCover;
  //20171006(finish)                 

  //20180314(start)
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
  //20180314(finish)

  //20180328(start)
  //20180329(start)
  // G4double           fBeampipe1OuterRadius;
  // G4double           fBeampipe1InnerRadius;

  G4double           fBeampipe1Innerdx1;
  G4double           fBeampipe1Innerdx2;
  G4double           fBeampipe1Innerdy1;
  G4double           fBeampipe1Innerdy2;
  G4double           fBeampipe1Outerdx1;
  G4double           fBeampipe1Outerdx2;
  G4double           fBeampipe1Outerdy1;
  G4double           fBeampipe1Outerdy2;
  //20180329(finish)
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
  //20180328(finish)

  G4double           fDetectorLength;
  G4double           fDetectorThickness;
  G4Material*        fDetectorMater;
  G4LogicalVolume*   fLogicDetector;

  //20171010(start)
  G4LogicalVolume*   fLogicCrystal;//in order to call from ConstructSDandField
  //20171010(finish)

  //20190121(start)
  G4LogicalVolume*   fLogicFluxOptFiber;
  //20190121(finish)

  //20171103(start)
  G4LogicalVolume*   fLogicMagnetic;//in order to call from ConstructSDandField//for the magnetic field in certain volume. can be used someday
  //20171103(finish)

  //20180412(start)
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
  //20180412(finish)

  //20180117(start)
  G4LogicalVolume*   fLogicPMT;
  G4LogicalVolume*   fLogicWrap;
  G4LogicalVolume*   fLogicPMTcover;
  G4LogicalVolume*   fLogicWrapFront;
  G4Material*        fPMTmater;
  G4Material*        fWrapMater;
  G4Material*        fPMTcoverMater;
  //20170117(finish)

  //20180222(start)
  G4Material* fFrameMater;
  //20180222(finish)

  //20171127(start)
  G4LogicalVolume*   fLogicFlux;//20180713, using it for pseudo sphere
  //20171127(finish)

  //20171114(start)
  G4LogicalVolume* logicEnv;
  G4double      env_size[3];
  G4Material*   env_material;
  G4Material*   coil_material;
  G4Material*   core_material;
  G4Material*   coilinsert_material;
  G4Material*   concrete_shield;
  //  beamline_material = nist->FindOrBuildMaterial("G4_Fe");
  G4Material*   radiator_material;
  SimpleField *fField;
G4Mag_UsualEqRhs *fEquation;
   G4MagIntegratorStepper* fStepper;
G4ChordFinder*          fChordFinder;
  //20171114(finish)

  //20180426(start)
  G4bool             field;
  //20180426(finish)

  //20181026(start)
  G4double      fFieldStr;
  //20181026(finish)

 /*20171006 deleted for new world box               
    G4double           fWorldLength;
    G4double           fWorldRadius;
  */
  //20171006 added new world box
  G4double           fWorld_X;
  G4double           fWorld_Y;
  G4double           fWorld_Z;
  //20171006

  //20181022(start)
  G4double           fNPS_distance;
  //20181022(finish)

  //20171006 added
  G4double             fMom_X;//mother volume#1(contains temp control box, detector, support frames)
  G4double             fMom_Y;
  G4double             fMom_Z;
  G4double             fMom_pos_X;//mother volume#1 position(to move it freely)
  G4double             fMom_pos_Y;
  G4double             fMom_pos_Z;

  //20171121(start)
  G4double             fMom_theta;
  //20171121(finish)

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

  //20180117(start)
  G4double             fWrapThickness;
  G4double             fPMTcoverThickness;
  //20180117(finish)

  G4double             fPMT_radius;
  G4double             fPMT_length;
  G4double             fPMT_pos_X;
  G4double             fPMT_pos_Y;
  G4double             fPMT_pos_Z;
  //20171006
  //20171121(start)
  //the actual size of the physical sweeping magnet
  G4double             fSweepingMagnet_X;
  G4double             fSweepingMagnet_Y;
  G4double             fSweepingMagnet_Z;
  //20181016(start)
  //20180412(start)
  // G4double             fSweepingMagnet_pos_X;
  // G4double             fSweepingMagnet_pos_Y;
  // //20180412(finish)
  // G4double             fSweepingMagnet_pos_Z;
  //20181016(finish)

  //20190905(start)
  G4double             fSweepingMagnet_angle;//magnet_center's angle with respect to the NPS: 4deg(default) or 5.5deg
  G4double             fSweepingMagnet_center_dist;//magnet_center's distance from the target : 1.6m (for WACS, some of the settings are 1.4m)
  //20190905(finish)

  //20180417(start)
  G4double             fSweepingMagnetSTPshift_X;
  G4double             fSweepingMagnetSTPshift_Y;
  G4double             fSweepingMagnetSTPshift_Z;
  G4double             fSweepingMagnetSTPcenter_X;
  G4double             fSweepingMagnetSTPcenter_Y;
  G4double             fSweepingMagnetSTPcenter_Z;
  //20180417(finish)
  //20181016(start)
  G4double             fSweepingMagnet_pos;
  G4double             fSweepingMagnet_arm_theta;
  //20181016(finish)
  G4double             fSweepingMagnet_theta;
  //20171121(finish)
  //20180928(start)
  G4double             fSweepingMagnetField_theta;
  //20180928(finish)
  //20180412(start)
  G4double             fSweepingMagnetShift_X;
  G4double             fSweepingMagnetShift_Y;
  G4double             fSweepingMagnetShift_Z;
  //20180412(finish)
  //20171130(start)
  //where the field is installed
  G4double             fSweepingMagnet_Mom_X;
  G4double             fSweepingMagnet_Mom_Y;
  G4double             fSweepingMagnet_Mom_Z;
  //20171130(finish)
  G4Material*          fWorldMater;     
  G4VPhysicalVolume*   fPhysiWorld;
                
  DetectorMessenger*   fDetectorMessenger;


  //20171009(start)
  static G4ThreadLocal B5MagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;
  //20171009(finish)

private:
    
  void               DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();     
  G4bool fCheckOverlaps; //20171009 added

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

