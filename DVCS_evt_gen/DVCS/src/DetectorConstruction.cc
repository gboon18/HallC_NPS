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
/// \file hadronic/Hadr03/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SDManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "B5HadCalorimeterSD.hh"
#include "G4AutoDelete.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"//coloring

#include "G4Trd.hh"
#include "G4GenericTrap.hh"

#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "SimpleField.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4ExplicitEuler.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"

#include "G4SubtractionSolid.hh"

#include "FluxSD.hh"

#include "ChambWindSD.hh"

#include "CrystalCoverSD.hh"
#include "CrystalFrontCoverSD.hh"
#include "PMTcoverSD.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

#include "dvcsGlobals.hh"
#include "TMath.h"

#include "G4Cons.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4double NPS_distance, G4double NPS_angle, G4bool field_input, G4double field_str)
:G4VUserDetectorConstruction(),
 fVacuumMater(0), fChamberMater(0), fWindowMater(0), fKaptonMater(0), lWorld(0),
 fLogicChamber(0), fLogicInnerChamber(0),
  fLogicInnerChamber2(0),
  fChamberWindowMater(0), fLogicChamberWindow(0),
  fLogicWindowFrame(0), fLogicWindowClamp(0),
 fTargetMater(0), fLogicTarget(0),
  fTargetCoverMater(0), fLogicTargetCover(0),
 fWindowInnerJoint1Mater(0), fWindowInnerJoint2Mater(0),
  fWindowOuterJoint1Mater(0), fWindowOuterJoint2Mater(0), fWindowOuterJoint2_1Mater(0), fWindowOuterJoint2_2Mater(0), fWindowOuterJoint2_3Mater(0),
  fLogicWindowInnerJoint1(0), fLogicWindowInnerJoint2(0),
  fLogicWindowOuterJoint1(0), fLogicWindowOuterJoint2(0), fLogicWindowOuterJoint2_1(0), fLogicWindowOuterJoint2_2(0), fLogicWindowOuterJoint2_3(0),
  fBeampipe1Mater(0), fBeampipe2Mater(0), fBeampipe3Mater(0), fLogicBeampipe1(0), fLogicBeampipe2(0), fLogicBeampipe3(0), fLogicBeampipe2FrontCover(0), fLogicBeampipe3FrontCover(0),
  fDetectorMater(0), fLogicDetector(0),
  fLogicCrystal(0),
  fLogicFluxOptFiber(0),
  fLogicMagnetic(0),
  fSweepingMagnet_1Mater(0),
  fSweepingMagnet_2Mater(0),
  fSweepingMagnet_3Mater(0),
  fSweepingMagnet_4Mater(0),
  fSweepingMagnet_5Mater(0),
  fSweepingMagnet_6_1Mater(0), fSweepingMagnet_6_234Mater(0),
  fLogicSweepingMagnet_1_1(0), fLogicSweepingMagnet_1_2(0), fLogicSweepingMagnet_1_3(0), fLogicSweepingMagnet_1_4_1(0), fLogicSweepingMagnet_1_4_2(0),
  fLogicSweepingMagnet_2_1_1(0), fLogicSweepingMagnet_2_1_2(0), fLogicSweepingMagnet_2_2(0),
  fLogicSweepingMagnet_3_1(0), fLogicSweepingMagnet_3_2(0), fLogicSweepingMagnet_3_3(0),
  fLogicSweepingMagnet_4_1(0), fLogicSweepingMagnet_4_2(0), fLogicSweepingMagnet_4_3(0),
  fLogicSweepingMagnet_5_1_1(0), fLogicSweepingMagnet_5_1_2(0), fLogicSweepingMagnet_5_2(0),
  fLogicSweepingMagnet_6_1(0), fLogicSweepingMagnet_6_2(0), fLogicSweepingMagnet_6_3(0), fLogicSweepingMagnet_6_4(0),
  fLogicPMT(0), fLogicWrap(0), fLogicPMTcover(0), fLogicWrapFront(0),
  fPMTmater(0), fWrapMater(0), fPMTcoverMater(0), 
  fFrameMater(0),
  fLogicFlux(0),
  fField(0),
  fEquation(0),
  fStepper(0),
  fChordFinder(0),
  field(field_input),
  fWorldMater(0), fPhysiWorld(0),
  fDetectorMessenger(0)
{
  fChamberOuterRadius = 22.5*inch,
    fChamberInnerRadius = 20.5*inch,
    fChamberHight       = 24.25*2*inch;//to be modified

  fWindowHight = 19*inch;
  fWindowSubHight = 15*inch;
  fWindowFrameThickness = 1.25*inch;//guessed
  fWindowFrameHight = 20*inch;//guessed
  fWindowClampThickness = 0.750*inch;
  fWindowClampHight = 20*inch;
  fWindowThickness = 0.020*inch;

  fTargetLength      = 150.*mm;//(129.759 - 29.759 + 50.)*mm; 
  fTargetRadius      = 0.5*50.*mm;//20.179*mm;
  fTargetCoverLength = (0.125 + 150.)*mm;//(5. + 129.759 - 29.759 + 50.)*mm; //(129.887 - 29.759 + 50.)*mm; 
  fTargetCoverRadius = (0.125 + 0.5*50.)*mm;//(20.179 + 5.)*mm;//20.32*mm;
  fTargetWindowThickness = 0.125*mm;//5.*mm;//0.128*mm;

  fWindowInnerJoint1InnerRadius   = 0.5*0.846*inch;
  fWindowInnerJoint1OuterRadius   = 0.5*1.469*inch;
  fWindowInnerJoint1Thickness     = 0.5*inch;
  fWindowInnerJoint2InnerRadius   = 0.5*1.068*inch;
  fWindowInnerJoint2OuterRadius   = 0.5*1.50*inch;
  fWindowInnerJoint2Thickness     = 0.109*inch;
  fWindowOuterJoint1InnerRadius   = 0.5*1.7*inch;//no need
  fWindowOuterJoint1OuterRadius   = 0.5*2.38*inch;//no need
  fWindowOuterJoint1Thickness     = 0.016*inch;//no need
  fWindowOuterJoint2InnerRadius   = 0.5*0.68*inch;
  fWindowOuterJoint2OuterRadius   = 0.5*1.062*inch;
  fWindowOuterJoint2Thickness     = (0.62 + 1.562 - 0.137 - 0.06)*inch;
  fWindowOuterJoint2_1InnerRadius = fWindowOuterJoint2OuterRadius;
  fWindowOuterJoint2_1OuterRadius = 0.5*1.5*inch;
  fWindowOuterJoint2_1Thickness   = 0.19*inch;
  fWindowOuterJoint2_1Position    = 0.62*inch;
  fWindowOuterJoint2_2InnerRadius = fWindowOuterJoint2InnerRadius;
  fWindowOuterJoint2_2dx          = 1.12*inch;
  fWindowOuterJoint2_2dy          = 1.12*inch ;
  fWindowOuterJoint2_2dz          = 0.137*inch;
  fWindowOuterJoint2_3InnerRadius   = 0.5*0.68*inch;//= fWindowOuterJoint2InnerRadius;
  fWindowOuterJoint2_3OuterRadius   = 0.5*0.738*inch;
  fWindowOuterJoint2_3Thickness     = 0.06*inch;

  fBeampipe1Innerdx1    = 18.915*mm;
  fBeampipe1Innerdx2    = 63.4*mm;
  fBeampipe1Innerdy1    = fBeampipe1Innerdx1;
  fBeampipe1Innerdy2    = fBeampipe1Innerdx2;
  fBeampipe1Outerdx1    = 25.265*mm;
  fBeampipe1Outerdx2    = 69.749*mm;
  fBeampipe1Outerdy1    = fBeampipe1Outerdx1;
  fBeampipe1Outerdy2    = fBeampipe1Outerdx2;
  fBeampipe1Length      = 1685.925*mm;
  fBeampipe2OuterRadius = 0.5*168.275*mm;
  fBeampipe2InnerRadius = 0.5*154.051*mm;
  fBeampipe2Length      = 129.633*inch;//changed from 2620.138*mm;<-this is the actual length of the pipe.
  //however, there is a gap between beampipe1 and 2. In order to fill the gap. The length of the beampipe2 currently is longer.
  fBeampipe2FrontCoverThickness = fBeampipe2OuterRadius - fBeampipe2InnerRadius; //made up
  fBeampipe3OuterRadius = 0.5*273.05*mm;
  fBeampipe3InnerRadius = 0.5*254.508*mm;
  fBeampipe3Length      = 5102.225*mm;
  fBeampipe3FrontCoverThickness = fBeampipe3OuterRadius - fBeampipe3InnerRadius; //made up

  fDetectorLength    = 5*cm; 
  fDetectorThickness = 2*cm;

  fWorld_X = 2*1.5*2.1*4.1*m;
  fWorld_Y = 2*1.5*2.1*4.1*m;
  fWorld_Z = 2*1.5*2.1*4.1*m;

  gap = 0.5*mm;//carbon frame thickness
  fFrame_length = 20.*mm;//carbon frame length;

  field = field_input;
  fFieldStr = field_str;

  fWrapThickness = 65*1e-3*mm;
  fPMTcoverThickness = 65*1e-3*mm;

  fNPS_distance = NPS_distance*mm;

  // fMom_X = 713.5*mm;//mother volume#1(contains temp control box, detector, support frames)
  fMom_pos_X = 0*mm;//mother volume#1 position(to move it freely)
  fMom_pos_Y = 0*mm;

  //Angle in radian
  fMom_theta = -NPS_angle;//[rad]

  // fTemp_X = 713*mm;//temp control box(mother volume#2)(contains detectors). this is inside the mother volume#1
  fTemp_pos_X = 0*m;//temp control box position inside mother volume#1
  fTemp_pos_Y = 0*m;
  fTemp_pos_Z = 0*m;

  fCrystal_X = 20.5*mm;//PbWO4
  fCrystal_Y = 20.5*mm;//inside mother volume#5 with PMT
  fCrystal_Z = 200.5*mm;
  
  fPMT_radius = 18.6*0.5*mm;
  fPMT_length = 88.*mm;

  fSweepingMagnet_X = 55.925*inch;
  fSweepingMagnet_Y = 60*inch;
  fSweepingMagnet_Z = 61.1*inch;

  fSweepingMagnet_angle = dvcsGlobals::SM_angle;//magnet_center's angle with respect to the NPS: 4deg(default) or 5.5deg
  fSweepingMagnet_center_dist = 1.6*m;//magnet_center's distance from the target : 1.6m (for WACS, some of the settings are 1.4m)

  fSweepingMagnetSTPshift_X = -1*(0.5*fSweepingMagnet_X - 15.*inch/*x of part1_3*/ - 21.*inch/*x of part1_2_1*/);//-1* is to use shift_x as a positive value
  fSweepingMagnetSTPshift_Y = 0.*inch;
  fSweepingMagnetSTPshift_Z = 0.5*fSweepingMagnet_Z - 5.2*inch/*z of part6_4*/ - 0.5*39.5*inch/*half of z of part1*/;

  fSweepingMagnetSTPcenter_X = (fSweepingMagnet_center_dist - fSweepingMagnetSTPshift_Z - fSweepingMagnetSTPshift_X*TMath::Tan(NPS_angle - fSweepingMagnet_angle))*TMath::Sin(NPS_angle - fSweepingMagnet_angle) + fSweepingMagnetSTPshift_X/TMath::Cos(NPS_angle - fSweepingMagnet_angle);
  fSweepingMagnetSTPcenter_Y = 0.*inch;
  fSweepingMagnetSTPcenter_Z = (fSweepingMagnet_center_dist - fSweepingMagnetSTPshift_Z - fSweepingMagnetSTPshift_X*TMath::Tan(NPS_angle - fSweepingMagnet_angle))*TMath::Cos(NPS_angle - fSweepingMagnet_angle); 

  fSweepingMagnet_pos = TMath::Sqrt(fSweepingMagnetSTPcenter_X*fSweepingMagnetSTPcenter_X + fSweepingMagnetSTPcenter_Y*fSweepingMagnetSTPcenter_Y + fSweepingMagnetSTPcenter_Z*fSweepingMagnetSTPcenter_Z);

  //This is the global angle of the magnet. If magnet is tilted, use fSweepingMagnet_theta(currently no tilting)
  fSweepingMagnet_arm_theta = -NPS_angle + fSweepingMagnet_angle - TMath::ATan2(fSweepingMagnetSTPshift_X, fSweepingMagnet_center_dist - fSweepingMagnetSTPshift_Z);

  //based on https://userweb.jlab.org/~medeiros/NPS/DETECTOR/STEP-detector_reference.stp//NPS is positioned 6 degrees.
  fSweepingMagnet_theta = -NPS_angle + fSweepingMagnet_angle;//to make -2degrees when NPS_angle is -6degrees. Sweeping Magnet moves along with NPS.
  fSweepingMagnetField_theta = NPS_angle - fSweepingMagnet_angle - .1*pi/180.;//to make 2.2degrees when NPS_angle is 6.3degrees. (When the angle between NPS and magnet is 4 deg.)

  fSweepingMagnetShift_X = 0.5*fSweepingMagnet_X - (0.5*43.8)*inch;
  fSweepingMagnetShift_Y = 0;
  fSweepingMagnetShift_Z = 0.5*fSweepingMagnet_Z - (0.5*39.5 + 5.2)*inch;

  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4double density, fractionmass;
  G4int ncomponents, natoms;

  G4NistManager* nist = G4NistManager::Instance();

  //World Material
  G4Element* N  = new G4Element("Nitrogen",  "N",   7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",    "O",   8, 16.00*g/mole);
  G4Element* Al  = new G4Element("Aluminium","Al", 13, 26.98*g/mole);

  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
				     kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20->AddElement(N, fractionmass=0.7);
  Air20->AddElement(O, fractionmass=0.3);
  fWorldMater = Air20;

  fVacuumMater = nist->FindOrBuildMaterial("G4_Galactic");//vacuum
  fChamberMater = nist->FindOrBuildMaterial("G4_Fe");
  fWindowMater = nist->FindOrBuildMaterial("G4_Al");//change it later for more detailed Al.
  fKaptonMater = nist->FindOrBuildMaterial("G4_KAPTON");

  int targ_type = dvcsGlobals::target_type;
  if(targ_type == 0){
    G4cout<<"DetectorConstruction, target type is liquid hydrogen"<<G4endl;
    fTargetMater = nist->FindOrBuildMaterial("G4_lH2");
  }
  if(targ_type == 1){
    G4cout<<"DetectorConstruction, target type is liquid deuterium"<<G4endl;
    fTargetMater =  new G4Material("LD2",  1, 2.01*g/mole, 0.167*g/cm3, kStateLiquid, 22*kelvin);
  }
  //Target Cover Material Al
  fTargetCoverMater =new G4Material("TargetCover", 2.70*g/cm3, ncomponents=1);
  fTargetCoverMater -> AddElement(Al, 1);

  //Detector Material PbWO4
  G4Isotope *Pb204 = new G4Isotope("204Pb", 82, 204, 203.97*g/mole);
  G4Isotope *Pb206 = new G4Isotope("204Pb", 82, 206, 205.97*g/mole);
  G4Isotope *Pb207 = new G4Isotope("204Pb", 82, 207, 206.98*g/mole);
  G4Isotope *Pb208 = new G4Isotope("204Pb", 82, 208, 207.98*g/mole);
  G4Element* Pb = new G4Element("Plomo"    ,"Pb", 4); //Introducimos el plomo  
  Pb->AddIsotope(Pb204, 1.400*perCent);
  Pb->AddIsotope(Pb206, 24.100*perCent);
  Pb->AddIsotope(Pb207, 22.100*perCent);
  Pb->AddIsotope(Pb208, 52.400*perCent);

  G4Isotope *W180 = new G4Isotope("180W", 74, 180, 179.95*g/mole);
  G4Isotope *W182 = new G4Isotope("182W", 74, 182, 181.95*g/mole);
  G4Isotope *W183 = new G4Isotope("183W", 74, 183, 182.95*g/mole);
  G4Isotope *W184 = new G4Isotope("184W", 74, 184, 183.95*g/mole);
  G4Isotope *W186 = new G4Isotope("186W", 74, 186, 185.95*g/mole);
  G4Element* W  = new G4Element("Wolframio","W" , 5); //Introducimos el Wolframio
  W->AddIsotope(W180, 0.120*perCent);
  W->AddIsotope(W182, 26.500*perCent);
  W->AddIsotope(W183, 14.310*perCent);
  W->AddIsotope(W184, 30.640*perCent);
  W->AddIsotope(W186, 28.430*perCent);

  //Material detector
  fDetectorMater = 
    new G4Material("PbWO4", density=8.3*g/cm3, ncomponents=3);
  fDetectorMater->AddElement(Pb, natoms=1);
  fDetectorMater->AddElement(W , natoms=1);
  fDetectorMater->AddElement(O , natoms=4);

  G4Element* Si = new G4Element("Silicon", "Si", 14, 28.09*g/mole);
  fPMTmater = new G4Material("SiO2", 2.648*g/cm3, ncomponents=2);
  fPMTmater -> AddElement(Si, fractionmass = 1./3.);
  fPMTmater -> AddElement(O, fractionmass = 2./3.);

  G4Element* C  = new G4Element("Carbon",    "C",   6, 12.01*g/mole);
  G4Element* H  = new G4Element("Hydrogen",  "H",   1, 1.008*g/mole);
  fWrapMater = new G4Material("VM2000", 1.38*g/cm3, ncomponents=3);
  fWrapMater->AddElement(C, fractionmass=10./22.);
  fWrapMater->AddElement(H, fractionmass=8./22.);
  fWrapMater->AddElement(O, fractionmass=4./22.);

  fPMTcoverMater = new G4Material("Al_cover", 2.70*g/cm3, ncomponents=1);
  fPMTcoverMater -> AddElement(Al, fractionmass = 1.);

  fFrameMater = new G4Material("Frame", 1.55*g/cm3, ncomponents=1);
  fFrameMater -> AddElement(C, fractionmass = 1.);

  G4Element* Cr = new G4Element("Chromium",    "Cr", 24, 51.9961*g/mole);
  G4Element* Cu = new G4Element("Copper",      "Cu", 29, 63.546*g/mole);
  G4Element* Fe = new G4Element("Iron",        "Fe", 26, 55.845*g/mole);
  G4Element* Mg = new G4Element("Magnesum",    "Mg", 12, 24.305*g/mole);
  G4Element* Mn = new G4Element("Manganese",   "Mn", 25, 54.938*g/mole);
  G4Element* Ti = new G4Element("Titanium",    "Ti", 22, 47.867*g/mole);
  G4Element* Zn = new G4Element("Zinc",        "Zn", 30, 65.38*g/mole);
  G4Element* P  = new G4Element("Phosphorous", "P", 15, 30.974*g/mole);
  G4Element* S  = new G4Element("Sulfur",      "S", 16, 32.06*g/mole);
  G4Element* Ni = new G4Element("Nickel",      "Ni", 28, 58.693*g/mole);

  G4Material* AL_7075_T6 = new G4Material("AL_7075_T6", 2.81*g/cm3, ncomponents = 9);//http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T6
  AL_7075_T6->AddElement(Al, fractionmass=(87.1 + 91.4)*0.5/100);
  AL_7075_T6->AddElement(Cr, fractionmass=(0.18 + 0.28)*0.5/100);
  AL_7075_T6->AddElement(Cu, fractionmass=(1.2  +    2)*0.5/100);
  AL_7075_T6->AddElement(Fe, fractionmass=(41./50.)*(5./(5.+3.+4.+2.))/100);//0.5/100);
  AL_7075_T6->AddElement(Mg, fractionmass=(2.1  +  2.9)*0.5/100);
  AL_7075_T6->AddElement(Mn, fractionmass=(41./50.)*(3./(5.+3.+4.+2.))/100);//0.3/100);
  AL_7075_T6->AddElement(Si, fractionmass=(41./50.)*(4./(5.+3.+4.+2.))/100);//0.4/100);
  AL_7075_T6->AddElement(Ti, fractionmass=(41./50.)*(2./(5.+3.+4.+2.))/100);//0.2/100);
  AL_7075_T6->AddElement(Zn, fractionmass=(5.1  +  6.1)*0.5/100);

  G4Material* AL_6061_T6 = new G4Material("AL_6061_T6", 2.7*g/cm3, ncomponents = 9);//http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma6061t6
  AL_6061_T6->AddElement(Al, fractionmass=(95.8 + 98.6)*0.5/100);
  AL_6061_T6->AddElement(Cr, fractionmass=(0.04 + 0.35)*0.5/100);
  AL_6061_T6->AddElement(Cu, fractionmass=(0.15 +  0.4)*0.5/100);
  AL_6061_T6->AddElement(Fe, fractionmass=0.73*(7./(7.+1.5+1.5+2.5))/100);
  AL_6061_T6->AddElement(Mg, fractionmass=(0.8  +  1.2)*0.5/100);
  AL_6061_T6->AddElement(Mn, fractionmass=0.73*(1.5/(7.+1.5+1.5+2.5))/100);
  AL_6061_T6->AddElement(Si, fractionmass=(0.4  +  0.8)*0.5/100);
  AL_6061_T6->AddElement(Ti, fractionmass=0.73*(1.5/(7.+1.5+1.5+2.5))/100);
  AL_6061_T6->AddElement(Zn, fractionmass=0.73*(2.5/(7.+1.5+1.5+2.5))/100);

  G4Material* AL_6063_T52 = new G4Material("AL_6063_T52", 2.7*g/cm3, ncomponents = 12);//http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T5
  //datat from 6063-T5//6063 depend greatly on the temper(T), heat treatment of the material
  AL_6063_T52->AddElement(Al, fractionmass=97.5/100);
  AL_6063_T52->AddElement(Cr, fractionmass=0.1/100);
  AL_6063_T52->AddElement(Cu, fractionmass=0.1/100);
  AL_6063_T52->AddElement(Fe, fractionmass=0.35/100);
  AL_6063_T52->AddElement(Mg, fractionmass=0.9/100);
  AL_6063_T52->AddElement(Mn, fractionmass=0.1/100);
  AL_6063_T52->AddElement(Si, fractionmass=0.6/100);
  AL_6063_T52->AddElement(Ti, fractionmass=0.1/100);
  AL_6063_T52->AddElement(Zn, fractionmass=0.1/100);
  //it is not 100%! 99.85%
  //these P, S, Si are to make AL_6063_T52 100%. No reference. I made it myself
  AL_6063_T52->AddElement(P,  fractionmass=0.05/100);
  AL_6063_T52->AddElement(S,  fractionmass=0.05/100);
  AL_6063_T52->AddElement(Ni, fractionmass=0.05/100);

  G4Material* STEEL_1010 = new G4Material("STEEL_1010", 7.87*g/cm3, ncomponents = 5);//https://www.azom.com/article.aspx?ArticleID=6539
  STEEL_1010->AddElement(Fe, fractionmass=(99.18 + 99.62)*0.5/100);
  STEEL_1010->AddElement(Mn, fractionmass=(0.30  +  0.60)*0.5/100);
  STEEL_1010->AddElement(C,  fractionmass=(0.08  +  0.13)*0.5/100);
  STEEL_1010->AddElement(S,  fractionmass=0.045*(5./9.)/100);
  STEEL_1010->AddElement(P,  fractionmass=0.045*(4./9.)/100);

  G4Material* SS_18_8 = new G4Material("SS_18_8", 8.*g/cm3, ncomponents = 8);//http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mq304a
  SS_18_8->AddElement(C,  fractionmass=1.5775*(0.08/(0.08+2+0.045+0.03+1))/100);
  SS_18_8->AddElement(Cr, fractionmass=(18     +   20)*0.5/100);
  SS_18_8->AddElement(Fe, fractionmass=(66.345 +   74)*0.5/100);
  SS_18_8->AddElement(Mn, fractionmass=1.5775*(2/(0.08+2+0.045+0.03+1))/100);
  SS_18_8->AddElement(Ni, fractionmass=(8      + 10.5)*0.5/100);
  SS_18_8->AddElement(P,  fractionmass=1.5775*(0.045/(0.08+2+0.045+0.03+1))/100);
  SS_18_8->AddElement(S,  fractionmass=1.5775*(0.03/(0.08+2+0.045+0.03+1))/100);
  SS_18_8->AddElement(Si, fractionmass=1.5775*(1/(0.08+2+0.045+0.03+1))/100);

  fWindowInnerJoint1Mater =    AL_7075_T6;
  fWindowInnerJoint2Mater =    AL_6061_T6;
  //  fWindowInnerJoint1Mater = no_use;
  fWindowOuterJoint2Mater =    AL_6061_T6;
  fWindowOuterJoint2_1Mater =  AL_6061_T6;
  fWindowOuterJoint2_2Mater =  AL_6061_T6;
  fWindowOuterJoint2_3Mater =  AL_6061_T6;
  fBeampipe2Mater =            AL_6061_T6;
  fBeampipe3Mater =            AL_6061_T6;
  fSweepingMagnet_6_234Mater = AL_6061_T6;
  fBeampipe1Mater =            AL_6063_T52;
  fSweepingMagnet_1Mater =     STEEL_1010;
  fSweepingMagnet_2Mater =     STEEL_1010;
  fSweepingMagnet_3Mater =     STEEL_1010;
  fSweepingMagnet_4Mater =     STEEL_1010;
  fSweepingMagnet_5Mater =     STEEL_1010;
  fSweepingMagnet_6_1Mater =   SS_18_8;

}

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  //
  //optical photons properties
  // 
  G4double hc = 1.23984193;
  
  const G4int n = 8;
  G4double OpticalPhotonWavelength[n] = //in micro meters
    { 0.400, 0.440, 0.480, 0.520, 0.560, 0.600, 0.640, 0.680};
  G4double OpticalPhotonEnergy[n] = {0.}; 
  for (G4int i = 0 ; i < n ; i++){
    OpticalPhotonEnergy[i] = (hc/OpticalPhotonWavelength[i])*eV;
  }

  //Refractive Index of Crystal
  G4double RefractiveIndexCrystal[n] = 
    { 2.35, 2.30, 2.27, 2.25, 2.23, 2.21, 2.201, 2.2}; 

  //Refractive Index of Air
  G4double RefractiveIndexAir[n] = 
    { 1., 1., 1., 1., 1., 1., 1., 1.};

  //Reflectivity. To be filled from down below
  G4double R[n] = {0.};

  //Theoretical Transmittance. To be filled from down below
  G4double Ts[n] = {0.};

  //Measured Transmittance(longitudinal)
  G4double T[n] = 
    { 0.33, 0.48, 0.62, 0.67, 0.68, 0.689, 0.69, 0.69};

  //Attenuation length
  G4double LAL[n] = {0.};


  for (G4int i = 0 ; i < n ; i++){
    R[i] = ((RefractiveIndexCrystal[i] - RefractiveIndexAir[i])*(RefractiveIndexCrystal[i] - RefractiveIndexAir[i]))/((RefractiveIndexCrystal[i] + RefractiveIndexAir[i])*(RefractiveIndexCrystal[i] + RefractiveIndexAir[i]));
    Ts[i] = (1 - R[i])/(1 + R[i]); 
    LAL[i] = fCrystal_Z/(log((T[i]*(1-Ts[i])*(1-Ts[i]))/(sqrt(4*Ts[i]*Ts[i]*Ts[i]*Ts[i] + T[i]*T[i]*(1-Ts[i]*Ts[i])*(1-Ts[i]*Ts[i]))-2*Ts[i]*Ts[i])));
    G4cout<<"n[i]"<<RefractiveIndexCrystal[i]<<"["<<i<<"], "<<"Wavelength[i] : "<<OpticalPhotonWavelength[i]<<"["<<i<<"], "<<"R[i] : "<<R[i]<<"["<<i<<"], "<<"T[i] : "<<T[i]<<"["<<i<<"], "<<"LAL[i] : "<<LAL[i]<<"["<<i<<"]"<<G4endl;  
  }

  G4double ScintilFast[n] =
    {10., 25., 45., 55., 40., 35., 20., 12.};

  G4MaterialPropertiesTable* CrystalOP = new G4MaterialPropertiesTable();
  CrystalOP->AddProperty("RINDEX",       OpticalPhotonEnergy, RefractiveIndexCrystal,n);
  CrystalOP->AddProperty("ABSLENGTH",    OpticalPhotonEnergy, LAL,     n);
  CrystalOP->AddProperty("FASTCOMPONENT",OpticalPhotonEnergy, ScintilFast,     n);
  //slow component is also filled with fast component
  CrystalOP->AddProperty("SLOWCOMPONENT",OpticalPhotonEnergy, ScintilFast,     n);
  CrystalOP->AddConstProperty("SCINTILLATIONYIELD",0.455*0.5*3.34*15*4.19673/MeV);//to get 15 p.e./MeV
  CrystalOP->AddConstProperty("RESOLUTIONSCALE",1.0);  
  CrystalOP->AddConstProperty("FASTTIMECONSTANT", 13.26*ns);
  CrystalOP->AddConstProperty("SLOWTIMECONSTANT", 412.2*ns);
  CrystalOP->AddConstProperty("YIELDRATIO",0.85);
  CrystalOP->DumpTable();

  fDetectorMater->SetMaterialPropertiesTable(CrystalOP);

  G4double RefractiveIndexPMT[n] =
    {1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50};

  G4double AbsorptionPMT[n] =
    { 1.*m, 1.*m, 1.*m, 1.*m, 1.*m, 1.*m, 1.*m, 1.*m};
  
  G4MaterialPropertiesTable* pmtOP = new G4MaterialPropertiesTable();
  pmtOP->AddProperty("RINDEX",       OpticalPhotonEnergy, RefractiveIndexPMT,n);
  pmtOP->AddProperty("ABSLENGTH",    OpticalPhotonEnergy, AbsorptionPMT,     n);
  pmtOP->DumpTable();

  fPMTmater->SetMaterialPropertiesTable(pmtOP);

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  //controlling the sizes  
  fSingle_X = fCrystal_X + 2*fWrapThickness + gap;
  fSingle_Y = fCrystal_Y + 2*fWrapThickness + gap;
  fSingle_Z = fCrystal_Z + fPMT_length + fPMTcoverThickness + fWrapThickness;

  //to make the NPS' volume as compact as possible
  fTemp_X = 30*fSingle_X;
  fTemp_Y = 36*fSingle_Y;
  fTemp_Z = fSingle_Z;
  fMom_X = fTemp_X + .1*mm;
  fMom_Y = fTemp_Y + .1*mm;
  fMom_Z = fTemp_Z + .1*mm;

  fCrystal_pos_X = 0*mm;
  fCrystal_pos_Y = 0*mm;
  fCrystal_pos_Z = -(0.5*fSingle_Z - (0.5*fCrystal_Z + fWrapThickness));
  fPMT_pos_X = 0*mm;//position of the PMT inside mother volume#5(single)
  fPMT_pos_Y = 0*mm;
  fPMT_pos_Z = 0.5*fSingle_Z - (0.5*fPMT_length + fPMTcoverThickness);

  fBulk_X = fSingle_X;
  fBulk_Y = fSingle_Y*36;//36 singles combined into one column (totally, 30 of them exists)
  fBulk_Z = fSingle_Z;

  fNPS_X = fBulk_X*30;//now 30 of the colums become one
  fNPS_Y = fBulk_Y;
  fNPS_Z = fBulk_Z;

  fMom_pos_Z = fNPS_distance*mm + 0.5*fSingle_Z;//the position of front face wrapper of the crystal as set distance of the calorimeter from the target.

  fCheckOverlaps = false;//true;//activate checing overlaps

  G4Box* sWorld =
    new G4Box("World_sol",0.5*fWorld_X, 0.5*fWorld_Y, 0.5*fWorld_Z); //its size
 			   
  lWorld =                         
    new G4LogicalVolume(sWorld,          //its solid
			fWorldMater,         //its material(air)
			"World_log");            //its name
                                   
  fPhysiWorld = 
    new G4PVPlacement(0,                     //no rotation
		      G4ThreeVector(),       //at (0,0,0)
		      lWorld,            //its logical volume
		      "World_pos",               //its name
		      0,                     //its mother  volume
		      false,                 //no boolean operation
		      0,                     //copy number
		      fCheckOverlaps);       // checking overlaps 

  //
  //Scattering Chamber
  //
  //HAll C
  G4RotationMatrix *xChambRot = new G4RotationMatrix;  // Rotates Y and Z axes only
  xChambRot->rotateX(90*degree);                     // Rotate 90 degrees

  //
  //Chamber
  //
  G4Tubs*
    sChamberOuter = new G4Tubs("ChamberOuter_sol",
			       fChamberInnerRadius,//inner radius
			       fChamberOuterRadius,//outer radius
			       0.5*fChamberHight,//hight
			       0.,//starting angle
			       twopi);//ending angle

  G4double fWindowStartTheta = (3.+10.+90.+180.)*pi/180.;//window beginning position
  G4double fWindowDeltaTheta = 124.*pi/180.;//window size 127deg
  G4double fWindowFrameStartTheta =  fWindowStartTheta-(4.5)*pi/180.;//to match the center of the window
  G4double fWindowFrameDeltaTheta = fWindowDeltaTheta+9*pi/180.;//total 133 deg.

  G4Tubs*
    sWindowSub = new G4Tubs("WindowSub_sol",
			    fChamberInnerRadius-1*cm
			    ,
			    fChamberOuterRadius+1*cm
			    ,
			    0.5*fWindowSubHight,
			    fWindowStartTheta,			    
			    fWindowDeltaTheta);

  G4RotationMatrix *zWindowRot = new G4RotationMatrix;  // Rotates Y and Z axes only
  zWindowRot->rotateZ(90*degree);                     // Rotate 90 degrees
  G4SubtractionSolid* sChamber_sub_Window = new G4SubtractionSolid("Chamber_sub_Window",
								   sChamberOuter,
								   sWindowSub,
								   zWindowRot,
								   G4ThreeVector());
									  
  fLogicChamber = new G4LogicalVolume(sChamber_sub_Window,//shape			
				      fChamberMater,//material
				      "ChamberOuter_log");

  new G4PVPlacement(xChambRot,
  		    G4ThreeVector(),
  		    fLogicChamber,
  		    "ChamberOuter_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  G4Tubs*
    sWindowFrame_before_sub = new G4Tubs("WindowFrame_before_sub_sol",
					 fChamberOuterRadius,
					 fChamberOuterRadius + fWindowFrameThickness,//1.25inch is the thickness
					 0.5*fWindowFrameHight,//20inch is the hight
					 fWindowFrameStartTheta,//to match the center of the window
					 fWindowFrameDeltaTheta);//total 133 deg.

  G4Tubs*
    sWindowFrame_sub = new G4Tubs("WindowFrame_sub_sol",
				  fChamberOuterRadius - 1*cm,
				  fChamberOuterRadius + fWindowFrameThickness + 1*cm,//1.25inch is the thickness of the window frame
				  0.5*fWindowSubHight,
				  fWindowStartTheta,
				  fWindowDeltaTheta);

  G4RotationMatrix *pseudoWindowRot = new G4RotationMatrix;  // 
  pseudoWindowRot->rotateZ(0*degree);                     // 

  G4SubtractionSolid* sWindowFrame = new G4SubtractionSolid("WindowFrame_sol",
							    sWindowFrame_before_sub,
							    sWindowFrame_sub,
							    pseudoWindowRot,
							    G4ThreeVector());

  fLogicWindowFrame = new G4LogicalVolume(sWindowFrame,
					  fChamberMater,
					  "WindowFrame_log");


  G4RotationMatrix *zxWindowRot = new G4RotationMatrix(0,-90*degree,-90*degree);
  new G4PVPlacement(zxWindowRot,
  		    G4ThreeVector(),
  		    fLogicWindowFrame,
  		    "WindowFrame_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);


  G4Tubs*
    sWindowClamp_before_sub = new G4Tubs("WindowClamp_before_sub_sol",
  					 fChamberOuterRadius + fWindowFrameThickness + fWindowThickness,//1.25inch is the thickness of the window frame
  					 fChamberOuterRadius + fWindowFrameThickness + fWindowThickness + fWindowClampThickness,//0.75inch is the thickness of the window clamp
  					 0.5*fWindowClampHight,//20inch is the hight
  					 fWindowFrameStartTheta,//to match the center of the window
  					 fWindowFrameDeltaTheta);//total 133 deg.

  G4Tubs*
    sWindowClamp_sub = new G4Tubs("WindowClamp_sub_sol",
				  fChamberOuterRadius + fWindowFrameThickness + fWindowThickness -1*cm,
				  fChamberOuterRadius + fWindowFrameThickness + fWindowThickness + fWindowClampThickness +1*cm,//1.25inch is the thickness
				  0.5*fWindowSubHight,
				  fWindowStartTheta,
				  fWindowDeltaTheta);

  G4SubtractionSolid* sWindowClamp = new G4SubtractionSolid("WindowClamp_sol",
  							    sWindowClamp_before_sub,
  							    sWindowClamp_sub,
  							    pseudoWindowRot,
  							    G4ThreeVector());


  fLogicWindowClamp = new G4LogicalVolume(sWindowClamp,
					  fChamberMater,
					  "WindowClamp_log");

  new G4PVPlacement(zxWindowRot,
  		    G4ThreeVector(),
  		    fLogicWindowClamp,
  		    "WindowClamp_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  //Vacuum inside the chamber
  G4Tubs*
    sChamberInner = new G4Tubs("ChamberInner_sol",
			       0.,//inner radius
			       fChamberInnerRadius,//outer radius
			       0.5*fChamberHight,//hight
			       0.,//starting angle
			       twopi);//ending angle

  fLogicInnerChamber = new G4LogicalVolume(sChamberInner,//shape
					   fVacuumMater,//material
					   "ChamberInner_log");

  new G4PVPlacement(xChambRot,
  		    G4ThreeVector(),
  		    fLogicInnerChamber,
  		    "ChamberInner_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  // Target Cover & Target
  //Target Cover //cylinder
  G4Tubs* 
    sTargetCover = new G4Tubs("TargetCover_sol",                                   //name
			      0., fTargetCoverRadius, 0.5*fTargetCoverLength, 0.,twopi); //dimensions


  fLogicTargetCover = new G4LogicalVolume(sTargetCover,           //shape
					  fTargetCoverMater,              //material
					  "TargetCover_log");                 //name

  G4RotationMatrix *xTargetRot = new G4RotationMatrix;  // Rotates Y and Z axes only
  xTargetRot->rotateX(-90*degree);                     // Rotate -90 degrees to reposition the target inside the chamber

  new G4PVPlacement(xTargetRot,
  		    G4ThreeVector(0.*mm, 0.*mm, 0.*mm),
  		    fLogicTargetCover,            //its logical volume
  		    "TargetCover_pos",               //its name
  		    fLogicInnerChamber,//its mother  volume
  		    false,                 //no boolean operation
  		    0,                     //copy number
  		    fCheckOverlaps);       // checking overlaps 

  //Target Cover //half sphere

  G4Sphere*
    sTargetCoverSphere = new G4Sphere("TargetCoverSphere_sol",
				      0.,//radius_min
				      fTargetCoverRadius,//radius_max
				      0.,//starting Phi angle
				      twopi,//ending Phi angle
				      0.,//starting Theta angle
				      0.5*pi);//ending Theta angle

  G4LogicalVolume* lTargetCoverSphere = new G4LogicalVolume(sTargetCoverSphere,
							    fTargetCoverMater,
							    "TargetCoverSphere_log");

  new G4PVPlacement(xTargetRot,
  		    G4ThreeVector(0.*mm, -0.5*fTargetCoverLength, 0.*mm ),
  		    lTargetCoverSphere,            //its logical volume
  		    "TargetCoverSphere_pos",               //its name
  		    fLogicInnerChamber,//its mother  volume
  		    false,                 //no boolean operation
  		    0,                     //copy number
  		    fCheckOverlaps);       // checking overlaps 

  //Target //cylinder
  G4Tubs* 
    sTarget = new G4Tubs("Target_sol",                                   //name
			 0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
				     fTargetMater,              //material
				     "Target_log");                 //name
  
  new G4PVPlacement(0,                     //no rotation
		    G4ThreeVector(0.*mm, 0.*mm, 0.5*fTargetWindowThickness),
		    fLogicTarget,            //its logical volume
		    "Target_pos",               //its name
		    fLogicTargetCover,                     //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    fCheckOverlaps);       // checking overlaps 

  //Target //half sphere

  
  G4Sphere*
    sTargetSphere = new G4Sphere("TargetSphere_sol",
				 0.,//radius_min
				 fTargetRadius,//radius_max
				 0.,//starting Phi angle
				 twopi,//ending Phi angle
				 0.,//starting Theta angle
				 0.5*pi);//ending Theta angle
  
  G4LogicalVolume* lTargetSphere = new G4LogicalVolume(sTargetSphere,
						       fTargetMater,
						       "TargetSphere_log");

  new G4PVPlacement(0,                     //no rotation
		    G4ThreeVector(0.*mm, 0.*mm, 0*mm),       //at (0,0,0)
		    lTargetSphere,            //its logical volume
		    "TargetSphere_pos",               //its name
		    lTargetCoverSphere,                     //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    fCheckOverlaps);       // checking overlaps 

  //
  //beamline
  //window joint
  //
  G4Tubs*
    sWindowInnerJoint1 = new G4Tubs("WindowInnerJoint1_sol",
				    fWindowOuterJoint2OuterRadius,//fWindowInnerJoint1InnerRadius,//inner radius//it is(will be) placed inside the vacuum
				    fWindowInnerJoint1OuterRadius,//outer radius
				    0.5*fWindowInnerJoint1Thickness,//hight. lenght
				    0.,//starting angle
				    twopi);//ending angle

  fLogicWindowInnerJoint1 = new G4LogicalVolume(sWindowInnerJoint1,
						fWindowInnerJoint1Mater,
						"WindowInnerJoint1_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness - fWindowThickness -fWindowInnerJoint2Thickness - 0.5*fWindowInnerJoint1Thickness),
  		    fLogicWindowInnerJoint1,
  		    "WindowInnerJoint1_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  G4Tubs*
    sWindowInnerJoint2 = new G4Tubs("WindowInnerJoint2_sol",
				    fWindowInnerJoint2InnerRadius,//inner radius//it is(will be) placed inside the vacuum
				    fWindowInnerJoint2OuterRadius,//outer radius
				    0.5*fWindowInnerJoint2Thickness,//hight. lenght
				    0.,//starting angle
				    twopi);//ending angle

  fLogicWindowInnerJoint2 = new G4LogicalVolume(sWindowInnerJoint2,
						fWindowInnerJoint2Mater,
						"WindowInnerJoint2_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness - fWindowThickness - 0.5*fWindowInnerJoint2Thickness),
  		    fLogicWindowInnerJoint2,
  		    "WindowInnerJoint2_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //no need
  // G4Tubs*
  //   sWindowOuterJoint1 = new G4Tubs("WindowOuterJoint1_sol",
  // 				    fWindowOuterJoint1InnerRadius,//inner radius
  // 				    fWindowOuterJoint1OuterRadius,//outer radius
  // 				    0.5*fWindowOuterJoint1Thickness,//hight. lenght
  // 				    0.,//starting angle
  // 				    twopi);//ending angle

  // fLogicWindowOuterJoint1 = new G4LogicalVolume(sWindowOuterJoint1,
  // 						fWindowOuterJoint1Mater,
  // 						"WindowOuterJoint1_log");


  // new G4PVPlacement(0,
  // 		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowThickness + 0.5*fWindowOuterJoint1Thickness),
  // 		    fLogicWindowOuterJoint1,
  // 		    "WindowOuterJoint1_pos",
  // 		    lWorld,
  // 		    false,
  // 		    0,
  // 		    fCheckOverlaps);

  G4Tubs*
    sWindowOuterJoint2 = new G4Tubs("WindowOuterJoint2_sol",
				    0.,//inner radius//will be filled with vacuum
				    fWindowOuterJoint2OuterRadius,//outer radius
				    0.5*fWindowOuterJoint2Thickness,//hight. lenght
				    0.,//starting angle
				    twopi);//ending angle

  fLogicWindowOuterJoint2 = new G4LogicalVolume(sWindowOuterJoint2,
						fWindowOuterJoint2Mater,
						"WindowOuterJoint2_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + 0.5*fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position),
  		    fLogicWindowOuterJoint2,
  		    "WindowOuterJoint2_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //vacuum inside the OuterJoint2
  //
  G4Tubs*
    sWindowOuterJoint2Vacuum = new G4Tubs("WindowOuterJoint2Vacuum_sol",
					  0.,//inner radius//will be filled with vacuum
					  fWindowOuterJoint2InnerRadius,//outer radius
					  0.5*fWindowOuterJoint2Thickness,//hight. lenght
					  0.,//starting angle
					  twopi);//ending angle

  G4LogicalVolume* fLogicWindowOuterJoint2Vacuum = new G4LogicalVolume(sWindowOuterJoint2Vacuum,
								       fVacuumMater,
								       "WindowOuterJoint2Vacuum_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,0),
  		    fLogicWindowOuterJoint2Vacuum,
  		    "WindowOuterJoint2Vacuum_pos",
  		    fLogicWindowOuterJoint2,
  		    false,
  		    0,
  		    fCheckOverlaps);

  G4Tubs*
    sWindowOuterJoint2_1 = new G4Tubs("WindowOuterJoint2_1_sol",
				      fWindowOuterJoint2_1InnerRadius,//inner radius//will be filled with vacuum
				      fWindowOuterJoint2_1OuterRadius,//outer radius
				      0.5*fWindowOuterJoint2_1Thickness,//hight. lenght
				      0.,//starting angle
				      twopi);//ending angle

  fLogicWindowOuterJoint2_1 = new G4LogicalVolume(sWindowOuterJoint2_1,
						  fWindowOuterJoint2_1Mater,
						  "WindowOuterJoint2_1_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + 0.5*fWindowOuterJoint2_1Thickness),
  		    fLogicWindowOuterJoint2_1,
  		    "WindowOuterJoint2_1_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  G4RotationMatrix *zPipeRot = new G4RotationMatrix;  // Rotates Y and Z axes only
  zPipeRot->rotateZ(45*degree);                     // Rotate 90 degrees

  G4Box*
    sWindowOuterJoint2_2 = new G4Box("WindowOuterJoint2_2_sol",
				     0.5*fWindowOuterJoint2_2dx,
				     0.5*fWindowOuterJoint2_2dy,
				     0.5*fWindowOuterJoint2_2dz);

  fLogicWindowOuterJoint2_2 = new G4LogicalVolume(sWindowOuterJoint2_2,
						  fWindowOuterJoint2_2Mater,
						  "WindowOuterJoint2_2_log");

  new G4PVPlacement(zPipeRot,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + 0.5*fWindowOuterJoint2_2dz),
  		    fLogicWindowOuterJoint2_2,
  		    "WindowOuterJoint2_2_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //vacuum inside the OuterJoint2_2
  //
  G4Tubs*
    sWindowOuterJoint2_2Vacuum = new G4Tubs("WindowOuterJoint2_2Vacuum_sol",
					    0.,//inner radius//will be filled with vacuum
					    fWindowOuterJoint2InnerRadius,//outer radius
					    0.5*fWindowOuterJoint2_2dz,//hight. lenght
					    0.,//starting angle
					    twopi);//ending angle

  G4LogicalVolume* fLogicWindowOuterJoint2_2Vacuum = new G4LogicalVolume(sWindowOuterJoint2_2Vacuum,
									 fVacuumMater,
									 "WindowOuterJoint2_2Vacuum_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,0),
  		    fLogicWindowOuterJoint2_2Vacuum,
  		    "WindowOuterJoint2_2Vacuum_pos",
  		    fLogicWindowOuterJoint2_2,
  		    false,
  		    0,
  		    fCheckOverlaps);

  G4Tubs*
    sWindowOuterJoint2_3 = new G4Tubs("WindowOuterJoint2_3_sol",
				      0.,//inner radius//will be filled with vacuum
				      fWindowOuterJoint2_3OuterRadius,//outer radius
				      0.5*fWindowOuterJoint2_3Thickness,//hight. lenght
				      0.,//starting angle
				      twopi);//ending angle

  fLogicWindowOuterJoint2_3 = new G4LogicalVolume(sWindowOuterJoint2_3,
						  fWindowOuterJoint2_3Mater,
						  "WindowOuterJoint2_3_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz + 0.5*fWindowOuterJoint2_3Thickness),
  		    fLogicWindowOuterJoint2_3,
  		    "WindowOuterJoint2_3_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //vacuum inside the OuterJoint2_3
  //
  G4Tubs*
    sWindowOuterJoint2_3Vacuum = new G4Tubs("WindowOuterJoint2_3Vacuum_sol",
					    0.,//inner radius//will be filled with vacuum
					    fWindowOuterJoint2_3InnerRadius,//outer radius
					    0.5*fWindowOuterJoint2_3Thickness,//hight. lenght
					    0.,//starting angle
					    twopi);//ending angle

  G4LogicalVolume* fLogicWindowOuterJoint2_3Vacuum = new G4LogicalVolume(sWindowOuterJoint2_3Vacuum,
									 fVacuumMater,
									 "WindowOuterJoint2_3Vacuum_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,0),
  		    fLogicWindowOuterJoint2_3Vacuum,
  		    "WindowOuterJoint2_3Vacuum_pos",
  		    fLogicWindowOuterJoint2_3,
  		    false,
  		    0,
  		    fCheckOverlaps);
  // G4Tubs*
  //   sWindowOuterJoint3_1 = new G4Tubs("WindowOuterJoint3_1_sol",
  // 				      fWindowOuterJoint3_1InnerRadius,//inner radius//will be filled with vacuum
  // 				      fWindowOuterJoint3_1OuterRadius,//outer radius
  // 				      0.5*fWindowOuterJoint3_1Thickness,//hight. lenght
  // 				      0.,//starting angle
  // 				      twopi);//ending angle

  // fLogicWindowOuterJoint3_1 = new G4LogicalVolume(sWindowOuterJoint3_1,
  // 						  fWindowOuterJoint3_1Mater,
  // 						  "WindowOuterJoint3_1_log");

  // new G4PVPlacement(0,
  // 		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness /* + fWindowThickness */ + fWindowOuterJoint1Thickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint3Thickness - 0.5*fWindowOuterJoint3_1Thickness),
  // 		    fLogicWindowOuterJoint3_1,
  // 		    "WindowOuterJoint3_1_pos",
  // 		    lWorld,
  // 		    false,
  // 		    0,
  // 		    fCheckOverlaps);

  //
  //Beampipe
  //
  //Beampipe1
  //
  G4Trd*
    sBeampipe1 = new G4Trd("Beampipe1_sol",
			   0.5*fBeampipe1Outerdx1,//half dx1
			   0.5*fBeampipe1Outerdx2,//half dx2
			   0.5*fBeampipe1Outerdy1,//half dy1
			   0.5*fBeampipe1Outerdy2,//half dy2
			   0.5*fBeampipe1Length);//half dz

  fLogicBeampipe1 = new G4LogicalVolume(sBeampipe1,
					fBeampipe1Mater,
					"Beampipe1_log");

  new G4PVPlacement(zPipeRot,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz  + 0.5*fBeampipe1Length),
  		    fLogicBeampipe1,
  		    "Beampipe1_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);
  //
  //vacuum inside the Beampipe1
  //
  G4Trd*
    sBeampipe1Vacuum = new G4Trd("Beampipe1Vacuum_sol",
				 0.5*fBeampipe1Innerdx1,//half dx1
				 0.5*fBeampipe1Innerdx2,//half dx2
				 0.5*fBeampipe1Innerdy1,//half dy1
				 0.5*fBeampipe1Innerdy2,//half dy2
				 0.5*fBeampipe1Length);//half dz

  G4LogicalVolume* fLogicBeampipe1Vacuum = new G4LogicalVolume(sBeampipe1Vacuum,
							       fVacuumMater,
							       "Beampipe1Vacuum_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,0),
  		    fLogicBeampipe1Vacuum,
  		    "Beampipe1Vacuum_pos",
  		    fLogicBeampipe1,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  //Beampipe2
  //
  G4Tubs*
    sBeampipe2 = new G4Tubs("Beampipe_sol",
			    0.,//fWindowOuterJoint3InnerRadius,//inner radius//will be filled with vacuum
			    fBeampipe2OuterRadius,//outer radius
			    0.5*fBeampipe2Length,//hight. lenght
			    0.,//starting angle
			    twopi);//ending angle

  fLogicBeampipe2 = new G4LogicalVolume(sBeampipe2,
					fBeampipe2Mater,
					"Beampipe2_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz + fBeampipe1Length + 0.5*fBeampipe2Length),
  		    fLogicBeampipe2,
  		    "Beampipe2_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);
  //
  //vacuum inside the Beampipe2
  //
  G4Tubs*
    sBeampipe2Vacuum = new G4Tubs("Beampipe2Vacuum_sol",
				  0.,//fWindowOuterJoint3InnerRadius,//inner radius//will be filled with vacuum
				  fBeampipe2InnerRadius,//outer radius
				  0.5*fBeampipe2Length,//hight. lenght
				  0.,//starting angle
				  twopi);//ending angle

  G4LogicalVolume* fLogicBeampipe2Vacuum = new G4LogicalVolume(sBeampipe2Vacuum,
							       fVacuumMater,
							       "Beampipe2Vacuum_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,0),
  		    fLogicBeampipe2Vacuum,
  		    "Beampipe2Vacuum_pos",
  		    fLogicBeampipe2,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  //Beampipe2 front cover
  //
  G4Tubs*
    sBeampipe2FrontCover_before_sub = new G4Tubs("Beampipe2FrontCover_before_sub_sol",
						 0.,//inner radius
						 fBeampipe2OuterRadius,//outer radius
						 0.5*fBeampipe2FrontCoverThickness,//length
						 0.,//starting angle
						 twopi);//ending angle


  G4SubtractionSolid* sBeampipe2FrontCover = new G4SubtractionSolid("Beampipe2FrontCover_sol",
								    sBeampipe2FrontCover_before_sub,
								    sBeampipe1,
								    zPipeRot,
								    G4ThreeVector(0, 0, -0.5*fBeampipe1Length + 0.5*fBeampipe2FrontCoverThickness));

  fLogicBeampipe2FrontCover = new G4LogicalVolume(sBeampipe2FrontCover,
						  fBeampipe2Mater,
						  "Beampipe2FrontCover_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz + fBeampipe1Length - 0.5*fBeampipe2FrontCoverThickness),
  		    fLogicBeampipe2FrontCover,
  		    "Beampipe2FrontCover_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);
  //
  //Beampipe3
  //
  G4Tubs*
    sBeampipe3 = new G4Tubs("Beampipe_sol",
			    0.,//fWindowOuterJoint3InnerRadius,//inner radius//will be filled with vacuum
			    fBeampipe3OuterRadius,//outer radius
			    0.5*fBeampipe3Length,//hight. lenght
			    0.,//starting angle
			    twopi);//ending angle

  fLogicBeampipe3 = new G4LogicalVolume(sBeampipe3,
					fBeampipe3Mater,
					"Beampipe3_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz + fBeampipe1Length + fBeampipe2Length + 0.5*fBeampipe3Length),
  		    fLogicBeampipe3,
  		    "Beampipe3_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);
  //
  //vacuum inside the Beampipe3
  //
  G4Tubs*
    sBeampipe3Vacuum = new G4Tubs("Beampipe3Vacuum_sol",
				  0.,//fWindowOuterJoint3InnerRadius,//inner radius//will be filled with vacuum
				  fBeampipe3InnerRadius,//outer radius
				  0.5*fBeampipe3Length,//hight. lenght
				  0.,//starting angle
				  twopi);//ending angle

  G4LogicalVolume* fLogicBeampipe3Vacuum = new G4LogicalVolume(sBeampipe3Vacuum,
							       fVacuumMater,
							       "Beampipe3Vacuum_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,0),
  		    fLogicBeampipe3Vacuum,
  		    "Beampipe3Vacuum_pos",
  		    fLogicBeampipe3,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  //Beampipe3 front cover
  //
  G4Tubs*
    sBeampipe3FrontCover = new G4Tubs("Beampipe3FrontCover_before_sub_sol",
				      fBeampipe2OuterRadius,//inner radius
				      fBeampipe3OuterRadius,//outer radius
				      0.5*fBeampipe3FrontCoverThickness,//length
				      0.,//starting angle
				      twopi);//ending angle

  fLogicBeampipe3FrontCover = new G4LogicalVolume(sBeampipe3FrontCover,
						  fBeampipe3Mater,
						  "Beampipe3FrontCover_log");

  new G4PVPlacement(0,
  		    G4ThreeVector(0,0,fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz + fBeampipe1Length + fBeampipe2Length - 0.5*fBeampipe3FrontCoverThickness),
  		    fLogicBeampipe3FrontCover,
  		    "Beampipe3FrontCover_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  //Al window for the chamber
  //
  G4Tubs*
    sWindow = new G4Tubs("Window_sol",
  			 fChamberOuterRadius + fWindowFrameThickness,//1.25inch is the thickness of the window frame
  			 fChamberOuterRadius + fWindowFrameThickness + fWindowThickness,
  			 0.5*fWindowHight,//actual size is 19*inch however, they are covered by the window clamp
			 fWindowFrameStartTheta,//to match the center of the window
			 fWindowFrameDeltaTheta);//total 133 deg.

  G4RotationMatrix *yRot = new G4RotationMatrix; 
  yRot->rotateY(-90*degree);                    

  G4SubtractionSolid* sWindow_sub_InnerJoint2 = new G4SubtractionSolid("Window_sub_InnerJoint2",
  								       sWindow,
  								       sWindowInnerJoint2,
  								       yRot,
  								       G4ThreeVector(fChamberOuterRadius + fWindowFrameThickness - fWindowThickness - 0.5*fWindowInnerJoint2Thickness, 0, 0));
  G4SubtractionSolid* sWindow_sub_InnerJoint2_sub_OuterJoint2 = new G4SubtractionSolid("Window_sub_InnerJoint2_sub_OutherJoint2",
  										       sWindow_sub_InnerJoint2,
  										       sWindowOuterJoint2,
  										       yRot,
  										       G4ThreeVector(fChamberOuterRadius + fWindowFrameThickness +0.5*fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position, 0, 0));

  G4SubtractionSolid* sWindow_sub_InnerJoint2_sub_OuterJoint2_sub_OuterJoint2_1 = new G4SubtractionSolid("Window_sub_InnerJoint2_sub_OutherJoint2",
													 sWindow_sub_InnerJoint2_sub_OuterJoint2,
													 sWindowOuterJoint2_1,
													 yRot,
													 G4ThreeVector(fChamberOuterRadius + fWindowFrameThickness + 0.5*fWindowOuterJoint2_1Thickness, 0, 0));

  fLogicChamberWindow = new G4LogicalVolume(sWindow_sub_InnerJoint2_sub_OuterJoint2_sub_OuterJoint2_1,
  					    fWindowMater,
  					    "Window_log");

  new G4PVPlacement(zxWindowRot,//xChambRot,
  		    G4ThreeVector(),
  		    fLogicChamberWindow,
  		    "Window_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //
  //extra detailed vaccum in the chamber
  //
  G4Tubs*
    sChamberWindowVacuum = new G4Tubs("ChamberWindowVacuum_sol",
				      fChamberInnerRadius,
				      fChamberOuterRadius+fWindowFrameThickness,
				      0.5*fWindowSubHight,
				      fWindowStartTheta,
				      fWindowDeltaTheta);


  G4SubtractionSolid* sChamberWindowVacuum_sub_InnerJoint2 = new G4SubtractionSolid("ChamberWindowVacuum_sub_InnerJoint2",
										    sChamberWindowVacuum,
										    sWindowInnerJoint2,
										    yRot,
										    G4ThreeVector(fChamberOuterRadius + fWindowFrameThickness - fWindowThickness - 0.5*fWindowInnerJoint2Thickness, 0, 0));
  G4SubtractionSolid* sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2 = new G4SubtractionSolid("ChamberWindowVacuum_sub_InnerJoint2_sub_OutherJoint2",
												    sChamberWindowVacuum_sub_InnerJoint2,
												    sWindowOuterJoint2,
												    yRot,
												    G4ThreeVector(fChamberOuterRadius + fWindowFrameThickness + 0.5*fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position, 0, 0));

  G4SubtractionSolid* sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2_sub_InnerJoint1 = new G4SubtractionSolid("ChamberWindowVacuum_sub_InnerJoint2_sub_OutherJoint2_sub_InnerJoint1",
														    sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2,
														    sWindowInnerJoint1,
														    yRot,
														    G4ThreeVector(fChamberOuterRadius + fWindowFrameThickness - fWindowThickness - fWindowInnerJoint2Thickness - 0.5*fWindowInnerJoint1Thickness, 0, 0));



  fLogicInnerChamber2 = new G4LogicalVolume(sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2_sub_InnerJoint1,//sChamberWindowVacuum,
					    fVacuumMater,
					    "ChamberWindowVacuum_log");

  new G4PVPlacement(zxWindowRot,
  		    G4ThreeVector(),
  		    fLogicInnerChamber2,
  		    "ChamberWindowVacuum_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  // Detector
  //
  //pseudo volume to check background everywhere

  G4Sphere*
    sFlux = new G4Sphere("Flux_sol",
			 fNPS_distance - 2*fWrapThickness,//radius_min//NPS distance from the target
			 fNPS_distance - fWrapThickness,//radius_max	
			 0.,//starting Phi angle
			 twopi,//ending Phi angle
			 0.,//starting Theta angle
			 pi);//ending Theta angle


  fLogicFlux = new G4LogicalVolume(sFlux,
				   fWorldMater,//air
				   "Flux_log");

  // new G4PVPlacement(0,
  // 		    G4ThreeVector(0, 0, 0),
  // 		    fLogicFlux,
  // 		    "Flux_pos",
  // 		    lWorld,
  // 		    false,
  // 		    0,
  // 		    fCheckOverlaps);

  //From Mark Jones email(jones@jlab.org) to Carlos forwarded to me(hosanko@ipno.inp2p3.fr)
  G4double fHMS_angle = 1.*dvcsGlobals::HMS_angle;//[rad]

  G4RotationMatrix* HMS_window_rot = new G4RotationMatrix();
  HMS_window_rot->rotateY(fHMS_angle);

  vector<G4TwoVector> vertices_HMS_1_opening;//must be clockwise ordered
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (-11.646)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575)*cm,    (-11.646/2.)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575)*cm,    (11.646/2.)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (11.646)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (-11.646)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575)*cm,    (-11.646/2.)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575)*cm,    (11.646/2.)*cm ) );
  vertices_HMS_1_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (11.646)*cm ) );

  vector<G4TwoVector> vertices_HMS_2_opening;//must be clockwise ordered
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (11.646)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575)*cm,     (11.646/2.)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575)*cm,     (-11.646/2.)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (-11.646)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (11.646)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575)*cm,     (11.646/2.)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575)*cm,     (-11.646/2.)*cm ) );
  vertices_HMS_2_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (-11.646)*cm ) );

  vector<G4TwoVector> vertices_HMS_3_opening;//must be clockwise ordered
  vertices_HMS_3_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (-11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (-11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (4.575/2.)*cm,  (-11.646)*cm ) );
  vertices_HMS_3_opening.push_back( G4TwoVector( (-4.575/2.)*cm, (-11.646)*cm ) );

  vector<G4TwoVector> vertices_HMS_1_exit;//must be clockwise ordered
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (-12.114)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759)*cm,    (-12.114/2.)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759)*cm,    (12.114/2.)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (12.114)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (-12.114)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759)*cm,    (-12.114/2.)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759)*cm,    (12.114/2.)*cm ) );
  vertices_HMS_1_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (12.114)*cm ) );

  vector<G4TwoVector> vertices_HMS_2_exit;//must be clockwise ordered
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (12.114)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759)*cm,     (12.114/2.)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759)*cm,     (-12.114/2.)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (-12.114)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (12.114)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759)*cm,     (12.114/2.)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759)*cm,     (-12.114/2.)*cm ) );
  vertices_HMS_2_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (-12.114)*cm ) );

  vector<G4TwoVector> vertices_HMS_3_exit;//must be clockwise ordered
  vertices_HMS_3_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (-12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (-12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (4.759/2.)*cm,  (-12.114)*cm ) );
  vertices_HMS_3_exit.push_back( G4TwoVector( (-4.759/2.)*cm, (-12.114)*cm ) );

  G4GenericTrap* sHMS_1_opening = new G4GenericTrap("HMS_1_opening_sol", 0.5*fWrapThickness*cm, vertices_HMS_1_opening);
  G4GenericTrap* sHMS_2_opening = new G4GenericTrap("HMS_2_opening_sol", 0.5*fWrapThickness*cm, vertices_HMS_2_opening);
  G4GenericTrap* sHMS_3_opening = new G4GenericTrap("HMS_3_opening_sol", 0.5*fWrapThickness*cm, vertices_HMS_3_opening);
  G4GenericTrap* sHMS_1_exit = new G4GenericTrap("HMS_1_exit_sol", 0.5*fWrapThickness*cm, vertices_HMS_1_exit);
  G4GenericTrap* sHMS_2_exit = new G4GenericTrap("HMS_2_exit_sol", 0.5*fWrapThickness*cm, vertices_HMS_2_exit);
  G4GenericTrap* sHMS_3_exit = new G4GenericTrap("HMS_3_exit_sol", 0.5*fWrapThickness*cm, vertices_HMS_3_exit);

  G4UnionSolid* sHMS_1_3_opening = new G4UnionSolid("HMS_1_3_opening_sol", sHMS_1_opening, sHMS_3_opening);
  G4UnionSolid* sHMS_opening = new G4UnionSolid("HMS_opening_sol", sHMS_1_3_opening, sHMS_2_opening);
  G4UnionSolid* sHMS_1_3_exit = new G4UnionSolid("HMS_1_3_exit_sol", sHMS_1_exit, sHMS_3_exit);
  G4UnionSolid* sHMS_exit = new G4UnionSolid("HMS_exit_sol", sHMS_1_3_exit, sHMS_2_exit);

  G4LogicalVolume *lHMS_opening = new G4LogicalVolume(sHMS_opening,
						      fWorldMater,
						      "HMS_opening_log");
  G4LogicalVolume *lHMS_exit = new G4LogicalVolume(sHMS_exit,
						   fWorldMater,
						   "HMS_exit_log");
  new G4PVPlacement(HMS_window_rot,
  		    G4ThreeVector(166.37*cm*sin(-fHMS_angle), 0, 166.37*cm*cos(-fHMS_angle)),
  		    lHMS_opening,
  		    "HMS_opening_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  new G4PVPlacement(HMS_window_rot,
  		    G4ThreeVector((166.37 + 6.3)*cm*sin(-fHMS_angle), 0, (166.37 + 6.3)*cm*cos(-fHMS_angle)),
  		    lHMS_exit,
  		    "HMS_exit_pos",
  		    lWorld,
  		    false,
  		    0,
  		    fCheckOverlaps);

  //Mother Volume #1 : to contain NPS and its temp control box. To move & rotate it freely.
  G4Box*
    sMother = new G4Box("Mother_sol", 0.5*fMom_X, 0.5*fMom_Y, 0.5*fMom_Z);

  G4LogicalVolume* lMother = new G4LogicalVolume(sMother,
						 fWorldMater,//air
						 "Mother_log");

  G4RotationMatrix *yMomRot = new G4RotationMatrix;  // Rotates X and Z axes only
  yMomRot->rotateY(fMom_theta);  //Angle in radian                     

  new G4PVPlacement(yMomRot,                     //no rotation
  		    G4ThreeVector(fMom_pos_Z*sin(-fMom_theta), fMom_pos_Y, fMom_pos_Z*cos(-fMom_theta)),
  		    lMother,            //its logical volume
  		    "Mother_pos",               //its name
  		    lWorld,                     //its mother  volume
  		    false,                 //no boolean operation
  		    0,                     //copy number
  		    fCheckOverlaps);       // checking overlaps 

  //Temp Control Box(Mother Volume #2)
  G4Box*
    sTemp = new G4Box("Temp_sol", 0.5*fTemp_X, 0.5*fTemp_Y, 0.5*fTemp_Z);

  G4LogicalVolume* lTemp = new G4LogicalVolume(sTemp,
					       fWorldMater,//air
					       "Temp_log");
  new G4PVPlacement(0,                     //no rotation
		    G4ThreeVector(fTemp_pos_X, fTemp_pos_Y, fTemp_pos_Z),       //at (0,0,0)
		    lTemp,            //its logical volume
		    "Temp_pos",               //its name
		    lMother,                     //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    fCheckOverlaps);       // checking overlaps 

  ////////////////////////////////////////////////////////////////////////

  //Each Crystals Logical Volume
  G4Box*
    sCrystal = new G4Box("Crystal_sol", 0.5*fCrystal_X, 0.5*fCrystal_Y, 0.5*fCrystal_Z);

  fLogicCrystal = new G4LogicalVolume(sCrystal,
				      fDetectorMater,//detectormater = crystmater
				      "Crystal_log");
  
  //Each Detectors(Crysital + PMT + gap) Logical Volume (Mother volume #5)
  G4Box*
    sSingle = new G4Box("Single_sol", 0.5*fSingle_X, 0.5*fSingle_Y, 0.5*fSingle_Z);

  G4LogicalVolume* lSingle = new G4LogicalVolume(sSingle,
						 fWorldMater,//air
						 "Single_log");
  //Each columns Logical Volume (Mother volume #4) (36 crystals in y-axis become one single bulk, 30 colums in total)
  G4Box*
    sBulk = new G4Box("Bulk_sol", 0.5*fBulk_X, 0.5*fBulk_Y, 0.5*fBulk_Z);

  G4LogicalVolume* lBulk = new G4LogicalVolume(sBulk,
					       fWorldMater,//air
					       "Bulk_log");
  
  //NPS made with 30 bulks (Mother volume #3)
  G4Box* sNPS = new G4Box("Detector_sol", 0.5*fNPS_X, 0.5*fNPS_Y, 0.5*fNPS_Z);
  G4LogicalVolume* lNPS = new G4LogicalVolume(sNPS,
					      fWorldMater,//air

					      "Detector_sol");

  //Replicas
  for(G4int ly = 0 ; ly < 30 ; ly++){
    for(G4int lx = 0 ; lx < 36 ; lx++){
      new G4PVPlacement(0,
			G4ThreeVector(14.5*fSingle_X - ly*fSingle_X, -17.5*fSingle_Y + lx*fSingle_Y,0),
			lSingle,
			"Single",
			lTemp,//place the NPS inside the Temp control box
			false,
			lx + ly*36,
			fCheckOverlaps);
    }
  }
  
  //Now, position the Crystals inside each Singles
  G4PVPlacement* fCrystalPos = new G4PVPlacement(0,                     //no rotation
						 G4ThreeVector(fCrystal_pos_X, fCrystal_pos_Y, fCrystal_pos_Z),       //at (0,0,0)
						 fLogicCrystal,            //its logical volume
						 "Crystal",               //its name
						 lSingle,                     //its mother  volume
						 false,                 //no boolean operation
						 0,                     //copy number
						 fCheckOverlaps);       // checking overlaps 

  //
  //crystal wrapper
  //
  G4Box*
    sWrap = new G4Box("Wrap_sol", 0.5*fCrystal_X + fWrapThickness, 0.5*fCrystal_Y + fWrapThickness, 0.5*fCrystal_Z + fWrapThickness);

  G4ThreeVector trans(0., 0., 0.5*(fCrystal_Z + fWrapThickness));
  G4Box*
    sWrap_end = new G4Box("Wrap_end_sol", 0.5*(fCrystal_X + 2*fWrapThickness), 0.5*(fCrystal_Y + 2*fWrapThickness), 0.5*fWrapThickness);

  G4RotationMatrix* Rot = new G4RotationMatrix;
  G4SubtractionSolid* sWrap_end_sol = new G4SubtractionSolid("Wrap_sub_sol_1", sWrap, sWrap_end, Rot, trans);

  G4Box*
    sWrap_inside = new G4Box("Wrap_inside_sol", 0.5*fCrystal_X, 0.5*fCrystal_Y, 0.5*fCrystal_Z);//to hollow out the crystal wrapper         
  G4SubtractionSolid* sWrap_inside_sol = new G4SubtractionSolid("Wrap_sub_sol", sWrap_end_sol, sWrap_inside);
  G4Box*
    sWrap_front = new G4Box("Wrap_front_sol", 0.5*(fCrystal_X + 2*fWrapThickness), 0.5*(fCrystal_Y + 2*fWrapThickness), 0.5*fWrapThickness);

  G4ThreeVector trans_front(0., 0., -0.5*(fCrystal_Z + fWrapThickness));
  G4SubtractionSolid* sWrap_front_sol = new G4SubtractionSolid("Wrap_sub_sol_1", sWrap_inside_sol, sWrap_front, Rot, trans_front);

  fLogicWrap = new G4LogicalVolume(sWrap_front_sol,
				   fWrapMater,
				   "Wrap_log");

  G4PVPlacement* fWrapPos = new G4PVPlacement(0,
					      G4ThreeVector(fCrystal_pos_X, fCrystal_pos_Y, fCrystal_pos_Z),
					      fLogicWrap,
					      "Wrap_pos",
					      lSingle,
					      false,
					      0,
					      fCheckOverlaps);

  fLogicWrapFront = new G4LogicalVolume(sWrap_front,
					fWrapMater,
					"WrapFront_log");

  G4PVPlacement* fWrapFrontPos = new G4PVPlacement(0,
						   G4ThreeVector(fCrystal_pos_X, fCrystal_pos_Y, fCrystal_pos_Z) + trans_front,
						   fLogicWrapFront,
						   "WrapFront_pos",
						   lSingle,
						   false,
						   0,
						   fCheckOverlaps);


  //
  //PMT
  //
  G4Tubs*
    sPMT = new G4Tubs("PMT_sol", 0.*mm, fPMT_radius, 0.5*fPMT_length, 0, twopi);
  fLogicPMT = new G4LogicalVolume(sPMT,
                                  fPMTmater,
                                  "PMT_log");
  G4PVPlacement* fPMTpos = new G4PVPlacement(0,
                                             G4ThreeVector(fPMT_pos_X, fPMT_pos_Y, fPMT_pos_Z),
                                             fLogicPMT,
                                             "PMT_pos",
                                             lSingle,
                                             false,
                                             0,
                                             fCheckOverlaps);

  //
  //PMTcover
  //
  G4ThreeVector trans2(0., 0., -0.5*(fPMT_length + fPMTcoverThickness));
  G4Tubs*
    sPMTcover = new G4Tubs("PMTcover_sol", 0.*mm, fPMT_radius + fPMTcoverThickness, 0.5*(fPMT_length + 2*fPMTcoverThickness), 0, twopi);
  G4Tubs*
    sPMTcover_front = new G4Tubs("PMTcover_front_sol", 0.*mm, fPMT_radius + fPMTcoverThickness, 0.5*fPMTcoverThickness, 0, twopi);  //to get rid of the cover at the front side of the crystal
  G4SubtractionSolid* sPMTcover_sol_1 = new G4SubtractionSolid("PMTcover_sub_sol_1", sPMTcover, sPMTcover_front, Rot, trans2);
  G4Tubs*
    sPMTcover_inside = new G4Tubs("PMTcover_inside_sol", 0.*mm, fPMT_radius, 0.5*fPMT_length, 0, twopi);

  G4SubtractionSolid* sPMTcover_sol = new G4SubtractionSolid("Wrap_sub_sol", sPMTcover_sol_1, sPMTcover_inside);

  fLogicPMTcover = new G4LogicalVolume(sPMTcover_sol,
				       fPMTcoverMater,
				       "PMTcover_log");



  G4PVPlacement* fPMTcoverPos = new G4PVPlacement(0,
						  G4ThreeVector(fPMT_pos_X, fPMT_pos_Y, fPMT_pos_Z),
						  fLogicPMTcover,
						  "PMTcover_pos",
						  lSingle,
						  false,
						  0,
						  fCheckOverlaps);
  //
  //NPS carbon frame
  //
  G4Box* sFrame_outer = new G4Box("Frame_outer_sol", 0.5*fSingle_X, 0.5*fSingle_Y, 0.5*fFrame_length);
  G4Box* sFrame_inner = new G4Box("Frame_inner_sol", 0.5*(fSingle_X - gap), 0.5*(fSingle_Y - gap), 0.5*fFrame_length);
  G4SubtractionSolid* sFrame = new G4SubtractionSolid("Frame_sol", sFrame_outer, sFrame_inner);

  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,
						fFrameMater,
						"Frame_log");

  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., -(0.5*fSingle_Z - (0.5*fFrame_length + fWrapThickness))),
		    lFrame,
		    "Frame_pos1",
		    lSingle,
		    false,
		    0,
		    fCheckOverlaps);
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., 0.5*fSingle_Z - (0.5*fFrame_length + fPMTcoverThickness + fPMT_length)),
                    lFrame,
                    "Frame_pos2",
                    lSingle,
                    false,
                    0,
                    fCheckOverlaps);

  //
  //opticalsurface
  //
  G4OpticalSurface* opWrapperSurface = new G4OpticalSurface("WrapperSurface");
  opWrapperSurface->SetType(dielectric_LUT);
  opWrapperSurface->SetFinish(polishedvm2000air);
  opWrapperSurface->SetModel(LUT);
  new G4LogicalBorderSurface("WrapperSurface",
			     fCrystalPos,fWrapPos,opWrapperSurface);
  G4MaterialPropertiesTable* opWS = new G4MaterialPropertiesTable();
  opWS->DumpTable();
  opWrapperSurface->SetMaterialPropertiesTable(opWS);

  G4OpticalSurface* opWrapperFrontSurface = new G4OpticalSurface("WrapperFrontSurface");
  opWrapperFrontSurface->SetType(dielectric_LUT);
  opWrapperFrontSurface->SetFinish(polishedvm2000air);
  opWrapperFrontSurface->SetModel(LUT);
  new G4LogicalBorderSurface("WrapperFrontSurface",
			     fCrystalPos,fWrapFrontPos,opWrapperFrontSurface);
  G4MaterialPropertiesTable* opWFS = new G4MaterialPropertiesTable();
  opWFS->DumpTable();
  opWrapperFrontSurface->SetMaterialPropertiesTable(opWFS);

  G4OpticalSurface* opPMTSurface = new G4OpticalSurface("PMTSurface");
  opPMTSurface->SetType(dielectric_dielectric);
  opPMTSurface->SetFinish(polished);
  opPMTSurface->SetModel(unified);
  new G4LogicalBorderSurface("PMTSurface",
			     fCrystalPos,fPMTpos,opPMTSurface);
  G4MaterialPropertiesTable* opPS = new G4MaterialPropertiesTable();
  opPS->DumpTable();
  opPMTSurface->SetMaterialPropertiesTable(opPS);

  G4OpticalSurface* opPMTcoverSurface = new G4OpticalSurface("PMTcoverSurface");
  opPMTcoverSurface->SetType(dielectric_metal);
  opPMTcoverSurface->SetFinish(polished);
  opPMTcoverSurface->SetModel(unified);
  const  G4int nEntriesPMTcover = 2;
  G4double PhotonEnergyPMTcover[nEntriesPMTcover] = { 2.30*eV, 3.26*eV}; 
  G4double reflectivityPMTcover[nEntriesPMTcover] = {0., 0.};
  G4double efficiencyPMTcover[nEntriesPMTcover] = {1., 1.};
  new G4LogicalBorderSurface("PMTcoverSurface",
			     fPMTpos,fPMTcoverPos,opPMTcoverSurface);
  G4MaterialPropertiesTable* opPcS = new G4MaterialPropertiesTable();
  opPcS -> AddProperty("REFLECTIVITY",PhotonEnergyPMTcover,reflectivityPMTcover,nEntriesPMTcover);
  opPcS -> AddProperty("EFFICIENCY",PhotonEnergyPMTcover,efficiencyPMTcover,nEntriesPMTcover);  
  opPcS->DumpTable();
  opPMTcoverSurface->SetMaterialPropertiesTable(opPcS);

  //These are all temporary/////////////////////////////////////////////////////////////////////////////////////////
  //Going to use to compare the number of particles hitting the NPS and number of particles reconstructed
  G4Box*
    sFluxOptFiber = new G4Box("FluxOptFiber_sol", 0.5*fMom_X, 0.5*fMom_Y, 0.5*fWrapThickness);

  fLogicFluxOptFiber = new G4LogicalVolume(sFluxOptFiber,
					   fWorldMater,
					   "FluxOptFiber_log");
  new G4PVPlacement(yMomRot,                     //no rotation
		    G4ThreeVector((fMom_pos_Z - 0.5*fWrapThickness - 0.5*fMom_Z)*sin(-fMom_theta), fMom_pos_Y, (fMom_pos_Z - 0.5*fWrapThickness - 0.5*fMom_Z)*cos(-fMom_theta)),
		    fLogicFluxOptFiber,
		    "FluxOptFiber_pos",
		    lWorld,
		    false,
		    0,
		    fCheckOverlaps);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //subtracting beam-pipe volume from sweeper-magnet volume
  G4RotationMatrix *SM_BP_SubRot = new G4RotationMatrix;
  SM_BP_SubRot->rotateY(-fSweepingMagnet_theta); //Angle in radian

  G4Box*
    sSweepingMagnetic_to_be_sub = new G4Box("SweepingMagnetic_to_be_sub_sol", 0.5*fSweepingMagnet_X, 0.5*fSweepingMagnet_Y, 0.5*fSweepingMagnet_Z);

  vector<G4TwoVector> vertices_beampipe1_sub;
  vertices_beampipe1_sub.push_back(G4TwoVector( 0., -25.265/(TMath::Sqrt(2))*mm ));
  vertices_beampipe1_sub.push_back(G4TwoVector( -25.265/(TMath::Sqrt(2))*mm, 0. ));
  vertices_beampipe1_sub.push_back(G4TwoVector( 0.,  25.265/(TMath::Sqrt(2))*mm ));
  vertices_beampipe1_sub.push_back(G4TwoVector(  25.265/(TMath::Sqrt(2))*mm, 0. ));
  vertices_beampipe1_sub.push_back(G4TwoVector( 0., -69.749/(TMath::Sqrt(2))*mm ));
  vertices_beampipe1_sub.push_back(G4TwoVector( -69.749/(TMath::Sqrt(2))*mm, 0. ));
  vertices_beampipe1_sub.push_back(G4TwoVector( 0.,  69.749/(TMath::Sqrt(2))*mm ));
  vertices_beampipe1_sub.push_back(G4TwoVector(  69.749/(TMath::Sqrt(2))*mm, 0. ));

  G4GenericTrap* sBeampipe1_sub = new G4GenericTrap("Beampipe1_sub_sol", 0.5*fBeampipe1Length, vertices_beampipe1_sub);

  G4double x_global = 0. - fSweepingMagnet_pos*sin(-fSweepingMagnet_arm_theta);//distance x between beam-pipe 1 and SM in global coordinate
  G4double z_global = (fChamberOuterRadius + fWindowFrameThickness + fWindowOuterJoint2Thickness - fWindowOuterJoint2_1Position + fWindowOuterJoint2_2dz + 0.5*fBeampipe1Length) - fSweepingMagnet_pos*cos(-fSweepingMagnet_arm_theta);//distance z between beam-pipe 1 and SM in global coordinate
  G4SubtractionSolid* sSweepingMagnetic = new G4SubtractionSolid("SweepingMagnetic_sol",
  								 sSweepingMagnetic_to_be_sub,
  								 sBeampipe1_sub,
  								 SM_BP_SubRot,
  								 G4ThreeVector(x_global*cos(fSweepingMagnet_theta) + z_global*sin(fSweepingMagnet_theta), 0. - fSweepingMagnetSTPcenter_Y, - x_global*sin(fSweepingMagnet_theta) + z_global*cos(fSweepingMagnet_theta) ));//beam-pipe 1 global position - sweeping magnet global position

  fLogicMagnetic = new G4LogicalVolume(sSweepingMagnetic,
				       fWorldMater,
				       "SweepingMagnetic_log");

  G4RotationMatrix *ySweepingMagnetRot = new G4RotationMatrix;  // Rotates X and Z axes only
  ySweepingMagnetRot->rotateY(fSweepingMagnet_theta); //Angle in radian                     

  //Put sweeping magnet ONLY when magnetic field is ON. When OFF there is no physical sweeping magent.
  if(field){
    new G4PVPlacement(ySweepingMagnetRot,
		      G4ThreeVector(fSweepingMagnet_pos*sin(-fSweepingMagnet_arm_theta), fSweepingMagnetSTPcenter_Y, fSweepingMagnet_pos*cos(-fSweepingMagnet_arm_theta)),
		      fLogicMagnetic,
		      "SweepingMagnetic_pos",
		      lWorld,
		      false,
		      0,
		      fCheckOverlaps);
  }

  //
  //Part 1
  //
  G4Box*
    sSweepingMagnet_1_1 = new G4Box("SweepingMagnet_1_1_sol", 0.5*7.8*inch, 0.5*27.6*inch, 0.5*39.5*inch);

  fLogicSweepingMagnet_1_1 = new G4LogicalVolume(sSweepingMagnet_1_1,
						 fSweepingMagnet_1Mater,
						 "SwepingMagnet_1_1_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*7.8)*inch + fSweepingMagnetShift_X, 16.2*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_1,
		    "SweepingMagnet_1_1_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*7.8)*inch + fSweepingMagnetShift_X, -16.2*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_1,
		    "SweepingMagnet_1_1_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Box*
    sSweepingMagnet_1_2 = new G4Box("SweepingMagnet_1_2_sol", 0.5*21*inch, 0.5*15*inch, 0.5*39.5*inch);

  fLogicSweepingMagnet_1_2 = new G4LogicalVolume(sSweepingMagnet_1_2,
						 fSweepingMagnet_1Mater,
						 "SwepingMagnet_1_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 - 0.5*21)*inch + fSweepingMagnetShift_X, 22.5*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_2,
		    "SweepingMagnet_1_2_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 - 0.5*21)*inch + fSweepingMagnetShift_X, -22.5*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_2,
		    "SweepingMagnet_1_2_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Box*
    sSweepingMagnet_1_3 = new G4Box("SweepingMagnet_1_3_sol", 0.5*15*inch, 0.5*60*inch, 0.5*39.5*inch);

  fLogicSweepingMagnet_1_3 = new G4LogicalVolume(sSweepingMagnet_1_3,
						 fSweepingMagnet_1Mater,
						 "SwepingMagnet_1_3_log");
  new G4PVPlacement(0,
		    G4ThreeVector(14.4*inch + fSweepingMagnetShift_X, 0*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_3,
		    "SweepingMagnet_1_3_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Box*
    sSweepingMagnet_1_4_1 = new G4Box("SweepingMagnet_1_4_1_sol", 0.5*4.6*inch, 0.5*30*inch, 0.5*39.5*inch);

  fLogicSweepingMagnet_1_4_1 = new G4LogicalVolume(sSweepingMagnet_1_4_1,
						   fSweepingMagnet_1Mater,
						   "SwepingMagnet_1_4_1_log");
  new G4PVPlacement(0,
		    G4ThreeVector(4.225*inch + fSweepingMagnetShift_X + 0.5*(5.35 - 4.6)*inch, 0*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_4_1,
		    "SweepingMagnet_1_4_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_1_4_2;
  vertices_1_4_2.push_back(G4TwoVector(-0.5*7.5905*inch, -0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector(-0.5*7.5905*inch,  0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector( 0.5*7.5905*inch,  0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector( 0.5*7.5905*inch, -0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector( (0.5*7.5905 - 0.10945)*inch, -0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector( (0.5*7.5905 - 0.10945)*inch,  0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector( 0.5*7.5905*inch,  0.5*14.0117*inch));
  vertices_1_4_2.push_back(G4TwoVector( 0.5*7.5905*inch, -0.5*14.0117*inch));

  G4GenericTrap*
    sSweepingMagnet_1_4_2 = new G4GenericTrap("SweepingMagnet_1_4_2_sol", 
					      0.5*39.5*inch,//dz
					      vertices_1_4_2);

  fLogicSweepingMagnet_1_4_2 = new G4LogicalVolume(sSweepingMagnet_1_4_2,
						   fSweepingMagnet_1Mater,
						   "SwepingMagnet_1_4_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(4.225*inch + fSweepingMagnetShift_X + 0.5*(5.35 - 4.6)*inch - 0.5*4.6*inch - 0.5*7.5905*inch, 0*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_1_4_2,
		    "SweepingMagnet_1_4_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  //
  //Part 2
  //
  vector<G4TwoVector> vertices_2_1_1;//must be clockwise ordered
  vertices_2_1_1.push_back(G4TwoVector( (0.5*3.14404 - 0.275987)*inch, -0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector( (0.5*3.14404 - 2.53596)*inch,  0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector( 0.5*3.14404*inch,  0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector( 0.5*3.14404*inch, -0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector( (0.5*3.14404 - 0.884073)*inch, -0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector(-0.5*3.14404*inch,  0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector( 0.5*3.14404*inch,  0.5*2.256*inch));
  vertices_2_1_1.push_back(G4TwoVector( 0.5*3.14404*inch, -0.5*2.256*inch));

  G4GenericTrap*
    sSweepingMagnet_2_1_1 = new G4GenericTrap("SweepingMagnet_2_1_1_sol", 
					      0.5*39.5*inch,//dz
					      vertices_2_1_1);

  fLogicSweepingMagnet_2_1_1 = new G4LogicalVolume(sSweepingMagnet_2_1_1,
						   fSweepingMagnet_2Mater,
						   "SwepingMagnet_2_1_1_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 + 0.5*3.14404)*inch + fSweepingMagnetShift_X, 1.272*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_2_1_1,
		    "SweepingMagnet_2_1_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_2_1_2;//must be clockwise ordered
  vertices_2_1_2.push_back(G4TwoVector( (0.5*3.14404 - 2.53596)*inch, -0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector( (0.5*3.14404 - 0.275987)*inch, 0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector( 0.5*3.14404*inch,  0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector( 0.5*3.14404*inch, -0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector(-0.5*3.14404*inch, -0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector( (0.5*3.14404 - 0.884073)*inch,  0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector( 0.5*3.14404*inch,  0.5*2.256*inch));
  vertices_2_1_2.push_back(G4TwoVector( 0.5*3.14404*inch, -0.5*2.256*inch));

  G4GenericTrap*
    sSweepingMagnet_2_1_2 = new G4GenericTrap("SweepingMagnet_2_1_2_sol", 
					      0.5*39.5*inch,//dz
					      vertices_2_1_2);

  fLogicSweepingMagnet_2_1_2 = new G4LogicalVolume(sSweepingMagnet_2_1_2,
						   fSweepingMagnet_2Mater,
						   "SwepingMagnet_2_1_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 + 0.5*3.14404)*inch + fSweepingMagnetShift_X, -1.272*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_2_1_2,
		    "SweepingMagnet_2_1_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_2_2;//must be clockwise ordered
  vertices_2_2.push_back(G4TwoVector( (0.5*0.884073 - 0.275987)*inch, -0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector( (0.5*0.884073 - 0.275987)*inch, 0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector( 0.5*0.884073*inch,  0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector( 0.5*0.884073*inch, -0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector(-0.5*0.884073*inch, -0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector(-0.5*0.884073*inch,  0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector( 0.5*0.884073*inch,  0.5*0.288*inch));
  vertices_2_2.push_back(G4TwoVector( 0.5*0.884073*inch, -0.5*0.288*inch));

  G4GenericTrap*
    sSweepingMagnet_2_2 = new G4GenericTrap("SweepingMagnet_2_2_sol", 
					    0.5*39.5*inch,//dz
					    vertices_2_2);

  fLogicSweepingMagnet_2_2 = new G4LogicalVolume(sSweepingMagnet_2_2,
						 fSweepingMagnet_2Mater,
						 "SwepingMagnet_2_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 + 0.5*0.884073)*inch + fSweepingMagnetShift_X, 0*inch + fSweepingMagnetShift_Y, 0 + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_2_2,
		    "SweepingMagnet_2_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  //
  //Part 3
  //
  vector<G4TwoVector> vertices_3_1;//must be clockwise ordered
  vertices_3_1.push_back(G4TwoVector( (0.5*0.25 - 12.125)*inch, -0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector( (0.5*0.25 - 12.125)*inch, 0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector( 0.5*0.25*inch,  0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector( 0.5*0.25*inch, -0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector(-0.5*0.25*inch, -0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector(-0.5*0.25*inch,  0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector( 0.5*0.25*inch,  0.5*1*inch));
  vertices_3_1.push_back(G4TwoVector( 0.5*0.25*inch, -0.5*1*inch));

  G4GenericTrap*
    sSweepingMagnet_3_1 = new G4GenericTrap("SweepingMagnet_3_1_sol", 
					    0.5*2.75*inch,//dz
					    vertices_3_1);

  fLogicSweepingMagnet_3_1 = new G4LogicalVolume(sSweepingMagnet_3_1,
						 fSweepingMagnet_3Mater,
						 "SwepingMagnet_3_1_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-22.025*inch + fSweepingMagnetShift_X, 2.9*inch + fSweepingMagnetShift_Y, (0.5*39.5 - 0.5*2.75)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_3_1,
		    "SweepingMagnet_3_1_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-22.025*inch + fSweepingMagnetShift_X, -2.9*inch + fSweepingMagnetShift_Y, (0.5*39.5 - 0.5*2.75)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_3_1,
		    "SweepingMagnet_3_1_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_3_2;//must be clockwise ordered
  vertices_3_2.push_back(G4TwoVector( (0.5*12.125 - 12.125)*inch, -0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector( (0.5*12.125 - 12.125)*inch, 0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector( 0.5*12.125*inch,  0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector( 0.5*12.125*inch, -0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector(-0.5*12.125*inch, -0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector(-0.5*12.125*inch,  0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector( 0.5*12.125*inch,  0.5*1*inch));
  vertices_3_2.push_back(G4TwoVector( 0.5*12.125*inch, -0.5*1*inch));

  G4GenericTrap*
    sSweepingMagnet_3_2 = new G4GenericTrap("SweepingMagnet_3_2_sol", 
					    0.5*28*inch,//dz
					    vertices_3_2);

  fLogicSweepingMagnet_3_2 = new G4LogicalVolume(sSweepingMagnet_3_2,
						 fSweepingMagnet_3Mater,
						 "SwepingMagnet_3_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 + 0.5*12.125)*inch + fSweepingMagnetShift_X, 2.9*inch + fSweepingMagnetShift_Y, (0.5*39.5 - 2.75 - 0.5*28)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_3_2,
		    "SweepingMagnet_3_2_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 + 0.5*12.125)*inch + fSweepingMagnetShift_X, -2.9*inch + fSweepingMagnetShift_Y, (0.5*39.5 - 2.75 - 0.5*28)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_3_2,
		    "SweepingMagnet_3_2_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_3_3;//must be clockwise ordered
  vertices_3_3.push_back(G4TwoVector( (0.5*12.125 - 6.125)*inch, -0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector( (0.5*12.125 - 6.125)*inch, 0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector( 0.5*12.125*inch,  0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector( 0.5*12.125*inch, -0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector(-0.5*12.125*inch, -0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector(-0.5*12.125*inch,  0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector( 0.5*12.125*inch,  0.5*1*inch));
  vertices_3_3.push_back(G4TwoVector( 0.5*12.125*inch, -0.5*1*inch));

  G4GenericTrap*
    sSweepingMagnet_3_3 = new G4GenericTrap("SweepingMagnet_3_3_sol", 
					    0.5*18*inch,//dz
					    vertices_3_3);

  fLogicSweepingMagnet_3_3 = new G4LogicalVolume(sSweepingMagnet_3_3,
						 fSweepingMagnet_3Mater,
						 "SwepingMagnet_3_3_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 + 0.5*12.125)*inch + fSweepingMagnetShift_X, 2.9*inch + fSweepingMagnetShift_Y, (0.5*39.5 - 2.75 - 28 - 0.5*18)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_3_3,
		    "SweepingMagnet_3_3_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 + 0.5*12.125)*inch + fSweepingMagnetShift_X, -2.9*inch + fSweepingMagnetShift_Y, (0.5*39.5 - 2.75 - 28 - 0.5*18)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_3_3,
		    "SweepingMagnet_3_3_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  //
  //Part 4
  //
  G4Box*
    sSweepingMagnet_4_1 = new G4Box("SweepingMagnet_4_1_sol", 0.5*7.8*inch, 0.5*26.8*inch, 0.5*5.20*inch);

  fLogicSweepingMagnet_4_1 = new G4LogicalVolume(sSweepingMagnet_4_1,
						 fSweepingMagnet_4Mater,
						 "SwepingMagnet_4_1_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*7.8)*inch + fSweepingMagnetShift_X, (0.5*2.4 + 0.5*26.8)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*5.20)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_4_1,
		    "SweepingMagnet_4_1_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*7.8)*inch + fSweepingMagnetShift_X, -(0.5*2.4 + 0.5*26.8)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*5.20)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_4_1,
		    "SweepingMagnet_4_1_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Box*
    sSweepingMagnet_4_2_before_sub = new G4Box("SweepingMagnet_4_2_before_sub_sol", 0.5*32*inch, 0.5*56*inch, 0.5*4*inch);
  G4Box*
    sSweepingMagnet_4_2_hole1 = new G4Box("SweepingMagnet_4_2_hole1_sol", (0.5*7.8 + 1/*extra cutting*/)*inch, 0.5*2.4*inch,  (0.5*4 + 1/*extra cutting*/)*inch);
  G4Box*
    sSweepingMagnet_4_2_hole2 = new G4Box("SweepingMagnet_4_2_hole2_sol", 0.5*7.4*inch, 0.5*16.8*inch, (0.5*4 + 1/*extra cutting*/)*inch);

  G4SubtractionSolid* sSweepingMagnet_4_2_sub = new G4SubtractionSolid("SweepingMagnet_4_2_sub1_sol",
								       sSweepingMagnet_4_2_before_sub,
								       sSweepingMagnet_4_2_hole1,
								       Rot,
								       G4ThreeVector(-(0.5*32 - 0.5*7.8)*inch, 0, 0));
  G4SubtractionSolid* sSweepingMagnet_4_2 = new G4SubtractionSolid("SweepingMagnet_4_2_sol",
								   sSweepingMagnet_4_2_sub,
								   sSweepingMagnet_4_2_hole2,
								   Rot,
								   G4ThreeVector(-(0.5*32 - 7.8 - 0.5*7.4)*inch, 0, 0));

  fLogicSweepingMagnet_4_2 = new G4LogicalVolume(sSweepingMagnet_4_2,
						 fSweepingMagnet_4Mater,
						 "SwepingMagnet_4_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32)*inch + fSweepingMagnetShift_X, 0*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 0.5*4)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_4_2,
		    "SweepingMagnet_4_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Box*
    sSweepingMagnet_4_3_before_sub = new G4Box("SweepingMagnet_4_3_before_sub_sol", 0.5*27*inch, 0.5*56*inch, 0.5*3*inch);
  G4Box*
    sSweepingMagnet_4_3_hole1 = new G4Box("SweepingMagnet_4_3_hole1_sol", (0.5*6 + 1/*extra cutting*/)*inch, 0.5*2.4*inch,  (0.5*3 + 1/*extra cutting*/)*inch);
  G4Box*
    sSweepingMagnet_4_3_hole2 = new G4Box("SweepingMagnet_4_3_hole2_sol", 0.5*7.4*inch, 0.5*15.225*inch, (0.5*3 + 1/*extra cutting*/)*inch);

  G4SubtractionSolid* sSweepingMagnet_4_3_sub = new G4SubtractionSolid("SweepingMagnet_4_3_sub1_sol",
								       sSweepingMagnet_4_3_before_sub,
								       sSweepingMagnet_4_3_hole1,
								       Rot,
								       G4ThreeVector(-(0.5*27 - 0.5*6)*inch, 0, 0));
  G4SubtractionSolid* sSweepingMagnet_4_3 = new G4SubtractionSolid("SweepingMagnet_4_3_sol",
								   sSweepingMagnet_4_3_sub,
								   sSweepingMagnet_4_3_hole2,
								   Rot,
								   G4ThreeVector(-(0.5*27 - 6 - 0.5*7.4)*inch, 0, 0));

  fLogicSweepingMagnet_4_3 = new G4LogicalVolume(sSweepingMagnet_4_3,
						 fSweepingMagnet_4Mater,
						 "SwepingMagnet_4_3_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 1.8 - 0.5*27)*inch + fSweepingMagnetShift_X, 0*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 3.2 - 0.5*3)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_4_3,
		    "SweepingMagnet_4_3_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  //
  //Part 5
  //
  vector<G4TwoVector> vertices_5_1_1;//must be clockwise ordered
  vertices_5_1_1.push_back(G4TwoVector( (0.5*1.33307 - 0.0381403)*inch, -0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector( (0.5*1.33307 - 1.096)*inch,  0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector( 0.5*1.33307*inch,  0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector( 0.5*1.33307*inch, -0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector( (0.5*1.33307 - 0.275217)*inch, -0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector(-0.5*1.33307*inch,  0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector( 0.5*1.33307*inch,  0.5*1.056*inch));
  vertices_5_1_1.push_back(G4TwoVector( 0.5*1.33307*inch, -0.5*1.056*inch));

  G4GenericTrap*
    sSweepingMagnet_5_1_1_before_sub = new G4GenericTrap("SweepingMagnet_5_1_1_before_sub_sol", 
							 0.5*15.4*inch,//dz
							 vertices_5_1_1);

  G4Box*
    sSweepingMagnet_5_sub = new G4Box("SweepingMagnet_5_sub_sol", 
				      0.5*1.33307*inch,
				      0.5*2.4*inch,
				      0.5*3.2*inch);

  G4SubtractionSolid* sSweepingMagnet_5_1_1 = new G4SubtractionSolid("SweepingMagnet_5_1_1_sol",
								     sSweepingMagnet_5_1_1_before_sub,
								     sSweepingMagnet_5_sub,
								     Rot,
								     G4ThreeVector(0, 0, (0.5*15.4 - 9.2 - 0.5*3.2)*inch));

  fLogicSweepingMagnet_5_1_1 = new G4LogicalVolume(sSweepingMagnet_5_1_1,
						   fSweepingMagnet_5Mater,
						   "SwepingMagnet_5_1_1_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 + 0.5*1.33307)*inch + fSweepingMagnetShift_X, (0.5*0.288 + 0.5*1.056)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*15.4)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_5_1_1,
		    "SweepingMagnet_5_1_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_5_1_2;//must be clockwise ordered
  vertices_5_1_2.push_back(G4TwoVector( (0.5*1.33307 - 1.096)*inch, -0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector( (0.5*1.33307 - 0.0381403)*inch, 0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector( 0.5*1.33307*inch,  0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector( 0.5*1.33307*inch, -0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector(-0.5*1.33307*inch, -0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector( (0.5*1.33307 - 0.275217)*inch,  0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector( 0.5*1.33307*inch,  0.5*1.056*inch));
  vertices_5_1_2.push_back(G4TwoVector( 0.5*1.33307*inch, -0.5*1.056*inch));

  G4GenericTrap*
    sSweepingMagnet_5_1_2_before_sub = new G4GenericTrap("SweepingMagnet_5_1_2_before_sub_sol", 
							 0.5*15.4*inch,//dz
							 vertices_5_1_2);

  G4SubtractionSolid* sSweepingMagnet_5_1_2 = new G4SubtractionSolid("SweepingMagnet_5_1_2_sol",
								     sSweepingMagnet_5_1_2_before_sub,
								     sSweepingMagnet_5_sub,
								     Rot,
								     G4ThreeVector(0, 0, (0.5*15.4 - 9.2 - 0.5*3.2)*inch));

  fLogicSweepingMagnet_5_1_2 = new G4LogicalVolume(sSweepingMagnet_5_1_2,
						   fSweepingMagnet_5Mater,
						   "SwepingMagnet_5_1_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 + 0.5*1.33307)*inch + fSweepingMagnetShift_X, -(0.5*0.288 + 0.5*1.056)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*15.4)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_5_1_2,
		    "SweepingMagnet_5_1_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  vector<G4TwoVector> vertices_5_2;//must be clockwise ordered
  vertices_5_2.push_back(G4TwoVector( (0.5*0.275217 - 0.0381403)*inch, -0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector( (0.5*0.275217 - 0.0381403)*inch, 0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector( 0.5*0.275217*inch,  0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector( 0.5*0.275217*inch, -0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector(-0.5*0.275217*inch, -0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector(-0.5*0.275217*inch,  0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector( 0.5*0.275217*inch,  0.5*0.288*inch));
  vertices_5_2.push_back(G4TwoVector( 0.5*0.275217*inch, -0.5*0.288*inch));

  G4GenericTrap*
    sSweepingMagnet_5_2_before_sub = new G4GenericTrap("SweepingMagnet_5_2_before_sub_sol", 
						       0.5*15.4*inch,//dz
						       vertices_5_2);

  G4SubtractionSolid* sSweepingMagnet_5_2 = new G4SubtractionSolid("SweepingMagnet_5_2_sol",
								   sSweepingMagnet_5_2_before_sub,
								   sSweepingMagnet_5_sub,
								   Rot,
								   G4ThreeVector(0, 0, (0.5*15.4 - 9.2 - 0.5*3.2)*inch));

  fLogicSweepingMagnet_5_2 = new G4LogicalVolume(sSweepingMagnet_5_2,
						 fSweepingMagnet_5Mater,
						 "SwepingMagnet_5_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 + 0.5*0.275217)*inch + fSweepingMagnetShift_X, 0*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*15.4)*inch + fSweepingMagnetShift_Z),
		    fLogicSweepingMagnet_5_2,
		    "SweepingMagnet_5_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  //
  //Part 6
  //
  G4Tubs* sSweepingMagnet_6_1 = new G4Tubs("SweepingMagnet_6_1_sol",
					   0,
					   0.5*1.5*inch,
					   0.5*1*inch,
					   0.,
					   twopi);

  fLogicSweepingMagnet_6_1 = new G4LogicalVolume(sSweepingMagnet_6_1,
						 fSweepingMagnet_6_1Mater,
						 "SweepingMagnet_6_1_log");

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 + 0.5*20.8)*inch + fSweepingMagnetShift_X, (0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 3.2 - 3 - 0.5*1)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_1,
		    "SweepingMagnet_6_1_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 - 0.5*20.8)*inch + fSweepingMagnetShift_X, (0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 3.2 - 3 - 0.5*1)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_1,
		    "SweepingMagnet_6_1_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 + 0.5*20.8)*inch + fSweepingMagnetShift_X, -(0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 3.2 - 3 - 0.5*1)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_1,
		    "SweepingMagnet_6_1_3_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 - 0.5*20.8)*inch + fSweepingMagnetShift_X, -(0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 3.2 - 3 - 0.5*1)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_1,
		    "SweepingMagnet_6_1_4_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Tubs* sSweepingMagnet_6_2 = new G4Tubs("SweepingMagnet_6_2_sol",
					   0,
					   0.5*6*inch,
					   0.5*3.2*inch,
					   0.,
					   twopi);

  fLogicSweepingMagnet_6_2 = new G4LogicalVolume(sSweepingMagnet_6_2,
						 fSweepingMagnet_6_234Mater,
						 "SweepingMagnet_6_2_log");
  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 + 0.5*20.8)*inch + fSweepingMagnetShift_X, (0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 0.5*3.2)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_2,
		    "SweepingMagnet_6_2_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 - 0.5*20.8)*inch + fSweepingMagnetShift_X, (0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 0.5*3.2)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_2,
		    "SweepingMagnet_6_2_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 + 0.5*20.8)*inch + fSweepingMagnetShift_X, -(0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 0.5*3.2)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_2,
		    "SweepingMagnet_6_2_3_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 - 0.5*20.8)*inch + fSweepingMagnetShift_X, -(0.5*56 - 5)*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 5.20 - 4 - 0.5*3.2)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_2,
		    "SweepingMagnet_6_2_4_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4Tubs* sSweepingMagnet_6_3 = new G4Tubs("SweepingMagnet_6_3_sol",
					   0,
					   0.5*7.5*inch,
					   0.5*5.2*inch,
					   0.,
					   twopi);

  fLogicSweepingMagnet_6_3 = new G4LogicalVolume(sSweepingMagnet_6_3,
						 fSweepingMagnet_6_234Mater,
						 "SweepingMagnet_6_3_log");

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 - 8.7)*inch + fSweepingMagnetShift_X, 18.62*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*5.20)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_3,
		    "SweepingMagnet_6_3_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 0.5*32 - 8.7)*inch + fSweepingMagnetShift_X, -18.63*inch + fSweepingMagnetShift_Y, (-0.5*39.5 - 0.5*5.20)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_4_2
		    fLogicSweepingMagnet_6_3,
		    "SweepingMagnet_6_3_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);


  G4Box *sSweepingMagnet_6_4_sub = new G4Box("SweepingMagnet_6_4_sub_sol",
					     0.5*7.5*inch,
					     0.5*1.75*inch,
					     0.5*5.2*inch);

  G4SubtractionSolid* sSweepingMagnet_6_4 = new G4SubtractionSolid("SweepingMagnet_6_4_sol",
								   sSweepingMagnet_6_3,
								   sSweepingMagnet_6_4_sub,
								   Rot,
								   G4ThreeVector(0, (-0.5*7.5 + 0.5*(0.5*7.5 - 2.0))*inch, 0));

  fLogicSweepingMagnet_6_4 = new G4LogicalVolume(sSweepingMagnet_6_4,
						 fSweepingMagnet_6_234Mater,
						 "SweepingMagnet_6_4_log");

  new G4PVPlacement(0,
		    G4ThreeVector(-(0.5*43.8 - 7.8 - 0.5*21 - (-0.5*21.25 + 17.5) )*inch + fSweepingMagnetShift_X, (22.5 + (-0.5*15 + 3.62) )*inch + fSweepingMagnetShift_Y, (0.5*39.5 + 0.5*5.2)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_1_2_1
		    fLogicSweepingMagnet_6_4,
		    "SweepingMagnet_6_4_1_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  G4RotationMatrix *zRot = new G4RotationMatrix;  // Rotates Y and Z axes only
  zRot->rotateZ(180*degree);                     // Rotate 90 degrees
  new G4PVPlacement(zRot,
		    G4ThreeVector(-(0.5*43.8 - 7.8 - 0.5*21 - (-0.5*21.25 + 17.5) )*inch + fSweepingMagnetShift_X, (-22.5 - (-0.5*15 + 3.62) )*inch + fSweepingMagnetShift_Y, (0.5*39.5 + 0.5*5.2)*inch + fSweepingMagnetShift_Z),//positioned relative to SweepingMagnet_1_2_2
		    fLogicSweepingMagnet_6_4,
		    "SweepingMagnet_6_4_2_pos",
		    fLogicMagnetic,
		    false,
		    0,
		    fCheckOverlaps);

  //Coloring
  G4VisAttributes* CrystalVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));//yellow
  fLogicCrystal->SetVisAttributes(CrystalVisAtt);

  G4VisAttributes* WrapVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.));//cyan
  fLogicWrap->SetVisAttributes(WrapVisAtt);

  G4VisAttributes* WrapFrontVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));//red
  fLogicWrapFront->SetVisAttributes(WrapFrontVisAtt);

  G4VisAttributes* PMTVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicPMT->SetVisAttributes(PMTVisAtt);
  G4VisAttributes* PMTcoverVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));//blue
  fLogicPMTcover->SetVisAttributes(PMTcoverVisAtt);

  G4VisAttributes* TempVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  lTemp->SetVisAttributes(TempVisAtt);
  G4VisAttributes* MomVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));//blue
  lMother->SetVisAttributes(MomVisAtt);

  G4VisAttributes* FrameVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));//magenta
  lFrame->SetVisAttributes(FrameVisAtt);
  
  G4VisAttributes* TargetVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));//blue
  fLogicTarget->SetVisAttributes(TargetVisAtt);
  G4VisAttributes* TargetCoverVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));//red
  fLogicTargetCover->SetVisAttributes(TargetCoverVisAtt);
  G4VisAttributes* TargetSphereVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));//blue
  lTargetSphere->SetVisAttributes(TargetSphereVisAtt);
  G4VisAttributes* TargetCoverSphereVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));//red
  lTargetCoverSphere->SetVisAttributes(TargetCoverSphereVisAtt);
  
  G4VisAttributes* InnerChamberVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
  fLogicInnerChamber->SetVisAttributes(InnerChamberVisAtt);

  G4VisAttributes* ChamberWindowVacuumVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
  fLogicInnerChamber2->SetVisAttributes(ChamberWindowVacuumVisAtt);

  G4VisAttributes* WindowOuterJoint2VacuumVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicWindowOuterJoint2Vacuum->SetVisAttributes(WindowOuterJoint2VacuumVisAtt);
  G4VisAttributes* WindowOuterJoint2_2VacuumVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicWindowOuterJoint2_2Vacuum->SetVisAttributes(WindowOuterJoint2_2VacuumVisAtt);
  G4VisAttributes* WindowOuterJoint2_3VacuumVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicWindowOuterJoint2_3Vacuum->SetVisAttributes(WindowOuterJoint2_3VacuumVisAtt);
  G4VisAttributes* Beampipe1VacuumVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicBeampipe1Vacuum->SetVisAttributes(Beampipe1VacuumVisAtt);
  G4VisAttributes* Beampipe2VacuumVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicBeampipe2Vacuum->SetVisAttributes(Beampipe2VacuumVisAtt);
  G4VisAttributes* Beampipe3VacuumVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  fLogicBeampipe3Vacuum->SetVisAttributes(Beampipe3VacuumVisAtt);

  G4VisAttributes* FluxVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));//blue
  fLogicFlux->SetVisAttributes(FluxVisAtt);

  PrintParameters();

  //always return the root volume
  //
  return fPhysiWorld;
}

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructSDandField()
    {
      // sensitive detectors -----------------------------------------------------
      G4SDManager* SDman = G4SDManager::GetSDMpointer();
      G4String SDname;
      G4VSensitiveDetector* PbWO4 
	= new B5HadCalorimeterSD(SDname="/HadCalorimeter");
      SDman->AddNewDetector(PbWO4);
      fLogicCrystal->SetSensitiveDetector(PbWO4);

      G4VSensitiveDetector* CrystalCover
	= new CrystalCoverSD(SDname="/CrystalCover");
      SDman->AddNewDetector(CrystalCover);
      fLogicWrap->SetSensitiveDetector(CrystalCover);

      G4VSensitiveDetector* CrystalFrontCover
	= new CrystalFrontCoverSD(SDname="/CrystalFrontCover");
      SDman->AddNewDetector(CrystalFrontCover);
      fLogicWrapFront->SetSensitiveDetector(CrystalFrontCover);

      G4VSensitiveDetector* PMTcover
	= new PMTcoverSD(SDname="/PMTcover");
      SDman->AddNewDetector(PMTcover);
      fLogicPMTcover->SetSensitiveDetector(PMTcover);

      //to measure the flux in the pseudo volume in NPS
      G4VSensitiveDetector* Flux 
	= new FluxSD(SDname="/Flux");
      SDman->AddNewDetector(Flux);
      fLogicFlux->SetSensitiveDetector(Flux);

      //to measure the flux in the pseudo volume in front of the NPS
      G4VSensitiveDetector* ChambWind 
	= new ChambWindSD(SDname="/ChambWind");
      SDman->AddNewDetector(ChambWind);
      fLogicChamberWindow->SetSensitiveDetector(ChambWind);
    
      G4FieldManager* localFieldMgr = new G4FieldManager();
      static G4bool fieldIsInitialized = field;//"true" : field on, "false" : field off
      if(fieldIsInitialized)    {
	fField = new SimpleField(fSweepingMagnetField_theta, fFieldStr);//Angle in radian
 
	fEquation = new G4Mag_UsualEqRhs (fField);
	fStepper = new G4ClassicalRK4 (fEquation);
	fChordFinder = new G4ChordFinder(fField,1e-4*m,fStepper);
	localFieldMgr->SetChordFinder(fChordFinder);
	localFieldMgr->SetDetectorField(fField);
	localFieldMgr->GetChordFinder()->SetDeltaChord(1e-4*m);
	localFieldMgr->SetDeltaIntersection(1e-4*m);
	localFieldMgr->SetDeltaOneStep(1e-4*m);
	lWorld->SetFieldManager(localFieldMgr,true);
	G4cout << "Magnetic field has been constructed " << 
	  "in DetectorConstruction::ConstructField()" << G4endl;
	fieldIsInitialized = true; 

      }
   
    }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::PrintParameters()
  {
    G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
	   << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
	   << " Material = " << fTargetMater->GetName();
    G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
	   << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")  
	   << " Material = " << fDetectorMater->GetName() << G4endl;          
    G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
  {
    // search the material by its name
    G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
    if (pttoMaterial) { 
      fTargetMater = pttoMaterial;
      if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    } else {
      G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
	     << materialChoice << " not found" << G4endl;
    }              
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
  {
    // search the material by its name
    G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
    if (pttoMaterial) { 
      fDetectorMater = pttoMaterial;
      if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    } else {
      G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
	     << materialChoice << " not found" << G4endl;
    }              
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::SetTargetRadius(G4double value)
  {
    fTargetRadius = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::SetTargetLength(G4double value)
  {
    fTargetLength = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::SetDetectorThickness(G4double value)
  {
    fDetectorThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::SetDetectorLength(G4double value)
  {
    fDetectorLength = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  void DetectorConstruction::SetDetectorGap(G4double value)
  {
    gap = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  void DetectorConstruction::SetMagneticField(G4bool value)
  {
    field = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4double DetectorConstruction::GetTargetLength()
  {
    return fTargetLength;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4double DetectorConstruction::GetTargetRadius()
  {
    return fTargetRadius;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4Material* DetectorConstruction::GetTargetMaterial()
  {
    return fTargetMater;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4LogicalVolume* DetectorConstruction::GetLogicTarget()
  {
    return fLogicTarget;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4double DetectorConstruction::GetDetectorLength()
  {
    return fDetectorLength;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4double DetectorConstruction::GetDetectorThickness()
  {
    return fDetectorThickness;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4Material* DetectorConstruction::GetDetectorMaterial()
  {
    return fDetectorMater;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4LogicalVolume* DetectorConstruction::GetLogicDetector()
  {
    return fLogicDetector;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
