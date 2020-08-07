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
/// \file analysis/AnaEx02/AnaEx02.cc
/// \brief Main program of the analysis/AnaEx02 example
//
//
// $Id: AnaEx02.cc 81444 2014-05-28 14:28:20Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "PhysicsList.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"

#include "B3StackingAction.hh"

#include "dvcsGlobals.hh"
#include "TDVCSDB.h"
#include "TCaloEvent.h"

#include <ctime>//for cout time

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int    dvcsGlobals::run_number = 0;
double dvcsGlobals::Ebeam = 0;
double dvcsGlobals::HMS_angle = 0;
double dvcsGlobals::HMS_momentum = 0;
double dvcsGlobals::Calo_distance = 0;
double dvcsGlobals::Calo_angle = 0;
int    dvcsGlobals::target_type = 0;
int    dvcsGlobals::target_gen_proc_type = 0;
double dvcsGlobals::target_density = 0;
double dvcsGlobals::target_offset = 0;
double dvcsGlobals::target_length = 0;
bool   dvcsGlobals::hit_HMS_CALO_flag = false;

bool   dvcsGlobals::SM_field_flag = true;
double dvcsGlobals::SM_field_str  = 1.;

double dvcsGlobals::SM_angle = 0.;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{ 
  
  // current date/time based on current system
  time_t now = time(0);

  TDVCSDB *db=new TDVCSDB("dvcs", "clrlpc", 3306, "munoz", "");

  G4int run_number = atoi(argv[1]);
  G4double *Calo_distance;
  G4double *Calo_angle;
  G4double *BeamEnergy;
  G4double *HMSangle;
  G4double *HMSmomentum;
  G4double *SM_angle;

  Calo_distance=db->GetEntry_d("CALO_geom_Dist", run_number);
  Calo_angle=db->GetEntry_d("CALO_geom_Yaw", run_number);
  BeamEnergy = db->GetEntry_d("BEAM_param_Energy", run_number);
  HMSangle = db->GetEntry_d("SIMU_param_HMSangle", run_number);
  HMSmomentum = db->GetEntry_d("SIMU_param_HMSmomentum", run_number);
  SM_angle = db->GetEntry_d("SM_geom_Yaw", run_number);

  dvcsGlobals::Ebeam = BeamEnergy[0];//[GeV]
  dvcsGlobals::HMS_angle = HMSangle[0];//[rad]
  dvcsGlobals::HMS_momentum = HMSmomentum[0];//[GeV]
  dvcsGlobals::Calo_distance = Calo_distance[0];//[cm]
  dvcsGlobals::Calo_angle = Calo_angle[0];//[rad]  

  dvcsGlobals::target_type = atoi(argv[3]); //"0" : hydrogen, "1" : deuterium
  dvcsGlobals::target_gen_proc_type = atoi(argv[4]);//"0" : proton, "1" : neutron, "2" : deutron

  dvcsGlobals::run_number = run_number;
  dvcsGlobals::SM_angle = SM_angle[0];

  int targ_type = dvcsGlobals::target_type;
  cout<<"targ_type = "<<targ_type<<endl;
  if( targ_type == 0 ) { dvcsGlobals::target_density = 0.0723; cout<<"Using Hydrogen Target"<<endl;} // Hydrogen
  else if( targ_type == 1 ) { dvcsGlobals::target_density = 0.167; cout<<"Using Deuterium Target"<<endl;}  // Deuterium
  int target_gen_proc_type = dvcsGlobals::target_gen_proc_type;
  if( target_gen_proc_type == 0 ) {cout<<"Process to be genarated is DVCS on proton"<<endl;}
  else if ( target_gen_proc_type == 1 ) {cout<<"Process to be genarated is DVCS on neutron"<<endl;}
  else if ( target_gen_proc_type == 2 ) {cout<<"Process to be genarated is DVCS on deutron"<<endl;}

  dvcsGlobals::target_length = 15.; // [cm]

  G4cout<<"calo angle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! : "<<dvcsGlobals::Calo_angle*TMath::RadToDeg()<<G4endl;
  G4cout<<"calo dist!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! : "<<dvcsGlobals::Calo_distance<<G4endl;
  G4cout<<"spectro angle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! : "<<dvcsGlobals::HMS_angle*TMath::RadToDeg()<<G4endl;

  delete Calo_distance;
  delete Calo_angle;
  delete BeamEnergy;
  delete HMSangle;
  delete HMSmomentum;

  delete db;

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  G4RunManager* runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  DetectorConstruction* detector = new DetectorConstruction(dvcsGlobals::Calo_distance*10./*cm to mm*/, dvcsGlobals::Calo_angle/*radian*/, dvcsGlobals::SM_field_flag, dvcsGlobals::SM_field_str);
  runManager->SetUserInitialization(detector);

  PhysicsList* phys = new PhysicsList;
  runManager->SetUserInitialization(phys);          

  // set an HistoManager
  //
  HistoManager*  histo = new HistoManager();
      
  G4String fileNamee;
  G4cout<<"Name of the output file?"<<G4endl;
  G4cin>>fileNamee;

  G4int index;
  long seed1, seed2;
  G4cout<<"index for Rancu seed table?"<<G4endl;
  G4cin>>index;
  G4cout<<"seed1?"<<G4endl;
  G4cin>>seed1;
  G4cout<<"seed2?"<<G4endl;
  G4cin>>seed2;

  RunAction* run_action = new RunAction(histo, fileNamee, index, seed1, seed2);  
  runManager->SetUserAction(run_action);

  EventAction* event_action = new EventAction(run_action, histo);
  runManager->SetUserAction(event_action);

  PrimaryGeneratorAction* gen_action = 
    new PrimaryGeneratorAction(event_action, seed1, seed2);
  runManager->SetUserAction(gen_action);

  SteppingAction* stepping_action =
    new SteppingAction(detector, event_action, histo);
  runManager->SetUserAction(stepping_action);

  //To kill track such as secondary particle(ex, "cerenkov")
  B3StackingAction* stack_action = new B3StackingAction;
  runManager->SetUserAction(stack_action);          

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc == 3 || argc == 11)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[2];
      UImanager->ApplyCommand(command+fileName);    
    }
  else
    {  // interactive mode : define visualization and UI terminal
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      delete ui;
#endif

#ifdef G4VIS_USE
      delete visManager;
#endif
    }

  // current date/time based on current system
  time_t then = time(0);
  time_t diff = then - now;
  // convert now to string form
  char* dt = ctime(&diff);
  // convert now to tm struct for UTC
  tm* gmtm = gmtime(&diff);
  dt = asctime(gmtm);
  G4cout << "It took :"<< dt <<"for the job to be done." << G4endl;

  // Job termination
  delete runManager;

  return 0;

}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
