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
/*
//20171013(start)
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
//20171013(finish)
*/
#include "G4RunManager.hh"
#include "G4UImanager.hh"
//20171016(start)
//20171123(start)
//#include "dvcsPhysicsList.hh"
#include "PhysicsList.hh"
//20171123(finish)
//#include "FTFP_BERT.hh"
//20171016(finish)

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"

//20180927(start)
#include "B3StackingAction.hh"
//20180927(finish)

//20190409(start)
#include "dvcsGlobals.hh"
#include "TDVCSDB.h"
#include "TCaloEvent.h"
//20190409(finish)

//20171019(start)
#include <ctime>//for cout time
//20171019(finish)

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
  //20200727(start)
  G4int pseudo_run_number = atoi(argv[3]);//0 : use hydrogen/proton target. 1 : use other targets.
  //20200727(finish)
  G4double *Calo_distance;
  G4double *Calo_angle;
  G4double *BeamEnergy;
  G4double *HMSangle;
  G4double *HMSmomentum;
  //20200727(start)
  G4double *SM_angle;
  //20200727(finish)

  // if(run_number <= 18){//DVCS run numbers
  Calo_distance=db->GetEntry_d("CALO_geom_Dist", run_number);
  Calo_angle=db->GetEntry_d("CALO_geom_Yaw", run_number);
  BeamEnergy = db->GetEntry_d("BEAM_param_Energy", run_number);
  HMSangle = db->GetEntry_d("SIMU_param_HMSangle", run_number);
  HMSmomentum = db->GetEntry_d("SIMU_param_HMSmomentum", run_number);
  //20200727(start)
  SM_angle = db->GetEntry_d("SM_geom_Yaw", run_number);
  //20200727(finish)
  dvcsGlobals::Ebeam = BeamEnergy[0];//[GeV]
  dvcsGlobals::HMS_angle = HMSangle[0];//[rad]
  dvcsGlobals::HMS_momentum = HMSmomentum[0];//[GeV]
  dvcsGlobals::Calo_distance = Calo_distance[0];//[cm]
  dvcsGlobals::Calo_angle = Calo_angle[0];//[rad]  
  dvcsGlobals::target_type = 0; //"0" : proton, "1" : deutron
  dvcsGlobals::target_gen_proc_type = 0;//"0" : proton, "1" : neutron, "2" : deutron
  // }
  if(pseudo_run_number == 1){
    //Getting fake information from the data base. 
    //No need to care, we are not going to use it.
    Calo_distance=db->GetEntry_d("CALO_geom_Dist", run_number);
    Calo_angle=db->GetEntry_d("CALO_geom_Yaw", run_number);
    BeamEnergy = db->GetEntry_d("BEAM_param_Energy", run_number);
    HMSangle = db->GetEntry_d("SIMU_param_HMSangle", run_number);
    HMSmomentum = db->GetEntry_d("SIMU_param_HMSmomentum", run_number);
    //
    dvcsGlobals::Calo_distance = atof(argv[4]);//[cm]
    dvcsGlobals::Calo_angle = atof(argv[5])*TMath::DegToRad();//[rad]
    dvcsGlobals::Ebeam = atof(argv[6]);//[GeV]
    dvcsGlobals::HMS_angle = atof(argv[7])*TMath::DegToRad();//[rad]
    dvcsGlobals::HMS_momentum = atof(argv[8]);//[GeV]
    dvcsGlobals::target_type = atoi(argv[9]); //"0" : hydrogen, "1" : deuterium
    dvcsGlobals::target_gen_proc_type = atoi(argv[10]);//"0" : proton, "1" : neutron, "2" : deutron
  }

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

  //20191017(start)
  delete Calo_distance;
  delete Calo_angle;
  delete BeamEnergy;
  delete HMSangle;
  delete HMSmomentum;

  delete db;
  //20191017(finish)

  //20190410(start)

  /*
  //20171013(start)
  // Construct the default run manager
  #ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  G4int nThreads = G4Threading::G4GetNumberOfCores();
  if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
  runManager->SetNumberOfThreads(nThreads);
  #else
  G4RunManager* runManager = new G4RunManager;
  #endif
  //20171013(finish)
  */

  //20180129(start)
  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  //20180129(finish)

  G4RunManager* runManager = new G4RunManager;
  // Set mandatory initialization classes
  //
  //20180118(start)
  //20181026(start)
  //20181026(finish)
  //  DetectorConstruction* detector = new DetectorConstruction(gap, Calo_distance, Calo_angle, field);
  //20190409(start)
  // DetectorConstruction* detector = new DetectorConstruction(gap, Calo_distance, Calo_angle, field, field_str);
  DetectorConstruction* detector = new DetectorConstruction(dvcsGlobals::Calo_distance*10./*cm to mm*/, dvcsGlobals::Calo_angle/*radian*/, dvcsGlobals::SM_field_flag, dvcsGlobals::SM_field_str);
  //20190409(finish)
  //20181026(finish)
  //  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  //20180118(finish)
  //
  //20171016(start)
  // runManager->SetUserInitialization(new FTFP_BERT);
  //20171123(start) 
  PhysicsList* phys = new PhysicsList;
  //dvcsPhysicsList* phys = new dvcsPhysicsList;
  runManager->SetUserInitialization(phys);          
  //20171123(finish)
  //20171016(finish)
  // set an HistoManager
  //
  HistoManager*  histo = new HistoManager();
      
  // Set user action classes
  //
  //20171013(start)
  /*
    PrimaryGeneratorAction* gen_action = 
    new PrimaryGeneratorAction(detector);
  */ 

  //20171201(start)
  G4String fileNamee;
  G4cout<<"Name of the output file?"<<G4endl;
  G4cin>>fileNamee;
  //20180129(start)
  G4int index;//20180720 added
  long seed1, seed2;
  G4cout<<"index for Rancu seed table?"<<G4endl;
  G4cin>>index;
  G4cout<<"seed1?"<<G4endl;
  G4cin>>seed1;
  G4cout<<"seed2?"<<G4endl;
  G4cin>>seed2;

  //20190418(start)
  //moved below EventAction
  // //20190410(start)
  // // PrimaryGeneratorAction* gen_action = 
  // //   new PrimaryGeneratorAction();
  // // runManager->SetUserAction(gen_action);
  // PrimaryGeneratorAction* gen_action = 
  //   new PrimaryGeneratorAction(seed1, seed2);
  // runManager->SetUserAction(gen_action);
  // //20190410(finish)
  // //
  // //20171013(finish)
  //20190418(finish)

  //  RunAction* run_action = new RunAction(histo, fileNamee);  
  RunAction* run_action = new RunAction(histo, fileNamee, index, seed1, seed2);  
  //20180129(finish)
  runManager->SetUserAction(run_action);
  //20171201(finish)

  //
  /*
    RunAction* run_action = new RunAction(histo);  
    runManager->SetUserAction(run_action);
  */
  //
  //20190408(start)
  // EventAction* event_action = new EventAction(run_action,histo);
  //20190409(start)
  // EventAction* event_action = new EventAction(run_action, gen_action, histo);
  EventAction* event_action = new EventAction(run_action, histo);
  //20190409(finish)
  //20190408(finish)
  runManager->SetUserAction(event_action);
  //
  //20190418(start)
  //moved from above to put EventAction
  PrimaryGeneratorAction* gen_action = 
    new PrimaryGeneratorAction(event_action, seed1, seed2);
  runManager->SetUserAction(gen_action);
  //20190418(finish)

  //
  SteppingAction* stepping_action =
    //20180622(start)
    // new SteppingAction(detector, event_action);
    new SteppingAction(detector, event_action, histo);
  //20180622(finish)  
  runManager->SetUserAction(stepping_action);

  //20180927(start)
  //To kill track such as secondary particle(ex, "cerenkov")
  B3StackingAction* stack_action = new B3StackingAction;
  runManager->SetUserAction(stack_action);          
  //20180927(finish)

  // Initialize G4 kernel
  //
  //20171019(start)
  //  runManager->Initialize();//in order to use "/testem/phys/addPhysics" command. Initialize in macro!
  //add physics process before initializing!!
  //20171019(finish)
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  //20190412(start)
  // if (argc!=1)   // batch mode
  if (argc == 3 || argc == 11)   // batch mode
    //20190412(finish)
    {
      G4String command = "/control/execute ";
      //20190412(start)
      // G4String fileName = argv[1];
      G4String fileName = argv[2];
      //20190412(finish)
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

  //20171019(start) 
  // current date/time based on current system
  time_t then = time(0);
  time_t diff = then - now;
  // convert now to string form
  char* dt = ctime(&diff);
  // convert now to tm struct for UTC
  tm* gmtm = gmtime(&diff);
  dt = asctime(gmtm);
  G4cout << "It took :"<< dt <<"for the job to be done." << G4endl;
  //20171019(finish) 

  // Job termination
  //20180622(start)
  //deleted
  //  delete histo;//20180622 deleted, why was I deleteing histo??
  //20180622(finish) 
  delete runManager;

  return 0;

}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
