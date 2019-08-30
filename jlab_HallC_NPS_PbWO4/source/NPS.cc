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
// Based on AnaEx02
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

#include <ctime>//for cout time

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{ 
  // current date/time based on current system
  time_t now = time(0);
  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  G4RunManager* runManager = new G4RunManager;
  // Set mandatory initialization classes
  //
  G4double gap;
  G4double NPS_distance;
  G4double NPS_angle;
  G4bool field;
  G4double field_str;

  G4cout<<"Size of the gap between the crystals(mm)?"<<G4endl;
  G4cin>>gap;
  G4cout<<"Distance of the NPS from the target?(mm)"<<G4endl;
  G4cin>>NPS_distance;
  G4cout<<"Anlge of the NPS from the beam line?(deg)"<<G4endl;
  G4cout<<"sweeping magnet angle : theta_NPS + 1.97"<<G4endl;
  G4cout<<"sweeping magnet field angle : theta_NPS - 4.1"<<G4endl;
  G4cin>>NPS_angle;
  G4cout<<"Sweeping Magnetic Field On? Off? (On : 1(true), Off : 0(false))"<<G4endl;
  G4cin>>field;
  G4cout<<"Field Strength?"<<G4endl;
  G4cin>>field_str;

  DetectorConstruction* detector = new DetectorConstruction(gap, NPS_distance, NPS_angle, field, field_str);
  runManager->SetUserInitialization(detector);

  PhysicsList* phys = new PhysicsList;
  runManager->SetUserInitialization(phys);          

  HistoManager*  histo = new HistoManager();

  PrimaryGeneratorAction* gen_action = 
    new PrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  G4String fileNamee;
  G4cout<<"Name of the output file?"<<G4endl;
  G4cin>>fileNamee;
  G4int index;//20180720 added
  long seed1, seed2;
  G4cout<<"index for Rancu seed table?"<<G4endl;
  G4cin>>index;
  G4cout<<"seed1?"<<G4endl;
  G4cin>>seed1;
  G4cout<<"seed2?"<<G4endl;
  G4cin>>seed2;

  RunAction* run_action = new RunAction(histo, fileNamee, index, seed1, seed2);  
  runManager->SetUserAction(run_action);

  EventAction* event_action = new EventAction(run_action,histo);
  runManager->SetUserAction(event_action);

  SteppingAction* stepping_action =
    new SteppingAction(detector, event_action, histo);
  runManager->SetUserAction(stepping_action);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
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

  time_t then = time(0);
  time_t diff = then - now;
  // convert now to string form
  char* dt = ctime(&diff);
  // convert now to tm struct for UTC
  tm* gmtm = gmtime(&diff);
  dt = asctime(gmtm);
  G4cout << "It took :"<< dt <<"for the job to be done." << G4endl;

  delete runManager;

  return 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
