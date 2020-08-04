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
/// \file electromagnetic/TestEm4/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
//#include "G4ParticleGun.hh"//20171006 deleted
#include "globals.hh"
//20190225(start)
//20171006(start)
// #include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"//20190225. Back to using this.
//20171006(finish)
//20190225(finish)

//20190212(start)
// #include "TGenGeo.h"
// #include "TGenBase.h"
#include "TGenDVCS.h"
#include <TObject.h>
//20190212(finish)

//20190408(start)
#include "TVector3.h"
//20190408(finish)

//20190417(start)
#include "TLorentzVector.h"
//20190417(finish)

//20190417(start)
class EventAction;
//20190417(finish)
class G4Event;
//20190225(start)
//20171006(start)
// class G4GeneralParticleSource;
class G4ParticleGun;//20190225. Back to using this.
//20171006(finish)
//20190225(finish)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
  //20190410(start)
    // PrimaryGeneratorAction();    
  //20190418(start)
  // PrimaryGeneratorAction(long, long);    
  PrimaryGeneratorAction(EventAction*, long, long);    
  //20190418(finish)
  //20190410(finish)
   ~PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);
//20190225(start)
  // G4GeneralParticleSource* GetParticleGun() {return fParticleGun;};//20171006 changed from G4ParticleGun
  G4ParticleGun* GetParticleGun() {return fParticleGun;}; //20190225. Back to using this.
//20190225(finish)

  private:

  //20190417(start)
  EventAction*    fEventAction;
  //20190417(finish)

  //20190419(start)
  G4double M_targ;
  //20190419(finish)

//20190225(start)
  // G4GeneralParticleSource*  fParticleGun; //20171006 changed from G4ParticleGun       //pointer a to G4 service class
  G4ParticleGun*  fParticleGun; //20190225. Back to using this.
//20190225(finish)

//20190212(start)
  TGenDVCS* gEv;
//20190212(finish)

//20190417(start)
  TLorentzVector* L_elec_init;
  TLorentzVector* L_elec_scat;
  TLorentzVector* L_phot_final;
  TLorentzVector* L_prot_final;

  G4int target_gen_proc_type;
//20190417(finish)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


