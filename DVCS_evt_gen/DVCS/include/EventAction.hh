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
/// \file analysis/shared/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 67226 2013-02-08 12:07:18Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//to compilable in cc-in2p3
#include <string>
#include <vector>

#include "G4ThreeVector.hh"

#include "TDVCSEvent.h"

class RunAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*, HistoManager*);
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event*);
  virtual void    EndOfEventAction(const G4Event*);

  void DefinePrimaries(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, G4double);
  void DefineWeights(G4double, G4double, G4double, G4double);
  void DefineHMS_elec(G4ThreeVector, G4double);
    
  void AddAbs(G4double de, G4double dl) {fEnergyAbs += de; fTrackLAbs += dl;};
  void AddGap(G4double de, G4double dl) {fEnergyGap += de; fTrackLGap += dl;};

  void GetFluxEnergy(G4double FluxEnergy) {fFluxEne = FluxEnergy;};

public:
  G4int GetEventNb();
    
private:
   RunAction*    fRunAct;
   HistoManager* fHistoManager;

  G4int fEvtNb;

  G4int fHCHCID;
  std::vector<G4double> fHadCalEdep;
  std::vector<G4int> fPID;

  std::vector<G4int> fOP_sc;
  std::vector<G4int> fOP_ce;

  G4int fCrystCoverHCID;// optical photon in crystal side cover.
  std::vector<G4int> fCrystCoverOP;
  G4int fCrystFrontCoverHCID;
  std::vector<G4int> fCrystFrontCoverOP;
  G4int fPMTcoverHCID;
  std::vector<G4int> fPMTcoverOP;

  G4double fFluxEne;

  // delete this, if you want
  // reuse for Pseudo Sphere
  G4int fFlxHCID;
  G4double fFluxEnergy;
  G4int fFluxPID;
  G4ThreeVector fPos;

   TDVCSEvent* dvcs_evt;
   TCaloEvent* calo_evt;

   G4bool    hit_HMS;

   TLorentzVector* L_calo_phot;

   G4double  M_targ;
   G4double  fPx[5];
   G4double  fPy[5];
   G4double  fPz[5];
   G4double  fE[5];
   G4double  fSmeared_vertex_z;
   G4double  fVertex_z;
   G4double  fPSF;
   G4double  fX_sum;
   G4double  fX_diff;
   G4double  fX_BH;

   G4double  fEnergyAbs, fEnergyGap;
   G4double  fTrackLAbs, fTrackLGap;
                     
   G4int     fPrintModulo;                             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif

    
