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

//20171023(start)
//to compilable in cc-in2p3
#include <string>
#include <vector>
//using namespace std;
//20171023(finish)

//20171127(start)
#include "G4ThreeVector.hh"
//20171127(finish)

//20190404(start)
#include "TDVCSEvent.h"
// #include "TCaloEvent.h"
//20190404(finish)

class RunAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  //20190408(start)
  //  EventAction(RunAction*, HistoManager*);
  //20190409(start)
  // EventAction(RunAction*, PrimaryGeneratorAction*, HistoManager*);
  EventAction(RunAction*, HistoManager*);
  //20190409(finish)
  //20190408(finish)
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event*);
  virtual void    EndOfEventAction(const G4Event*);

  //20190417(start)
  void DefinePrimaries(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, G4double);
  void DefineWeights(G4double, G4double, G4double, G4double);
  void DefineHMS_elec(G4ThreeVector, G4double);
  //20190417(finish)

    
  void AddAbs(G4double de, G4double dl) {fEnergyAbs += de; fTrackLAbs += dl;};
  void AddGap(G4double de, G4double dl) {fEnergyGap += de; fTrackLGap += dl;};

  //20180623(start)
  void GetFluxEnergy(G4double FluxEnergy) {fFluxEne = FluxEnergy;};
  //20180623(finish)

  //20190408(start)
public:
  G4int GetEventNb();
  //20190408(finish)
    
private:
   RunAction*    fRunAct;
   HistoManager* fHistoManager;

  //20190408(start)
  G4int fEvtNb;
  //20190408(finish)

  //20171013(start)
  G4int fHCHCID;
  std::vector<G4double> fHadCalEdep;
  //20171013(finish)
  //20171124(start)
  std::vector<G4int> fPID;
  //20171124(finish)      

  //20180117(start)
  std::vector<G4int> fOP_sc;
  std::vector<G4int> fOP_ce;

  G4int fCrystCoverHCID;//optical photon in crystal side cover.
  std::vector<G4int> fCrystCoverOP;
  G4int fCrystFrontCoverHCID;
  std::vector<G4int> fCrystFrontCoverOP;
  G4int fPMTcoverHCID;
  std::vector<G4int> fPMTcoverOP;
  //20180117(finish)
  //20180425(start)

  //20180623(start)
  G4double fFluxEne;
  //20180623(finish)

  //delete
  //20180713(start) reuse for Pseudo Sphere
  //20171127(start)
  G4int fFlxHCID;
  G4double fFluxEnergy;
  G4int fFluxPID;
  G4ThreeVector fPos;
  //20180713(finish)
  /*
  //20171127(finish)
  //20171129(start)
  G4int fChambWindHCID;
  G4double fChambWindEnergy;
  G4int fChambWindPID;
  G4ThreeVector fChambWindPos;
  //20171129(finish)
  */
  //20180425(finish)

   TDVCSEvent* dvcs_evt;
   TCaloEvent* calo_evt;

   //20190417(start)
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
   //20190417(finish)

   G4double  fEnergyAbs, fEnergyGap;
   G4double  fTrackLAbs, fTrackLGap;
                     
   G4int     fPrintModulo;                             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif

    
