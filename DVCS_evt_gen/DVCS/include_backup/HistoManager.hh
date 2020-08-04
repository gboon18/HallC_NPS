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
/// \file analysis/AnaEx02/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 92322 2015-08-27 14:54:05Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

//20171127(start)
#include "G4ThreeVector.hh"
//20171127(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TFile;
class TTree;
class TH1D;

const G4int MaxHisto = 4;
const G4int MaxNtuple = 1080;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
public:
  HistoManager();
  ~HistoManager();
   
  //20171201(start)
  //  void Book();
  void Book(G4String);
  //20171201(finish)

  void Save();

  //20171013(start)
  //20171017(start) 
  //  void FillNtuple2(G4int id, G4double edep);
  void FillNtuple();
  //20190417(start)
  void FillNtupleDVCS();
  //20190417(finish)
  //20190408(start)
  void FillNtupleClust();
  //20190408(finish)
  //20180622(start)
  void FillNtuple2();
  //20180622(finish)
  //20180713(start)
  //20190124(start)
  //  void FillNtuple3(G4int);
  void FillNtupleFluxSphere(G4int);//20190124. Changed its name
  //20190124(finish)
  //20180713(finish)
  //20190124(start)
  void FillNtupleFluxOptFiber(G4int);
  //20190124(finish)
  void SetEnergyandPID(G4int , G4int, G4double, G4int, G4int, G4int, G4int, G4int);//20171017/20180117

  //20190417(start)
  void SetEvtNb(G4int);
  void SetCaloData(G4double, G4double, G4double, G4int);
  void SetClustW2(G4double);
  void SetInitElecRec(G4double, G4double, G4double);
  void SetInitElecGen(G4double, G4double, G4double);
  void SetScatElecRec(G4double, G4double, G4double);
  void SetScatElecGen(G4double, G4double, G4double);
  void SetPhotRec(G4double, G4double, G4double);
  void SetPhotGen(G4double, G4double, G4double);
  void SetVertexz(G4double, G4double);
  void SetKineRec(G4double, G4double, G4double, G4double);
  void SetKineGen(G4double, G4double, G4double, G4double);
  void SetWeights(G4double, G4double, G4double, G4double);
  void SetRvalRec(G4double);
  void SetRvalGen(G4double);
  //20190417(finish)

  //20190408(start)
  void SetClusterPositionandMomentum(G4int, G4double, G4double, G4double, G4double, G4double, G4double);
  //20190408(finish)
  //20180622(start)
  void SetFluxEnergyandPID(G4int , G4int, G4double);
  //20180622(finish)
  //20180713(start)
  void SetFluxSphereEnergyandPID(G4int , G4double, G4ThreeVector);
  //20180713(finish)
  //20190121(start)
  //20190408(start)
  // void SetFluxOptFiberEnergyandPID(G4int , G4double, G4ThreeVector, G4ThreeVector);
  //20190424(start)
  // void SetFluxOptFiberEnergyandPID(G4int, G4int , G4double, G4ThreeVector, G4ThreeVector);
  void SetFluxOptFiberEnergyandPID(G4int, G4int , G4double, G4ThreeVector, G4ThreeVector, G4ThreeVector);
  //20190424(finish)
  //20190408(finish)
  //20190121(finish)
  //  void SetFluxEnergyandPID(G4int, G4double, G4ThreeVector);//20171127
  //  void SetFluxEnergyandPID(G4int, G4double, G4ThreeVector, G4int, G4double, G4ThreeVector);//20171129
  //20171017(finish)
  //20171013(finish)

  //20171020(start)
  void SetOP(G4int , G4int );//20171017
  //20171020(finish)
  void PrintStatistic();
        
private:
  TFile*   fRootFile;
  TTree*   fNtuple;
  //20190417(start)
  TTree*   fNtuple_DVCS;
  //20190417(finish)
  //20190408(start)
  TTree*   fNtuple_Clust;
  //20190408(finish)
  //20180622(start)
  TTree*   fNtuple2;
  //20180622(finish)
  //20180713(start)
  TTree*   fNtuple_gamma;
  TTree*   fNtuple_pos;
  TTree*   fNtuple_elec;
  TTree*   fNtuple_unknown;
  //20180713(finish)
  //20190124(start)
  TTree*   fNtuple_OptFiber_gamma;
  TTree*   fNtuple_OptFiber_pos;
  TTree*   fNtuple_OptFiber_elec;
  TTree*   fNtuple_OptFiber_unknown;
  //20190124(finish)
  //20171127(start)
  //  TTree*   fNtuple2;//20180425 deleted
  //20171127(finish)
  //  G4double fEdep[MaxNtuple+1];
  G4double fEdep[MaxNtuple];
  //20171124(start)
  G4int    fPID[MaxNtuple];
  //20171124(finish)
  //20180117(start)
  G4int    fOP_sc[MaxNtuple];
  G4int    fOP_ce[MaxNtuple];
  G4int    fOP_cover[MaxNtuple];
  G4int    fOP_frontcover[MaxNtuple];
  G4int    fOP_pmtcover[MaxNtuple];
  //20180117(finish)
  //20180425(start)

  //20180622(start)
  G4double fFluxEne[MaxNtuple];
  G4int fFluxPID[MaxNtuple];
  //20180622(finish)

  //20190417(start)
  G4int    fEvtNb;
  G4double fClust_Ene;
  G4double fClust_X;
  G4double fClust_Y;
  G4double fClust_W2;
  G4int    fClust_size;
  G4double fPSF;
  G4double fRIE_Px;
  G4double fRIE_Py;
  G4double fRIE_Pz;
  G4double fGIE_Px;
  G4double fGIE_Py;
  G4double fGIE_Pz;
  G4double fRSE_Px;
  G4double fRSE_Py;
  G4double fRSE_Pz;
  G4double fGSE_Px;
  G4double fGSE_Py;
  G4double fGSE_Pz;
  G4double fRP_Px;
  G4double fRP_Py;
  G4double fRP_Pz;
  G4double fGP_Px;
  G4double fGP_Py;
  G4double fGP_Pz;
  G4double fRV_Z;
  G4double fGV_Z;
  G4double fRt;
  G4double fGt;
  G4double fRxB;
  G4double fGxB;
  G4double fRQ2;
  G4double fGQ2;
  G4double fRphi;
  G4double fGphi;
  G4double fX_sum;
  G4double fX_diff;
  G4double fX_BH;
  G4double fRr_val;
  G4double fGr_val;
  //20190417(finish)

  //20190408(start)
  G4int     fEvtNb_;
  G4double  fClust_X_;
  G4double  fClust_Y_;
  G4double  fClust_Ene_;
  G4double  fClust_Px_;
  G4double  fClust_Py_;
  G4double  fClust_Pz_;
  //20190408(finish)

  //delete
  //20180713(start) reuse for Pseudo Sphere   
 
  //20171127(start)
  //20180719(start)(modified)
  G4double fFluxEnergy_gamma;
  G4double fTheta_gamma; 
  G4double fPhi_gamma; 
  G4double fFluxEnergy_pos;
  G4double fTheta_pos; 
  G4double fPhi_pos; 
  G4double fFluxEnergy_elec;
  G4double fTheta_elec; 
  G4double fPhi_elec; 
  G4double fFluxEnergy_unknown;
  G4int    fFluxPID_unknown;
  G4double fTheta_unknown; 
  G4double fPhi_unknown; 
  //20180719(finish)(modified)
  //20180713(finish)

  //20190124(start)
  //20190408(start)
  G4int    fFluxEvtNb;
  //20190408(finish)
  //20190424(start)
  G4double fFluxMomX_OptFiber_gamma;
  G4double fFluxMomY_OptFiber_gamma;
  G4double fFluxMomZ_OptFiber_gamma;
  G4double fFluxMomX_OptFiber_pos;
  G4double fFluxMomY_OptFiber_pos;
  G4double fFluxMomZ_OptFiber_pos;
  G4double fFluxMomX_OptFiber_elec;
  G4double fFluxMomY_OptFiber_elec;
  G4double fFluxMomZ_OptFiber_elec;
  G4double fFluxMomX_OptFiber_unknown;
  G4double fFluxMomY_OptFiber_unknown;
  G4double fFluxMomZ_OptFiber_unknown;
  //20190424(finish)
  G4double fFluxEnergy_OptFiber_gamma;
  G4double fGlobalX_OptFiber_gamma;
  G4double fGlobalY_OptFiber_gamma;
  G4double fGlobalZ_OptFiber_gamma;
  G4double fLocalX_OptFiber_gamma;
  G4double fLocalY_OptFiber_gamma;
  G4double fLocalZ_OptFiber_gamma;
  G4double fFluxEnergy_OptFiber_pos;
  G4double fGlobalX_OptFiber_pos;
  G4double fGlobalY_OptFiber_pos;
  G4double fGlobalZ_OptFiber_pos;
  G4double fLocalX_OptFiber_pos;
  G4double fLocalY_OptFiber_pos;
  G4double fLocalZ_OptFiber_pos;
  G4double fFluxEnergy_OptFiber_elec;
  G4double fGlobalX_OptFiber_elec;
  G4double fGlobalY_OptFiber_elec;
  G4double fGlobalZ_OptFiber_elec;
  G4double fLocalX_OptFiber_elec;
  G4double fLocalY_OptFiber_elec;
  G4double fLocalZ_OptFiber_elec;
  G4double fFluxEnergy_OptFiber_unknown;
  G4int    fFluxPID_OptFiber_unknown;
  G4double fGlobalX_OptFiber_unknown;
  G4double fGlobalY_OptFiber_unknown;
  G4double fGlobalZ_OptFiber_unknown;
  G4double fLocalX_OptFiber_unknown;
  G4double fLocalY_OptFiber_unknown;
  G4double fLocalZ_OptFiber_unknown;

  //20190124(finish)

  /* 
  //20171127(finish)
  //20171129(start)
  G4double fChambWindEnergy;
  G4int    fChambWindPID;
  G4double fChambWindPos_x; 
  G4double fChambWindPos_y; 
  G4double fChambWindPos_z; 
 //20171129(finish)
 */
  //20180425(finish)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
