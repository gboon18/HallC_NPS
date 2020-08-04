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
/// \file analysis/AnaEx02/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 92374 2015-08-31 08:52:09Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//20171127(start)
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//20171127(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fRootFile(0),
   //20171127(start)
   fNtuple(0),
   //20180425(start)

   //20190417(start)
   fNtuple_DVCS(0),
   //20190417(start)

   //20190408(start)
   fNtuple_Clust(0),
   //20190408(finish)

   //20180622(start)
   fNtuple2(0),
   //20180622(finish)

   //delete
   //20180713(start) reuse for Pseudo Sphere   
   //20180719(start)(modified)
   fNtuple_gamma(0),
   fNtuple_pos(0),
   fNtuple_elec(0),
   fNtuple_unknown(0),
   //20180719(finish)(modified)
   //20180713(finish)

   //20190124(start)
   fNtuple_OptFiber_gamma(0),
   fNtuple_OptFiber_pos(0),
   fNtuple_OptFiber_elec(0),
   fNtuple_OptFiber_unknown(0),
//20190124(finish)

   //20190417(start)
   fEvtNb(-1),
  fClust_Ene(-999.),
   fClust_X(-999.),
  fClust_Y(-999.),
  fClust_W2(-999.),
   fClust_size(-1),
   fPSF(-999.),
   fRIE_Px(-999.),
   fRIE_Py(-999.),
  fRIE_Pz(-999.),
  fGIE_Px(-999.),
  fGIE_Py(-999.),
  fGIE_Pz(-999.),
  fRSE_Px(-999.),
  fRSE_Py(-999.),
  fRSE_Pz(-999.),
  fGSE_Px(-999.),
  fGSE_Py(-999.),
  fGSE_Pz(-999.),
  fRP_Px(-999.),
  fRP_Py(-999.),
  fRP_Pz(-999.),
  fGP_Px(-999.),
  fGP_Py(-999.),
  fGP_Pz(-999.),
  fRV_Z(-999.),
  fGV_Z(-999.),
  fRt(-999.),
  fGt(-999.),
  fRxB(-999.),
  fGxB(-999.),
  fRQ2(-999.),
  fGQ2(-999.),
  fRphi(-999.),
  fGphi(-999.),
  fX_sum(-999.),
  fX_diff(-999.),
  fX_BH(-999.),
  fRr_val(-999.),
  fGr_val(-999.),
//20190417(finish)


//20190408(start)
  fEvtNb_(-1),
  fClust_X_(-999.),
  fClust_Y_(-999.),
  fClust_Ene_(-999.),
  fClust_Px_(-999.),
  fClust_Py_(-999.),
  fClust_Pz_(-999.),
//20190408(finish)


//delete
//20180713(start) reuse for Pseudo Sphere   
//20180719(start)(modified)
  fFluxEnergy_gamma(-999.),
  fTheta_gamma(-999.), 
  fPhi_gamma(-999.), 
  fFluxEnergy_pos(-999.),
  fTheta_pos(-999.), 
  fPhi_pos(-999.), 
  fFluxEnergy_elec(-999.),
  fTheta_elec(-999.), 
  fPhi_elec(-999.), 
  fFluxEnergy_unknown(-999.),
  fFluxPID_unknown(-999),
  fTheta_unknown(-999.), 
  fPhi_unknown(-999.),
//20180719(finish)(modified)
//20180713(finish)

//20190124(start)
//20190408(start)
  fFluxEvtNb(-1),
//20190408(finish)
//20190424(start)
  fFluxMomX_OptFiber_gamma(-999.),
  fFluxMomY_OptFiber_gamma(-999.),
  fFluxMomZ_OptFiber_gamma(-999.),
  fFluxMomX_OptFiber_pos(-999.),
  fFluxMomY_OptFiber_pos(-999.),
  fFluxMomZ_OptFiber_pos(-999.),
  fFluxMomX_OptFiber_elec(-999.),
  fFluxMomY_OptFiber_elec(-999.),
  fFluxMomZ_OptFiber_elec(-999.),
  fFluxMomX_OptFiber_unknown(-999.),
  fFluxMomY_OptFiber_unknown(-999.),
  fFluxMomZ_OptFiber_unknown(-999.),
//20190424(finish)
  fFluxEnergy_OptFiber_gamma(-999.),
  fGlobalX_OptFiber_gamma(-999.),
  fGlobalY_OptFiber_gamma(-999.),
  fGlobalZ_OptFiber_gamma(-999.),
  fLocalX_OptFiber_gamma(-999.),
  fLocalY_OptFiber_gamma(-999.),
  fLocalZ_OptFiber_gamma(-999.),
  fFluxEnergy_OptFiber_pos(-999.),
  fGlobalX_OptFiber_pos(-999.),
  fGlobalY_OptFiber_pos(-999.),
  fGlobalZ_OptFiber_pos(-999.),
  fLocalX_OptFiber_pos(-999.),
  fLocalY_OptFiber_pos(-999.),
  fLocalZ_OptFiber_pos(-999.),
  fFluxEnergy_OptFiber_elec(-999.),
  fGlobalX_OptFiber_elec(-999.),
  fGlobalY_OptFiber_elec(-999.),
  fGlobalZ_OptFiber_elec(-999.),
  fLocalX_OptFiber_elec(-999.),
  fLocalY_OptFiber_elec(-999.),
  fLocalZ_OptFiber_elec(-999.),
  fFluxEnergy_OptFiber_unknown(-999.),
  fFluxPID_OptFiber_unknown(-999),
  fGlobalX_OptFiber_unknown(-999.),
  fGlobalY_OptFiber_unknown(-999.),
  fGlobalZ_OptFiber_unknown(-999.),
  fLocalX_OptFiber_unknown(-999.),
  fLocalY_OptFiber_unknown(-999.),
  fLocalZ_OptFiber_unknown(-999.)
  //20190124(finish)

  /*
  //20171127(finish)
  //20171129(start)
  fChambWindEnergy(-999.),
  fChambWindPID(-999),
  fChambWindPos_x(-999.), 
  fChambWindPos_y(-999.), 
  fChambWindPos_z(-999.) 
  //20171129(finish)
  */
  //20180425(finish)
{
  fEdep[MaxNtuple] = {0.};
  fPID[MaxNtuple] = {0};
  fOP_sc[MaxNtuple] = {0};
  fOP_ce[MaxNtuple] = {0};
  fOP_cover[MaxNtuple] = {0};
  fOP_frontcover[MaxNtuple] = {0};
  fOP_pmtcover[MaxNtuple] = {0};

  //20180622(start)
  fFluxEne[MaxNtuple] = {0.};
  fFluxPID[MaxNtuple] = {0};
  //20180622(finish)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  if (fRootFile) delete fRootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//20171201(start)
//void HistoManager::Book()
void HistoManager::Book(G4String fileName)
//20171201(finish)
{ 
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  // G4String fileName = "AnaEx02_20171130_mapmagon_2.root";
  fRootFile = new TFile(fileName,"RECREATE");
  if (! fRootFile) {
    G4cout << " HistoManager::Book :" 
           << " problem creating the ROOT TFile "
           << G4endl;
    return;
  }
  
  fNtuple = new TTree("t","Energy deposition and OP in crystas");
  //20190423(start)
  //temporary
  fNtuple->Branch("evtNb", &fEvtNb, "Event Number/I");
  //20190423(finish)
  fNtuple->Branch("edep", fEdep, "energy_deposition[1080]/D");//20171017(changed from"energy_deposition[1080]/F"<--does not work!!!)
  // fNtuple->Branch("pid", fPID, "PID[1080]/I");

  //20180425(start)
  //temporarily disabled, no op physics
  // //20180117(start)
  // fNtuple->Branch("sc", fOP_sc, "scintillated OP[1080]/I");
  // fNtuple->Branch("ce", fOP_ce, "cerenkov OP[1080]/I");
  // fNtuple->Branch("op_cover", fOP_cover, "OP on the side of the crystal wrapper[1080]/I");
  // fNtuple->Branch("op_frontcover", fOP_frontcover, "OP on the front side of the crystal wrapper[1080]/I");
  // fNtuple->Branch("op_pc", fOP_pmtcover, "OP arrived at the pmt cover[1080]/I");
  // //20180117(finish)
  //20180425(finish)

  //20190417(start)
  //DVCS
  fNtuple_DVCS = new TTree("t_dvcs", "DVCS events generated and reconstructed");
  fNtuple_DVCS->Branch("evtNb", &fEvtNb, "Event Number/I");
  fNtuple_DVCS->Branch("clust_ene", &fClust_Ene, "Cluster Energy/D");
  fNtuple_DVCS->Branch("clust_x", &fClust_X, "Cluster X/D");
  fNtuple_DVCS->Branch("clust_y", &fClust_Y, "Cluster Y/D");
  fNtuple_DVCS->Branch("clust_W2", &fClust_W2, "Invariant mass from Clustering/D");
  fNtuple_DVCS->Branch("clust_size", &fClust_size, "Cluster size/I");
  fNtuple_DVCS->Branch("psf", &fPSF, "Phase Space Factor from DVCS gen/D");
  fNtuple_DVCS->Branch("RIE_px", &fRIE_Px, "Initial electron Px from Geant4, Actually, it is beam energy/D");
  fNtuple_DVCS->Branch("RIE_py", &fRIE_Py, "Initial electron Py from Geant4, Actually, it is beam energy/D");
  fNtuple_DVCS->Branch("RIE_pz", &fRIE_Pz, "Initial electron Pz from Geant4, Actually, it is beam energy/D");
  fNtuple_DVCS->Branch("GIE_px", &fGIE_Px, "Initial electron Px from DVCS gen/D");
  fNtuple_DVCS->Branch("GIE_py", &fGIE_Py, "Initial electron Py from DVCS gen/D");
  fNtuple_DVCS->Branch("GIE_pz", &fGIE_Pz, "Initial electron Pz from DVCS gen/D");
  fNtuple_DVCS->Branch("RSE_px", &fRSE_Px, "Scattered electron Px detected in HMS window in Geant4/D");
  fNtuple_DVCS->Branch("RSE_py", &fRSE_Py, "Scattered electron Py detected in HMS window in Geant4/D");
  fNtuple_DVCS->Branch("RSE_pz", &fRSE_Pz, "Scattered electron Pz detected in HMS window in Geant4/D");
  fNtuple_DVCS->Branch("GSE_px", &fGSE_Px, "Scattered electron Px from DVCS gen/D");
  fNtuple_DVCS->Branch("GSE_py", &fGSE_Py, "Scattered electron Py from DVCS gen/D");
  fNtuple_DVCS->Branch("GSE_pz", &fGSE_Pz, "Scattered electron Pz from DVCS gen/D");
  fNtuple_DVCS->Branch("RP_px", &fRP_Px, "Real photon Px reconstructed from clustering/D");
  fNtuple_DVCS->Branch("RP_py", &fRP_Py, "Real photon Py reconstructed from clustering/D");
  fNtuple_DVCS->Branch("RP_pz", &fRP_Pz, "Real photon Pz reconstructed from clustering/D");
  fNtuple_DVCS->Branch("GP_px", &fGP_Px, "Real photon Px from DVCS gen/D");
  fNtuple_DVCS->Branch("GP_py", &fGP_Py, "Real photon Py from DVCS gen/D");
  fNtuple_DVCS->Branch("GP_pz", &fGP_Pz, "Real photon Pz from DVCS gen/D");
  fNtuple_DVCS->Branch("RV_z", &fRV_Z, "Vertex position Z smeared from Geant4/D");//Considereing HMS resolution.
  fNtuple_DVCS->Branch("GV_z", &fGV_Z, "Vertex position Z from DVCS gen/D");
  fNtuple_DVCS->Branch("Rt", &fRt, "t from Geant4/D");
  fNtuple_DVCS->Branch("Gt", &fGt, "t from DVCS gen/D");
  fNtuple_DVCS->Branch("RxB", &fRxB, "xB from Geant4/D");
  fNtuple_DVCS->Branch("GxB", &fGxB, "xB from DVCS gen/D");
  fNtuple_DVCS->Branch("RQ2", &fRQ2, "Q2 from Geant4/D");
  fNtuple_DVCS->Branch("GQ2", &fGQ2, "Q2 from DVCS gen/D");
  fNtuple_DVCS->Branch("Rphi", &fRphi, "phi from Geant4/D");
  fNtuple_DVCS->Branch("Gphi", &fGphi, "phi from DVCS gen/D");
  fNtuple_DVCS->Branch("X_sum", &fX_sum, "XSecSum(0) from DVCS gen/D");
  fNtuple_DVCS->Branch("X_diff", &fX_diff, "XSecDif() from DVCS gen/D");
  fNtuple_DVCS->Branch("X_BH", &fX_BH, "XSecSum(1) from DVCS gen/D");
  fNtuple_DVCS->Branch("Rr_val", &fRr_val, "RFunction using values from Geant4/D");
  fNtuple_DVCS->Branch("Gr_val", &fGr_val, "RFunction using values from DVCS gen/D");
  //20190417(finish)

  //20190626(start)
  //These are not necessary for now.
  // //20190408(start)
  // //Clustering
  // fNtuple_Clust = new TTree("t_clust", "Cluster");
  // fNtuple_Clust->Branch("evtNb", &fEvtNb_, "Event Number/I");
  // fNtuple_Clust->Branch("x", &fClust_X_, "Cluster's x position on NPS/D");
  // fNtuple_Clust->Branch("y", &fClust_Y_, "Cluster's y position on NPS/D");
  // fNtuple_Clust->Branch("Ene", &fClust_Ene_, "Reconstructed particle's Energy/D");
  // fNtuple_Clust->Branch("Px", &fClust_Px_, "Reconstructed particle's Px/D");
  // fNtuple_Clust->Branch("Py", &fClust_Py_, "Reconstructed particle's Py/D");
  // fNtuple_Clust->Branch("Pz", &fClust_Pz_, "Reconstructed particle's Pz/D");
  // //20190408(finish)

  // //20180622(start)
  // fNtuple2 = new TTree("t2","Checking the flux at the Crystal front face of the NPS");
  // fNtuple2->Branch("ene", fFluxEne, "Energy of the flux in Crystal's front face[1080]/D");
  // fNtuple2->Branch("pid", fFluxPID, "PID of the flux in Crystal's front face[1080]/I");
  // //20180622(finish)

  // //20180425(start)
  // //delete
  // //20180713(start) reuse for Pseudo Sphere   
  // //20171127(start)
  // //20180719(start)(modified)
  // fNtuple_gamma = new TTree("t_gamma","Flux, 4m from the target, gamma");
  // fNtuple_gamma->Branch("fluxEnergy", &fFluxEnergy_gamma, "Energy of the flux in psuedo volume/D");
  // fNtuple_gamma->Branch("theta", &fTheta_gamma, "theta/D");
  // fNtuple_gamma->Branch("phi", &fPhi_gamma, "phi/D");

  // fNtuple_pos = new TTree("t_pos","Flux, 4m from the target, positron");
  // fNtuple_pos->Branch("fluxEnergy", &fFluxEnergy_pos, "Energy of the flux in psuedo volume/D");
  // fNtuple_pos->Branch("theta", &fTheta_pos, "theta/D");
  // fNtuple_pos->Branch("phi", &fPhi_pos, "phi/D");

  // fNtuple_elec = new TTree("t_elec","Flux, 4m from the target, electron");
  // fNtuple_elec->Branch("fluxEnergy", &fFluxEnergy_elec, "Energy of the flux in psuedo volume/D");
  // fNtuple_elec->Branch("theta", &fTheta_elec, "theta/D");
  // fNtuple_elec->Branch("phi", &fPhi_elec, "phi/D");

  // fNtuple_unknown = new TTree("t_unknown","Flux, 4m from the target, unknown");
  // fNtuple_unknown->Branch("fluxEnergy", &fFluxEnergy_unknown, "Energy of the flux in psuedo volume/D");
  // fNtuple_unknown->Branch("fluxPIS", &fFluxPID_unknown, "PID/I");
  // fNtuple_unknown->Branch("theta", &fTheta_unknown, "theta/D");
  // fNtuple_unknown->Branch("phi", &fPhi_unknown, "phi/D");
  // //20180719(finish)(modified)
  // //20180713(finish)

  // //20190124(start)
  // fNtuple_OptFiber_gamma = new TTree("t_OptFiber_gamma","Flux on Optical Fiber, mm unit, gamma");
  // //20190408(start)
  // fNtuple_OptFiber_gamma->Branch("evtNb", &fFluxEvtNb, "Event Number/I");
  // //20190408(finish)
  // fNtuple_OptFiber_gamma->Branch("fluxEnergy", &fFluxEnergy_OptFiber_gamma, "Energy of the flux in psuedo volume/D");
  // //20190424(start)
  // fNtuple_OptFiber_gamma->Branch("px", &fFluxMomX_OptFiber_gamma, "Momentum x of the flux in psuedo volume/D");
  // fNtuple_OptFiber_gamma->Branch("py", &fFluxMomY_OptFiber_gamma, "Momentum y of the flux in psuedo volume/D");
  // fNtuple_OptFiber_gamma->Branch("pz", &fFluxMomZ_OptFiber_gamma, "Momentum z of the flux in psuedo volume/D");
  // //20190424(finish)
  // fNtuple_OptFiber_gamma->Branch("x_global", &fGlobalX_OptFiber_gamma, "global_x/D");
  // fNtuple_OptFiber_gamma->Branch("y_global", &fGlobalY_OptFiber_gamma, "global_y/D");
  // fNtuple_OptFiber_gamma->Branch("z_global", &fGlobalZ_OptFiber_gamma, "global_z/D");
  // fNtuple_OptFiber_gamma->Branch("x_local", &fLocalX_OptFiber_gamma, "local_x/D");
  // fNtuple_OptFiber_gamma->Branch("y_local", &fLocalY_OptFiber_gamma, "local_y/D");
  // fNtuple_OptFiber_gamma->Branch("z_local", &fLocalZ_OptFiber_gamma, "local_z/D");

  // fNtuple_OptFiber_pos = new TTree("t_OptFiber_pos","Flux on Optical Fiber, mm unit, positron");
  // //20190408(start)
  // fNtuple_OptFiber_pos->Branch("evtNb", &fFluxEvtNb, "Event Number/I");
  // //20190408(finish)
  // fNtuple_OptFiber_pos->Branch("fluxEnergy", &fFluxEnergy_OptFiber_pos, "Energy of the flux in psuedo volume/D");
  // //20190424(start)
  // fNtuple_OptFiber_pos->Branch("px", &fFluxMomX_OptFiber_pos, "Momentum x of the flux in psuedo volume/D");
  // fNtuple_OptFiber_pos->Branch("py", &fFluxMomY_OptFiber_pos, "Momentum y of the flux in psuedo volume/D");
  // fNtuple_OptFiber_pos->Branch("pz", &fFluxMomZ_OptFiber_pos, "Momentum z of the flux in psuedo volume/D");
  // //20190424(finish)
  // fNtuple_OptFiber_pos->Branch("x_global", &fGlobalX_OptFiber_pos, "global_x/D");
  // fNtuple_OptFiber_pos->Branch("y_global", &fGlobalY_OptFiber_pos, "global_y/D");
  // fNtuple_OptFiber_pos->Branch("z_global", &fGlobalZ_OptFiber_pos, "global_z/D");
  // fNtuple_OptFiber_pos->Branch("x_local", &fLocalX_OptFiber_pos, "local_x/D");
  // fNtuple_OptFiber_pos->Branch("y_local", &fLocalY_OptFiber_pos, "local_y/D");
  // fNtuple_OptFiber_pos->Branch("z_local", &fLocalZ_OptFiber_pos, "local_z/D");

  // fNtuple_OptFiber_elec = new TTree("t_OptFiber_elec","Flux on Optical Fiber, mm unit, electron");
  // //20190408(start)
  // fNtuple_OptFiber_elec->Branch("evtNb", &fFluxEvtNb, "Event Number/I");
  // //20190408(finish)
  // fNtuple_OptFiber_elec->Branch("fluxEnergy", &fFluxEnergy_OptFiber_elec, "Energy of the flux in psuedo volume/D");
  // //20190424(start)
  // fNtuple_OptFiber_elec->Branch("px", &fFluxMomX_OptFiber_elec, "Momentum x of the flux in psuedo volume/D");
  // fNtuple_OptFiber_elec->Branch("py", &fFluxMomY_OptFiber_elec, "Momentum y of the flux in psuedo volume/D");
  // fNtuple_OptFiber_elec->Branch("pz", &fFluxMomZ_OptFiber_elec, "Momentum z of the flux in psuedo volume/D");
  // //20190424(finish)
  // fNtuple_OptFiber_elec->Branch("x_global", &fGlobalX_OptFiber_elec, "global_x/D");
  // fNtuple_OptFiber_elec->Branch("y_global", &fGlobalY_OptFiber_elec, "global_y/D");
  // fNtuple_OptFiber_elec->Branch("z_global", &fGlobalZ_OptFiber_elec, "global_z/D");
  // fNtuple_OptFiber_elec->Branch("x_local", &fLocalX_OptFiber_elec, "local_x/D");
  // fNtuple_OptFiber_elec->Branch("y_local", &fLocalY_OptFiber_elec, "local_y/D");
  // fNtuple_OptFiber_elec->Branch("z_local", &fLocalZ_OptFiber_elec, "local_z/D");

  // fNtuple_OptFiber_unknown = new TTree("t_OptFiber_unknown","Flux on Optical Fiber, mm unit, unknown");
  // //20190408(start)
  // fNtuple_OptFiber_unknown->Branch("evtNb", &fFluxEvtNb, "Event Number/I");
  // //20190408(finish)
  // fNtuple_OptFiber_unknown->Branch("fluxEnergy", &fFluxEnergy_OptFiber_unknown, "Energy of the flux in psuedo volume/D");
  // //20190424(start)
  // fNtuple_OptFiber_unknown->Branch("px", &fFluxMomX_OptFiber_unknown, "Momentum x of the flux in psuedo volume/D");
  // fNtuple_OptFiber_unknown->Branch("py", &fFluxMomY_OptFiber_unknown, "Momentum y of the flux in psuedo volume/D");
  // fNtuple_OptFiber_unknown->Branch("pz", &fFluxMomZ_OptFiber_unknown, "Momentum z of the flux in psuedo volume/D");
  // //20190424(finish)
  // fNtuple_OptFiber_unknown->Branch("fluxPID", &fFluxPID_OptFiber_unknown, "PID/I");
  // fNtuple_OptFiber_unknown->Branch("x_global", &fGlobalX_OptFiber_unknown, "global_x/D");
  // fNtuple_OptFiber_unknown->Branch("y_global", &fGlobalY_OptFiber_unknown, "global_y/D");
  // fNtuple_OptFiber_unknown->Branch("z_global", &fGlobalZ_OptFiber_unknown, "global_z/D");
  // fNtuple_OptFiber_unknown->Branch("x_local", &fLocalX_OptFiber_unknown, "local_x/D");
  // fNtuple_OptFiber_unknown->Branch("y_local", &fLocalY_OptFiber_unknown, "local_y/D");
  // fNtuple_OptFiber_unknown->Branch("z_local", &fLocalZ_OptFiber_unknown, "local_z/D");
  // //20190124(finish)


  // /*
  // //20171127(finish)
  // //20171129(start)
  // fNtuple2->Branch("ChambWindEnergy", &fChambWindEnergy, "Energy of the flux in chamber window/D");
  // fNtuple2->Branch("ChambWindPID", &fChambWindPID, "PID of the flux in flux in chamber window/I");
  // fNtuple2->Branch("chambx", &fChambWindPos_x, "Position(x) of the flux in flux in chamber window/D");
  // fNtuple2->Branch("chamby", &fChambWindPos_y, "Position(y) of the flux in flux in chamber window/D");
  // fNtuple2->Branch("chambz", &fChambWindPos_z, "Position(z) of the flux in flux in chamber window/D");
  // //20171129(finish)
  // */
  // //20180425(finish)
  //20190626(finish)
  G4cout << "\n----> Output file is open in " << fileName << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{ 
  if (! fRootFile) return;
  fRootFile->Write();       // Writing the histograms to the file
  fRootFile->Close();       // and closing the tree (and the file)
  
  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::FillNtuple()
{
  //  G4cout<<"Getting Called, FillNtuple"<<G4endl;
  //  G4cout<<"filling"<<G4endl;
  fNtuple->Fill();//20171127(temp)
  //  fNtuple2->Fill();
}
//20171013(finish)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20190417(start)
void HistoManager::FillNtupleDVCS()
{
  fNtuple_DVCS->Fill();
}
//20190417(finish)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20190408(start)
void HistoManager::FillNtupleClust()
{
  //20190629(start)
  // fNtuple_Clust->Fill();
  //20190629(finish)
}
//20190408(finish)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20180622(start)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::FillNtuple2()
{
  //   G4cout<<"filling"<<G4endl;
  //20190629(start)
  // fNtuple2->Fill();
  //20190629(finish)
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20180622(finish)
//20180713(start)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20190124(start)
// void HistoManager::FillNtuple3(G4int FluxPID)
void HistoManager::FillNtupleFluxSphere(G4int FluxPID)//Changed its name
//20190124(finish)
{
  //   G4cout<<"filling"<<G4endl;
  //20190629(start)
  // if(FluxPID == 22)       fNtuple_gamma->Fill();
  // else if(FluxPID == -11) fNtuple_pos->Fill();
  // else if(FluxPID ==  11) fNtuple_elec->Fill();
  // else                    fNtuple_unknown->Fill();
  //20190629(finish)
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20180713(finish)
//20190124(start)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::FillNtupleFluxOptFiber(G4int FluxPID)
{
  //20190629(start)
  // if(FluxPID == 22)       fNtuple_OptFiber_gamma->Fill();
  // else if(FluxPID == -11) fNtuple_OptFiber_pos->Fill();
  // else if(FluxPID ==  11) fNtuple_OptFiber_elec->Fill();
  // else                    fNtuple_OptFiber_unknown->Fill();
  //20190629(finish)
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20190124(finish)

//20171017(start)
void HistoManager::SetEnergyandPID(G4int id, G4int PID, G4double edep, G4int sc, G4int ce, G4int opc, G4int opfc, G4int oppc)
{
  //  G4cout<<"Getting Called, SetEnergyandPID, "<<id<<G4endl;
  //20190423(start)
  // if(edep) //20190423 Do not use this!!. Setting the edep should be done for all 1080 crystals.
              //If not, the edep from last events with edep>9 will be stored again.
  fEdep[id] = edep;
  //20190423(finish)
  //  if(fEdep[id])G4cout<<"Getting fEdep[id] : "<<fEdep[id]<<"["<<id<<"]"<<G4endl;
  //  if(edep) G4cout<<id<<", "<<fEdep[id]<<G4endl;

  //20180425(start)
  /*
  //temporarily disabled. no op physics
  fPID[id] = PID;

  //20180117(start)
  fOP_sc[id] = sc;
  fOP_ce[id] = ce;
  fOP_cover[id] = opc;
  fOP_frontcover[id] = opfc;
  fOP_pmtcover[id] = oppc;
  //20180117(finish)
  */
  //20180425(finish)
}
//20171017(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//20190417(start)
void HistoManager::SetEvtNb(G4int evtNb)
{
  fEvtNb = evtNb;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetCaloData(G4double clust_ene, G4double clust_x, G4double clust_y, G4int clust_size)
{
  fClust_Ene = clust_ene;
  fClust_X = clust_x;
  fClust_Y = clust_y;
  fClust_size = clust_size;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetClustW2(G4double W2)
{
  fClust_W2 = W2;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetInitElecRec(G4double px, G4double py, G4double pz)
{
  fRIE_Px = px;
  fRIE_Py = py;
  fRIE_Pz = pz;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetInitElecGen(G4double px, G4double py, G4double pz)
{
  fGIE_Px = px;
  fGIE_Py = py;
  fGIE_Pz = pz;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetScatElecRec(G4double px, G4double py, G4double pz)
{
  fRSE_Px = px;
  fRSE_Py = py;
  fRSE_Pz = pz;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetScatElecGen(G4double px, G4double py, G4double pz)
{
  fGSE_Px = px;
  fGSE_Py = py;
  fGSE_Pz = pz;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetPhotRec(G4double px, G4double py, G4double pz)
{
  fRP_Px = px;
  fRP_Py = py;
  fRP_Pz = pz;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetPhotGen(G4double px, G4double py, G4double pz)
{
  fGP_Px = px;
  fGP_Py = py;
  fGP_Pz = pz;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetVertexz(G4double rv_z, G4double gv_z)
{
  fRV_Z = rv_z;
  fGV_Z = gv_z;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetKineRec(G4double t, G4double xB, G4double Q2, G4double phi)
{
  fRt = t;
  fRxB = xB;
  fRQ2 = Q2;
  fRphi = phi;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetKineGen(G4double t, G4double xB, G4double Q2, G4double phi)
{
  fGt = t;
  fGxB = xB;
  fGQ2 = Q2;
  fGphi = phi;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetWeights(G4double PSF, G4double X_sum, G4double X_diff, G4double X_BH)
{
  fPSF = PSF;
  fX_sum = X_sum;
  fX_diff = X_diff;
  fX_BH = X_BH;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetRvalRec(G4double r_val)
{
  fRr_val = r_val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetRvalGen(G4double r_val)
{
  fGr_val = r_val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20190417(finish)











//20190408(start)
void HistoManager::SetClusterPositionandMomentum(G4int evtNb, G4double x, G4double y, G4double Energy, G4double Px, G4double Py, G4double Pz)
{
  // G4cout<<"histo x : "<<x<<",\t histo y : "<<y<<G4endl;
  fEvtNb_ = evtNb;
  fClust_X_ = x;
  fClust_Y_ = y;
  fClust_Ene_ = Energy;
  fClust_Px_ = Px;
  fClust_Py_ = Py;
  fClust_Pz_ = Pz;
  // G4cout<<"gotit x : "<<fClust_X<<",\t gotit y : "<<fClust_Y<<G4endl;
}
//20190408(finish)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//20180622(start)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::SetFluxEnergyandPID(G4int id, G4int PID, G4double ene)
{

  //20180704(start)//need to initialize!!
  for(int i = 0 ; i < MaxNtuple ; i++){fFluxPID[i] = 0; fFluxEne[i] = 0.;}
  //20180704(finish)
  fFluxPID[id] = PID;
  fFluxEne[id] = ene;
  // G4cout<<"HistoManager, "<<fFluxEne[id]<<"["<<id<<"]"<<G4endl;

}
//20180622(finish)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//20180713(start)
void HistoManager::SetFluxSphereEnergyandPID(G4int FluxPID, G4double energy, G4ThreeVector position)
{
  // G4cout<<"!!!!!!!!!!!!!!!!!!!!!!!"<<FluxPID<<G4endl;
  if(FluxPID == 22){
    fFluxEnergy_gamma = energy;
    fTheta_gamma = (180./pi)*(acos(position.z()/sqrt(position.x()*position.x() + position.y()*position.y() + position.z()*position.z())));
    fPhi_gamma = (180./pi)*(atan2(position.y(), position.x()));
  }
  else if(FluxPID == -11){
    fFluxEnergy_pos = energy;
    fTheta_pos = (180./pi)*(acos(position.z()/sqrt(position.x()*position.x() + position.y()*position.y() + position.z()*position.z())));
    fPhi_pos = (180./pi)*(atan2(position.y(), position.x()));
  }
  else if(FluxPID == 11){
    fFluxEnergy_elec = energy;
    fTheta_elec = (180./pi)*(acos(position.z()/sqrt(position.x()*position.x() + position.y()*position.y() + position.z()*position.z())));
    fPhi_elec = (180./pi)*(atan2(position.y(), position.x()));
  }
  else{
    fFluxEnergy_unknown = energy;
    fFluxPID_unknown = FluxPID;
    fTheta_unknown = (180./pi)*(acos(position.z()/sqrt(position.x()*position.x() + position.y()*position.y() + position.z()*position.z())));
    fPhi_unknown = (180./pi)*(atan2(position.y(), position.x()));
  }
}
//20180713(finish)
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //20190121(start)
//20190424(start)
// void HistoManager::SetFluxOptFiberEnergyandPID(G4int evtNb, G4int FluxPID, G4double energy, G4ThreeVector position, G4ThreeVector localPosition)
void HistoManager::SetFluxOptFiberEnergyandPID(G4int evtNb, G4int FluxPID, G4double energy, G4ThreeVector mom, G4ThreeVector position, G4ThreeVector localPosition)
//20190424(finish)
{
  fFluxEvtNb = evtNb;
  // G4cout<<"!!!!!!!!!!!!!!!!!!!!!!!"<<FluxPID<<G4endl;
  if(FluxPID == 22){
    fFluxEnergy_OptFiber_gamma = energy;
    //20190424(start)
    fFluxMomX_OptFiber_gamma = mom.x();
    fFluxMomY_OptFiber_gamma = mom.y();
    fFluxMomZ_OptFiber_gamma = mom.z();
    //20190424(finish)
    fGlobalX_OptFiber_gamma = position.x();
    fGlobalY_OptFiber_gamma = position.y();
    fGlobalZ_OptFiber_gamma = position.z();
    fLocalX_OptFiber_gamma = localPosition.x();
    fLocalY_OptFiber_gamma = localPosition.y();
    fLocalZ_OptFiber_gamma = localPosition.z();
  }
  else if(FluxPID == -11){
    fFluxEnergy_OptFiber_pos = energy;
    //20190424(start)
    fFluxMomX_OptFiber_pos = mom.x();
    fFluxMomY_OptFiber_pos = mom.y();
    fFluxMomZ_OptFiber_pos = mom.z();
    //20190424(finish)
    fGlobalX_OptFiber_pos = position.x();
    fGlobalY_OptFiber_pos = position.y();
    fGlobalZ_OptFiber_pos = position.z();
    fLocalX_OptFiber_pos = localPosition.x();
    fLocalY_OptFiber_pos = localPosition.y();
    fLocalZ_OptFiber_pos = localPosition.z();
  }
  else if(FluxPID == 11){
    fFluxEnergy_OptFiber_elec = energy;
    //20190424(start)
    fFluxMomX_OptFiber_elec = mom.x();
    fFluxMomY_OptFiber_elec = mom.y();
    fFluxMomZ_OptFiber_elec = mom.z();
    //20190424(finish)
    fGlobalX_OptFiber_elec = position.x();
    fGlobalY_OptFiber_elec = position.y();
    fGlobalZ_OptFiber_elec = position.z();
    fLocalX_OptFiber_elec = localPosition.x();
    fLocalY_OptFiber_elec = localPosition.y();
    fLocalZ_OptFiber_elec = localPosition.z();
  }
  else{
    fFluxEnergy_OptFiber_unknown = energy;
    //20190424(start)
    fFluxMomX_OptFiber_unknown = mom.x();
    fFluxMomY_OptFiber_unknown = mom.y();
    fFluxMomZ_OptFiber_unknown = mom.z();
    //20190424(finish)
    fFluxPID_OptFiber_unknown = FluxPID;
    fGlobalX_OptFiber_unknown = position.x();
    fGlobalY_OptFiber_unknown = position.y();
    fGlobalZ_OptFiber_unknown = position.z();
    fLocalX_OptFiber_unknown = localPosition.x();
    fLocalY_OptFiber_unknown = localPosition.y();
    fLocalZ_OptFiber_unknown = localPosition.z();
  }
}
  //20190121(finish)


  //20180425(start)
  //delete(finish)
  /*
  //20171127(start)
  //void HistoManager::SetFluxEnergyandPID(G4int FluxPID, G4double energy, G4ThreeVector position)
  //20171129(start)
  void HistoManager::SetFluxEnergyandPID(G4int FluxPID, G4double energy, G4ThreeVector position, G4int ChambWindPID, G4double ChambWindenergy, G4ThreeVector ChambWindposition)
  //20171129(finish)
  {

  G4double fMom_theta = 8.93*pi/180.;
  G4double fCrystal_Z = 20.5*mm; 
  G4double fSingle_Z = fCrystal_Z;//+fPMT_Z;// later day
  G4double fBulk_Z = fSingle_Z;
  G4double fNPS_Z = fBulk_Z;
  G4double fMom_Z = 794*mm;
  G4double fMom_pos_Z = 4000*mm - 0.5*fMom_Z + 0.5*fCrystal_Z;//front face of the NPS is 4m from the target
  G4double fFlux_Z = 65*1e-3*mm;
  G4double fPseudo_pos_Z = (fMom_pos_Z-0.5*(fFlux_Z + fNPS_Z))*cos(fMom_theta);
  G4double fPseudo_pos_X = (fMom_pos_Z-0.5*(fFlux_Z + fNPS_Z))*sin(fMom_theta);
  */
  // G4double fRel_pos_Z_1 = position.z() - fPseudo_pos_Z;//*cos(fMom_theta);
  // G4double fRel_pos_X_1 = position.x() - fPseudo_pos_X;//*sin(fMom_theta);
  /*
    G4double fRel_pos_Z = fRel_pos_Z_1*cos(fMom_theta) + fRel_pos_X_1*sin(fMom_theta);
    G4double fRel_pos_X = fRel_pos_X_1*cos(fMom_theta) - fRel_pos_Z_1*sin(fMom_theta);
    //  if(energy!=0){
    fFluxEnergy = energy;
    fFluxPID = FluxPID;
    fPos_x = fRel_pos_X;
    fPos_y = position.y();
    fPos_z = fRel_pos_Z;
    //  }

    //20171129(start)
    fChambWindEnergy = ChambWindenergy;
    fChambWindPID = ChambWindPID;
    fChambWindPos_x = ChambWindposition.x();
    fChambWindPos_y = ChambWindposition.y();
    fChambWindPos_z = ChambWindposition.z();
    //20171129(finish)


    }
    //20171127(finish)
    */
  //20180425(finish)
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
