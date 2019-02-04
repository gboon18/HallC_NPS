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
// $Id: B5HadCalorimeterHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file B5HadCalorimeterHit.hh
/// \brief Definition of the B5HadCalorimeterHit class

#ifndef CrystalCoverHit_h
#define CrystalCoverHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Hadron Calorimeter hit
///
/// It records:
/// - the cell column ID and row ID
/// - the energy deposit 
/// - the cell position and rotation

class CrystalCoverHit : public G4VHit
{
public:
  CrystalCoverHit();
  CrystalCoverHit(G4int iCol,G4int iRow);
  CrystalCoverHit(const CrystalCoverHit &right);
  virtual ~CrystalCoverHit();

  const CrystalCoverHit& operator=(const CrystalCoverHit &right);
  int operator==(const CrystalCoverHit &right) const;
    
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
    
  virtual void Draw();
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;
  virtual void Print();
    
  void SetColumnID(G4int z) { fColumnID = z; }
  G4int GetColumnID() const { return fColumnID; }

  void SetRowID(G4int z) { fRowID = z; }
  G4int GetRowID() const { return fRowID; }

  void SetEdep(G4double de) { fEdep = de; }
  void AddEdep(G4double de) { fEdep += de; }
  G4double GetEdep() const { return fEdep; }

  void SetPos(G4ThreeVector xyz) { fPos = xyz; }
  G4ThreeVector GetPos() const { return fPos; }

  void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
  G4RotationMatrix GetRot() const { return fRot; }

  void SetOPString (G4String OP) {fOP_String = OP;};
  G4String GetOPString() const { return fOP_String;};
  void AddOPInt (G4int    OP) {fOP_int += OP;};
  G4int GetOPInt() const { return fOP_int;};

  void SetOPOrigin (G4int OP) {fOP_origin = OP;};
  G4int GetOPOrigin() const { return fOP_origin;};

private:
  G4int fColumnID;
  G4int fRowID;
  G4double fEdep;
  G4ThreeVector fPos;
  G4RotationMatrix fRot;

  G4String fOP_String;
  G4int fOP_int;
  G4int fOP_origin;
};

typedef G4THitsCollection<CrystalCoverHit> CrystalCoverHitsCollection;

extern G4ThreadLocal G4Allocator<CrystalCoverHit>* CrystalCoverHitAllocator;

inline void* CrystalCoverHit::operator new(size_t)
{
    if (!CrystalCoverHitAllocator)
        CrystalCoverHitAllocator = new G4Allocator<CrystalCoverHit>;
    return (void*)CrystalCoverHitAllocator->MallocSingle();
}

inline void CrystalCoverHit::operator delete(void* aHit)
{
    CrystalCoverHitAllocator->FreeSingle((CrystalCoverHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
