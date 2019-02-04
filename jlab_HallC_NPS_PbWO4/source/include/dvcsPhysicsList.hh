#ifndef DVCS_PHYSICS_LIST
#define DVCS_PHYSICS_LIST

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class dvcsPhysicsList: public G4VUserPhysicsList
{
public:
  dvcsPhysicsList();
  ~dvcsPhysicsList();

protected:
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

  
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();

  
  void ConstructGeneral();
  void ConstructEM();
  void AddStepMax();
  
};

#endif
