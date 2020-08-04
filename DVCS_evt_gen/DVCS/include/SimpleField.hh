#ifndef SimpleField_H
#define SimpleField_H 1

/**
  * @file
  * @brief Defines class SimpleField
*/
#include "globals.hh"
#include "G4MagneticField.hh"
 
#include <fstream>
#include <vector>

using namespace std;

class SimpleField 

#ifndef STANDALONE
: public G4MagneticField
#endif

{
public:
  SimpleField(G4double, G4double);
  ~SimpleField();
  void GetFieldValue( const  double Point[4],
		      double *Bfield ) const;
private:
  G4double Bvalue;
  G4double zMax;
  G4double str;

  G4double tmp;
  int ix, iy, iz;
  float xval,yval,zval,bx,by,bz;

  vector< vector< vector< double > > > fXField_1cm;
  vector< vector< vector< double > > > fYField_1cm;
  vector< vector< vector< double > > > fZField_1cm;
  G4int fNx_1cm,fNy_1cm,fNz_1cm; 
  G4double fMinix_1cm, fMaxix_1cm, fMiniy_1cm, fMaxiy_1cm, fMiniz_1cm, fMaxiz_1cm;
  G4double fDx_1cm, fDy_1cm, fDz_1cm;
  G4float fGradient1_1cm, fGradient2_1cm, fGradient3_1cm, fGradient4_1cm;

  G4double fSM_theta;
};

#endif
