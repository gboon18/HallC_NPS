/**
  * @file
  * @brief implements class SimpleField
*/

// GN 2015 from above model

#include "SimpleField.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicalConstants.hh"//pi = 3.14159.....

#include <iomanip> //to use "setw(#)"

SimpleField::SimpleField(G4double SM_theta, G4double field_str)
{
  fSM_theta = SM_theta;//SM_theta : already in radian.
  str = field_str;

  // first things first, let's read the data from the file
  // hardwired for now, later we'll see.
  const char * filename_1cm ="/work/hallc/nps/hosan/SM_magnetic_field/bogdan/1cm/SAM-Mc.table";
  double lenUnit= cm;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      3D Magnetic field from TOSCA MODEL (BW) "
	 << "\n-----------------------------------------------------------";
  G4cout << "\n ---> " "Reading the field grid from " << filename_1cm << " ... " << endl; 

  ifstream file_1cm( filename_1cm ); // Open the file for reading.
  // Read table dimensions
  file_1cm >> fNz_1cm >> fNy_1cm >> fNx_1cm >> tmp;
  G4cout << "  [ Grid size in x,y,z: " 
	 << fNx_1cm << " " << fNy_1cm << " " << fNz_1cm << " ] "
	 << G4endl;
  //
  // Set up storage space for table
  fXField_1cm.resize( fNx_1cm );
  fYField_1cm.resize( fNx_1cm );
  fZField_1cm.resize( fNx_1cm );
  for (ix=0; ix<fNx_1cm; ix++) {
    fXField_1cm[ix].resize(fNy_1cm);
    fYField_1cm[ix].resize(fNy_1cm);
    fZField_1cm[ix].resize(fNy_1cm);
    for (iy=0; iy<fNy_1cm; iy++) {
      fXField_1cm[ix][iy].resize(fNz_1cm);
      fYField_1cm[ix][iy].resize(fNz_1cm);
      fZField_1cm[ix][iy].resize(fNz_1cm);
    }
  }
  // Read in the data
  // double xval,yval,zval,bx,by,bz;
  // double permeability; // Not used in this example.
  for (ix=0; ix<fNx_1cm; ix++) {
    for (iy=0; iy<fNy_1cm; iy++) {
      for (iz=0; iz<fNz_1cm; iz++) {
        file_1cm >> xval >> yval >> zval >> bx >> by >> bz;
	if (iz<5 && iy==0 && ix==0) {
	  G4cout << xval << " "<<yval << " "<<zval<<G4endl;
	}
        if ( ix==0 && iy==0 && iz==0 ) {
          fMinix_1cm = xval * lenUnit;
          fMiniy_1cm = yval * lenUnit;
          fMiniz_1cm = zval * lenUnit;
        }
        fXField_1cm[ix][iy][iz] = bx *gauss;
        fYField_1cm[ix][iy][iz] = by *gauss;
        fZField_1cm[ix][iy][iz] = bz *gauss;
      }
    }
  }
  fMaxix_1cm = xval * lenUnit;
  fMaxiy_1cm = yval * lenUnit;
  fMaxiz_1cm = zval * lenUnit;
  G4cout << "\n ---> ... done reading " << G4endl;
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << fMinix_1cm/cm << " " << fMiniy_1cm/cm << " " << fMiniz_1cm/cm << " cm "
	 <<"\n ---> Min position field Bx,By,Bz: "
	 << fXField_1cm[0][0][0] << " " << fYField_1cm[0][0][0] << " " << fZField_1cm[0][0][0] << ""
	 << "\n ---> Max values x,y,z: " 
	 << fMaxix_1cm/cm << " " << fMaxiy_1cm/cm << " " << fMaxiz_1cm/cm << " cm " 
	 <<"\n ---> Max position field Bx,By,Bz: "
	 << fXField_1cm[fNx_1cm - 1][fNy_1cm - 1][fNz_1cm - 1] << " " << fYField_1cm[fNx_1cm - 1][fNy_1cm - 1][fNz_1cm - 1] << " " << fZField_1cm[fNx_1cm - 1][fNy_1cm - 1][fNz_1cm - 1] << endl;

  fDx_1cm = fMaxix_1cm - fMinix_1cm;
  fDy_1cm = fMaxiy_1cm - fMiniy_1cm;
  fDz_1cm = fMaxiz_1cm - fMiniz_1cm;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << fDx_1cm/cm << " " << fDy_1cm/cm << " " << fDz_1cm/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << endl;
  //close the file
  file_1cm.close();

  G4cout<< "----> Using "<<str<<" for field strength."<<endl;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n  Done Initializing Simple Magnetic Field! "
	 << "\n-----------------------------------------------------------"<<G4endl;;
  return;
}
 
SimpleField::~SimpleField()
{;}

void SimpleField::GetFieldValue(const double point[4],double *Bfield) const
{

  //******************************************************************

  // Interpolate field MAP (a la TabulatedField3D!)
  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  Bfield[3] = 0.0;
  Bfield[4] = 0.0;
  Bfield[5] = 0.0;
 
  // Changing the coordinate from target coordinate to magnetic field coordinate
  // comment : you need to put the unit in mm since the default unit of Geant4 is mm. not 157*cm, just 1570.
  double x = ( (point[0] - 1570*sin(fSM_theta))*cos(fSM_theta) - (point[2] - 1570*cos(fSM_theta))*sin(fSM_theta) );
  double y = point[1];
  double z = ( (point[0] - 1570*sin(fSM_theta))*sin(fSM_theta) + (point[2] - 1570*cos(fSM_theta))*cos(fSM_theta) );

  // Check that the point is within the defined region 
  if   (
  	x>=fMinix_1cm && x<=fMaxix_1cm &&
  	y>=fMiniy_1cm && y<=fMaxiy_1cm &&
  	z>=fMiniz_1cm && z<=fMaxiz_1cm  )
    {
      double xfraction = (x - fMinix_1cm) / fDx_1cm;
      double yfraction = (y - fMiniy_1cm) / fDy_1cm; 
      double zfraction = (z - fMiniz_1cm) / fDz_1cm;

      // Need addresses of these to pass to modf below.
      // modf uses its second argument as an OUTPUT argument.
      double xdindex, ydindex, zdindex;
    
      // Position of the point within the cuboid defined by the
      // nearest surrounding tabulated points
      double xlocal = ( std::modf(xfraction*(fNx_1cm-1), &xdindex));
      double ylocal = ( std::modf(yfraction*(fNy_1cm-1), &ydindex));
      double zlocal = ( std::modf(zfraction*(fNz_1cm-1), &zdindex));
    
      // The indices of the nearest tabulated point whose coordinates
      // are all less than those of the given point
      int xindex = static_cast<int>(xdindex);
      int yindex = static_cast<int>(ydindex);
      int zindex = static_cast<int>(zdindex);

      // Interpolated field
      Bfield[0] =
      	((fXField_1cm[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      	  fXField_1cm[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      	  fXField_1cm[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      	  fXField_1cm[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      	  fXField_1cm[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      	  fXField_1cm[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      	  fXField_1cm[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      	  fXField_1cm[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)
      	 + Bfield[0])*str;
      
      Bfield[1] =
      	((fYField_1cm[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      	  fYField_1cm[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      	  fYField_1cm[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      	  fYField_1cm[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      	  fYField_1cm[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      	  fYField_1cm[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      	  fYField_1cm[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      	  fYField_1cm[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal) 
      	 + Bfield[1])*str;

      Bfield[2] =
      	((fZField_1cm[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      	  fZField_1cm[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      	  fZField_1cm[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      	  fZField_1cm[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      	  fZField_1cm[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      	  fZField_1cm[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      	  fZField_1cm[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      	  fZField_1cm[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal)
      	 + Bfield[2])*str;

      //rotate the B-field
      Bfield[0] = ( (Bfield[0])*cos(-fSM_theta) - (Bfield[2])*sin(-fSM_theta) );
      Bfield[1] = Bfield[1];
      Bfield[2] = ( (Bfield[0])*sin(-fSM_theta) + (Bfield[2])*cos(-fSM_theta) );
    }
  //end  MAP
}
