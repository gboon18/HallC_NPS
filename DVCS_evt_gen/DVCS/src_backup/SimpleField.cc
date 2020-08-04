/**
  * @file
  * @brief implements class SimpleField
*/

// GN 2015 from above model

#include "SimpleField.hh"
#include "G4SystemOfUnits.hh"

//20171121(start)
#include "G4PhysicalConstants.hh"//pi = 3.14159.....
//20171121(finish)

//20180919(start)
#include <iomanip> //to use "setw(#)"
//20180919(finish)

//20180928(start)
//SimpleField::SimpleField()
SimpleField::SimpleField(G4double SM_theta, G4double field_str)
//20180928(finish)
{
  //20180928(start)
  fSM_theta = SM_theta;//SM_theta : already in radian.
  //20180928(finish)
  //20181026(start)
  str = field_str;
  //20181026(finish)

  // first things first, let's read the data from the file
  // hardwired for now, later we'll see.
  //  const char * filename ="../DumpMagnet/Dump-field.table.txt";
  //20180927(start)
  const char * filename_1cm ="/pbs/throng/clas/hosanko/NPS/NPS_Sweeping_Magnet/manetic_field_table/bogdan/1cm/SAM-Mc.table";
  //20180927(finish)
  double lenUnit= cm;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      3D Magnetic field from TOSCA MODEL (BW) "
	 << "\n-----------------------------------------------------------";
  G4cout << "\n ---> " "Reading the field grid from " << filename_1cm << " ... " << endl; 

  //20180927(start)
  ifstream file_1cm( filename_1cm ); // Open the file for reading.
  if(file_1cm)G4cout<<"YEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAAAAAAAAAAAAAH"<<G4endl;
  if(!file_1cm)G4cout<<"FUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCK"<<G4endl;
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
  //  double permeability; // Not used in this example.
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
  //  G4cout<<tesla<<endl;
  G4cout << "\n ---> ... done reading " << G4endl;
  // G4cout << " Read values of field from file " << filename_1cm << endl; 
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
  //20180927(finish)


  //20181026(start)(disable reading field_str from "field_str.txt".
  //Coming from the input in the initialization of the simulation.
  //  ifstream file2( "../DVCS/DumpMagnet/field_str.txt" ); // Open the file for reading.
  // Read field strength
  //  file2 >> str;
  G4cout<< "----> Using "<<str<<" for field strength."<<endl;
  //close the file
  //  file2.close();     
  //20181026(finish)
  /////
  G4cout << "\n-----------------------------------------------------------"
	 << "\n  Done Initializing Simple Magnetic Field! "
	 << "\n-----------------------------------------------------------"<<G4endl;;
  return;
}
 
SimpleField::~SimpleField()
{;}

void SimpleField::GetFieldValue(const double point[4],double *Bfield) const
{

  // G4double coef, G0;
  //G0 = 0;
  //coef=1; // for protons
  //coef=2; // for alphas

  //******************************************************************

  // Interpolate field MAP (a la TabulatedField3D!)
  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  Bfield[3] = 0.0;
  Bfield[4] = 0.0;
  Bfield[5] = 0.0;

  //20171121(start)  

  //20180425(start)
  // double x = point[0];
  // double y = point[1];
  // double z = point[2]; 
 
  //changing the coordinate from target coordinate to magnetic field coordinate
  //20180426 comment : you need to put the unit in mm since the default unit of Geant4 is mm. not 157*cm, just 1570.
  //20180928(start)
  // double x = ( (point[0] - 1570*sin((2.2*pi)/180.))*cos((2.2*pi)/180.) - (point[2] - 1570*cos((2.2*pi)/180.))*sin((2.2*pi)/180.) );
  // double y = point[1];
  // double z = ( (point[0] - 1570*sin((2.2*pi)/180.))*sin((2.2*pi)/180.) + (point[2] - 1570*cos((2.2*pi)/180.))*cos((2.2*pi)/180.) );
  // G4cout<<"fSM_theta : "<<fSM_theta*180./pi<<G4endl;
  //fSM_theta : already in radian.
  double x = ( (point[0] - 1570*sin(fSM_theta))*cos(fSM_theta) - (point[2] - 1570*cos(fSM_theta))*sin(fSM_theta) );
  double y = point[1];
  double z = ( (point[0] - 1570*sin(fSM_theta))*sin(fSM_theta) + (point[2] - 1570*cos(fSM_theta))*cos(fSM_theta) );
  //20180928(finish)
  //20180425(finish)
  
  // ofstream test("text.txt");
  // test<<
  //   "\n\n!!!!!!!!!!!!!!!"<<position_X<<", "<<position_Z<<""
  // 			 <<"\npoints : "<<point[0]<<", "<<point[1]<<", "<<point[2]<<""
  // 			 <<"\nx-, y-, z- : "<<x_prime<<", "<<y<<", "<<z_prime<<""
  // 			 <<"\nx, y, z : "<<x<<", "<<y<<", "<<z;
  // test.close();
  
  // G4cout<<"!!!!!!!!!!!!!!!"<<position_X<<", "<<position_Z<<G4endl;
  // G4cout<<"points : "<<point[0]<<", "<<point[1]<<", "<<point[2]<<G4endl;
  // G4cout<<"x-, y-, z- : "<<x_prime<<", "<<y<<", "<<z_prime<<G4endl;
  // G4cout<<"x, y, z : "<<x<<", "<<y<<", "<<z<<G4endl;

  //20171121(finish)

  //  G4int quad;
  //G4double gradient[5];
  /*  
      gradient[0]=fGradient1*(gauss/cm)/coef;
      gradient[1]=fGradient2*(gauss/cm)/coef;
      gradient[2]=fGradient3*(gauss/cm)/coef; 
      gradient[3]=fGradient4*(gauss/cm)/coef;
      gradient[4]=-fGradient3*(gauss/cm)/coef;
  */
  // Check that the point is within the defined region 

  // ofstream outfile("/afs/in2p3.fr/home/h/hosanko/public/playground/geant4.10.03.p01_wMT/NPS/30102018/build/field_e-_3GeV_field_center.txt", ios::app
  // 		   );
  // outfile<<fixed<<setprecision(10)<<"global x (cm)"<<"\t"<<"global y (cm)"<<"\t"<<"global z (cm)"<<"\t"<<"local x (cm)"<<"\t"<<"local y (cm)"<<"\t"<<"local z (cm)"<<"\t"<<"Bx (gaus)"<<"\t"<<"By (gaus)"<<"\t"<<"Bz (gaus)"<<"\t"<<"Bx_table (gaus)"<<"\t"<<"By_table (gaus)"<<"\t"<<"Bz_table (gaus)"<<endl;
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

      //20180924(start)
      //rotate the B-field
      //20180928(start)
      // Bfield[0] = ( (Bfield[0])*cos((-2.2*pi)/180.) - (Bfield[2])*sin((-2.2*pi)/180.) );
      // Bfield[1] = Bfield[1];
      // Bfield[2] = ( (Bfield[0])*sin((-2.2*pi)/180.) + (Bfield[2])*cos((-2.2*pi)/180.) );
      //fSM_theta : already in radian.
      Bfield[0] = ( (Bfield[0])*cos(-fSM_theta) - (Bfield[2])*sin(-fSM_theta) );
      Bfield[1] = Bfield[1];
      Bfield[2] = ( (Bfield[0])*sin(-fSM_theta) + (Bfield[2])*cos(-fSM_theta) );
      //20180928(finish)
      //20180924(finish)

      //20180919(start)
      // if(point[0] > -1 && point[0] < 1 && point[1] > -1 && point[1] < 1){
      // if(x > -10 && x < 10 && y > -10 && y < 10){
      // outfile<<fixed<<setprecision(10)<<point[0]*1e-1<<"\t"<<point[1]*1e-1<<"\t"<<point[2]*1e-1<<"\t"<<x*1e-1<<"\t"<<y*1e-1<<"\t"<<z*1e-1<<"\t"<<Bfield[0]*1e7<<"\t"<<Bfield[1]*1e7<<"\t"<<Bfield[2]*1e7<<"\t"<<fXField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<fYField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<fZField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<endl;
      // }
      // cout<<fixed<<setprecision(10)<<point[0]*1e-1<<"\t"<<point[1]*1e-1<<"\t"<<point[2]*1e-1<<"\t"<<x*1e-1<<"\t"<<y*1e-1<<"\t"<<z*1e-1<<"\t"<<Bfield[0]*1e7<<"\t"<<Bfield[1]*1e7<<"\t"<<Bfield[2]*1e7<<"\t"<<fXField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<fYField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<fZField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<endl;
      // <<" fField_1cm = "<<fXField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<fYField_1cm[xindex][yindex][zindex]*1e7<<"\t"<<fZField_1cm[xindex][yindex][zindex]*1e7<<"\n"
      // <<"/////////////////////////////////////////////////////////////////////////////////////"<<G4endl;

      //20180919(finish)
    }
  //end  MAP
}
