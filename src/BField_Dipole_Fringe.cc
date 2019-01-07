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
//
// $Id: BField_Dipole_Fringe.cc 69786 2013-05-15 09:38:51Z gcosmo $
//
// -------------------------------------------------------------------

#include "BField_Dipole_Fringe.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>
#include <math.h>
#define DEBUG_BFIELD_DIPOLE_FRINGE 1

extern"C" {
  void dfringe_(float* x,    float* y,  float* z,
		float* fact, float* a,
		float* fx,   float* fy, float* fz);
}

static G4RotationMatrix IdentityMatrix; 

BField_Dipole_Fringe::BField_Dipole_Fringe(G4double pBend)
{

   fBend = pBend ;
   fOrigin      = G4ThreeVector( 0.0, 0.0, 0.0) ;
   fpMatrix      = &IdentityMatrix;

}

/////////////////////////////////////////////////////////////////////////

BField_Dipole_Fringe::BField_Dipole_Fringe(G4double pBend, G4ThreeVector
pOrigin, G4RotationMatrix* pMatrix)
{
   fBend    = pBend ;
   fOrigin      = pOrigin ;
   fpMatrix      = pMatrix ;




}

/////////////////////////////////////////////////////////////////////////

BField_Dipole_Fringe::~BField_Dipole_Fringe()
{
  delete fpMatrix;
}

////////////////////////////////////////////////////////////////////////
//  Allow displaced origin and rotation 
//  Extensions by Björn Riese (GSI)

void BField_Dipole_Fringe::GetFieldValue( const G4double y[7],
				 G4double B[3]  ) const  
{
  using namespace std;
  if (0)
  {
     for (float z_tmp=0.; z_tmp < 60; z_tmp+=0.1)
     {
       float snake_xx_tmp=0.;
       float snake_yy_tmp=0.1;
       float snake_zz_tmp=0.;
       float grd1_tmp = 1.;
       float width_tmp = 250.;
       snake_yy_tmp=z_tmp*10.;
       float bx_tmp=0, by_tmp=0, bz_tmp=0;
//       short int n_tmp=3;

       dfringe_(&snake_xx_tmp, &snake_zz_tmp, &snake_yy_tmp,
              &grd1_tmp, &width_tmp,
              &bx_tmp, &bz_tmp, &by_tmp);

       std::cout<<"Fillldd("<<z_tmp<<",  "<<bz_tmp<<";"<<std::endl;
     }
  }
  //This old way doesn't work, I believe.
  //G4ThreeVector r_global= G4ThreeVector
  //(y[0], 
  //y[1],
  //y[2]);
  
  //G4ThreeVector r_local = G4ThreeVector
  //(fpMatrix->colX() * r_global - fOrigin.x(),
  //fpMatrix->colY() * r_global - fOrigin.y(),
  //fpMatrix->colZ() * r_global - fOrigin.z());

  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z());

  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

  //indexed dipole
  //B = ( r / r_0 ) ^ -n  * B_0
  G4double B_0 = -fBend;
  //G4double r   = sqrt( pow(r_local.z(), 2) + pow(r_local.y(), 2) );
  G4double r_0 = 8.4 * m;
  //G4double n   =-1.25;
  G4double n   =-1.263;//SNAKE
  //G4double indexed_field = -pow( r / r_0, -n ) * B_0;
  // y is -x and x is y
  //G4double indexed_field_x = ( 1. - n * ( r - r_0 ) / r_0  - 0.5 * n * r_local.x() * r_local.x()  / r_0 / r_0 ) * B_0;
  //G4double indexed_field_x = ( 1. - n * ( r - r_0 ) / r_0 ) * B_0;
  //G4double indexed_field_y = ( 0. + n * r_local.x() / r_0 ) * B_0;
  //G4double chi = 0.;
  //if( r_local.z() < 0.){
  //chi = -atan( r_local.z() / abs(r_local.y()) );
  //}else{
  //chi = atan( r_local.z() / abs(r_local.y()) );
  //}
  //G4double indexed_field_x = ( 1. + n * ( r_0 - r_local.y() ) / r_0 + 0.5 * n * r_local.x() * r_local.x() / r_0 / r_0 ) * B_0;
  //G4double indexed_field_y = ( 0. + n *         r_local.x()   / r_0                                                   ) * B_0;// * cos( chi );
  //G4double indexed_field_z = 0.;//( 0. - n * r_local.x() / r_0 ) * B_0 * sin( chi );
  /*
  G4cout << "r, r_0, x, y, z (local)" << G4endl;
  G4cout << r << " " << r_0 << " " << r_local.x() << " " << r_local.y() << " " << r_local.z() << G4endl;
  G4cout << "chi (rad), chi (deg), Bx, By, Bz (local)" << G4endl;
  G4cout << chi << " " << chi * 180. / 3.141592654 << " " << indexed_field_x << " " << indexed_field_y << " " << indexed_field_z << G4endl;
  */
  float indexed_fringe_field_x = 0.;
  float indexed_fringe_field_y = 0.;
  float indexed_fringe_field_z = 0.;

  //G4cout << snake_x << " " << snake_y << " " << snake_z << G4endl;


  double pi=acos(-1);
  float dipole_width = 0.25 * m;
  //float dipole_width = 0.4;
  //G4cout << B_0 << " ";
  float fact         = (float)B_0 / -tesla;
  if (0)
  {
    G4cout<<"dipole fringe input G4 fact="<<fact<<"  a="<<dipole_width<<G4endl;
    float snake_x = - r_local.y() - r_0;
    //float snake_y = - r_local.z() + sqrt( 3. ) / r_local.x();
    float snake_y = - r_local.z() + ( r_local.y() + r_0 ) * tan( pi / 6 );
    float snake_z =   r_local.x();
//    dfringe_( &snake_x               , &snake_y                , &snake_z,
//              &fact                  , &dipole_width           ,
//              &indexed_fringe_field_x, &indexed_fringe_field_y, &indexed_fringe_field_z);

    double indexed_fringe_field_x_tmp[3];

    dfringecpp( snake_x               , snake_y                , snake_z,
                fact                  , dipole_width           ,
                indexed_fringe_field_x_tmp);
    indexed_fringe_field_x = indexed_fringe_field_x_tmp[0];
    indexed_fringe_field_y = indexed_fringe_field_x_tmp[1];
    indexed_fringe_field_z = indexed_fringe_field_x_tmp[2];
//    cout<<"old f:"<<setw(20)<<indexed_fringe_field_x    <<setw(20)<<indexed_fringe_field_y    <<setw(20)<<indexed_fringe_field_z<<endl;
//    cout<<"new c:"<<setw(20)<<indexed_fringe_field_x_tmp[0]<<setw(20)<<indexed_fringe_field_x_tmp[1]<<setw(20)<<indexed_fringe_field_x_tmp[2]<<endl<<endl;
  }

  if( r_local.z() < -1.5 * m || r_local.z() > 8.0 * m){
    //do nothing
  }else if( r_local.z() < 1.0 * m && r_local.z() > -1.5 * m ){//conservatively large since there is a 30 deg tilt which puts events in fringe past z = 0 (den)
    //G4cout << tesla << " " << fact << G4endl;
    float snake_x = - r_local.y() - r_0;
    //float snake_y = - r_local.z() + sqrt( 3. ) / r_local.x();
    float snake_y = - r_local.z() + ( r_local.y() + r_0 ) * tan( pi / 6 );
    float snake_z =   r_local.x();

/*
    dfringe_( &snake_x               , &snake_y                , &snake_z,
              &fact                  , &dipole_width           ,
              &indexed_fringe_field_x, &indexed_fringe_field_y, &indexed_fringe_field_z);
*/
    double indexed_fringe_field_x_tmp[3];

    dfringecpp( snake_x               , snake_y                , snake_z,
                fact                  , dipole_width           ,
                indexed_fringe_field_x_tmp);

    indexed_fringe_field_x = indexed_fringe_field_x_tmp[0];
    indexed_fringe_field_y = indexed_fringe_field_x_tmp[1];
    indexed_fringe_field_z = indexed_fringe_field_x_tmp[2];
//    cout<<"old f:"<<setw(20)<<indexed_fringe_field_x    <<setw(20)<<indexed_fringe_field_y    <<setw(20)<<indexed_fringe_field_z<<endl;
//    cout<<"new c:"<<setw(20)<<indexed_fringe_field_x_tmp[0]<<setw(20)<<indexed_fringe_field_x_tmp[1]<<setw(20)<<indexed_fringe_field_x_tmp[2]<<endl<<endl;

  }else{
    float tilt    = pi / 4.;
    float snake_x =   r_local.z() * sin( tilt ) - r_local.y() * cos( tilt ) - r_0;                     //-y
    //float snake_y = - r_local.z() + sqrt( 3. ) / r_local.x();
    float snake_y =   r_local.z() * cos( tilt ) + r_local.y() * sin( tilt ) - snake_x * tan( pi / 6 ); //+z
    float snake_z = - r_local.x();                                                                     //-x

    //G4cout << "Snake: " << snake_x << " " << snake_y << " " << snake_z << G4endl;
//    dfringe_( &snake_x               , &snake_y                , &snake_z,
//              &fact                  , &dipole_width           ,
//              &indexed_fringe_field_x, &indexed_fringe_field_y, &indexed_fringe_field_z);

    double indexed_fringe_field_x_tmp[3];

    dfringecpp( snake_x               , snake_y                , snake_z,
                fact                  , dipole_width           ,
                indexed_fringe_field_x_tmp);

    indexed_fringe_field_x = indexed_fringe_field_x_tmp[0];
    indexed_fringe_field_y = indexed_fringe_field_x_tmp[1];
    indexed_fringe_field_z = indexed_fringe_field_x_tmp[2];

//    cout<<"old f:"<<setw(20)<<indexed_fringe_field_x    <<setw(20)<<indexed_fringe_field_y    <<setw(20)<<indexed_fringe_field_z<<endl;
//    cout<<"new c:"<<setw(20)<<indexed_fringe_field_x_tmp[0]<<setw(20)<<indexed_fringe_field_x_tmp[1]<<setw(20)<<indexed_fringe_field_x_tmp[2]<<endl<<endl;
  }

  double rad_displ = sqrt(r_local.z()*r_local.z() + r_local.y()*r_local.y())-r_0;
  
  double  ndx_ratio=1./(1.+n*rad_displ/r_0);
  
//  if ((fabs(y[1]/cm) < 20.) && (y[2]/cm>970/sqrt(2.)) && (y[2]/cm<1170/sqrt(2.)) )  G4cout << "coord: "<<y[0]/cm << "  " << y[1]/cm << "  " << y[2]/cm <<"   local:"<<r_local.x()/cm<<"  "<<r_local.y()/cm<<"  "<<r_local.z()/cm<<"     displacement="<<rad_displ/m<<"    ratio="<<ndx_ratio<<G4endl;
  
  G4ThreeVector B_local = G4ThreeVector
    ( - indexed_fringe_field_z * tesla*ndx_ratio,//snake
      - indexed_fringe_field_x * tesla*ndx_ratio,//snake
        indexed_fringe_field_y * tesla*ndx_ratio//snake
     );
  
  
  
  G4ThreeVector B_global = G4ThreeVector
    (fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local);

  B[0] = -1.* B_global.x() ;
  B[1] = -1.* B_global.y() ;
  B[2] = -1.* B_global.z() ;
  
  //G4cout << B_0 << G4endl;
  //G4cout << r_local.x() << " " << r_local.y() << " " << r_local.z() << G4endl;
  //G4cout << B[0] << " " << B[1] << " " << B[2] << G4endl;

  //Print what the field would be, but don't use it: CAREFUL!!
  //B[0] = 0;
  //B[1] = 0;
  //B[2] = 0;
#ifdef DEBUG_BFIELD_DIPOLE_FRINGE

  //G4cout << deg << G4endl;
//vardan  G4cout << "Dipole fringe: global, local, local field: " << G4endl;
  //G4cout << y[0] << " " << y[1] << " " << y[2] << " " << G4endl;
//vardan  G4cout << r_global.x() << " " << r_global.y() << " " << r_global.z() << " " << G4endl;
  //G4cout << "Local coordinates:";
//vardan  G4cout << r_local.x() << " " << r_local.y() << " " << r_local.z() << " " << G4endl;
  //G4cout << sqrt(r_local.x() * r_local.x() + r_local.y() * r_local.y() + r_local.z() * r_local.z()) << G4endl;
  //G4cout << "Bend strength: " << fBend << G4endl;
  //G4cout << "Origin coordinates:"
  //	 << fOrigin.x()  << " " << fOrigin.y()  << " " << fOrigin.z()  << G4endl;
  //G4cout << "DOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOP!" << G4endl;
  //G4cout << "Rotation matrix:" << G4endl;
  //G4cout << fpMatrix->xx() << " " << fpMatrix->xy() << " " << fpMatrix->xz() << G4endl;
  //G4cout << fpMatrix->yx() << " " << fpMatrix->yy() << " " << fpMatrix->yz() << G4endl;
  //G4cout << fpMatrix->zx() << " " << fpMatrix->zy() << " " << fpMatrix->zz() << G4endl;
  //G4cout << "Transpose matrix:" << G4endl;
  //G4cout << fpMatrix->colX().x() << " " << fpMatrix->colX().y() << " "<< fpMatrix->colX().z() << G4endl;
  //G4cout << fpMatrix->colY().x() << " " << fpMatrix->colY().y() << " "<< fpMatrix->colY().z() << G4endl;
  //G4cout << fpMatrix->colZ().x() << " " << fpMatrix->colZ().y() << " "<< fpMatrix->colZ().z() << G4endl;
  //G4cout << "Rotation matrix:" << G4endl;
  //G4cout << fpMatrix->rowX().x() << " " << fpMatrix->rowX().y() << " "<< fpMatrix->rowX().z() << G4endl;
  //G4cout << fpMatrix->rowY().x() << " " << fpMatrix->rowY().y() << " "<< fpMatrix->rowY().z() << G4endl;
  //G4cout << fpMatrix->rowZ().x() << " " << fpMatrix->rowZ().y() << " "<< fpMatrix->rowZ().z() << G4endl;
  //G4cout << "Local field:" << G4endl;
//vardan  G4cout << B_local .x() << " " << B_local .y() << " " << B_local .z() << " " << tan(r_local.y() / r_local.x()) << G4endl;
  //G4cout << "Global field:" << G4endl;
  //G4cout << B_global.x() << " " << B_global.y() << " " << B_global.z() << G4endl;
#endif

}


void BField_Dipole_Fringe::dfringecpp(double x,double y,double z,double fact,double a,double fb[3]) const
{
      double delta;
      double s;
      double b[5][5];
//      b(-2:2,-2:2)
      double c0[6] = {0.2383,1.7395,-0.4768,0.5288,-0.1299,0.0222};
      delta=80.;
//c      print *, x,y,z,fact,a,fx,fy,fz
//c
//c   calculates  entrance fringe field of dipole magnet
//c   in the same manner as raytrace, for use in snake. jjl 2/17/87
//c   n.b. can only be used for bending in the x-y plane.
//c
//c          fact = central field
//c             a = gap size (radius)
//c         x,y,z = coordinates relative to magnet entrance
//c                  (set entrance=0 in the relative ref frame)
//c            fz = returned field value
//c         c0(i) = fringe field coefficients
//c
//c            fz = fact/(1+exp(s))
//c             s = c0(0) + c0(1)*(y/a) + c0(2)*((y/a)**2) ....etc.
//c****************************************************************
//c
//c  set up grid of b's for expansion (see raytrace manual p. 11-12)
//c
      for (int j =-2; j<=2; j++)
      for (int k =-2; k<=2; k++)
      {
        if((abs(j)+abs(k)) <  3)
        {
          s=c0[0];
          for (int i =1; i<=5; i++)
          {
            s=s+(c0[i]*(pow(((y+(j*delta))/a),i)));
          }
        }
        b[j+2][k+2]=fact/(1+exp(s));
      }
//c  calculate fields
//c        1         2         3         4         5         6         7
      fb[2]=b[0+2][0+2]
        -(((z*z)/(delta*delta))*
           (((2./3.)*(b[1+2][0+2]+b[-1+2][0+2]+b[0+2][1+2]+b[0+2][-1+2]-(4.*b[0+2][0+2])))
          -((1./24.)*(b[2+2][0+2]+b[-2+2][0+2]+b[0+2][2+2]+b[0+2][-2+2]-(4.*b[0+2][0+2]))))
        +(((z*z*z*z)/(delta*delta*delta*delta))*
           (((-1./6.)*(b[1+2][0+2]+b[-1+2][0+2]+b[0+2][1+2]+b[0+2][-1+2]-(4.*b[0+2][0+2])))
           +((1./24.)*(b[2+2][0+2]+b[-2+2][0+2]+b[0+2][2+2]+b[0+2][-2+2]-(4.*b[0+2][0+2])))
           +((1./12.)*(b[1+2][1+2]+b[-1+2][1+2]+b[1+2][-1+2]+b[-1+2][-1+2]))
           -((1./6.0)*(b[1+2][0+2]+b[-1+2][0+2]+b[0+2][1+2]+b[0+2][-1+2]-(2.*b[0+2][0+2]))))));
//c        1         2         3         4         5         6         7
      fb[1]=((z/delta)*
             (((2./3.)*(b[1+2][0+2]-b[-1+2][0+2]))
             -((1./12.)*(b[2+2][0+2]-b[-2+2][0+2]))))
         +(((z*z*z)/(delta*delta*delta))*
             (((1./6.)*(b[1+2][0+2]-b[-1+2][0+2]))
            -((1./12.)*(b[2+2][0+2]-b[-2+2][0+2]))
            -((1./12.)*(b[1+2][1+2]+b[1+2][-1+2]-b[-1+2][1+2]-b[-1+2][-1+2]
                 -(2.*b[1+2][0+2])   +  (2.*b[-1+2][0+2])))));
//c        1         2         3         4         5         6         7
      fb[0]=((z/delta)*
             (((2./3.)*(b[0+2][1+2]-b[0+2][-1+2]))
            -((1./12.)*(b[0+2][2+2]-b[0+2][-2+2]))))
        +(((z*z*z)/(delta*delta*delta))*
             (((1./6.)*(b[0+2][1+2]-b[0+2][-1+2]))
            -((1./12.)*(b[0+2][2+2]-b[0+2][-2+2]))
            -((1./12.)*(b[1+2][1+2]+b[-1+2][1+2]-b[1+2][-1+2]-b[-1+2][-1+2]
                -(2.*b[0+2][1+2])+(2.*b[0+2][-1+2])))));
//c      print *, fx, fy, fz
}
