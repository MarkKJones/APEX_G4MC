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
// $Id: BFIELD_FRINGE_Q.cc 69786 2016-07-27 09:38:51Z gcosmo $
//
// -------------------------------------------------------------------

/*
snake
        xdum=xyz(1)
        ydum=-xyz(2)
        zdum=xyz(3)
        write(6,*)'1     fact1=',fact1,'grad1=',qrad1
        call mpoles(1,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
        f(2)=-f(2)
        return
*/


#include "UsageManager.hh"
#include "BField_Quad_Snake.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>
#include <iostream>
#include <iomanip>
#define DEBUG_BFIELD_QUAD_SNAKE 1
/*
extern"C" {
  void mpoles1_(short int* n, float* x, float* y, float* z,
                 float* fact1, float* qrad1,
                 float* bx,float* by,float* bz);
}
*/
extern UsageManager* gConfig;

static G4RotationMatrix IdentityMatrix;

using namespace std;

BField_Quad_Snake::BField_Quad_Snake(G4double Length, G4double pGradient1, G4double qrad1, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4int QuadNumber)
{
/*

  double x=140., y, z=0.1, ngradn=1., rad=150., bb[3]={0};

  for (double yy=1000; yy>-900; yy-=10.)
  {
    mpolescpp(1, x, z, yy, ngradn,rad, bb);
    cout<<z<<setw(15)<<x<<setw(15)<<yy<<setw(15)<<bb[0]<<setw(15)<<bb[1]<<setw(15)<<bb[2]<<endl;
  }
*/


   gConfig->GetParameter("Q1Sos",nQ1Sos);


   fGradient1   = pGradient1 ;
   fOrigin      = pOrigin ;
   fpMatrix     = pMatrix ;
   fRadius1     = qrad1;
   flength  = Length;
   pQuadNumber = QuadNumber;

}

BField_Quad_Snake::~BField_Quad_Snake()
{
  delete fpMatrix;
}

void BField_Quad_Snake::GetFieldValue( const G4double y[7],
				 G4double Bsn[3]  ) const  
{


  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z());

  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
//    (14.*cm,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

   float bx=0, by=0, bz=0;
   float grd1=fGradient1*1000.;
   float rad1=fRadius1;
   short int n;
   
   float snake_xx, snake_yy, snake_zz;


   G4ThreeVector B_local;
   B_local = G4ThreeVector( 0.0, 0.0, 0.0);

   float r_tmp = sqrt(r_local.y()*r_local.y() + r_local.x()*r_local.x());
   float dz=fabs(r_local.z());

   if (r_tmp < fRadius1)
   if (r_tmp > 0.0)
   if ( fabs(dz) <0.5*flength)
   {
     if (r_local.z()<=0)
     {
       snake_xx =   r_local.y();
       snake_zz =   r_local.x();
       snake_yy =   -1.*(0.5*flength+r_local.z());
       n=1;
     }

     if (r_local.z()>0)
     {
       snake_xx =   r_local.y();
       snake_yy =   r_local.z()-(0.5*flength);
       snake_zz =   r_local.x();
       n=3;
     }

//     G4cout<<"    fringe input "<<n<<"  x="<<snake_xx<<"  y="<<snake_yy<<"  z="<<snake_zz<<"  grad="<<grd1<<"   rad="<<rad1<<G4endl;

/*
     mpoles1_(&n, &snake_xx, &snake_zz, &snake_yy,
              &grd1, &rad1,
              &bx, &bz, &by);
*/
//       mpolescpp(1, x, z, yy, ngradn,rad, bx, by, bz);
     double b_new[3]={0};
     mpolescpp(n, snake_xx, snake_zz, snake_yy, grd1, rad1, b_new);

     bx=b_new[0];
     by=b_new[2];
     bz=b_new[1];

//     cout<< snake_xx<<setw(15)<< snake_zz<<setw(15)<< snake_yy<<setw(15)<< grd1<<setw(15)<< rad1<<endl;
//     cout<<"compare old"<<setw(15)<<bx<<setw(15)<<by<<setw(15)<<bz<<endl;
//     cout<<"compare new"<<setw(15)<<b_new[0]<<setw(15)<<b_new[2]<<setw(15)<<b_new[1]<<endl;
     if (n==1)
     if ( (bx>-1000.) && (by>-1000.) && (bz>-1000.) && (bx<1000.) && (by<1000.) && (bz<1000.) )
     {
       B_local = G4ThreeVector
        ( bz * tesla,//snake
          bx * tesla,//snake
          -1.*by * tesla//snake
        );
     }

     if (n==3)
     if ( (bx>-1000.) && (by>-1000.) && (bz>-1000.) && (bx<1000.) && (by<1000.) && (bz<1000.) )
     {
       B_local = G4ThreeVector
        ( bz * tesla,//snake
          bx * tesla,//snake
          by * tesla//snake
        );
     }
//     G4cout<<"snake:      local[cm]: "<<r_local.x()/cm<<", "<<r_local.y()/cm<<", "<<r_local.z()/cm<<"      b_local:"<<B_local.x()/tesla<<","<<B_local.y()/tesla<<","<<B_local.z()/tesla<<G4endl;
//     G4cout<<"h->Fill("<<r_local.z()/cm<<", "<<B_local.y()/tesla<<");"<<G4endl;
//     G4cout<<"exp:"<<y[0]/cm<<", "<<y[1]/cm<<", "<<y[2]/cm<<"      local[mm]: "<<r_local.x()/cm<<", "<<r_local.y()/cm<<", "<<r_local.z()/cm<<"      b_local:"<<B_local.x()/tesla<<","<<B_local.y()/tesla<<","<<B_local.z()/tesla<<G4endl;
   }

  G4ThreeVector B_global = G4ThreeVector
    (fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local);

  Bsn[0] = -1*B_global.x() ;
  Bsn[1] = -1*B_global.y() ;
  Bsn[2] = -1*B_global.z() ;
  
#ifdef DEBUG_BField_Quad_Snake

#endif

}


void BField_Quad_Snake::mpolescpp(int in, double nx, double ny, double nz, double ngradn,double nrad,double nb[3]) const
{

//      subroutine mpoles1(in,nx,ny,nz,ngradn,nrad,nbx,nby,nbz)
//c stolen from raytrace jjl 10/8/86
//c****                                                                   
//c**** calculation of multipole(poles) field components                  
//c****                                                                   
//c****                                                                   
//c****                                                                   
//c**** 2 - quadrupole  (grad1)                                           
//c**** 3 - hexapole    (grad2)                                           
//c**** 4 - octapole    (grad3)                                           
//c**** 5 - decapole    (grad4)                                           
//c**** 6 - dodecapole  (grad5)                                           
//c****                                                                   
//c****                                                                   
//      implicit real*8(a-h,o-z)                                          
//cvardan      real nx,ny,nz,ngrad(5),nrad,nbx,nby,nbz
/*
      double c0=0.1039;
      double c1=6.27108;
      double c2=-1.51247;
      double c3=3.59946;
      double c4=-2.1323;
      double c5=1.7230;
*/
      double ngrad[5];

      if (nrad == 0.)
      {
        cout<<" error in mpoles1,  nrad= 0."<<endl;
//        call exit(0)
      }
      ngrad[0] = ngradn;
      ngrad[1] = 0;
      ngrad[2] = 0;
      ngrad[3] = 0;
      ngrad[4] = 0;

      double grad1 = -ngrad[0]/nrad;
      double grad2 =  ngrad[1]/nrad*nrad;
      double grad3 = -ngrad[2]/nrad*nrad*nrad;
      double grad4 =  ngrad[3]/nrad*nrad*nrad*nrad;
      double grad5 = -ngrad[4]/nrad*nrad*nrad*nrad*nrad;

      double rad=nrad;
      double d = 2. * rad;
      double frh  = 1.;
      double fro  = 1.;
      double frd  = 1.;
      double frdd = 1.;
      double dh  = frh *d;
      double dov  = fro *d;
      double dd  = frd *d;
      double ddd = frdd*d;
      double x = nx;
      double y = ny;
      double z = nz;
      double x2 = x*x;
      double x3 = x2*x;
      double x4 = x3*x;
      double x5 = x4*x;
      double x6 = x5*x;
      double x7 = x6*x;
      double y2 = y*y;
      double y3 = y2*y;
      double y4 = y3*y;
      double y5 = y4*y;
      double y6 = y5*y;
      double y7 = y6*y;

      double s = z/d;
      double re , g1 , g2 , g3 , g4 , g5 , g6 ;
      double gg[7]={0};
      bplscpp( 2, d, s, gg );
      re= gg[0];
      g1= gg[1];
      g2= gg[2];
      g3= gg[3];
      g4= gg[4];
      g5= gg[5];
      g6= gg[6];
      if(0) cout<<grad2<<grad3<<grad4<<grad5<<dh<<dov<<dd<<ddd<<endl;
//      cout<<"ger@   "<<d<<setw(15)<<s<<setw(15)<<re<<setw(15)<<g1<<setw(15)<<g2<<setw(15)<<g3<<setw(15)<<g4<<setw(15)<<g5<<setw(15)<<g6<<endl;
//      cout<<"---->"<<setw(15)<<tt[0]<<endl;
      double b2x = grad1*( re*y - (g2/12.)*(3.*x2*y + y3) +
                   (g4/384.)*(5.*x4*y + 6.*x2*y3 + y5 ) -
                   (g6/23040.)*(7.*x6*y + 15.*x4*y3 + 9.*x2*y5 + y7)  );
      double b2y = grad1*( re*x - (g2/12.)*(x3 + 3.*x*y2) +
                   (g4/384.)*(x5 + 6.*x3*y2 + 5.*x*y4 ) -
                   (g6/23040.)*(x7 + 9.*x5*y2 + 15.*x3*y4 + 7.*x*y6) );
      double b2z = grad1*( g1*x*y - (g3/12.)*(x3*y + x*y3 ) +
                   (g5/384.)*(x5*y +2.*x3*y3 + x*y5)  );

      nb[0]=b2x;
      nb[1]=b2y;
      nb[2]=b2z;
//      cout<<nbx<<"    "<<nby<<"    "<<nbz<<endl;
}

void BField_Quad_Snake::bplscpp ( int igp, double d, double s, double gg[6]) const
{

      double c0=0.1039;
      double c1=6.27108;
      double c2=-1.51247;
      double c3=3.59946;
      double c4=-2.1323;
      double c5=1.7230;
/*
      double c0=0.1122;
      double c1=6.2671;
      double c2=-1.4982;
      double c3=3.5882;
      double c4=-2.1209;
      double c5=1.723;
*/
      if ((pQuadNumber==1) && (nQ1Sos) )
      {
//        std::cout<<"quad sos:"<<nQ1Sos<<std::endl;
        c0=0.2486;
        c1=5.3626;
        c2=-2.412;
        c3=0.9874;
        c4=0.0;
        c5=0.0;
      }

      double re, g1, g2, g3, g4, g5, g6;
      double s2 = s*s;
      double s3 = s2*s;
      double s4 = s2*s2;
      double s5 = s4*s;
      double cs = c0 + c1*s + c2*s2 + c3*s3 + c4*s4 + c5*s5;
      double cp1 =(c1 + 2.*c2*s + 3.*c3*s2 + 4.*c4*s3 + 5.*c5*s4) / d;
      double cp2 = (2.*c2 + 6.*c3*s + 12.*c4*s2 + 20.*c5*s3  ) / (d*d);
      double cp3 = ( 6.*c3 + 24.*c4*s + 60.*c5*s2 ) / (d*d*d);
      double cp4 = ( 24.*c4 + 120.*c5*s ) / (d*d*d*d);

      double cp5 = 120.*c5/(d*d*d*d*d);



      if( fabs(cs) > 70. )  cs = 70. * fabs(cs)/cs;
      double e = exp(cs);
//      cout<<cs<<endl;
      re = 1./(1. + e);
      double ere = e*re;
      double ere1= ere*re;
      double ere2= ere*ere1;
      double ere3= ere*ere2;
      double ere4= ere*ere3;

      double ere5= ere*ere4;
      double ere6= ere*ere5;


      double cp12 = cp1*cp1;
      double cp13 = cp1*cp12;
      double cp14 = cp12*cp12;
      double cp22 = cp2*cp2;
//c****
      double cp15 = cp12*cp13;
      double cp16 = cp13*cp13;
      double cp23 = cp2*cp22;
      double cp32 = cp3*cp3;
//c****
//c****
      if( igp == 6 ) return;
      g1 = -cp1*ere1;
//c****
//c****
      if( igp == 5 ) return;
      if( igp != 4 )
      {
        g2 =-1.*( cp2+cp12   )*ere1    + 2.*cp12 * ere2;
        g3 =-1.*(cp3 + 3.*cp1*cp2 + cp13  ) * ere1 + 6.*(cp1*cp2 + cp13)*ere2 - 6.*cp13*ere3;
        if( igp == 3 ) return;
      };
      g4 = -(cp4 + 4.*cp1*cp3 + 3.*cp22 + 6.*cp12*cp2 + cp14)*ere1  +
           (8.*cp1*cp3 + 36.*cp12*cp2 + 6.*cp22 + 14.*cp14)*ere2    -
           36.*(cp12*cp2 + cp14)*ere3 + 24.*cp14*ere4;

      if( igp != 2 ) return;
      g5 = (-cp5 - 5.*cp1*cp14 - 10.*cp2*cp3 - 10.*cp12*cp3 -
          15.*cp1*cp22 - 10.*cp13*cp2 - cp15)*ere1 +
          (10.*cp1*cp4 +20.*cp2*cp3 +60.*cp12*cp3 + 90.*cp1*cp22 +
          140.*cp13*cp2 +30.*cp15)*ere2 + (-60.*cp12*cp3 -
          90.*cp1*cp22 - 360.*cp13*cp2 - 150.*cp15)*ere3 +
          (240.*cp13*cp2 +240.*cp15)*ere4 + (-120.*cp15)*ere5;
      g6 = (-6.*cp1*cp5 - 15.*cp2*cp4 - 15.*cp12*cp4 - 10.*cp32 -
          60.*cp1*cp2*cp3 - 20.*cp13*cp3 - 15.*cp23 - 45.*cp12*cp22 -
          15.*cp14*cp2 - cp16)*ere1 + (12.*cp1*cp5 + 30.*cp2*cp4 +
          90.*cp12*cp4 +20.*cp32 + 360.*cp1*cp2*cp3 +280.*cp13*cp3 +
          90.*cp23 + 630.*cp12*cp22 + 450.*cp14*cp2 + 62.*cp16)*ere2 +
          (-90.*cp12*cp4 - 360.*cp1*cp2*cp3 -720.*cp13*cp3 -90.*cp23 -
          1620.*cp12*cp22 -2250.*cp14*cp2 - 540.*cp16)*ere3 +
          (480.*cp13*cp3 + 1080.*cp12*cp22 + 3600.*cp14*cp2 +
          1560.*cp16)*ere4 + (-1800.*cp14*cp2 - 1800.*cp16)*ere5 +
          720.*cp16*ere6;
//        g1 = 1.;
          gg[0]=re;
          gg[1]=g1;
          gg[2]=g2;
          gg[3]=g3;
          gg[4]=g4;
          gg[5]=g5;
          gg[6]=g6;
//      cout<<"                             ger@  tazuc"<<setw(15)<<g1<<setw(15)<<g2<<setw(15)<<g3<<setw(15)<<g4<<setw(15)<<g5<<setw(15)<<g6<<endl;

//      return;
}
