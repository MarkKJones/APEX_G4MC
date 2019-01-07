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
// $Id: BFIELD_FRINGE_QQ.cc 69786 2013-05-15 09:38:51Z gcosmo $
//
// -------------------------------------------------------------------

#include "BField_Fringe_QQ.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>
#define DEBUG_BFIELD_FRINGE_QQ 1

extern"C" {
  void qqfringe_(double* x, double* y, double* z, double* dy,
                 double* bx,double* by,double* bz,
                 double* fact1,double* fact2,
                 double* qrad1,double* qrad2);
}

static G4RotationMatrix IdentityMatrix; 

BField_Fringe_QQ::BField_Fringe_QQ(G4double pGradient1, G4double pGradient2, G4double qrad1, G4double qrad2, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix)
{
   fGradient1   = pGradient1 ;
   fGradient2   = pGradient2 ;
   fOrigin      = pOrigin ;
   fpMatrix     = pMatrix ;
   fRadius1      = qrad1;
   fRadius2      = qrad2;
   dist_z = 117.2*cm;

}

BField_Fringe_QQ::~BField_Fringe_QQ()
{
  delete fpMatrix;
}

void BField_Fringe_QQ::GetFieldValue( const G4double y[7],
				 G4double B[3]  ) const  
{

  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z());

  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

   double snake_x = - r_local.y();
   double snake_y = - r_local.z();
   double snake_z =   r_local.x();
   
   double bx=0, by=0, bz=0;
   double dist=dist_z;
   double grd1=fGradient1;
   double grd2=fGradient2;
   double rad1=fRadius1;
   double rad2=fRadius2;

  qqfringe_(&snake_x, &snake_y, &snake_z, &dist,
                 &bx, &by, &bz,
                 &grd1,&grd2,
                 &rad1, &rad2);

  G4ThreeVector B_local = G4ThreeVector
    ( bz * tesla,//snake
      bx * tesla,//snake
      by * tesla//snake
     );
  
  G4ThreeVector B_global = G4ThreeVector
    (fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local);

  B[0] = B_global.x() ;
  B[1] = B_global.y() ;
  B[2] = B_global.z() ;
  
#ifdef DEBUG_BFIELD_FRINGE_QQ
  G4cout << "Dipole fringe: global, local, local field: " << G4endl;
  G4cout << r_global.x() << " " << r_global.y() << " " << r_global.z() << " " << G4endl;
  G4cout << r_local.x() << " " << r_local.y() << " " << r_local.z() << " " << G4endl;
  G4cout << B_local .x() << " " << B_local .y() << " " << B_local .z() << " " << tan(r_local.y() / r_local.x()) << G4endl;
#endif

}
