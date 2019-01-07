#ifndef BFIELD_FRINGE_QQ_HH
#define BFIELD_FRINGE_QQ_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class BField_Fringe_QQ : public G4MagneticField
{
  public: // with description
    BField_Fringe_QQ(G4double pGradient1, 
		    G4double pGradient2, 
		    G4double qrad1, 
		    G4double qrad2, 
		    G4ThreeVector pOrigin, 
		    G4RotationMatrix* pMatrix);
   ~BField_Fringe_QQ();

    void GetFieldValue(const G4double yTrack[],
                             G4double B[]     ) const;
  private:
   double fGradient1;
   G4double fGradient2;
   G4ThreeVector fOrigin;
   G4RotationMatrix* fpMatrix;
   G4double fRadius1;
   G4double fRadius2;
   G4double dist_z;
   
};
#endif

