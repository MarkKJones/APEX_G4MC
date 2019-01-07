#ifndef BFIELD_QUAD_SNAKE_HH
#define BFIELD_QUAD_SNAKE_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class BField_Quad_Snake : public G4MagneticField
{
  public: // with description
    BField_Quad_Snake(G4double Length,
                    G4double pGradient1, 
		    G4double qrad1, 
		    G4ThreeVector pOrigin, 
		    G4RotationMatrix* pMatrix,
		    G4int QuadNumber);
   ~BField_Quad_Snake();
    void bplscpp ( int , double , double , double *) const;
    //void bplscpp ( int , double , double , double , double , double , double , double , double, double );
    void mpolescpp(int , double , double , double , double ,double ,double *) const;

    void GetFieldValue(const G4double yTrack[],
                             G4double B[]     ) const;
  private:
   float fGradient1;
   bool onceonly;
   G4ThreeVector fOrigin;
   G4RotationMatrix* fpMatrix;
   float fRadius1;
   G4double flength;
   G4int pQuadNumber;
   int nQ1Sos;

//   G4int Q_Num;
//   double re, g1, g2, g3, g4, g5, g6;
};
#endif

