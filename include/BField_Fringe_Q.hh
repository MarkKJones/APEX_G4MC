#ifndef BFIELD_FRINGE_Q_HH
#define BFIELD_FRINGE_Q_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class BField_Fringe_Q : public G4MagneticField
{
  public: // with description
    BField_Fringe_Q(G4int num,
                    G4double pGradient1, 
		    G4double qrad1, 
		    G4ThreeVector pOrigin, 
		    G4RotationMatrix* pMatrix,
		    G4int QuadNumber);
   ~BField_Fringe_Q();

    void bplscpp ( int , double , double , double *) const;
    void mpolescpp(int , double , double , double , double ,double ,double *) const;

    void GetFieldValue(const G4double yTrack[],
                             G4double B[]     ) const;
  private:
   float fGradient1;
   bool onceonly;
   G4ThreeVector fOrigin;
   G4RotationMatrix* fpMatrix;
   G4int pQuadNumber;
   float fRadius1;
   G4int fnum;
   G4int nQ1Sos;

};
#endif

