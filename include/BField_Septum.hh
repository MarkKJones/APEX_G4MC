// BField_Septum.h: interface for the BField_Septum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BFIELD_Septum_H)
#define BFIELD_Septum_H
#include <math.h>
#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Vector/ThreeVector.h>
using namespace CLHEP;

class BField_Septum
{
public:
        BField_Septum(double pMomentumL=0,double pMomentumR=0,
          const char *mapfile="g2p_septumfield.dat");
        virtual ~BField_Septum();
        void GetBField(double Pos[3],double B[3]);

private:
        void ReadMap(const char *filename);

private:
        //mBField[indexX][indexY][indexZ][0] x
        //mBField[indexX][indexY][indexZ][1] y
        //mBField[indexX][indexY][indexZ][2] z
        //mBField[indexX][indexY][indexZ][3] Bx
        //mBField[indexX][indexY][indexZ][4] By
        //mBField[indexX][indexY][indexZ][5] Bz
        double ****mBField;
        double mLHRSMomentum;
        double mRHRSMomentum;
//        float Btxt[101][51][201][3] ;
        float ****Btxt;
        float ****BCoord;
        int    nx;
        int    ny;
        int    nz;

        //parameters from ini file
//        CLHEP::HepRotation *mRotF2L;

};
typedef BField_Septum HRSSeptumField;

#endif // !defined(BFIELD_Septum_H)
