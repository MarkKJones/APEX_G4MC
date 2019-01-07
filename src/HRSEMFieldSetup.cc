// ********************************************************************
//
// $Id: HRSEMField.hh,v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//
//   User Field class Setup implementation.
//
#include "HRSEMFieldSetup.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "BField_Dipole.hh"
#include "BField_Quad.hh"
#include "BField_Quad_Snake.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "UsageManager.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ios.hh"
//#include "QuadFringe.hh"

#include <iostream>
using namespace std;

extern UsageManager* gConfig;

//////////////////////////////////////////////////////////////////////////
HRSEMFieldSetup* HRSEMFieldSetup::fHRSEMFieldSetup=0;
HRSEMFieldSetup* HRSEMFieldSetup::GetHRSEMFieldSetup()
{ 
	if(!fHRSEMFieldSetup)  
	{
		G4cout<<"HRSEMFieldSetup is not initialized yet...exit...\n";
		exit(-99);
	}
	return fHRSEMFieldSetup; 
}

//////////////////////////////////////////////////////////////////////////
//
HRSEMFieldSetup::HRSEMFieldSetup()
: fChordFinder(0), fStepper(0), fIntgrDriver(0)
{
//  gConfig->GetArgument("KAPPA1",KAPPA1);
//  gConfig->GetArgument("KAPPA2",KAPPA2);
//  gConfig->GetArgument("KAPPA3",KAPPA3);
//  gConfig->GetArgument("DipField",DipField);

  gConfig->GetArgument("LHRSMomentum",mLHRSMomentum);
  gConfig->GetArgument("RHRSMomentum",mRHRSMomentum);
  mLHRSMomentum /= 1000.;
  mRHRSMomentum /= 1000.;
  gConfig->GetArgument("SnakeModel",mSnakeModel);
  gConfig->GetParameter("LHRSAngle",mLHRSAngle);
  mLHRSAngle*=deg;
  gConfig->GetParameter("RHRSAngle",mRHRSAngle);
  mRHRSAngle*=deg;
  gConfig->GetParameter("LSeptumAngle",mLSeptumAngle);
  gConfig->GetParameter("RSeptumAngle",mRSeptumAngle);
  gConfig->GetParameter("Q1Sos",nQ1Sos);

  G4cout << "HRS angles: " << mLHRSAngle << " " << mRHRSAngle << G4endl;
  G4cout << "HRS momentums: " << mLHRSMomentum << "  " << mRHRSMomentum << G4endl;
  KAPPA1=0.243838 * tesla * mLHRSMomentum/0.83756;;
  KAPPA2=-0.1934 * tesla * mLHRSMomentum/0.83756;;
  KAPPA3=0.17892 * tesla * mLHRSMomentum/0.83756;;
  DipField=-0.4192 * tesla * mLHRSMomentum / 1.063;
  if (nQ1Sos)
  {
      double a_sup = 0.1492;                    //Bore radius, in meters
      double l_sup = 0.9413;                    //Length of quad, in meters
      double a_sos = 0.12827;
      double l_sos = 0.70;
      KAPPA1  *= l_sup / l_sos * a_sos / a_sup;
  }

  cout<<"setup     KAPPA1="<<KAPPA1<<"     KAPPA2="<<KAPPA2<<"     KAPPA3="<<KAPPA3<<"     DipField="<<DipField<<endl;

/*
  KAPPA1=-1.*KAPPA1;
  KAPPA2=-1.*KAPPA2;
  KAPPA3=-1.*KAPPA3;
*/
  //G4cout << "Quad fringe?" << G4endl;
  //QuadFringe* fringe = new QuadFringe();
  //G4cout << "Quad fringe!" << G4endl;
  fHRSEMFieldSetup=this;

  //global EM field
  fEMfield = new HRSEMField();
  messenger = new HRSEMFieldSetupMessenger(this) ;
  fEquation = new G4EqMagElectricField(fEMfield);
  fMinStep  = 0.001*mm ; // minimal step of 1 miron, default is 0.01 mm, Nickie finely
  fStepperType = 4 ;     // ClassicalRK4 -- the default stepper
  
  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  
  fChordFinder = 0;   //will be set in UpdateField()
  UpdateField();
  
  //G4double snakemagnumber = 1. / -4.77577734;//unitless, STD - correct, I believe
  //G4double snakemagnumber = 1. / -4.77577734 * ( 4.00 / 1.063 ) ;//unitless, PREX tune B - correct, I believe
  G4double snakemagnumber = -4.77577734 / 1.063;//Nickie's calculation//this is just 4.77577734 / 1.063
  //snakemagnumber *= 0.83756;
  snakemagnumber *= 0.83762;
  //snakemagnumber *= 0.7;//This is tuned with hrstrans #1
  //snakemagnumber *= 1.33;//This is tuned with hrstrans when hrstrans is tuned to JLR
  //snakemagnumber *= 1.24;//This is my new tune, with snake d.dat problems fixed #2
  //snakemagnumber *= 0.965; #3
  
  if( mSnakeModel == 53 || mSnakeModel == 55){
    //snakemagnumber *= 1.063 / 2.2;
    snakemagnumber *= 1.063;    
  }
  
  G4int    quads_on = 1;
//  G4double KAPPA1 =  0.;
//  G4double KAPPA2 =  0.;
//  G4double KAPPA3 =  0.;
  //BELOW ARE SNAKE VALUES, NOT NIM VALUES                                                                                  
  double pTarget =             0.0  * cm;
  double pQ1en   = pTarget + 159.03  * cm;
  double pQ1ex   = pQ1en   +  94.13 * cm;
  double pQ2en   = pQ1ex   + 117.2  * cm;
  double pQ2ex   = pQ2en   + 182.66 * cm;
  if (nQ1Sos)
  {
    pQ1en   = pTarget + 171.1*cm;
    pQ1ex   = pQ1en   +  70.0 * cm;
  }
//  double pQ1en   = pTarget + 160.0  * cm;
//  double pQ1ex   = pQ1en   +  94.13 * cm;
//  double pQ2en   = pQ1ex   + 115.58 * cm;
//  double pQ2ex   = pQ2en   + 182.66 * cm;


  //ABOVE ARE SNAKE VALUES, NOT NIM VALUES  
  //Local field  FZB2, Q1
  //double pHallCenter2Q1Face=1.69*m;//NIM
  double pHallCenter2Q1Face=pQ1en;//SNAKE
  //double pQ1Length= 80*cm;//NIM
  double pQ1Length= nQ1Sos ? 70.    * cm : 94.13*cm;//SNAKE
  double pQ1Radius= nQ1Sos ? 12.827 * cm : 0.1492 * m;//SNAKE
  double q1shift = nQ1Sos ? 0.0 * m : 0.0 * m;
  double pQ1Pos_Z=(pHallCenter2Q1Face+94.13*cm/2.0+q1shift);//NIM
  if (nQ1Sos) pQ1Pos_Z=(pHallCenter2Q1Face+70.0*cm/2.0+q1shift);//NIM

  //double pHallCenter2Q2Face=3.74*m;//NIM
  double pHallCenter2Q2Face=pQ2en;//SNAKE
  //double pQ2Length=180*cm;//NIM
  double pQ2Length=182.66*cm;//SNAKE
  double pQ2Radius=  0.3 * m;//SNAKE
  double pQ2Pos_Z=(pHallCenter2Q2Face+pQ2Length/2.0);
  //double pQ3Length=180*cm;//NIM
  double pQ3Length=182.68*cm;//SNAKE
  double pQ3Radius=  pQ2Radius;//SNAKE

  if( quads_on == 1 && ( mSnakeModel == 49 || mSnakeModel > 50 ) ){
//    KAPPA1 =  0.243838  * tesla * mLHRSMomentum/0.83756;// / snakemagnumber / .1492 / m; //STD tune
//    KAPPA2 = -0.1934    * tesla * mLHRSMomentum/0.83756;;// / snakemagnumber / .300  / m; //STD tune
//    KAPPA3 = -0.17892   * tesla * mLHRSMomentum/0.83756;;// / snakemagnumber / .300  / m; //STD tune
//    KAPPA1 =  0.2445  * tesla * mLHRSMomentum/0.83756;// / snakemagnumber / .1492 / m; //STD tune
//    KAPPA2 = -0.1939  * tesla * mLHRSMomentum/0.83756;;// / snakemagnumber / .300  / m; //STD tune
//    KAPPA3 = -0.1794  * tesla * mLHRSMomentum/0.83756;;// / snakemagnumber / .300  / m; //STD tune
//   Q1 := -0.243838 ;//cosy
//   Q2 :=  0.1934 ; //cosy
//   Q3 :=  0.17892 ;//cosy
//    KAPPA1 = -0.2445  * tesla * 1000*mLHRSMomentum/ snakemagnumber / .1492 / m; //STD tune
//    KAPPA2 =  0.1939  * tesla * 1000*mLHRSMomentum/ snakemagnumber / .300  / m; //STD tune
//    KAPPA3 =  0.1794  * tesla * 1000*mLHRSMomentum/ snakemagnumber / .300  / m; //STD tune

    G4cout << "gradients: " << KAPPA1/tesla << " " << KAPPA2/tesla << " " << KAPPA3/tesla << "    q1rad=" << pQ1Radius/cm << "  length=" << pQ1Length/cm << "       q1en" << pQ1en/cm << "       q2en" << pQ2en/cm << G4endl;

  }
  //G4double dipoleField = -0.4192 * tesla; //using dipolert2 from SNAKE, thin vb, this one is reconciled with SNAKE
//  G4double dipoleField = -0.4205 * tesla; //matches nicely with DATA (not snake comparison) 
////    G4double dipoleField = -0.4192 * tesla;
    G4double dipoleField = DipField;
//  }
  G4cout << " dip field: " << dipoleField/tesla << G4endl;
  
  G4RotationMatrix* LROTATED = new G4RotationMatrix;
  G4RotationMatrix* RROTATED = new G4RotationMatrix;
  LROTATED->rotateY( mLHRSAngle );
  RROTATED->rotateY( mRHRSAngle );
  double pDipoleRCenterY=8.4  * m;
  double pDipoleRCenterZ=9.961 * m;//SNAKE

  G4ThreeVector     LORIGIND(pDipoleRCenterZ * sin( mLHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mLHRSAngle ));
  G4ThreeVector     RORIGIND(pDipoleRCenterZ * sin( mRHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mRHRSAngle ));

//fringe left

  G4ThreeVector     LORIGINQ1_Fr_En(pQ1en * sin(mLHRSAngle), 0., pQ1en * cos(mLHRSAngle));
  G4RotationMatrix* LROTATEQ1_Fr_En = new G4RotationMatrix;
  LROTATEQ1_Fr_En->rotateY( mLHRSAngle);

  G4ThreeVector     LORIGINQ1_Fr_Ex(pQ1ex * sin(mLHRSAngle), 0., pQ1ex * cos(mLHRSAngle));
  G4RotationMatrix* LROTATEQ1_Fr_Ex = new G4RotationMatrix;
  LROTATEQ1_Fr_Ex->rotateY( mLHRSAngle);

  G4ThreeVector     LORIGINQ2_Fr_En(pQ2en * sin(mLHRSAngle), 0., pQ2en * cos(mLHRSAngle));
  G4RotationMatrix* LROTATEQ2_Fr_En = new G4RotationMatrix;
  LROTATEQ2_Fr_En->rotateY( mLHRSAngle);

  G4ThreeVector     LORIGINQ2_Fr_Ex(pQ2ex * sin(mLHRSAngle), 0., pQ2ex * cos(mLHRSAngle));
  G4RotationMatrix* LROTATEQ2_Fr_Ex = new G4RotationMatrix;
  LROTATEQ2_Fr_Ex->rotateY( mLHRSAngle);


  double pLQ3_Fr_EnX = (17.0267042 * m  ) *  sin( mLHRSAngle );//SNAKE
  double pLQ3_Fr_EnZ = (17.0267042 * m  ) *  cos( mLHRSAngle );//SNAKE
  double pLQ3_Fr_EnY = ( 3.58637   * m  );//SNAKE	
  G4ThreeVector     LORIGINQ3_Fr_En(pLQ3_Fr_EnX, pLQ3_Fr_EnY, pLQ3_Fr_EnZ);//try this one
  G4RotationMatrix* LROTATEQ3_Fr_En = new G4RotationMatrix;
  LROTATEQ3_Fr_En->rotateX(-45.0 * deg);
  LROTATEQ3_Fr_En->rotateY( mLHRSAngle);

  double pLQ3_Fr_ExX = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  sin( mLHRSAngle );//SNAKE
  double pLQ3_Fr_ExZ = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  cos( mLHRSAngle );//SNAKE
  double pLQ3_Fr_ExY = ( 3.58637   * m  + pQ3Length / sqrt(2.) );//SNAKE	
  G4ThreeVector     LORIGINQ3_Fr_Ex(pLQ3_Fr_ExX, pLQ3_Fr_ExY, pLQ3_Fr_ExZ);//try this one
  G4RotationMatrix* LROTATEQ3_Fr_Ex = new G4RotationMatrix;
  LROTATEQ3_Fr_Ex->rotateX(-45.0 * deg);
  LROTATEQ3_Fr_Ex->rotateY( mLHRSAngle);


//fringe right

  G4ThreeVector     RORIGINQ1_Fr_En(pQ1en * sin(mRHRSAngle), 0., pQ1en * cos(mRHRSAngle));
  G4RotationMatrix* RROTATEQ1_Fr_En = new G4RotationMatrix;
  RROTATEQ1_Fr_En->rotateY( mRHRSAngle);

  G4ThreeVector     RORIGINQ1_Fr_Ex(pQ1ex * sin(mRHRSAngle), 0., pQ1ex * cos(mRHRSAngle));
  G4RotationMatrix* RROTATEQ1_Fr_Ex = new G4RotationMatrix;
  RROTATEQ1_Fr_Ex->rotateY( mRHRSAngle);

  G4ThreeVector     RORIGINQ2_Fr_En(pQ2en * sin(mRHRSAngle), 0., pQ2en * cos(mRHRSAngle));
  G4RotationMatrix* RROTATEQ2_Fr_En = new G4RotationMatrix;
  RROTATEQ2_Fr_En->rotateY( mRHRSAngle);

  G4ThreeVector     RORIGINQ2_Fr_Ex(pQ2ex * sin(mRHRSAngle), 0., pQ2ex * cos(mRHRSAngle));
  G4RotationMatrix* RROTATEQ2_Fr_Ex = new G4RotationMatrix;
  RROTATEQ2_Fr_Ex->rotateY( mRHRSAngle);


  double pRQ3_Fr_EnX = (17.0267042 * m  ) *  sin( mRHRSAngle );//SNAKE
  double pRQ3_Fr_EnZ = (17.0267042 * m  ) *  cos( mRHRSAngle );//SNAKE
  double pRQ3_Fr_EnY = ( 3.58637   * m  );//SNAKE	
  G4ThreeVector     RORIGINQ3_Fr_En(pRQ3_Fr_EnX, pRQ3_Fr_EnY, pRQ3_Fr_EnZ);//try this one
  G4RotationMatrix* RROTATEQ3_Fr_En = new G4RotationMatrix;
  RROTATEQ3_Fr_En->rotateX(-45.0 * deg);
  RROTATEQ3_Fr_En->rotateY( mRHRSAngle);

  double pRQ3_Fr_ExX = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  sin( mRHRSAngle );//SNAKE
  double pRQ3_Fr_ExZ = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  cos( mRHRSAngle );//SNAKE
  double pRQ3_Fr_ExY = ( 3.58637   * m  + pQ3Length / sqrt(2.) );//SNAKE	
  G4ThreeVector     RORIGINQ3_Fr_Ex(pRQ3_Fr_ExX, pRQ3_Fr_ExY, pRQ3_Fr_ExZ);//try this one
  G4RotationMatrix* RROTATEQ3_Fr_Ex = new G4RotationMatrix;
  RROTATEQ3_Fr_Ex->rotateX(-45.0 * deg);
  RROTATEQ3_Fr_Ex->rotateY( mRHRSAngle);




  G4ThreeVector     LORIGINQ1(pQ1Pos_Z * sin(mLHRSAngle), 0., pQ1Pos_Z * cos(mLHRSAngle));
  G4ThreeVector     RORIGINQ1(pQ1Pos_Z * sin(mRHRSAngle), 0., pQ1Pos_Z * cos(mRHRSAngle));
  G4RotationMatrix* LROTATEQ1 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ1 = new G4RotationMatrix;
  LROTATEQ1->rotateY( mLHRSAngle);
  RROTATEQ1->rotateY( mRHRSAngle);
  
  G4ThreeVector     LORIGINQ2(pQ2Pos_Z * sin(mLHRSAngle), 0., pQ2Pos_Z * cos(mLHRSAngle));
  G4ThreeVector     RORIGINQ2(pQ2Pos_Z * sin(mRHRSAngle), 0., pQ2Pos_Z * cos(mRHRSAngle));
  G4RotationMatrix* LROTATEQ2 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ2 = new G4RotationMatrix;
  
  LROTATEQ2->rotateY( mLHRSAngle);
  RROTATEQ2->rotateY( mRHRSAngle);
  
  //double pFPR       = 8.4  * m;//radius of curvature of dipole
  //double pFPA       = 9.96 * m;//distance from pivot to entrance of dipole//NIM
  //double pFPA       = pDipoleRCenterZ ;//distance from pivot to entrance of dipole//SNAKE
  //double pFPH       = pFPR * tan ( 22.5 * deg );//height
  //double pFPCenterX = ( pFPA + pFPH + ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2) ) *-sin( mRHRSAngle ); //not including half of Q3//NIM
  //double pFPCenterZ = ( pFPA + pFPH + ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2) ) * cos( mRHRSAngle ); //not including half of Q3//NIM
  //double pFPCenterY = ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2); //including half of Q3//NIM
/*
  double pLQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( mLHRSAngle );//SNAKE
  double pLQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( mLHRSAngle );//SNAKE
  double pLQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
  double pRQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( mRHRSAngle );//SNAKE
  double pRQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( mRHRSAngle );//SNAKE
  double pRQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
*/
  double pLQ3CenterX = (17.0267042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( mLHRSAngle );//SNAKE
  double pLQ3CenterZ = (17.0267042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( mLHRSAngle );//SNAKE
  double pLQ3CenterY = ( 3.58637   * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
  double pRQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( mRHRSAngle );//SNAKE
  double pRQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( mRHRSAngle );//SNAKE
  double pRQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
  
  //cout << pQ3CenterX << " " << pFPCenterX << endl;
  //cout << pQ3CenterY << " " << pFPCenterY << endl;
  //cout << pQ3CenterZ << " " << pFPCenterZ << endl;
  G4ThreeVector     LORIGINQ3(pLQ3CenterX, pLQ3CenterY, pLQ3CenterZ);//try this one
  G4ThreeVector     RORIGINQ3(pRQ3CenterX, pRQ3CenterY, pRQ3CenterZ);//try this one
  G4RotationMatrix* LROTATEQ3 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ3 = new G4RotationMatrix;
  
  LROTATEQ3->rotateX(-45.0 * deg);
  LROTATEQ3->rotateY( mLHRSAngle);
  RROTATEQ3->rotateX(-45.0 * deg);
  RROTATEQ3->rotateY( mRHRSAngle);
  
/*
  //fMagFieldFZBL1 = new BField_Quad(KAPPA1, LORIGINQ1, LROTATEQ1, pQ1Length, pQ1Radius, 1);
  fMagFieldFZBL1 = new BField_Quad_Snake(pQ1Length, KAPPA1, pQ1Radius, LORIGINQ1, LROTATEQ1, 1);
  fEquationFZBL1 = new G4Mag_UsualEqRhs(fMagFieldFZBL1);
  fStepperFZBL1  = new G4ClassicalRK4(fEquationFZBL1);
  fLocalFieldManagerFZBL1 = new G4FieldManager();
  fChordFinderFZBL1 = 0;
  UpdateFieldFZBL1();
//  fMagFieldFZBR1 = new BField_Quad(KAPPA1, RORIGINQ1, RROTATEQ1, pQ1Length, pQ1Radius, 1);
  fMagFieldFZBR1 = new BField_Quad_Snake(pQ1Length, -1.*KAPPA1, pQ1Radius, RORIGINQ1, RROTATEQ1, 1);
  fEquationFZBR1 = new G4Mag_UsualEqRhs(fMagFieldFZBR1);	
  fStepperFZBR1  = new G4ClassicalRK4(fEquationFZBR1);
  fLocalFieldManagerFZBR1 = new G4FieldManager();
  fChordFinderFZBR1 = 0;
  UpdateFieldFZBR1();
*/
/*
  //Local field  FZB2, Q2
//  fMagFieldFZBL2 = new BField_Quad(KAPPA2, LORIGINQ2, LROTATEQ2, pQ2Length, pQ2Radius, 2);
  fMagFieldFZBL2 = new BField_Quad_Snake(pQ2Length, KAPPA2, pQ2Radius, LORIGINQ2, LROTATEQ2, 2);
  fEquationFZBL2 = new G4Mag_UsualEqRhs(fMagFieldFZBL2);	
  fStepperFZBL2  = new G4ClassicalRK4(fEquationFZBL2);
  fLocalFieldManagerFZBL2 = new G4FieldManager();
  fChordFinderFZBL2 = 0;
  UpdateFieldFZBL2();
//  fMagFieldFZBR2 = new BField_Quad(KAPPA2, RORIGINQ2, RROTATEQ2, pQ2Length, pQ2Radius, 2);
  fMagFieldFZBR2 = new BField_Quad_Snake(pQ2Length, -1.*KAPPA2, pQ2Radius, RORIGINQ2, RROTATEQ2, 2);
  fEquationFZBR2 = new G4Mag_UsualEqRhs(fMagFieldFZBR2);	
  fStepperFZBR2  = new G4ClassicalRK4(fEquationFZBR2);
  fLocalFieldManagerFZBR2 = new G4FieldManager();
  fChordFinderFZBR2 = 0;
  UpdateFieldFZBR2();
*/

  fMagFieldFZBL3 = new BField_Dipole( dipoleField, LORIGIND, LROTATED );
  fEquationFZBL3 = new G4Mag_UsualEqRhs(fMagFieldFZBL3);	
  fLocalFieldManagerFZBL3 = new G4FieldManager();
  fChordFinderFZBL3 = 0;
  UpdateFieldFZBL3();
  fMagFieldFZBR3 = new BField_Dipole( -1.*dipoleField, RORIGIND, RROTATED );
  fEquationFZBR3 = new G4Mag_UsualEqRhs(fMagFieldFZBR3);	
  fLocalFieldManagerFZBR3 = new G4FieldManager();
  fChordFinderFZBR3 = 0;
  UpdateFieldFZBR3();

/*
  //Local field  FZB4, Q3
//  fMagFieldFZBL4 = new BField_Quad(KAPPA3, LORIGINQ3, LROTATEQ3, pQ3Length, pQ3Radius, 3);
  fMagFieldFZBL4 = new BField_Quad_Snake(pQ3Length, KAPPA3, pQ3Radius, LORIGINQ3, LROTATEQ3, 3);
  //fMagFieldFZBL4 = new BField_Quad(KAPPA3, ORIGINACTUAL, ROTATETEST);
  fEquationFZBL4 = new G4Mag_UsualEqRhs(fMagFieldFZBL4);	
  fStepperFZBL4  = new G4ClassicalRK4(fEquationFZBL4);
  fLocalFieldManagerFZBL4 = new G4FieldManager();
  fChordFinderFZBL4 = 0;
  UpdateFieldFZBL4();
  //Local field  FZB4, Q3
  fMagFieldFZBR4 = new BField_Quad_Snake(pQ3Length, -1.*KAPPA3, pQ3Radius, RORIGINQ3, RROTATEQ3, 3);
//  fMagFieldFZBR4 = new BField_Quad(KAPPA3, RORIGINQ3, RROTATEQ3, pQ3Length, pQ3Radius, 3);
  //fMagFieldFZBR4 = new BField_Quad(KAPPA3, ORIGINACTUAL, ROTATETEST);
  fEquationFZBR4 = new G4Mag_UsualEqRhs(fMagFieldFZBR4);	
  fStepperFZBR4  = new G4ClassicalRK4(fEquationFZBR4);
  fLocalFieldManagerFZBR4 = new G4FieldManager();
  fChordFinderFZBR4 = 0;
  UpdateFieldFZBR4();
*/
/*
  //Fringe field  FEnBQ2, Q2
  fMagFieldFEnBQ2 = new BField_Fringe_Q(1, KAPPA2, pQ2Radius, LORIGINQ2_Fr_En, LROTATEQ2_Fr_En, 2);
  fEquationFEnBQ2 = new G4Mag_UsualEqRhs(fMagFieldFEnBQ2);	
  fStepperFEnBQ2  = new G4ClassicalRK4(fEquationFEnBQ2);
  fLocalFieldManagerFEnBQ2 = new G4FieldManager();
  fChordFinderFEnBQ2 = 0;
  UpdateFieldFEnBQ2();
  //Fringe field  FExBQ2, Q2
  fMagFieldFExBQ2 = new BField_Fringe_Q(3, KAPPA2, pQ2Radius, LORIGINQ2_Fr_Ex, LROTATEQ2_Fr_Ex, 2);
  fEquationFExBQ2 = new G4Mag_UsualEqRhs(fMagFieldFExBQ2);
  fStepperFExBQ2  = new G4ClassicalRK4(fEquationFExBQ2);
  fLocalFieldManagerFExBQ2 = new G4FieldManager();
  fChordFinderFExBQ2 = 0;
  UpdateFieldFExBQ2();
*/




    
}

/////////////////////////////////////////////////////////////////////////////////
//

HRSEMFieldSetup::~HRSEMFieldSetup()
{
	if(fChordFinder) delete fChordFinder;
	if(fStepper)     delete fStepper;
	if(fEquation)    delete fEquation;
	if(fEMfield)     delete fEMfield;
	if(messenger)    delete messenger;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void HRSEMFieldSetup::UpdateField()
{
  //SetStepper();
  fStepper = new G4ClassicalRK4( fEquation, 8 );
//  fStepper = new G4HelixSimpleRunge( fEquation );
  G4cout<<"HRSEMFieldSetup:: The minimal rung-kutha step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

  fFieldManager->SetDetectorField(fEMfield);
//  fFieldManager->CreateChordFinder(fEMfield);
/*
   
   fIntgrDriver = new G4MagInt_Driver(fMinStep, fStepper, fStepper->GetNumberOfVariables() );
   fChordFinder = new G4ChordFinder(fIntgrDriver);
   fFieldManager->SetChordFinder( fChordFinder );
*/
  double fMinStep     = 0.001*1.*mm ; // minimal step of 1. mm
  if(fChordFinder) delete fChordFinder;
  fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldManager->SetChordFinder( fChordFinder );
/*
        G4FieldManager *globalFieldManager;
        G4TransportationManager *transportMgr= G4TransportationManager::GetTransportationManager();
        globalFieldManager = transportMgr->GetFieldManager();
        G4double minEps= 1.0e-5;  //   Minimum & value for smallest steps
        G4double maxEps= 1.0e-4;  //   Maximum & value for largest steps
        globalFieldManager->SetMinimumEpsilonStep( minEps );  
        globalFieldManager->SetMaximumEpsilonStep( maxEps );
        globalFieldManager->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
        G4cout << "EpsilonStep: set min= " << minEps << " max= " << maxEps << G4endl;
*/
}

/*
/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZBL1()
{
  G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL1->SetDetectorField(fMagFieldFZBL1);

	if(fChordFinderFZBL1) delete fChordFinderFZBL1;
	fIntgrDriverFZBL1 = new G4MagInt_Driver(fMinStep,fStepperFZBL1,fStepperFZBL1->GetNumberOfVariables());
	fChordFinderFZBL1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL1, fMinStep, fStepperFZBL1);
	fLocalFieldManagerFZBL1->SetChordFinder( fChordFinderFZBL1 );
	
}
void HRSEMFieldSetup::UpdateFieldFZBR1()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR1->SetDetectorField(fMagFieldFZBR1);

	if(fChordFinderFZBR1) delete fChordFinderFZBR1;
	fIntgrDriverFZBR1 = new G4MagInt_Driver(fMinStep,fStepperFZBR1,fStepperFZBR1->GetNumberOfVariables());
	fChordFinderFZBR1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR1, fMinStep, fStepperFZBR1);
	fLocalFieldManagerFZBR1->SetChordFinder( fChordFinderFZBR1 );
	
}
*/
/*

/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZBL2()
{
  G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL2->SetDetectorField(fMagFieldFZBL2);

	if(fChordFinderFZBL2) delete fChordFinderFZBL2;
	fIntgrDriverFZBL2 = new G4MagInt_Driver(fMinStep,fStepperFZBL2,fStepperFZBL2->GetNumberOfVariables());
	fChordFinderFZBL2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL2, fMinStep, fStepperFZBL2);
	fLocalFieldManagerFZBL2->SetChordFinder( fChordFinderFZBL2 );
}

void HRSEMFieldSetup::UpdateFieldFZBR2()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR2->SetDetectorField(fMagFieldFZBR2);

	if(fChordFinderFZBR2) delete fChordFinderFZBR2;
	fIntgrDriverFZBR2 = new G4MagInt_Driver(fMinStep,fStepperFZBR2,fStepperFZBR2->GetNumberOfVariables());
	fChordFinderFZBR2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR2, fMinStep, fStepperFZBR2);
	fLocalFieldManagerFZBR2->SetChordFinder( fChordFinderFZBR2 );
}
*/


void HRSEMFieldSetup::UpdateFieldFZBL3()
{
  fStepperFZBL3 = new G4ClassicalRK4( fEquationFZBL3, 8 );
  //G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;

  G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL3 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fLocalFieldManagerFZBL3->SetDetectorField(fMagFieldFZBL3);
  
  if(fChordFinderFZBL3) delete fChordFinderFZBL3;
  fChordFinderFZBL3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL3, fMinStep, fStepperFZBL3);
  fLocalFieldManagerFZBL3->SetChordFinder( fChordFinderFZBL3 );

  fLocalFieldManagerFZBL3->SetMinimumEpsilonStep( 1.0e-5 );
  fLocalFieldManagerFZBL3->SetMaximumEpsilonStep( 1.0e-2 );
  fLocalFieldManagerFZBL3->SetDeltaOneStep( 0.01*0.1 * mm );  // 0.5 micrometer
}
void HRSEMFieldSetup::UpdateFieldFZBR3()
{
  fStepperFZBR3 = new G4ClassicalRK4( fEquationFZBR3, 8 );
  //G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;

  G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR3 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fLocalFieldManagerFZBR3->SetDetectorField(fMagFieldFZBR3);
  
  if(fChordFinderFZBR3) delete fChordFinderFZBR3;
  fChordFinderFZBR3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR3, fMinStep, fStepperFZBR3);
  fLocalFieldManagerFZBR3->SetChordFinder( fChordFinderFZBR3 );

//  fLocalFieldManagerFZBR3->SetMinimumEpsilonStep( 1.0e-5 );
//  fLocalFieldManagerFZBR3->SetMaximumEpsilonStep( 1.0e-4 );
//  fLocalFieldManagerFZBR3->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
  fLocalFieldManagerFZBR3->SetMinimumEpsilonStep( 1.0e-5 );
  fLocalFieldManagerFZBR3->SetMaximumEpsilonStep( 1.0e-2 );
  fLocalFieldManagerFZBR3->SetDeltaOneStep( 0.1 * mm );  // 0.5 micrometer
}

/*
/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZBL4()
{
  G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL4 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL4->SetDetectorField(fMagFieldFZBL4);

	if(fChordFinderFZBL4) delete fChordFinderFZBL4;
	fIntgrDriverFZBL4 = new G4MagInt_Driver(fMinStep,fStepperFZBL4,fStepperFZBL4->GetNumberOfVariables());
	fChordFinderFZBL4 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL4, fMinStep, fStepperFZBL4);
	fLocalFieldManagerFZBL4->SetChordFinder( fChordFinderFZBL4 );
}
void HRSEMFieldSetup::UpdateFieldFZBR4()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR4 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR4->SetDetectorField(fMagFieldFZBR4);

	if(fChordFinderFZBR4) delete fChordFinderFZBR4;
	fIntgrDriverFZBR4 = new G4MagInt_Driver(fMinStep,fStepperFZBR4,fStepperFZBR4->GetNumberOfVariables());
	fChordFinderFZBR4 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR4, fMinStep, fStepperFZBR4);
	fLocalFieldManagerFZBR4->SetChordFinder( fChordFinderFZBR4 );
}
*/
/*
void HRSEMFieldSetup::UpdateFieldFEnBQ2()
{
  G4cout<<"HRSEMFieldSetup:: The minimal step for FEnBQ2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFEnBQ2->SetDetectorField(fMagFieldFEnBQ2);

	if(fChordFinderFEnBQ2) delete fChordFinderFEnBQ2;
	fIntgrDriverFEnBQ2 = new G4MagInt_Driver(fMinStep,fStepperFEnBQ2,fStepperFEnBQ2->GetNumberOfVariables());
	fChordFinderFEnBQ2 = new G4ChordFinder((G4MagneticField*) fMagFieldFEnBQ2, fMinStep, fStepperFEnBQ2);
	fLocalFieldManagerFEnBQ2->SetChordFinder( fChordFinderFEnBQ2 );
}
void HRSEMFieldSetup::UpdateFieldFExBQ2()
{
  G4cout<<"HRSEMFieldSetup:: The minimal step for FExBQ2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFExBQ2->SetDetectorField(fMagFieldFExBQ2);

	if(fChordFinderFExBQ2) delete fChordFinderFExBQ2;
	fIntgrDriverFExBQ2 = new G4MagInt_Driver(fMinStep,fStepperFExBQ2,fStepperFExBQ2->GetNumberOfVariables());
	fChordFinderFExBQ2 = new G4ChordFinder((G4MagneticField*) fMagFieldFExBQ2, fMinStep, fStepperFExBQ2);
	fLocalFieldManagerFExBQ2->SetChordFinder( fChordFinderFExBQ2 );
}
*/
/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//


void HRSEMFieldSetup::SetStepper()
{
	G4int nvar = 8;

	if(fStepper) delete fStepper;

	switch ( fStepperType )
	{
	case 0:
		fStepper = new G4ExplicitEuler( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4ExplicitEuler is calledS"<<G4endl;
		break;
	case 1:
		fStepper = new G4ImplicitEuler( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4ImplicitEuler is called"<<G4endl;
		break;
	case 2:
		fStepper = new G4SimpleRunge( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4SimpleRunge is called"<<G4endl;
		break;
	case 3:
		fStepper = new G4SimpleHeum( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4SimpleHeum is called"<<G4endl;
		break;
	case 4:
		fStepper = new G4ClassicalRK4( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;
		break;
	case 5:
		fStepper = new G4CashKarpRKF45( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4CashKarpRKF45 is called"<<G4endl;
		break;

	//The following not working for electric field
	case 6:
		fStepper = 0 ; // new G4RKG3_Stepper( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
		break;
	case 7:
		fStepper = 0 ; // new G4HelixExplicitEuler( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 8:
		fStepper = 0 ; //  new G4HelixImplicitEuler( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 9:
		fStepper = 0 ; //  new G4HelixSimpleRunge( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
		break;
	default: fStepper = 0;
	}
}


///////////////////////////////////////////////////////////////////////////////
