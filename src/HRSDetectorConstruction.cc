// ********************************************************************
//
// $Id: HRSDetectorConstruction.cc,v 1.02, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include <stdio.h>
#include <math.h>
#include <fstream>

#include "HRSMaterial.hh"
#include "HRSDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4UniformMagField.hh"
#include "G4SubtractionSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "BField_Dipole.hh"
#include "BField_Dipole_Fringe.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polycone.hh"
#include "G4AssemblyVolume.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "HRSEMFieldSetup.hh"
#include "HRSStdSD.hh"
#include "UsageManager.hh"

#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
//To verify some geometry, I use this flag to place some unit
//just for debugging
//#define G4DEBUG_GEOMETRY 0

//The G4RotationMatrix rotates the whole coordinate system, 
//Looking from top, it always rotate clockwise  

extern UsageManager* gConfig;	


/////////////////////////////////////////////////////////////////////
HRSDetectorConstruction::HRSDetectorConstruction()
{
	this->GetConfig();
	mEMFieldSetup = 0;

	//construct the material manager, this line should behind GetConfig()
	//since it need to access the buffer of gConfig
	mMaterialManager=HRSMaterial::GetHRSMaterialManager();
        double LarmStepLimit=2.000 * mm;
        LarmStepLimits = new G4UserLimits(LarmStepLimit);
}

HRSDetectorConstruction::~HRSDetectorConstruction()
{
	if(mEMFieldSetup) delete mEMFieldSetup;
}

/////////////////////////////////////////////////////////////////////
void HRSDetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aPhysVol,G4int depth,
												  ostream &fout)
{
	//define the aLogVol to speed up
	G4LogicalVolume* aLogVol=aPhysVol->GetLogicalVolume();
	fout << "#";
	for(int isp=0;isp<depth;isp++) { fout << "--"; }
	fout << aPhysVol->GetName() << "[" << aPhysVol->GetCopyNo() << "] "
		<< aLogVol->GetName() << " "
		<< aLogVol->GetNoDaughters() << " "
		<< aLogVol->GetMaterial()->GetName();
	if(aLogVol->GetSensitiveDetector())
	{
		fout << " " << aLogVol->GetSensitiveDetector()->GetFullPathName();
	}
	fout << endl;
	for(int i=0;i<aLogVol->GetNoDaughters();i++)
	{ 
		DumpGeometricalTree(aLogVol->GetDaughter(i),depth+1,fout); 
	}
}

/////////////////////////////////////////////////////////////////////
void HRSDetectorConstruction::GenerateGeometryMacro(G4VPhysicalVolume* aPhysVol,G4int depth)
{
	ofstream fout;
	if(depth==0) 
	{
		fout.open("geometry.mac",ios_base::out);
		//DumpGeometricalTree(aPhysVol,depth,fout);
	}
	else fout.open("geometry.mac",ios_base::app);


	G4LogicalVolume* aLogVol = aPhysVol->GetLogicalVolume();
	string logName = aLogVol->GetName().c_str();

	//some phys volumn share the same logical volumn
	//each log vol needs to print cmds for just one time
	if(mIsLogVolPrinted.size() > 0)
	{
		if(mIsLogVolPrinted.find(logName) != mIsLogVolPrinted.end()) return;
	}
	//add this log vol into the map
	mIsLogVolPrinted[logName] = 1;

	int visibility = (aLogVol->GetVisAttributes())?1:0;
	if(visibility) visibility=aLogVol->GetVisAttributes()->IsVisible()?1:0;
	G4Colour aColor;
	if(visibility) aColor = aLogVol->GetVisAttributes()->GetColor();
	int nDaughter = aLogVol->GetNoDaughters();

	fout << "#";
	for(int isp=0;isp<depth;isp++) { fout << "--"; }
	fout << aPhysVol->GetName() << "[" << aPhysVol->GetCopyNo() << "] "
		<< " " << logName 
		<< " " << aLogVol->GetMaterial()->GetName()
		<< " Ndaughter=" << nDaughter;
	if(visibility) fout<< " G4Color"<< aColor;
	if(aLogVol->GetSensitiveDetector())
		fout << " SD=" << aLogVol->GetSensitiveDetector()->GetFullPathName();
	fout << endl;

	char space[100];
	sprintf(space," ");
	for(int isp=0;isp<depth;isp++) { sprintf(space,"%s  ",space); }

	//now print the commands
	fout<<space<<"/vis/geometry/set/visibility "<<logName<<" 0 "
		<<(visibility?"true":"false")<<" \n";
	fout<<space<<"/vis/geometry/set/forceSolid "<<logName<<" 0 "<<" true \n";
	fout<<space<<"#/vis/geometry/set/forceWireframe "<<logName<<" 0 "<<" false \n";
	if(nDaughter)
		fout<<space<<"#/vis/geometry/set/daughtersInvisible "<<logName<<" 0 "<<" false \n";
	if(visibility)
		fout<<space<<"#/vis/geometry/set/colour "<<logName<<" 0 "
		<<aColor.GetRed()<<" "
		<<aColor.GetGreen()<<" "
		<<aColor.GetBlue()<<" "
		<<aColor.GetAlpha()<<"\n";

	fout.close();

	for(int i=0;i<nDaughter;i++)
	{ 
		GenerateGeometryMacro(aLogVol->GetDaughter(i),depth+1); 
	}

}
/////////////////////////////////////////////////////////////////////
void HRSDetectorConstruction::GetConfig()
{
	gConfig->GetParameter("HallX",mHallX);
	mHallX*=mm;
	gConfig->GetParameter("HallY",mHallY);
	mHallY*=mm;
	gConfig->GetParameter("HallZ",mHallZ);
	mHallZ*=mm;

	gConfig->GetParameter("FieldX",mFieldX);
	mFieldX*=mm;
	gConfig->GetParameter("FieldY",mFieldY);
	mFieldY*=mm;
	gConfig->GetParameter("FieldZ",mFieldZ);
	mFieldZ*=mm;

	mBMaterialType=2;  
	gConfig->GetParameter("BMaterialType",mBMaterialType);

	mSetupAirAxis=0;
	gConfig->GetParameter("SetupAirAxis",mSetupAirAxis);

	gConfig->GetParameter("SetupVirtualDetector",mSetupVirtualDetector);
	gConfig->GetParameter("VirtualDetectorWidth",mVirtualDetectorWidth);
	mVirtualDetectorWidth*=mm;
	gConfig->GetParameter("VirtualDetectorHeight",mVirtualDetectorHeight);
	mVirtualDetectorHeight*=mm;
	gConfig->GetParameter("VirtualDetectorThick",mVirtualDetectorThick);
	mVirtualDetectorThick*=mm;
	gConfig->GetParameter("VDRotYAngle",mVDRotYAngle);
	mVDRotYAngle*=deg;
	gConfig->GetParameter("VDRotXAngle",mVDRotXAngle);
	mVDRotXAngle*=deg;
	gConfig->GetParameter("Pivot2VDFace",mPivot2VDFace);
	mPivot2VDFace*=mm;

	mVDPhysVolName=gConfig->GetParameter("VDPhysVolName");

	/////////////////////////////////////////////////////////////////////////
	
	gConfig->GetParameter("SetupRadiator",mSetupRadiator);
	gConfig->GetParameter("SetupRadiatorVD",mSetupRadiatorVD);

	gConfig->GetParameter("RadiatorType",mRadiatorType);
	gConfig->GetParameter("RaditorThickness",mRaditorThickness);
	mRaditorThickness*=mm;
	gConfig->GetParameter("Radiator2Pivot",mRadiator2Pivot);
	mRadiator2Pivot*=mm;

	/////////////////////////////////////////////////////////////////////////

	mSetupCREXGeometry=0;
	gConfig->GetParameter("SetupCREXGeometry",mSetupCREXGeometry);

	mSetupAPEXGeometry=0;
	gConfig->GetParameter("SetupAPEXGeometry",mSetupAPEXGeometry);

	mSetupRTPCGeometry=0;
	gConfig->GetParameter("SetupRTPCGeometry",mSetupRTPCGeometry);

	//G4cout<<"\n****Load detector config parameters done!***"<<G4endl;

	///////////////////////////////////////////////////////////////////////////
	//global variables

	gConfig->GetParameter("SetupLHRS",mSetupLHRS);
	gConfig->GetParameter("SetupRHRS",mSetupRHRS);
	gConfig->GetParameter("SetupLSieveSlit",mSetupLSieveSlit);
	gConfig->GetParameter("SetupRSieveSlit",mSetupRSieveSlit);

        gConfig->GetParameter("sieve_thickness",sieve_thickness);
        gConfig->GetParameter("Rad_tail_foil_postn", Rad_tail_foil_postn);
        gConfig->GetParameter("SeptumOn",Septum_On);
	gConfig->GetParameter("Q1Sos",nQ1Sos);
	gConfig->GetParameter("LFocalPlaneAngle",fp_angle);
	fp_angle*=deg;
	gConfig->GetParameter("LHRSAngle",mLHRSAngle);
	mLHRSAngle*=deg;
	gConfig->GetParameter("LSeptumAngle",mLSeptumAngle);
	mLSeptumAngle*=deg;
	gConfig->GetParameter("RHRSAngle",mRHRSAngle);
	mRHRSAngle*=deg;
	gConfig->GetParameter("RSeptumAngle",mRSeptumAngle);
	mRSeptumAngle*=deg;

	gConfig->GetParameter("ScatChamberXOffset",mScatChamberXOffset);
	mScatChamberXOffset*=mm;
	gConfig->GetParameter("ScatChamberYOffset",mScatChamberYOffset);
	mScatChamberYOffset*=mm;
	gConfig->GetParameter("ScatChamberZOffset",mScatChamberZOffset);
	mScatChamberZOffset*=mm;

	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	gConfig->GetParameter("TargetXOffset",mTargetXOffset);
	mTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",mTargetYOffset);
	mTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",mTargetZOffset);
	mTargetZOffset*=mm;


	return ;
}


/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* HRSDetectorConstruction::Construct()
{
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	//G4VSensitiveDetector* virtualDetectorSD=new HRSStdSD(SDname="virtualDetector");
	G4VSensitiveDetector* virtualBoundarySD=new HRSStdSD(SDname="virtualBoundary");

	// sensitive detectors   
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	/////////////////////////////////////////////////////////////////////////////

	//set "user limits" for drawing smooth curve, only supported by Jixie's model
	//G4UserLimits* uHallStepLimits = new G4UserLimits(1000.*mm);//Jixie default
	G4UserLimits* uHallStepLimits = new G4UserLimits(.001 * mm);//Nickie default
	//double pBStepLimit=100;//Jixie default
	double pBStepLimit=0.001 * mm;//Nickie's edit
	gConfig->GetArgument("BStepLimit",pBStepLimit);
	pBStepLimit*=mm;
	pBStepLimit=2.037 * mm;
	G4cout << "B limits: " << pBStepLimit << G4endl;
	G4UserLimits* uBStepLimits = new G4UserLimits(pBStepLimit);

	///////////////////////////////////////////////////////////////////////////////
	// Magnetic field ----------------------------------------------------------
	mEMFieldSetup = new HRSEMFieldSetup();  //setup the field, 
	mEMFieldSetup->UpdateField();
	
	
//	G4FieldManager* dipoleFieldManager = mEMFieldSetup->GetFieldManagerFZBL3();
//	GetFieldManagerFZBL3
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//Shrink the hall size if bigbite is not built, this will mke it run fast
	//note that this is optimized for g2p only
	//These 2 from Detector_G2P.ini
	//int pSetupChicane=1,pSetupThirdArm=1;
	//if (!mSetupG2PGeometry) 
	//{
	//pSetupChicane=pSetupThirdArm=0;
	//}

	///////////////////////////////////////////////////////////////////////////

	// start to construct the geometries -----------------------------------------

	/////////////////////////
	// experimental hall (hall volume)
	/////////////////////////
	G4VSolid* hallSolid = new G4Box("hallBox",mHallX/2.0,mHallY/2.0,mHallZ/2.0);
	G4LogicalVolume* hallLogical = new G4LogicalVolume(hallSolid,
		mMaterialManager->vacuum,"hallLogical",0,0,uHallStepLimits);
	hallLogical->SetVisAttributes(HallVisAtt);
	mHallPhysVol = new G4PVPlacement(0,G4ThreeVector(),
		hallLogical,"hallPhys",0,0,0);
	hallLogical->SetUserLimits(uBStepLimits);
	cout<<"hall step is "<<pBStepLimit<<endl;


	/////////////////////////
	// air axes
	/////////////////////////
	//show the axis to help to understand the geometry in in the visualization
	//This part is also used to verify that the G4PVPlacement() will do the rotation at the origin 
	//then put the LogicalVolume at the given position 
	if(mSetupAirAxis)
	{
		G4RotationMatrix* pRotX270deg = new G4RotationMatrix();
		pRotX270deg->rotateX(270.*deg);
		G4RotationMatrix* pRotY270deg = new G4RotationMatrix();
		pRotY270deg->rotateY(270.*deg);

		G4VSolid* xAxisSolid = new G4Tubs("xAxisTubs",0.,5.0*mm,500.0*mm,0.,360.*deg);
		G4LogicalVolume* xAxisLogical = new G4LogicalVolume(xAxisSolid,
			mMaterialManager->vacuum,"xAxisLogical",0,0,LarmStepLimits);
		new G4PVPlacement(pRotY270deg,G4ThreeVector(650*mm,0,0),
			xAxisLogical,"xAxisPhys",hallLogical,0,0);
		xAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

		G4VSolid* yAxisSolid = new G4Tubs("yAxisTubs",0.,5.0*mm,500.0*mm,0.,360.*deg);
		G4LogicalVolume* yAxisLogical = new G4LogicalVolume(yAxisSolid,
			mMaterialManager->vacuum,"yAxisLogical",0,0,LarmStepLimits);
		new G4PVPlacement(pRotX270deg,G4ThreeVector(0,650*mm,0),
			yAxisLogical,"yAxisPhys",hallLogical,0,0);
		yAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

		G4VSolid* zAxisSolid = new G4Tubs("zAxisTubs",0.,5.0*mm,500.0*mm,0.,360.*deg);
		G4LogicalVolume* zAxisLogical = new G4LogicalVolume(zAxisSolid,
			mMaterialManager->vacuum,"zAxisLogical",0,0,LarmStepLimits);
		new G4PVPlacement(0,G4ThreeVector(0,0,650*mm),
			zAxisLogical,"zAxisPhys",hallLogical,0,0);
		zAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
	}



	/////////////////////////
	// magnetic field region
	/////////////////////////
	G4Material* pBMaterial=0;
	if(mBMaterialType==0) pBMaterial=mMaterialManager->vacuum;
	else if(mBMaterialType==2) pBMaterial=mMaterialManager->heliumGas;
	else pBMaterial=mMaterialManager->air;

	G4VSolid* magneticSolid = new G4Box("magneticBox",mFieldX/2.0,mFieldY/2.0,mFieldZ/2.0);
	pBMaterial=mMaterialManager->vacuum;
	G4LogicalVolume* magneticLogical = new G4LogicalVolume(magneticSolid,
		pBMaterial,"magneticLogical",0,0,LarmStepLimits);
	G4VisAttributes* MagneticVisAttTest = new G4VisAttributes(G4Colour(1., 0., 1.));
	magneticLogical->SetUserLimits(uBStepLimits);
	magneticLogical->SetVisAttributes(MagneticVisAtt);
	//magneticLogical->SetVisAttributes(MagneticVisAttTest);
	//magneticLogical->SetVisAttributes(G4VisAttributes::Invisible);

	//the field region is the mother volumn of a lot of dauthers, it should not rotate
	new G4PVPlacement(0,G4ThreeVector(),magneticLogical,"magneticPhys",hallLogical,0,0);


	/////////////////////////
	//the radiator
	/////////////////////////
//	if(mSetupRadiator) ConstructRadiator(magneticLogical);




	/////////////////////////
	//CREX geometry
	/////////////////////////	
	//scattering chamber, target, 
	//need to add sieve, septum, HRSVB, 
/*
	if(mSetupCREXGeometry) 
	{
		CREXDetectorConstruction* theCREX = new CREXDetectorConstruction(magneticLogical); 
		theCREX->Construct();
		//update the parameters
		this->GetConfig(); 
	}
*/
	/////////////////////////
	// Sieve slit, septum window and HRS
	/////////////////////////
	if(mSetupLHRS || mSetupRHRS)  this->ConstructHRS(magneticLogical);



	/////////////////////////
	// Cylinder Virtual Boundary
	/////////////////////////
	//-VB or -VirtualBoundary <SetupVirtualBoundary(0)> [VBRin(537)] [VBRout(540)] [VBHeight(2000)] [VBRotAxis(1)]
	//   [VBRotAngle(90)] [VBPosX(0)] [VBPosY(0)] [VBPosZ(0)]:
	//description: A switch to tell the program to place the cylinder virtual boundary. This 
	//solid is a cylinder made of air. One can specify the solid dimention, rotation and center position. 
	//Note that VBRotAxis can be 0(NO rotation) or 1 (X), 2(Y) and 3(Z) axis. The VBRotAngle is in deg. Only
	//one rotation is supported. The center position VBPosX,VBPosY and VBPosZ are all in mm and relative to 
	//the pivot, where is the center of the target in most settings.

	//one can use it for fast reason ...... 

	int pSetupVirtualBoundary=0;
	gConfig->GetArgument("SetupVirtualBoundary",pSetupVirtualBoundary); 
	if(pSetupVirtualBoundary)
	{
		double pVBRin=537.0,pVBRout=538.0,pVBHeight=2000.0;
		gConfig->GetArgument("VBRin",pVBRin); 
		gConfig->GetArgument("VBRout",pVBRout); 
		gConfig->GetArgument("VBHeight",pVBHeight); 
		pVBRin*=mm,pVBRout*=mm,pVBHeight*=mm;

		G4VSolid* virtualBoundary0Solid = new G4Tubs("virtualBoundary0Tubs",
			pVBRin,(pVBRin+pVBRout)/2,pVBHeight/2.0,-240*deg,300*deg);
		G4VSolid* virtualBoundary1Solid = new G4Tubs("virtualBoundary1Tubs",
			(pVBRin+pVBRout)/2,pVBRout,pVBHeight/2.0,-240*deg,300*deg);
		G4LogicalVolume* virtualBoundary0Logical = new G4LogicalVolume(virtualBoundary0Solid,
			mMaterialManager->air,"virtualBoundary0Logical",0,0,LarmStepLimits);
		G4LogicalVolume* virtualBoundary1Logical = new G4LogicalVolume(virtualBoundary1Solid,
			mMaterialManager->air,"virtualBoundary1Logical",0,0,LarmStepLimits);

		virtualBoundary0Logical->SetVisAttributes(LightYellowVisAtt);  
		virtualBoundary1Logical->SetVisAttributes(LightYellowVisAtt); 
		SDman->AddNewDetector(virtualBoundarySD);
		virtualBoundary0Logical->SetSensitiveDetector(virtualBoundarySD);

		double pVBPosX=0,pVBPosY=0,pVBPosZ=0;
		gConfig->GetArgument("VBPosX",pVBPosX); 
		gConfig->GetArgument("VBPosY",pVBPosY); 
		gConfig->GetArgument("VBPosZ",pVBPosZ); 
		pVBPosX=pVBPosX*mm+mPivotXOffset;
		pVBPosY=pVBPosX*mm+mPivotYOffset;
		pVBPosZ=pVBPosZ*mm+mPivotZOffset;

		G4RotationMatrix *pRotVB=new G4RotationMatrix();
		int pVBRotAxis=0;
		double pVBRotAngle=0; 
		gConfig->GetArgument("VBRotAxis",pVBRotAxis); 
		gConfig->GetArgument("VBRotAngle",pVBRotAngle); pVBRotAngle*=deg;
		if(pVBRotAxis && fabs(pVBRotAngle)>1.0E-05)
		{
			if(pVBRotAxis==1) pRotVB->rotateX(pVBRotAngle); 
			else if(pVBRotAxis==2) pRotVB->rotateY(pVBRotAngle); 
			else if(pVBRotAxis==3) pRotVB->rotateZ(pVBRotAngle); 
		}

		new G4PVPlacement(pRotVB,G4ThreeVector(pVBPosX,pVBPosY,pVBPosZ),			
			virtualBoundary0Logical,"virtualDetectorPhys",magneticLogical,0,0);
		new G4PVPlacement(pRotVB,G4ThreeVector(pVBPosX,pVBPosY,pVBPosZ),			
			virtualBoundary1Logical,"virtualBoundaryPhys",magneticLogical,0,0);
	}



        if (0)
        {
//          G4double target_width=0.5*0.0185*mm;
          G4double target_width=0.35*mm;
          G4double target_x=2.5*mm;
          G4double target_y=80.*mm;
          G4VSolid* targ1=new G4Box("Targ foil1",target_x/2.,target_y/2.,target_width/2.);
          G4LogicalVolume* targ1Logical = new G4LogicalVolume(targ1,mMaterialManager->tungsten,"targ1Logical",0,0,0);//LarmStepLimits);
          targ1Logical->SetVisAttributes((G4Colour(0.998, 0.008, 0.000)));

          G4LogicalVolume* targ2Logical = new G4LogicalVolume(targ1,mMaterialManager->tungsten,"targ2Logical",0,0,LarmStepLimits);
          targ2Logical->SetVisAttributes((G4Colour(0.00998, 0.99008, 0.000)));

          G4VSolid* targtrans = new G4Box("Targ trans",target_x/2.,target_y/2.,0.001*target_width/2.);
          G4LogicalVolume* targtransLogical = new G4LogicalVolume(targtrans,mMaterialManager->vacuum,"targtransLogical",0,0,LarmStepLimits);
          targtransLogical->SetVisAttributes((G4Colour(0.00998, 0.99008, 0.000)));

//          new G4PVPlacement(0,G4ThreeVector(0,0,0),targ1Logical,"targ1",hallLogical,0,0,0);
//          new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-90.0*cm),targ1Logical,"targ1",hallLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm-0.5*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
//          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);
        }

        // targets 1-10 for e_loss study
        if (0)
        {
//          G4double target_width=0.5*1.185*mm;
          G4double target_width=0.5*0.0185*mm;
//          G4double target_width=0.5*0.0098*mm;
//          G4double target_width=0.35*mm;
          G4double target_x=2.5*mm;
          G4double target_y=80.*mm;
          G4VSolid* targ1=new G4Box("Targ foil1",target_x/2.,target_y/2.,target_width/2.);
          G4LogicalVolume* targ1Logical = new G4LogicalVolume(targ1,mMaterialManager->tungsten,"targ1Logical",0,0,0);//LarmStepLimits);
          targ1Logical->SetVisAttributes((G4Colour(0.998, 0.008, 0.000)));

          G4LogicalVolume* targ2Logical = new G4LogicalVolume(targ1,mMaterialManager->tungsten,"targ2Logical",0,0,LarmStepLimits);
          targ2Logical->SetVisAttributes((G4Colour(0.00998, 0.99008, 0.000)));

          G4VSolid* targtrans = new G4Box("Targ trans",target_x/2.,target_y/2.,0.001*target_width/2.);
          G4LogicalVolume* targtransLogical = new G4LogicalVolume(targtrans,mMaterialManager->vacuum,"targtransLogical",0,0,LarmStepLimits);
          targtransLogical->SetVisAttributes((G4Colour(0.00998, 0.99008, 0.000)));

//          new G4PVPlacement(0,G4ThreeVector(0,0,0),targ1Logical,"targ1",hallLogical,0,0,0);
//          new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-90.0*cm),targ1Logical,"targ1",hallLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-150.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-150.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-150.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-145.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-145.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-145.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-140.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-140.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-140.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-135.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-135.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-135.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-130.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-130.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-130.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-125.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-125.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-125.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-120.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-120.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-120.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-115.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-115.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-115.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

          new G4PVPlacement(0,G4ThreeVector(0,0,-110.*cm-0.005*mm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-110.*cm+0.005*mm),targ2Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-110.*cm),     targtransLogical,"targtrans",magneticLogical,0,0,0);

//          new G4PVPlacement(0,G4ThreeVector(0,0,-149.*cm-100.*0.005*mm),targ1Logical,"targ",hallLogical,0,0,0);
//          new G4PVPlacement(0,G4ThreeVector(0,0,-149.*cm),targtransLogical,"targ",hallLogical,0,0,0);
//          new G4PVPlacement(0,G4ThreeVector(0,0,-149.*cm+100.*0.005*mm),targ2Logical,"targ",hallLogical,0,0,0);
        }


        if (0) //wab
        {
          G4double target_width=0.0185*mm; //4.4 GeV
//          G4double target_width=0.0098*mm; //2.2 GeV

          G4double target_x=2.5*mm;
          G4double target_y=80.*mm;
          G4VSolid* targ1=new G4Box("Targ foil1",target_x/2.,target_y/2.,target_width/2.);
          G4LogicalVolume* targ1Logical = new G4LogicalVolume(targ1,mMaterialManager->tungsten,"targ1Logical",0,0,0);//LarmStepLimits);
          targ1Logical->SetVisAttributes((G4Colour(0.998, 0.008, 0.000)));

          new G4PVPlacement(0,G4ThreeVector(0,0,-155.*cm),targ1Logical,"targ1",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-150.*cm),targ1Logical,"targ2",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-145.*cm),targ1Logical,"targ3",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-140.*cm),targ1Logical,"targ4",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-135.*cm),targ1Logical,"targ5",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-130.*cm),targ1Logical,"targ6",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-125.*cm),targ1Logical,"targ7",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-120.*cm),targ1Logical,"targ8",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-115.*cm),targ1Logical,"targ9",magneticLogical,0,0,0);
          new G4PVPlacement(0,G4ThreeVector(0,0,-110.*cm),targ1Logical,"targ10",magneticLogical,0,0,0);
        }



        cout<<" septum@ ka te che? "<<Septum_On<<endl;
        //septum vac box
        if (Septum_On)
        {
            double inch           = 2.54*cm;
            double y_en = 2.44*inch;
            double y_ex = 4.7 *inch;
            double zlength        = 173.939*cm;
            double z_sep_cen      = zlength*tan(5.*deg)/tan(12.5*deg);
            double z_real_tar     = z_sep_cen-zlength;
            double z_sept_en_min1 = z_real_tar + (17.31*inch+22.50*inch);
            double z_sept_en_max1 = z_sept_en_min1 - 3.15*inch*tan(5.*deg);
            double x_cen_sep_en   = 3.966*inch;
            double x_width_sep_en = 3.15*inch;
            double xmin_sep_en1   = x_cen_sep_en -x_width_sep_en*0.5*cos(5.*deg);
            double xmax_sep_en1   = x_cen_sep_en +x_width_sep_en*0.5*cos(5.*deg);
            double ang_en_min_1   = 5.*deg;
            double ang_en_max_1   = 5.*deg;
            double ang_en_min_2   = 6.6*deg;
            double ang_en_max_2   = 10.8*deg;
            double length_max_1   = 50.19* cm; // 19.76*inch
            double length_min_1   = 52.1 * cm; // 19.76*inch
            double length_min_2   = 59.82* cm; // 23.55*inch
            double length_max_2   = 60.3 * cm; // 23.74*inch

            double xmin_sep_ex1   = xmin_sep_en1 + length_min_1 * sin(ang_en_min_1);
            double xmax_sep_ex1   = xmax_sep_en1 + length_max_1 * sin(ang_en_max_1);
            double z_sept_ex_min1 = z_sept_en_min1 + length_min_1 * cos(ang_en_min_1);
            double z_sept_ex_max1 = z_sept_en_max1 + length_max_1 * cos(ang_en_max_1);

            double xmin_sep_en2   = xmin_sep_ex1;
            double xmax_sep_en2   = xmax_sep_ex1;
            double z_sept_en_min2 = z_sept_ex_min1;
            double z_sept_en_max2 = z_sept_ex_max1;

            double xmin_sep_ex2   = xmin_sep_en2 + length_min_2 * sin(ang_en_min_2);
            double xmax_sep_ex2   = xmax_sep_en2 + length_max_2 * sin(ang_en_max_2);
            double z_sept_ex_min2 = z_sept_en_min2 + length_min_2 * cos(ang_en_min_2);
            double z_sept_ex_max2 = z_sept_en_max2 + length_max_2 * cos(ang_en_max_2);

            double ymin_sep_ex1 = y_en + (z_sept_ex_min1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
            double ymax_sep_ex1 = y_en + (z_sept_ex_max1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
            double ymin_sep_en1 = y_en + (z_sept_en_min1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
            double ymax_sep_en1 = y_en + (z_sept_en_max1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);


            double ymin_sep_ex2 = y_en + (z_sept_ex_min2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
            double ymax_sep_ex2 = y_en + (z_sept_ex_max2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
            double ymin_sep_en2 = y_en + (z_sept_en_min2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
            double ymax_sep_en2 = y_en + (z_sept_en_max2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);


          G4double  pDz_1=0.31*inch/2.;
          G4double  pDz_2=0.25*inch/2.;
//          G4double  pDz_1=2.31*inch/2.;
//          G4double  pDz_2=2.25*inch/2.;
          vector<G4TwoVector> vertices_en1_min;
          vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
          vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
          vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
          vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );
          vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
          vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
          vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
          vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );

          G4VSolid* testTrap = new G4GenericTrap("vac Box 1", pDz_1, vertices_en1_min);
          G4LogicalVolume* trrap = new G4LogicalVolume(testTrap, mMaterialManager->tungsten,"vac Box 1",0,0,LarmStepLimits);
          trrap->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *pRotY90deg_en1_min=new G4RotationMatrix();
          pRotY90deg_en1_min->rotateY(90*deg-ang_en_max_1);
//          new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(-1.*pDz_1,0,0),
//minchev          new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(-1.*pDz_1+xmin_sep_en1,0,z_sept_en_min1), trrap,"vac Box en1 min",hallLogical,0,0,0);
          new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(-1.*pDz_1+xmin_sep_en1-0.02*mm,0,z_sept_en_min1), 
                            trrap,"vac Box en1 min",magneticLogical,0,0,0);


          vector<G4TwoVector> vertices_en1_max;
          vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
          vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
          vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
          vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );
          vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
          vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
          vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
          vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );



          G4VSolid* testTrap_en1_max = new G4GenericTrap("testTrap_en1_max", pDz_1, vertices_en1_max);
          G4LogicalVolume* trrap_en1_max = new G4LogicalVolume(testTrap_en1_max, mMaterialManager->tungsten,"trrap_en1_max",0,0,LarmStepLimits);
          trrap_en1_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *pRotY90deg_en1_max=new G4RotationMatrix();
          pRotY90deg_en1_max->rotateY(90.*deg-ang_en_min_1);
//          new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1,0,0),
//minchev          new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1+xmax_sep_en1,0,z_sept_en_max1), trrap_en1_max,"trrap_en1_max",hallLogical,0,0,0);
          new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1+xmax_sep_en1+0.02*mm,0,z_sept_en_max1), 
                            trrap_en1_max,"trrap_en1_max",magneticLogical,0,0,0);

          vector<G4TwoVector> vertices_en2_max;
          vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
          vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
          vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
          vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );
          vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
          vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
          vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
          vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );

          G4VSolid* testTrap_en2_max = new G4GenericTrap("testTrap_en2_max", pDz_1, vertices_en2_max);
          G4LogicalVolume* trrap_en2_max = new G4LogicalVolume(testTrap_en2_max, mMaterialManager->tungsten,"trrap_en2_max",0,0,LarmStepLimits);
          trrap_en2_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *pRotY90deg_en2_max=new G4RotationMatrix();
          pRotY90deg_en2_max->rotateY(90.*deg-ang_en_max_2);
//          new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1,0,0),
//          new G4PVPlacement(pRotY90deg_en2_max,G4ThreeVector(pDz_1+xmax_sep_en2,0,z_sept_en_max2), trrap_en2_max,"trrap_en2_max",hallLogical,0,0,0);
          new G4PVPlacement(pRotY90deg_en2_max,G4ThreeVector(pDz_1+xmax_sep_en2+0.02*mm,0,z_sept_en_max2), 
                            trrap_en2_max,"trrap_en2_max",magneticLogical,0,0,0);


          vector<G4TwoVector> vertices_en2_min;
          vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
          vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
          vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
          vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );
          vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
          vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
          vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
          vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );

          G4VSolid* testTrap_en2_min = new G4GenericTrap("testTrap_en2_min", pDz_1, vertices_en2_min);
          G4LogicalVolume* trrap_en2_min = new G4LogicalVolume(testTrap_en2_min, mMaterialManager->tungsten,"trrap_en2_min",0,0,LarmStepLimits);
          trrap_en2_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *pRotY90deg_en2_min=new G4RotationMatrix();
          pRotY90deg_en2_min->rotateY(90.*deg-ang_en_min_2);
//          new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(pDz_1,0,0),
//          new G4PVPlacement(pRotY90deg_en2_min,G4ThreeVector(-1.*pDz_1+xmin_sep_en2,0,z_sept_en_min2), trrap_en2_min,"trrap_en2_min",hallLogical,0,0,0);
          new G4PVPlacement(pRotY90deg_en2_min,G4ThreeVector(-1.*pDz_1+xmin_sep_en2-0.02*mm,0,z_sept_en_min2), 
                            trrap_en2_min,"trrap_en2_min",magneticLogical,0,0,0);

          vector<G4TwoVector> vertices_cov1_up;
          vertices_cov1_up.push_back( G4TwoVector(xmin_sep_en1-1.95*pDz_1, z_sept_en_min1) );
          vertices_cov1_up.push_back( G4TwoVector(xmax_sep_en1+1.95*pDz_1, z_sept_en_max1) );
          vertices_cov1_up.push_back( G4TwoVector(xmax_sep_ex1+1.95*pDz_1, z_sept_ex_max1) );
          vertices_cov1_up.push_back( G4TwoVector(xmin_sep_ex1-1.95*pDz_1, z_sept_ex_min1) );

          vertices_cov1_up.push_back( G4TwoVector(xmin_sep_en1-1.95*pDz_1, z_sept_en_min1) );
          vertices_cov1_up.push_back( G4TwoVector(xmax_sep_en1+1.95*pDz_1, z_sept_en_max1) );
          vertices_cov1_up.push_back( G4TwoVector(xmax_sep_ex1+1.95*pDz_1, z_sept_ex_max1) );
          vertices_cov1_up.push_back( G4TwoVector(xmin_sep_ex1-1.95*pDz_1, z_sept_ex_min1) );

          G4VSolid* testTrap_cov_1_up = new G4GenericTrap("testTrap_cov_1_up", pDz_2, vertices_cov1_up);
          G4LogicalVolume* trrap_cov_1_up = new G4LogicalVolume(testTrap_cov_1_up, mMaterialManager->tungsten,"trrap_cov_1_up",0,0,LarmStepLimits);
          trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
//          trrap_cov_1_up->SetVisAttributes(G4VisAttributes::Invisible);

          G4RotationMatrix *pRotX90deg_cov_1_up=new G4RotationMatrix();
          pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
//          new G4PVPlacement(pRotX90deg_cov_1_up,G4ThreeVector(0,ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0), trrap_cov_1_up,"cov1_up",hallLogical,0,0,0);
          new G4PVPlacement(pRotX90deg_cov_1_up,G4ThreeVector(0,0.02*mm+ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0), 
                            trrap_cov_1_up,"cov1_up",magneticLogical,0,0,0);

          G4VSolid* testTrap_cov_1_bot = testTrap_cov_1_up;
          G4LogicalVolume* trrap_cov_1_bot = new G4LogicalVolume(testTrap_cov_1_bot, mMaterialManager->tungsten,"trrap_cov_1_bot",0,0,LarmStepLimits);
          trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          G4RotationMatrix *pRotX90deg_cov_1_bot=new G4RotationMatrix();
          pRotX90deg_cov_1_bot->rotateX(-90.*deg-atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
          new G4PVPlacement(pRotX90deg_cov_1_bot,G4ThreeVector(0,-0.02*mm-1.*ymax_sep_en1-pDz_2-fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0), 
                            trrap_cov_1_bot,"cov1_bot",magneticLogical,0,0,0);
//          trrap_cov_1_bot->SetVisAttributes(G4VisAttributes::Invisible);
          trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));




          vector<G4TwoVector> vertices_cov2_up;
          vertices_cov2_up.push_back( G4TwoVector(xmin_sep_en2-1.95*pDz_1, z_sept_en_min2) );
          vertices_cov2_up.push_back( G4TwoVector(xmax_sep_en2+1.95*pDz_1, z_sept_en_max2) );
          vertices_cov2_up.push_back( G4TwoVector(xmax_sep_ex2+1.95*pDz_1, z_sept_ex_max2) );
          vertices_cov2_up.push_back( G4TwoVector(xmin_sep_ex2-1.95*pDz_1, z_sept_ex_min2) );

          vertices_cov2_up.push_back( G4TwoVector(xmin_sep_en2-1.95*pDz_1, z_sept_en_min2) );
          vertices_cov2_up.push_back( G4TwoVector(xmax_sep_en2+1.95*pDz_1, z_sept_en_max2) );
          vertices_cov2_up.push_back( G4TwoVector(xmax_sep_ex2+1.95*pDz_1, z_sept_ex_max2) );
          vertices_cov2_up.push_back( G4TwoVector(xmin_sep_ex2-1.95*pDz_1, z_sept_ex_min2) );

          G4VSolid* testTrap_cov_2_up = new G4GenericTrap("testTrap_cov_2_up", pDz_2, vertices_cov2_up);
          G4LogicalVolume* trrap_cov_2_up = new G4LogicalVolume(testTrap_cov_2_up, mMaterialManager->tungsten,"trrap_cov_2_up",0,0,LarmStepLimits);
          trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *pRotX90deg_cov_2_up=new G4RotationMatrix();
          pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
          new G4PVPlacement(pRotX90deg_cov_2_up,G4ThreeVector(0, 0.02*mm+ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0), 
                            trrap_cov_2_up,"cov2_up",magneticLogical,0,0,0);

          G4VSolid* testTrap_cov_2_bot = testTrap_cov_2_up;
          G4LogicalVolume* trrap_cov_2_bot = new G4LogicalVolume(testTrap_cov_2_bot, mMaterialManager->tungsten,"trrap_cov_2_bot",0,0,LarmStepLimits);
          trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          G4RotationMatrix *pRotX90deg_cov_2_bot=new G4RotationMatrix();
          pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
          new G4PVPlacement(pRotX90deg_cov_2_bot,G4ThreeVector(0,-0.02*mm-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0), 
                            trrap_cov_2_bot,"cov2_bot",magneticLogical,0,0,0);

/*
//          G4VSolid* vac_box_L;//= testTrap;
//          testTrap_en1_max
//          testTrap_en2_max
//          testTrap_en2_min
//          testTrap_cov_1_up
//          testTrap_cov_1_bot
//          testTrap_cov_2_up
//          testTrap_cov_2_bot
          G4UnionSolid* vac_box_L = new G4UnionSolid("Left Vac Box",
              testTrap, testTrap_en1_max, 0, G4ThreeVector(0.,0.,0.));

          G4UnionSolid* vac_box_L1 = new G4UnionSolid("Left Vac Box",
              vac_box_L, testTrap_en2_max, 0, G4ThreeVector(0.,0.,0.));

          G4UnionSolid* vac_box_L2 = new G4UnionSolid("Left Vac Box",
              vac_box_L1, testTrap_en2_min, 0, G4ThreeVector(0.,0.,0.));

          G4UnionSolid* vac_box_L3 = new G4UnionSolid("Left Vac Box",
              vac_box_L2, testTrap_cov_1_up, 0, G4ThreeVector(0.,0.,0.));

          G4UnionSolid* vac_box_L4 = new G4UnionSolid("Left Vac Box",
              vac_box_L3, testTrap_cov_1_bot, 0, G4ThreeVector(0.,0.,0.));

          G4UnionSolid* vac_box_L5 = new G4UnionSolid("Left Vac Box",
              vac_box_L4, testTrap_cov_2_up, 0, G4ThreeVector(0.,0.,0.));

          G4UnionSolid* vac_box_L6 = new G4UnionSolid("Left Vac Box",
              vac_box_L5, testTrap_cov_2_bot, 0, G4ThreeVector(0.,0.,0.));


          G4LogicalVolume* vac_box_R6 = new G4LogicalVolume(vac_box_L6, mMaterialManager->tungsten,"Right Vac Box",0,0,LarmStepLimits);
          vac_box_R6->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          G4RotationMatrix *vac_box_R_rot=new G4RotationMatrix();
          vac_box_R_rot->rotateX(180.*deg);
          new G4PVPlacement(vac_box_R_rot,G4ThreeVector(0,0,0),
               vac_box_R6,"Right Vac Box",hallLogical,0,0,0);

*/








          vector<G4TwoVector> r_vertices_en1_min;
          r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
          r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
          r_vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
          r_vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );
          r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
          r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
          r_vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
          r_vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );

          G4VSolid* r_testTrap = new G4GenericTrap("vac Box 1", pDz_1, r_vertices_en1_min);
          G4LogicalVolume* r_trrap = new G4LogicalVolume(r_testTrap, mMaterialManager->tungsten,"vac Box 1",0,0,LarmStepLimits);
          r_trrap->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *r_pRotY90deg_en1_min=new G4RotationMatrix();
          r_pRotY90deg_en1_min->rotateY(90*deg+ang_en_max_1);
//          new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(-1.*pDz_1,0,0),
          new G4PVPlacement(r_pRotY90deg_en1_min,G4ThreeVector(1.*pDz_1-xmin_sep_en1,0,z_sept_en_min1), 
                            r_trrap,"vac Box en1 min",magneticLogical,0,0,0);


          vector<G4TwoVector> r_vertices_en1_max;
          r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
          r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
          r_vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
          r_vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );
          r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
          r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
          r_vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
          r_vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );



          G4VSolid* r_testTrap_en1_max = new G4GenericTrap("testTrap_en1_max", pDz_1, r_vertices_en1_max);
          G4LogicalVolume* r_trrap_en1_max = new G4LogicalVolume(r_testTrap_en1_max, mMaterialManager->tungsten,"trrap_en1_max",0,0,LarmStepLimits);
          r_trrap_en1_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *r_pRotY90deg_en1_max=new G4RotationMatrix();
          r_pRotY90deg_en1_max->rotateY(90.*deg+ang_en_min_1);
//          new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1,0,0),
          new G4PVPlacement(r_pRotY90deg_en1_max,G4ThreeVector(-1.*(pDz_1+xmax_sep_en1),0,z_sept_en_max1),
               r_trrap_en1_max,"trrap_en1_max",magneticLogical,0,0,0);


          vector<G4TwoVector> r_vertices_en2_max;
          r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
          r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
          r_vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
          r_vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );
          r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
          r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
          r_vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
          r_vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );

          G4VSolid* r_testTrap_en2_max = new G4GenericTrap("testTrap_en2_max", pDz_1, vertices_en2_max);
          G4LogicalVolume* r_trrap_en2_max = new G4LogicalVolume(r_testTrap_en2_max, mMaterialManager->tungsten,"trrap_en2_max",0,0,LarmStepLimits);
          r_trrap_en2_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *r_pRotY90deg_en2_max=new G4RotationMatrix();
          r_pRotY90deg_en2_max->rotateY(90.*deg+ang_en_max_2);
//          new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1,0,0),
          new G4PVPlacement(r_pRotY90deg_en2_max,G4ThreeVector(-1.*(pDz_1+xmax_sep_en2),0,z_sept_en_max2),
               r_trrap_en2_max,"trrap_en2_max",magneticLogical,0,0,0);


          vector<G4TwoVector> r_vertices_en2_min;
          r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
          r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
          r_vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
          r_vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );
          r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
          r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
          r_vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
          r_vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );

          G4VSolid* r_testTrap_en2_min = new G4GenericTrap("testTrap_en2_min", pDz_1, vertices_en2_min);
          G4LogicalVolume* r_trrap_en2_min = new G4LogicalVolume(r_testTrap_en2_min, mMaterialManager->tungsten,"trrap_en2_min",0,0,LarmStepLimits);
          r_trrap_en2_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *r_pRotY90deg_en2_min=new G4RotationMatrix();
          r_pRotY90deg_en2_min->rotateY(90.*deg+ang_en_min_2);
//          new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(pDz_1,0,0),
          new G4PVPlacement(r_pRotY90deg_en2_min,G4ThreeVector(pDz_1-xmin_sep_en2,0,z_sept_en_min2),
               r_trrap_en2_min,"trrap_en2_min",magneticLogical,0,0,0);


          vector<G4TwoVector> r_vertices_cov1_up;
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_en1-1.95*pDz_1), z_sept_en_min1) );
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_en1+1.95*pDz_1), z_sept_en_max1) );
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_ex1+1.95*pDz_1), z_sept_ex_max1) );
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_ex1-1.95*pDz_1), z_sept_ex_min1) );

          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_en1-1.95*pDz_1), z_sept_en_min1) );
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_en1+1.95*pDz_1), z_sept_en_max1) );
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_ex1+1.95*pDz_1), z_sept_ex_max1) );
          r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_ex1-1.95*pDz_1), z_sept_ex_min1) );

          G4VSolid* r_testTrap_cov_1_up = new G4GenericTrap("testTrap_cov_1_up", pDz_2, r_vertices_cov1_up);
          G4LogicalVolume* r_trrap_cov_1_up = new G4LogicalVolume(r_testTrap_cov_1_up, mMaterialManager->tungsten,"trrap_cov_1_up",0,0,LarmStepLimits);
          r_trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
//          trrap_cov_1_up->SetVisAttributes(G4VisAttributes::Invisible);

          G4RotationMatrix *r_pRotX90deg_cov_1_up=new G4RotationMatrix();
          r_pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
          new G4PVPlacement(r_pRotX90deg_cov_1_up,G4ThreeVector(0,ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0),
               r_trrap_cov_1_up,"cov1_up",magneticLogical,0,0,0);

          G4VSolid* r_testTrap_cov_1_bot = r_testTrap_cov_1_up;
          G4LogicalVolume* r_trrap_cov_1_bot = new G4LogicalVolume(r_testTrap_cov_1_bot, mMaterialManager->tungsten,"trrap_cov_1_bot",0,0,LarmStepLimits);
          r_trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          G4RotationMatrix *r_pRotX90deg_cov_1_bot=new G4RotationMatrix();
          r_pRotX90deg_cov_1_bot->rotateX(-90.*deg-atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
          new G4PVPlacement(r_pRotX90deg_cov_1_bot,G4ThreeVector(0,-1.*ymax_sep_en1-pDz_2-fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0),
               r_trrap_cov_1_bot,"cov1_bot",magneticLogical,0,0,0);
//          trrap_cov_1_bot->SetVisAttributes(G4VisAttributes::Invisible);
          r_trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));



          vector<G4TwoVector> r_vertices_cov2_up;
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_en2-1.95*pDz_1), z_sept_en_min2) );
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_en2+1.95*pDz_1), z_sept_en_max2) );
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_ex2+1.95*pDz_1), z_sept_ex_max2) );
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_ex2-1.95*pDz_1), z_sept_ex_min2) );

          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_en2-1.95*pDz_1), z_sept_en_min2) );
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_en2+1.95*pDz_1), z_sept_en_max2) );
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_ex2+1.95*pDz_1), z_sept_ex_max2) );
          r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_ex2-1.95*pDz_1), z_sept_ex_min2) );

          G4VSolid* r_testTrap_cov_2_up = new G4GenericTrap("testTrap_cov_2_up", pDz_2, r_vertices_cov2_up);
          G4LogicalVolume* r_trrap_cov_2_up = new G4LogicalVolume(r_testTrap_cov_2_up, mMaterialManager->tungsten,"trrap_cov_2_up",0,0,LarmStepLimits);
          r_trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));

          G4RotationMatrix *r_pRotX90deg_cov_2_up=new G4RotationMatrix();
          r_pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
          new G4PVPlacement(r_pRotX90deg_cov_2_up,G4ThreeVector(0,     ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0),
               r_trrap_cov_2_up,"cov2_up",magneticLogical,0,0,0);

          G4VSolid* r_testTrap_cov_2_bot = r_testTrap_cov_2_up;
          G4LogicalVolume* r_trrap_cov_2_bot = new G4LogicalVolume(r_testTrap_cov_2_bot, mMaterialManager->tungsten,"trrap_cov_2_bot",0,0,LarmStepLimits);
          r_trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          G4RotationMatrix *r_pRotX90deg_cov_2_bot=new G4RotationMatrix();
          r_pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
          new G4PVPlacement(r_pRotX90deg_cov_2_bot,G4ThreeVector(0,-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0),
               r_trrap_cov_2_bot,"cov2_bot",magneticLogical,0,0,0);

        } // end of septum vacuum box



       //sieve setup. Both Sieve Slits are built in this block
       if (mSetupLSieveSlit)
       {
//         if (0)
         {
//           G4double target_width=0.5*1.185*mm;
           G4double target_width=0.200*mm;
//           G4double target_width=0.5*0.0098*mm;
//           G4double target_width=0.35*mm;
           G4double target_x=25.*cm;
           G4double target_y=25.*cm;
           G4VSolid* targ1=new G4Box("rad_tail_solid",target_x/2.,target_y/2.,target_width/2.);
           G4LogicalVolume* RadTailLogical = new G4LogicalVolume(targ1,mMaterialManager->tungsten,"RadTailLogical",0,0,0);//LarmStepLimits);
           RadTailLogical->SetVisAttributes((G4Colour(0.998, 0.008, 0.000)));

           G4VSolid* tail_det = new G4Box("Targ trans",target_x/2.,target_y/2.,0.001*target_width/2.);
           G4LogicalVolume* RadTailDetLogical = new G4LogicalVolume(tail_det,mMaterialManager->vacuum,"RadTailDetLogical",0,0,LarmStepLimits);
           RadTailDetLogical->SetVisAttributes((G4Colour(0.00998, 0.99008, 0.000)));

           new G4PVPlacement(0,G4ThreeVector(0,0, Rad_tail_foil_postn*cm),RadTailLogical ,"RadTail",magneticLogical,0,0,0);
           new G4PVPlacement(0,G4ThreeVector(0,0, (Rad_tail_foil_postn+1.)*cm),       RadTailDetLogical ,"RadTailDet",magneticLogical,0,0,0);
           cout<<"Rad_tail_foil_position="<<Rad_tail_foil_postn<<endl;
         }





         //sieve geometry is 8.38 x 8.38 x 1 inches, gap distance is 0.492 inch (horizontal) and 0.984 inch (vertical), gap diameter is 0.157e-3 inch, two gaps with 0.236e-3 inch
         double inch=2.54*cm;

         double pSieveSlitX=3.25*inch; //
         double pSieveSlitY=4.25*inch; //
         double pSieveSlitZ=sieve_thickness*inch;
         cout<<"sieve thickness is "<<sieve_thickness<<" inches"<<endl;

         double sieve_distance = 31.23*inch;

         double sieve_pos_z=-105.30001365*cm+(sieve_distance + pSieveSlitZ/2.)*cos(5.81*deg);
         double sieve_pos_x=(sieve_distance + pSieveSlitZ/2.)*sin(5.81*deg);
/*

         double startphi=0.0*deg;
         double deltaphi=360.0*deg;

//         double pSieveSlitHoleR=2.*mm;           //radius of small hole 0.055/2 inch
//         double pSieveSlitLargeHoleR=3.*mm;      //radius of large hole 0.106/2 inch
         double pSieveSlitHoleR=0.5*0.055*inch;           //radius of small hole 0.055/2 inch
         double pSieveSlitLargeHoleR=0.5*0.106*inch;      //radius of large hole 0.106/2 inch
         double pSieveSlitMediumHoleR=1.*0.5*0.075*inch;           //radius of small hole 0.055/2 inch
         double pSieveSlitHoldL=pSieveSlitZ+0.1*mm;  //need to make it longer to avoid round off in the subtraction

         //the whole position relative to the slit center 
         double pSieveSlitHolePosH[15]={0.295, 0.485, 0.675, 0.865, 1.055, 1.245, 1.435, 1.625, 1.815, 2.005, 2.195, 2.385, 2.575, 2.765, 2.955};
         double pSieveSlitHolePosV[9] = {0.285, 0.745, 1.205, 1.665, 2.125, 2.585, 3.045, 3.505, 3.965};
         for(int ii=0;ii<15;ii++)
         {
            pSieveSlitHolePosH[ii] = (pSieveSlitHolePosH[ii]-1.625)*inch;
         }
         for(int ii=0;ii<9;ii++)
         {
            pSieveSlitHolePosV[ii] = (pSieveSlitHolePosV[ii]-2.125)*inch;
         }
         //now start to build box then subtract 63 holes 
         G4VSolid* sieveSlitWholeSolid=new G4Box("sieveSlitWholeBox",pSieveSlitX/2.0,
              pSieveSlitY/2.0,pSieveSlitZ/2.0);
         startphi=0.0*deg;deltaphi=360.0*deg;

         G4VSolid* sieveSlitHoleSolid=new G4Tubs("sieveSlitHoleTubs",0,pSieveSlitHoleR,
              pSieveSlitHoldL/2.0,startphi,deltaphi); 

         G4VSolid* sieveSlitLargeHoleSolid=new G4Tubs("sieveSlitLargeHoleTubs",0,
              pSieveSlitLargeHoleR,pSieveSlitHoldL/2.0,startphi,deltaphi); 

         G4VSolid* sieveSlitMediumHoleSolid=new G4Tubs("sieveSlitMediumHoleTubs",0,
              pSieveSlitMediumHoleR,pSieveSlitHoldL/2.0,startphi,deltaphi); 

         G4SubtractionSolid* sieveSlitSolid=(G4SubtractionSolid*)sieveSlitWholeSolid;
         char strName[100];
         for(int ih=0;ih<15;ih++)
         {
           for(int iv=0;iv<9;iv++)
           {
             sprintf(strName,"sieveSlitHole_H%d_V%d",ih,iv);
             if((ih==7 && iv==4) || (ih==3 && iv==2)) 
             {
             //now dig large holes in the block
                sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
                    sieveSlitLargeHoleSolid,0,
                    G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv],0));
             }
             else
             {
                //now dig small holes in the block
                sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
                    sieveSlitHoleSolid,0,
                    G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv], 0));
             }
             if ((iv>=7) || (iv<2))
             {
               sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
                    sieveSlitMediumHoleSolid,0,
                    G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv], 0.3*inch));
             }
           }
         }


         for(int ih=0;ih<14-2;ih++)
         {
           for(int iv=0;iv<8;iv++)
           {
             sprintf(strName,"sieveSlitHole_H%d_V%d",ih,iv);
             if ( !( ( ih == 6 ) && (iv>0) && (iv<7) ) )
             {
                //now dig small holes in the block
                sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
                    sieveSlitHoleSolid,0,
                    G4ThreeVector(0.5*(pSieveSlitHolePosH[ih]+pSieveSlitHolePosH[ih+1]), 0.5*(pSieveSlitHolePosV[iv]+pSieveSlitHolePosV[iv+1]), 0));
             }
             if ((iv==7) || (iv==0))
             {
               sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
                    sieveSlitMediumHoleSolid,0,
                    G4ThreeVector(0.5*(pSieveSlitHolePosH[ih]+pSieveSlitHolePosH[ih+1]), 0.5*(pSieveSlitHolePosV[iv]+pSieveSlitHolePosV[iv+1]), 0.3*inch));
             }
           }
         }

         sieveSlitSolid->SetName("sieveSlitSolid");

         G4LogicalVolume* sieveSlitLogical = new G4LogicalVolume(sieveSlitSolid,
              mMaterialManager->tungsten,"sieveSlitLogical",0,0,LarmStepLimits);
         sieveSlitLogical->SetVisAttributes(IronVisAtt);

//         new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-sieve_pos_z,0),
//                         sieveSlitLogical,"SievePhys",LHRSContainerLogical,0,0,0);


*/

         G4RotationMatrix *R_RotY90deg_sieve_slit=new G4RotationMatrix();
         R_RotY90deg_sieve_slit->rotateY(-5.*deg);
//         new G4PVPlacement(R_RotY90deg_sieve_slit,G4ThreeVector(sieve_pos_x, 0, sieve_pos_z),
//               sieveSlitLogical,"RSievePhys", magneticLogical, 0, 0, 0);

         G4RotationMatrix *L_RotY90deg_sieve_slit=new G4RotationMatrix();
         L_RotY90deg_sieve_slit->rotateY(5.*deg);
//         new G4PVPlacement(L_RotY90deg_sieve_slit,G4ThreeVector(-sieve_pos_x, 0, sieve_pos_z),
//               sieveSlitLogical,"LSievePhys", magneticLogical, 0, 0, 0);


         G4VSolid* sieveBackPlaneSolid=new G4Box("sieveBackPlaneSolid", 0.75*pSieveSlitX, 0.75*pSieveSlitY,0.001*mm);
         G4LogicalVolume* sieveSlitBackLogical = new G4LogicalVolume(sieveBackPlaneSolid, mMaterialManager->vacuum, "sieveBackPlaneLogical",0,0,LarmStepLimits);
         new G4PVPlacement(R_RotY90deg_sieve_slit,G4ThreeVector(sieve_pos_x, 0, sieve_pos_z+1.5*cm), sieveSlitBackLogical,"RSvSlBack", magneticLogical, 0, 0, 0);
         new G4PVPlacement(L_RotY90deg_sieve_slit,G4ThreeVector(-sieve_pos_x, 0, sieve_pos_z+1.5*cm), sieveSlitBackLogical,"LSvSlBack", magneticLogical, 0, 0, 0);
         sieveSlitBackLogical->SetVisAttributes(G4VisAttributes::Invisible);
       }




#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY>=0)
	{
		DumpGeometricalTree(mHallPhysVol);
	}
#endif
	//create the geometry macro to manually modify the geometry in visulization
	GenerateGeometryMacro(mHallPhysVol,0);
	return mHallPhysVol;
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* HRSDetectorConstruction::ConstructHRS(G4LogicalVolume* motherLogical)
{
	int mSnakeModel;
	gConfig->GetArgument("SnakeModel",mSnakeModel);

	int IS_NIM = 0;//if it's 0, then we go by SNAKE, if it is 1, we go by NIM.
	//As it turns out, I think we want to go by SNAKE, so that is why it is zero right now.
	G4VPhysicalVolume* theHRSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
	pRotX45deg->rotateX(-45*deg); 
	G4RotationMatrix *pRotX30deg=new G4RotationMatrix();
	pRotX30deg->rotateX(60*deg); 
	G4RotationMatrix *pRotX105deg=new G4RotationMatrix();
	pRotX105deg->rotateX(165*deg); 
	G4RotationMatrix *pRotXLHRSdeg=new G4RotationMatrix();
	pRotXLHRSdeg->rotateY(mLHRSAngle); 
	G4RotationMatrix *pRotXRHRSdeg=new G4RotationMatrix();
	pRotXRHRSdeg->rotateY(mRHRSAngle); 

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDname;
	G4VSensitiveDetector* Q1WindowSD=new HRSStdSD(SDname="Q1Window");
	//Build the HRS, positions were taken fron the NIM paper
	//The Apertures are defined as:
	//Q1 exit is a circle of radius 0.15m
	//Q3 entrance and exit are circles of radius 0.30m
	//The dipole entrance and exit are trapezoids:
	//-0.40 m<x<0.40 m (+x is down) 
	//y=+-(0.125*(1-(1.25*x/8.40)) (smallest gap at positive x, x in m)

	/////////////////////////
	// HRS QQDQ Containner
	/////////////////////////
	//Build a container using polycone, covering 270+/-12 degrees, 1.46m to 20.76m
	//3.05m below the beam line is the ground, need to subtract everything beneath that
	//looks like this:
	/*
	HRS container:covering 24 degrees in X-Z plane, 10.05 m height.
	Stuff inside this containner will also use the pivot as the origin, but do not 
	need to worry about the rotation of the HRS.
	//                         7.0m above beam line 
	//                       -----------------------------| 20.76m from pivot,
	//                      /                             |
	//                     /                        Q3    |
	//                    /                               |
	//                   /                       E        |
	//                  /                     L           |
	//          --------                   O              |
	//         /                        P                 |
	//    -----                     I                     |
	//----|----- Q1 --- Q2 ---- D     --------------------|------ beam line -----
	//    |                                               |
	//    ------------------------------------------------|
	//    1.46m from povot, 3.05 m below beam line

	*/

	//double pHRSContainerRin=1.46*m,pHRSContainerRout=20.76*m;//JixieMode
	double pHRSContainerRin=0.46*m,pHRSContainerRout=25*m;//Vardan mode for the sieve
	//double pBeamLine2Ground;
	//pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {pHRSContainerRin,pHRSContainerRin,2.5*m,
		3.7*m,9.0*m,pHRSContainerRout-3.0*m,pHRSContainerRout};
	double rOuter_HRSContainer[] = {pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,
		pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,pHRSContainerRout};
	//double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
	//2.0*m,2.0*m,7.0*m,7.0*m};
	double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
					2.0*m,2.0*m,11.0*m,11.0*m};
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",258.0*deg,24.0*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	////build the container using tube
	//double pHRSContainerHeight=14.0*m;
	//G4VSolid* HRSContainerTub = new G4Tubs("HRSContainerTub",pHRSContainerRin,
	//	pHRSContainerRout,pHRSContainerHeight/2.0,258.0*deg,24.0*deg);
	////ground
	//double pGoundHeight=10*m;
	//G4VSolid* GroundTub = new G4Tubs("GroundTub",0,30*m,pGoundHeight/2,0*deg,360.0*deg);
	////tube subtract the ground
	//G4SubtractionSolid* HRSContainerSolid = new G4SubtractionSolid("HRSContainer",
	//		HRSContainerTub,GroundTub,
	//		0,G4ThreeVector(0,0,pBeamLine2Ground-pGoundHeight/2));

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"LHRSContainerLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"RHRSContainerLogical",0,0,LarmStepLimits);

//	G4VisAttributes* MagFieldVisAtt = new G4VisAttributes(G4VisAttributes::Invisible);	
	G4VisAttributes* MagFieldVisAtt = new G4VisAttributes(G4Colour(1., 0., 1.));

	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt);
	//RHRSContainerLogical->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-mLHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(-mRHRSAngle);  

	cout<<" -----> The Setup is ----------> "<<mSetupLHRS<<endl;
	if(mSetupLHRS>=2)
	{
		new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0,0,0),
			LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	}
	if(mSetupRHRS>=2)
	{
		new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0),
			RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);
	}


	//////////////////////////
	//Q1 Entrance Window //
	//////////////////////////
	bool pSetupQ1EntranceWindow=true;
	//G2P and CREX have their own VB defined at a different location
	if( mSetupCREXGeometry)  pSetupQ1EntranceWindow=false;
	if( mSetupAPEXGeometry)  pSetupQ1EntranceWindow=false;
	if(pSetupQ1EntranceWindow) 
	{

		//this part is trying to place a virtual boundary at the Q1 entrance

		//Place both left and right VB for HRS, which is pHRSContainerRin+4*mm away from the 
		//hall center(1.462m). This aperture is a round disk of 29.8 cm diameter and 2 mm thick
		//The real Q1 vacumn entrance to hall center is 1.312m, 

		double pHRSQ1WinThick = 2*mm;
		G4VSolid* HRSQ1WinSolid = new G4Tubs("HRSQ1WinTub",0.0,14.9*cm,
			pHRSQ1WinThick/2,0.0,360.0*deg);
		G4LogicalVolume* HRSQ1WinLogical = new G4LogicalVolume(HRSQ1WinSolid,
			mMaterialManager->mylar,"HRSQ1WinLogical",0,0,LarmStepLimits);
		SDman->AddNewDetector(Q1WindowSD);
		HRSQ1WinLogical->SetSensitiveDetector(Q1WindowSD);
		HRSQ1WinLogical->SetVisAttributes(LightYellowVisAtt); 

		//since the container has been rotated by 90 deg about x axis,
		//y'=z  z'=-y ==> I have to put this offset as -y
		double pHallCenter2Q1Win=pHRSContainerRin+4*mm;  //place it at the first 1.464 m
		gConfig->SetParameter("Pivot2LHRSVBFace",pHallCenter2Q1Win-mPivotZOffset*cos(mLHRSAngle));
		gConfig->SetParameter("Pivot2RHRSVBFace",pHallCenter2Q1Win-mPivotZOffset*cos(mRHRSAngle)); 

		if(mSetupLHRS)
		{
			new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pHallCenter2Q1Win,0),
				HRSQ1WinLogical,"virtualBoundaryPhys_LHRSQ1Win",LHRSContainerLogical,0,0,0);
		}
		if(mSetupRHRS)
		{
			new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pHallCenter2Q1Win,0),
				HRSQ1WinLogical,"virtualBoundaryPhys_RHRSQ1Win",RHRSContainerLogical,0,0,0);
		}
	}
	//BELOW ARE SNAKE VALUES, NOT NIM VALUES
	double pTarget =             0.0  * cm;
	double pQ1en   = pTarget + 159.03 * cm;
	double pQ1ex   = pQ1en   +  94.13 * cm;
	double pQ2en   = pQ1ex   + 117.2  * cm;
	double pQ2ex   = pQ2en   + 182.66 * cm;
	int    sos;//   = 1;
	sos=nQ1Sos;
	if (sos) 
	{
	  pQ1en=(pQ1ex+pQ1en)/2.-35.*cm;  //171.095 cm
	  pQ1ex=pQ1en   +  70.00 * cm;
	}
//	if( mSnakeModel >= 52){
//	  sos = 1;
//	}
	//ABOVE ARE SNAKE VALUES, NOT NIM VALUES


	/////////////////////////
	// HRS Q1              // Nickie is also adding the Q1 collimator here
	/////////////////////////
	double pHallCenter2LQ1Face;//=1.69*m;//NIM
	double pHallCenter2RQ1Face;//=1.69*m;//NIM
	//pHallCenter2LQ1Face = ( IS_NIM == 1 ) ? 169.0 * cm : pQ1en;
	//pHallCenter2RQ1Face = ( IS_NIM == 1 ) ? 169.0 * cm : pQ1en;
	pHallCenter2LQ1Face = pQ1en;
	pHallCenter2RQ1Face = pQ1en;
	//double pHallCenter2LQ1FaceMag=1.600*m;//SNAKE
	//double pHallCenter2RQ1FaceMag=1.600*m;//SNAKE
	//double fringe_extension = 1.25;
	//double fringe_extension = 1.2;
	double fringe_extension = 1.00;

	//flag to trigger sos quad instead of old quad for the proper runs
	//K. Allada and B. Schmookler:
	//Magnetic Length = 70 cm
	//Radius to Pole Tip = 12.827  cm
	
	sos=nQ1Sos;
	double q1shift = sos ? 0.0 * m : 0.0 * m;
	double pQ1Rin  = sos ? 12.827 * cm : 15.0*cm;
	double pQ1Rout = sos ? 35.0   * cm : 35.0*cm;//for now, keep the outer same for either
	double pQ1Length;//=80*cm;
	//pQ1Length = ( IS_NIM == 1 ) ? 80.0 * cm : pQ1ex - pQ1en;
	pQ1Length = sos ?  70. * cm : pQ1ex - pQ1en;
	double pQ1PosAct = pQ1ex - pQ1en;
	double pQ1LengthMag=pQ1Length * fringe_extension;
	//double pQ1LengthMag=94.1*cm;
	//double pQ1Length=(1.698M-1.36*m)*2+0.8*m;


        //By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
        G4VSolid* Q1Solid = new G4Tubs("Q1Tub",pQ1Rin,pQ1Rout,pQ1Length/2.0,0.0,360.0*deg);

        //build 2 copy since there are differennt fields involved in it
        G4LogicalVolume* LQ1Logical = new G4LogicalVolume(Q1Solid,
                mMaterialManager->siliconsteel,"LQ1Logical",0,0,LarmStepLimits);
        G4LogicalVolume* RQ1Logical = new G4LogicalVolume(Q1Solid,
                mMaterialManager->siliconsteel,"RQ1Logical",0,0,LarmStepLimits);

        LQ1Logical->SetVisAttributes(IronVisAtt);
        RQ1Logical->SetVisAttributes(IronVisAtt); 
//      LQ1Logical->SetVisAttributes(LightYellowVisAtt); 
//      RQ1Logical->SetVisAttributes(LightYellowVisAtt); 



        G4VSolid* Q1FrontWin= new G4Tubs("Q1Front",0 ,2.*pQ1Rout, 0.001*mm,0.0,360.0*deg);
        G4LogicalVolume* LQ1FrontLogical = new G4LogicalVolume(Q1FrontWin, mMaterialManager->vacuum,"LQ1FrontLogical",0,0,LarmStepLimits);
        new G4PVPlacement(pRotX90deg, G4ThreeVector(0, -1*(pQ1en-0.002*mm), 0), LQ1FrontLogical,"LQ1FrontPhys", LHRSContainerLogical, 0, 0, 0);

        G4LogicalVolume* RQ1FrontLogical = new G4LogicalVolume(Q1FrontWin, mMaterialManager->vacuum,"RQ1FrontLogical",0,0,LarmStepLimits);
        new G4PVPlacement(pRotX90deg, G4ThreeVector(0, -1*(pQ1en-0.002*mm), 0), RQ1FrontLogical,"RQ1FrontPhys", RHRSContainerLogical, 0, 0, 0);

//        LQ1FrontLogical->SetVisAttributes(G4Colour::Blue());
//        RQ1FrontLogical->SetVisAttributes(G4Colour::Blue());
        LQ1FrontLogical->SetVisAttributes(G4VisAttributes::Invisible);
        RQ1FrontLogical->SetVisAttributes(G4VisAttributes::Invisible);


        //dipole entr rectangles
        {
          double dx1, dx2, dy1, dy2, dz;
          dx1=32.33*cm;
          dx2=15.*cm;
          dy1=32.33*cm;
          dy2=50.*cm;
          dz=819.4887*cm-641.8*cm;


          G4VSolid* dip_en_1 = new G4Trd("dip_en_1",dx1,dx2,dy1,dy2,dz/2.+1.*cm); //just a bit longer to have a clean Subtraction
          G4VSolid* dip_en_out = new G4Trd("dip_en_out",dx1+5.*cm,dx2+5.*cm,dy1+5.*cm,dy2+5.*cm,dz/2.);
          G4RotationMatrix *r_dip_anRotX90deg=new G4RotationMatrix();
          r_dip_anRotX90deg->rotateX(-90.*deg);
          G4SubtractionSolid* DipoleEntrance = new G4SubtractionSolid("DipoleEntrance", dip_en_out,dip_en_1,
            0,G4ThreeVector(0,0,0));
          G4LogicalVolume* dip_en_logical_new = new G4LogicalVolume(DipoleEntrance, mMaterialManager->siliconsteel,"dip_en_logical_new",0,0,LarmStepLimits);
          dip_en_logical_new->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          new G4PVPlacement( r_dip_anRotX90deg,G4ThreeVector(0, -1.*(819.4887*cm+641.8*cm)/2., 0) ,
               dip_en_logical_new,"dip_en_new",LHRSContainerLogical,0,0,0);
          new G4PVPlacement( r_dip_anRotX90deg,G4ThreeVector(0, -1.*(819.4887*cm+641.8*cm)/2., 0) ,
               dip_en_logical_new,"dip_en_new",RHRSContainerLogical,0,0,0);
        }


        //q3 exit rectangles
        {
          double dx1=355.6*mm;
          double dx2u=985.34*mm;
          double dx2d=997.664/sqrt(2.)*mm;
          double dy1=171.450;
          double dz_q3_ex=3176.43/2.*mm;
          double y01=5463.582*mm;
          double y02=7210.738*mm;
          double z01=18903.966*mm;
          double z02=20651.122*mm;
          double z0_new=(z01+z02+997.664*mm/2.)/2.;
          double y0_new=(y01+y02+997.664*mm/2.)/2.;

          vector<G4TwoVector> vertices_q3_exit;
          vertices_q3_exit.push_back( G4TwoVector(-1.*dx1, -1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector(-1.*dx1,  1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector( 1.*dx1,  1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector( 1.*dx1, -1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector(-1.*dx2u, -1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector(-1.*dx2u,  1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector( 1.*dx2d,  1.0*dy1) );
          vertices_q3_exit.push_back( G4TwoVector( 1.*dx2d, -1.0*dy1) );

          vector<G4TwoVector> vertices_q3_exit_out;
          vertices_q3_exit_out.push_back( G4TwoVector(-1.*dx1-25.*cm, -1.0*dy1-25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector(-1.*dx1-25.*cm,  1.0*dy1+25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector( 1.*dx1+25.*cm,  1.0*dy1+25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector( 1.*dx1+25.*cm, -1.0*dy1-25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector(-1.*dx2u-25.*cm, -1.0*dy1-25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector(-1.*dx2u-25.*cm,  1.0*dy1+25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector( 1.*dx2d+25.*cm,  1.0*dy1+25.*cm) );
          vertices_q3_exit_out.push_back( G4TwoVector( 1.*dx2d+25.*cm, -1.0*dy1-25.*cm) );

          G4RotationMatrix *pRotY45deg_q3ex=new G4RotationMatrix();
          G4RotationMatrix *pRotY45deg_q3ex_all=new G4RotationMatrix();
          pRotY45deg_q3ex->rotateY(45.*deg);
          pRotY45deg_q3ex_all->rotateZ(90.*deg);
          pRotY45deg_q3ex_all->rotateY(-45.*deg);

          G4VSolid* q3_ex_sol_in = new G4GenericTrap("q3_ex_sol_in", dz_q3_ex, vertices_q3_exit);
          G4VSolid* q3_ex_sol_out = new G4GenericTrap("q3_ex_sol_out", dz_q3_ex-0.01*cm, vertices_q3_exit_out);
          double bs_tmp=100.*cm;
          double bs_tmp2=20.*cm;
          G4VSolid* q3_ex_dummy1 = new G4Box("q3_ex_dummy1", 500.*cm+bs_tmp/2., bs_tmp/2., bs_tmp2+bs_tmp/2.);

          G4SubtractionSolid* q3_ex_sol1  = new G4SubtractionSolid("q3_ex_sol1", q3_ex_sol_out,q3_ex_sol_in,0,G4ThreeVector(0,0,0));
          G4SubtractionSolid* q3_ex_sol  = new G4SubtractionSolid("q3_ex_sol2", q3_ex_sol1,q3_ex_dummy1,pRotY45deg_q3ex,G4ThreeVector(0,0,sqrt(2.)*bs_tmp2+sqrt(2.)*0.5*bs_tmp+(2470.85*mm-dz_q3_ex)));
          G4LogicalVolume* q3_ex_logic = new G4LogicalVolume(q3_ex_sol, mMaterialManager->siliconsteel,"q3_ex_logic",0,0,LarmStepLimits);
          q3_ex_logic->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
          new G4PVPlacement(pRotY45deg_q3ex_all,G4ThreeVector(0,-1.*z0_new,y0_new), q3_ex_logic,"q3_ex_logic left",LHRSContainerLogical,0,0,0);
          new G4PVPlacement(pRotY45deg_q3ex_all,G4ThreeVector(0,-1.*z0_new,y0_new), q3_ex_logic,"q3_ex_logic left",RHRSContainerLogical,0,0,0);

        }





	if(mSetupLHRS>=2){
	  //put it in the container, which also center at the hall center
	  //therefore only the z_at_lab position need to be considered
	  double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1PosAct/2.0) + q1shift;
	  //double pLQ1Pos_Z=pQ1PosAct + q1shift;
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			    LQ1Logical,"LQ1Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=2){
	  double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1PosAct/2.0) + q1shift;
	  //double pRQ1Pos_Z=pQ1PosAct + q1shift;
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			    RQ1Logical,"RQ1Phys",RHRSContainerLogical,0,0,0);
	}

	//vac
	double pQ1vacRin  = pQ1Rin;
	double pQ1vacRout = pQ1Rin + 5. * cm;
	double pQ1vacLength = 0.2 * mm; 
	double pQ1vacCenterZ = pHallCenter2LQ1Face + pQ1PosAct + pQ1vacLength / 2.;
	double pQ1vacCenterY = 0;
	G4VSolid* Q1vacSolid = new G4Tubs("Q1vacTub",pQ1vacRin,pQ1vacRout,pQ1vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ1vacLogical = new G4LogicalVolume(Q1vacSolid,
		mMaterialManager->siliconsteel,"LQ1vacLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ1vacLogical = new G4LogicalVolume(Q1vacSolid,
		mMaterialManager->siliconsteel,"RQ1vacLogical",0,0,LarmStepLimits);

	
	LQ1vacLogical->SetVisAttributes(G4Colour::Green()); 
	RQ1vacLogical->SetVisAttributes(G4Colour::Green()); 
	//LQ1vacLogical->SetVisAttributes(IronVisAtt); 
	//RQ1vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ1vacInContainer=new G4RotationMatrix();
	pRotQ1vacInContainer->rotateX(90*deg); 
	
	if(mSetupLHRS>=4)
	{
	        cout<<"mSetupLHRS="<<mSetupLHRS<<"    pQ1vacCenterZ="<<pQ1vacCenterZ<<"    pQ1vacCentery="<<pQ1vacCenterY<<endl;
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ1vacInContainer,G4ThreeVector(0, -pQ1vacCenterZ,pQ1vacCenterY),
			LQ1vacLogical,"LQ1vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ1vacInContainer,G4ThreeVector(0, -pQ1vacCenterZ,pQ1vacCenterY),
			RQ1vacLogical,"RQ1vacPhys",RHRSContainerLogical,0,0,0);
	}
	

	/////////////////////////
	// HRS Q2 
	/////////////////////////
	double pHallCenter2LQ2Face;//=3.74*m;//NIM
	double pHallCenter2RQ2Face;//=3.74*m;//NIM
	//pHallCenter2LQ2Face = ( IS_NIM == 1 ) ? 374.0 * cm : pQ2en;
	//pHallCenter2RQ2Face = ( IS_NIM == 1 ) ? 374.0 * cm : pQ2en;
	pHallCenter2LQ2Face = pQ2en;
	pHallCenter2RQ2Face = pQ2en;
	//double pHallCenter2LQ2FaceMag=3.696*m;//SNAKE
	//double pHallCenter2RQ2FaceMag=3.696*m;//SNAKE
	double pQ2Rin=30.0*cm;
	double pQ2Rout=75.0*cm;
	double pQ2Length;//=180*cm;//NIM
	//pQ2Length = ( IS_NIM == 1 ) ? 180.0 * cm : pQ2ex - pQ2en;
	pQ2Length = pQ2ex - pQ2en;
	double pQ2LengthMag=pQ2Length * fringe_extension;
	//double pQ2LengthMag=182.6*cm;//SNAKE
//??delete
	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q2Solid = new G4Tubs("Q2Tub",pQ2Rin,pQ2Rout,pQ2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"LQ2Logical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"RQ2Logical",0,0,LarmStepLimits);

	LQ2Logical->SetVisAttributes(IronVisAtt); 
	RQ2Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupLHRS>=3)
	{
		//put it in the container, which also center at the hall center
		double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			LQ2Logical,"LQ2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3)
	{
		double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			RQ2Logical,"RQ2Phys",RHRSContainerLogical,0,0,0);
	}


      if (nQ1Sos==1)  //new quad
      {
        double da=1124.2*mm;  double ra=198.4/2.*mm;
        double db=1197.9*mm; double  rb=203.7/2.*mm;
        double dc=1199.9*mm; double  rc=247.9/2.*mm;
        double dc1=1710.0*mm;
        double dd0=2420.0*mm; double rd0 = 14.9225*cm;
        double dd=2600.7*mm; double  rd=590.6/2.*mm;
        double de=3129.2*mm;
        G4VSolid* vacRa_Solid = new G4Tubs("vacRaSol", ra, (5.*cm+ra), (db - da)/2.0,    0.0,360.0*deg);
        G4VSolid* vacRb_Solid = new G4Tubs("vacRbSol", rb, (5.*cm+rb), (dc - db)/2.0,    0.0,360.0*deg);
        G4VSolid* vacRc1_Solid = new G4Tubs("vacRcSol", rc, (5.*cm+rc), (dc1 - dc)/2.0,    0.0,360.0*deg);
        G4VSolid* vacRc2_Solid = new G4Tubs("vacRcSol", rd0, (5.*cm+rd0), (dd - dd0)/2.0,    0.0,360.0*deg);
        G4VSolid* vacRd_Solid = new G4Tubs("vacRdSol", rd, (5.*cm+rd), (de - dd)/2.0,    0.0,360.0*deg);
        G4LogicalVolume* vacRa_Logical = new G4LogicalVolume(vacRa_Solid, mMaterialManager->siliconsteel,"vacRa_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRb_Logical = new G4LogicalVolume(vacRb_Solid, mMaterialManager->siliconsteel,"vacRb_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRc1_Logical = new G4LogicalVolume(vacRc1_Solid, mMaterialManager->siliconsteel,"vacRc1_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRc2_Logical = new G4LogicalVolume(vacRc2_Solid, mMaterialManager->siliconsteel,"vacRc2_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRd_Logical = new G4LogicalVolume(vacRd_Solid, mMaterialManager->siliconsteel,"vacRd_Logical",0,0,LarmStepLimits);
        G4RotationMatrix *RRotvacContainer=new G4RotationMatrix();
        RRotvacContainer->rotateX(90*deg);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(db + da)/2.,  0), vacRa_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + db)/2.,  0), vacRb_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc1 + dc)/2.,  0), vacRc1_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dd0 + dd)/2.,  0), vacRc2_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(de + dd)/2.,  0), vacRd_Logical, "vacPhys",RHRSContainerLogical,0,0,0);

        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(db + da)/2.,  0), vacRa_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + db)/2.,  0), vacRb_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc1 + dc)/2.,  0), vacRc1_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dd + dd0)/2.,  0), vacRc2_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(de + dd)/2.,  0), vacRd_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        vacRa_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRb_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRc1_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRc2_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRd_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
      }


      if (nQ1Sos==0) //old quad
      {
        double da = 110.911;
        double db = 118.9;
        double dc = 120.864; //12.522
        double dd = 135.064; //15.0cm
        double de = 159.03;
        double de1 = 253.16;
        double df = 300.464;
        double dg = 314.464;
        double dh = 320.36;
        //vacuum tunel from Barak
//        G4VSolid* vacRa_Solid = new G4Tubs("vacRaSol", 9.92*cm,   (5.+9.92)*cm,   3.*cm/2.0,                  0.0,360.0*deg);
//        G4VSolid* vacRb_Solid = new G4Tubs("vacRaSol", 10.185*cm, (5.+10.185)*cm, (db - da)*cm/2.0,    0.0,360.0*deg);
//        G4VSolid* vacRd_Solid = new G4Tubs("vacRaSol", 12.395*cm, (5.+12.395)*cm, (dc - db)*cm/2.0,  0.0,360.0*deg);
        //G4VSolid* vacRe_Solid = new G4Tubs("vacRaSol", 12.395*cm, (5.+12.395)*cm, (241.05  - 171.05)*cm/2.0,  0.0,360.0*deg);
        G4VSolid* vacRd_Solid = new G4Tubs("vacRaSol", 12.522*cm, (5.+12.522)*cm, (dd - dc)*cm/2.0,  0.0,360.0*deg);
        G4VSolid* vacRf_Solid = new G4Tubs("vacRaSol", 15.0*cm, (5.+15.0)*cm, (de - dd)*cm/2.0,  0.0,360.0*deg);
        G4VSolid* vacRg_Solid = new G4Tubs("vacRaSol", 14.922*cm, (5.+14.922)*cm, (df - de1)*cm/2.0,  0.0,360.0*deg);
        G4VSolid* vacRh_Solid = new G4Tubs("vacRaSol", 20.955*cm, (5.+20.955)*cm, (dg - df)*cm/2.0, 0.0,360.0*deg);
        G4VSolid* vacRi_Solid = new G4Tubs("vacRaSol", 30.22*cm,  (5.+30.22)*cm,  (dh - dg)*cm/2.0, 0.0,360.0*deg);
//        G4LogicalVolume* vacRa_Logical = new G4LogicalVolume(vacRa_Solid, mMaterialManager->siliconsteel,"vacRa_Logical",0,0,LarmStepLimits);
//        G4LogicalVolume* vacRb_Logical = new G4LogicalVolume(vacRb_Solid, mMaterialManager->siliconsteel,"vacRb_Logical",0,0,LarmStepLimits);

        G4LogicalVolume* vacRd_Logical = new G4LogicalVolume(vacRd_Solid, mMaterialManager->siliconsteel,"vacRd_Logical",0,0,LarmStepLimits);
//        G4LogicalVolume* vacLa_Logical = new G4LogicalVolume(vacRa_Solid, mMaterialManager->siliconsteel,"vacLa_Logical",0,0,LarmStepLimits);
//        G4LogicalVolume* vacLb_Logical = new G4LogicalVolume(vacRb_Solid, mMaterialManager->siliconsteel,"vacLb_Logical",0,0,LarmStepLimits);
//        G4LogicalVolume* vacLd_Logical = new G4LogicalVolume(vacRd_Solid, mMaterialManager->siliconsteel,"vacLd_Logical",0,0,LarmStepLimits);
        //G4LogicalVolume* vacRe_Logical = new G4LogicalVolume(vacRe_Solid, mMaterialManager->siliconsteel,"vacRe_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRf_Logical = new G4LogicalVolume(vacRf_Solid, mMaterialManager->siliconsteel,"vacRf_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRg_Logical = new G4LogicalVolume(vacRg_Solid, mMaterialManager->siliconsteel,"vacRg_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRh_Logical = new G4LogicalVolume(vacRh_Solid, mMaterialManager->siliconsteel,"vacRh_Logical",0,0,LarmStepLimits);
        G4LogicalVolume* vacRi_Logical = new G4LogicalVolume(vacRi_Solid, mMaterialManager->siliconsteel,"vacRi_Logical",0,0,LarmStepLimits);
        G4RotationMatrix *RRotvacContainer=new G4RotationMatrix();
        RRotvacContainer->rotateX(90*deg);
/*
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*112.48*cm, 0), vacRa_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(db + da)/2.*cm,  0), vacRb_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + db)*cm/2.0, 0), vacRd_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*112.48*cm, 0), vacRa_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(da + db)/2.*cm,  0), vacRb_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + db)*cm/2.0, 0), vacRd_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        //new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, coord, 0), vacRe_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(de + dd)*cm/2.0, 0), vacRf_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(df + de)*cm/2.0, 0), vacRg_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dg + df)*cm/2.0, 0), vacRh_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dh + dg)*cm/2.0, 0), vacRi_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
*/
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*112.48*cm, 0), vacRa_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(db + da)/2.*cm,  0), vacRb_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + db)*cm/2.0, 0), vacRd_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*112.48*cm, 0), vacRa_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(da + db)/2.*cm,  0), vacRb_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + dd)*cm/2.0, 0), vacRd_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        //new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, coord, 0), vacRe_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(de + dd)*cm/2.0, 0), vacRf_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(df + de1)*cm/2.0, 0), vacRg_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dg + df)*cm/2.0, 0), vacRh_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dh + dg)*cm/2.0, 0), vacRi_Logical, "vacPhys",RHRSContainerLogical,0,0,0);

//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*112.48*cm, 0), vacRa_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(db + da)/2.*cm,  0), vacRb_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + db)*cm/2.0, 0), vacRd_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*112.48*cm, 0), vacRa_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
//        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(da + db)/2.*cm,  0), vacRb_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dc + dd)*cm/2.0, 0), vacRd_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        //new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, coord, 0), vacRe_Logical, "vacPhys",RHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(de + dd)*cm/2.0, 0), vacRf_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(df + de1)*cm/2.0, 0), vacRg_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dg + df)*cm/2.0, 0), vacRh_Logical, "vacPhys",LHRSContainerLogical,0,0,0);
        new G4PVPlacement(RRotvacContainer, G4ThreeVector(0, -1.*(dh + dg)*cm/2.0, 0), vacRi_Logical, "vacPhys",LHRSContainerLogical,0,0,0);


/*
        vacRa_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacRb_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacRd_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacLa_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacLb_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacLd_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        //vacRe_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacRf_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacRg_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacRh_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
        vacRi_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
*/
//        vacRa_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
//        vacRb_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
//        vacLa_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
//        vacLb_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
//        vacLd_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRd_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        //vacRe_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRf_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRg_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRh_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
        vacRi_Logical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

        //at 641.8 rectangular (x-sec) material with circle tunnel  d = 60.0146
        //at 819.4887 rectangular (x-sec) 30.00 cm x 100.00 cm
        //at the Dip exit 30. cm width + 60.655 diam circle.
        //dip exit Q3 entrance (1890.3) rectangular
        //z=2065.1122   y = 721.0738  ==> titanium 0.125 mm detector exit
      }
      



        //vac

        double pQ2vacRin  = pQ2Rin;
        double pQ2vacRout = pQ2Rin + 5. * cm;
        double pQ2vacLength = ( 1166.06 - 675. ) * mm;
        pQ2vacLength = ( 370.-313.) * cm;
        double pQ2vacCenterZ = pHallCenter2LQ2Face - pQ2vacLength / 2.;
        double pQ2vacCenterY = 0;
        G4VSolid* Q2vacSolid = new G4Tubs("Q2vacTub",pQ2vacRin,pQ2vacRout,pQ2vacLength/2.0,0.0,360.0*deg);

//??delete
	G4LogicalVolume* LQ2vacLogical = new G4LogicalVolume(Q2vacSolid,
		mMaterialManager->siliconsteel,"LQ2vacLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ2vacLogical = new G4LogicalVolume(Q2vacSolid,
		mMaterialManager->siliconsteel,"RQ2vacLogical",0,0,LarmStepLimits);

	LQ2vacLogical->SetVisAttributes(G4Colour::Green()); 
	RQ2vacLogical->SetVisAttributes(G4Colour::Green()); 
//	LQ2vacLogical->SetVisAttributes(IronVisAtt); 
//	RQ2vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ2vacInContainer=new G4RotationMatrix();
	pRotQ2vacInContainer->rotateX(90*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ2vacInContainer,G4ThreeVector(0, -pQ2vacCenterZ,pQ2vacCenterY),
			LQ2vacLogical,"LQ2vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ2vacInContainer,G4ThreeVector(0, -pQ2vacCenterZ,pQ2vacCenterY),
			RQ2vacLogical,"RQ2vacPhys",RHRSContainerLogical,0,0,0);
	}
	

	//vac2
	double pQ2vac2Rin     = pQ2Rin;
	double pQ2vac2Rout    = pQ2Rin + 1. * mm;
	double pQ2vac2Length  = 561.7 * mm;
	double pQ2vac2CenterZ = pHallCenter2RQ2Face + pQ2Length + pQ2vac2Length / 2.;
	double pQ2vac2CenterY = 0;
	G4VSolid* Q2vac2Solid = new G4Tubs("Q2vac2Tub",pQ2vac2Rin,pQ2vac2Rout,pQ2vac2Length/2.0,0.0,360.0*deg);
//??delete
	G4LogicalVolume* LQ2vac2Logical = new G4LogicalVolume(Q2vac2Solid,
		mMaterialManager->siliconsteel,"LQ2vac2Logical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ2vac2Logical = new G4LogicalVolume(Q2vac2Solid,
		mMaterialManager->siliconsteel,"RQ2vac2Logical",0,0,LarmStepLimits);

	LQ2vac2Logical->SetVisAttributes(G4Colour::Green()); 
	RQ2vac2Logical->SetVisAttributes(G4Colour::Green()); 
//	LQ2vac2Logical->SetVisAttributes(IronVisAtt); 
//	RQ2vac2Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ2vac2InContainer=new G4RotationMatrix();
	pRotQ2vac2InContainer->rotateX(90*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ2vac2InContainer,G4ThreeVector(0, -pQ2vac2CenterZ,pQ2vac2CenterY),
			LQ2vac2Logical,"LQ2vac2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ2vac2InContainer,G4ThreeVector(0, -pQ2vac2CenterZ,pQ2vac2CenterY),
			RQ2vac2Logical,"RQ2vac2Phys",RHRSContainerLogical,0,0,0);
	}

	/////////////////////////
	// HRS Dipole 
	/////////////////////////
	//The dipole is built as a disc subtraced subtract another disc to get the side 
	//then subtract by the tunnel disc
	double pDipoleBendAngle=45.*deg, pDipoleFaceAngle=30.*deg;
	double pDipoleR=8.4*m;
	double pDipoleRprime=pDipoleR*sin(pDipoleBendAngle/2.)/
		sin(180.*deg-pDipoleBendAngle/2.-pDipoleFaceAngle);
	//double pDipoleRprime=4.0518*m;

	double pDipoleRCenterY=pDipoleR;
	double pDipoleRCenterZ;//=9.96*m;
	//pDipoleRCenterZ = ( IS_NIM == 1 ) ? 9.96 * m : 9.961 * m;//snake
	pDipoleRCenterZ = 9.961 * m;//snake
	//pDipoleRCenterZ = ( IS_NIM == 1 ) ? 9.96 * m : 9.9547 * m;//seamus
	G4cout << pDipoleRCenterZ << G4endl;
	//double pDipoleRprimeCenterY=pDipoleRprime*cos(pDipoleFaceAngle);			// =2.865*m;
	//double pDipoleRprimeCenterZ=9.96*m+pDipoleRprime*sin(pDipoleFaceAngle);	//12.825*m;

	//the center of Rprime relative to R 
	double pRprime2R_Y=pDipoleRprime*cos(pDipoleFaceAngle)-pDipoleR;	// =-5.535*m;
	double pRprime2R_Z=pDipoleRprime*sin(pDipoleFaceAngle);				// =1.865*m;

	//the original disc
	G4VSolid* DipoleWholeTub = new G4Tubs("DipoleWholeTub",
		pDipoleR-0.8*m,pDipoleR+0.8*m,0.4*m,
		172.*deg,pDipoleBendAngle+16.*deg);
	//the disc to be subtracted, musu be thicker 
	G4VSolid* DipolePrimeTub = new G4Tubs("DipolePrimeTub",
		0,pDipoleR,0.5*m,
		180.*deg+pDipoleFaceAngle+pDipoleBendAngle,360.*deg-pDipoleFaceAngle*2.-pDipoleBendAngle);
	//subtract the small tube to form the shape of the sides
	G4SubtractionSolid* DipoleWithSides = new G4SubtractionSolid("DipoleWithSides",
		DipoleWholeTub,DipolePrimeTub,
		0,G4ThreeVector(pRprime2R_Y,-pRprime2R_Z,0));

	//the tunnel disc, I use a rectangle shape here
	//G4VSolid* DipoleTunnelTub = new G4Tubs("DipoleTunnelTub",
	//	pDipoleR-0.4*m,pDipoleR+0.4*m,0.125*m,
	//	170*deg,pDipoleBendAngle+20*deg);
	//The shape of the dipole tunnel is a trapzoid, x=+/-0.4; y=+/-(0.125*(1-(1.25*x/8.40)) Jixie
	//To build this shape, I have to use polycone Jixie
	// -5.22008 < x < -4.98099 J.L.R.
	//  -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784 J.L.R.
	double dy=0.125*(1.25*0.40/8.40)*m;
	double dy_up=(13.43-12.5)*cm;
	double dy_dn=(12.5-11.477)*cm;
	const int kNPlane_DipoleTunnel=4;
	double rInner_DipoleTunnel[] = {pDipoleR-0.433*m,pDipoleR-0.433*m,pDipoleR-0.433*m,pDipoleR-0.433*m};
	double rOuter_DipoleTunnel[] = {pDipoleR-0.433*m,pDipoleR+0.4763*m,pDipoleR+0.4763*m,pDipoleR-0.433*m};
	double zPlane_DipoleTunnel[] = {-0.125*m-dy_up,-0.125*m+dy_dn,0.125*m-dy_dn,0.125*m+dy_up};
	G4Polycone* DipoleTunnelCone = new G4Polycone("DipoleTunnelCone",
		    170.0*deg,pDipoleBendAngle+20.*deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);//Jixie's original
	
	//subtract the tunnel disc
	G4SubtractionSolid* DipoleSolid = new G4SubtractionSolid("DipoleSolid",
		DipoleWithSides,DipoleTunnelCone);

//	double pBStepLimit=0.00001 * mm;//vardan's edit
//	G4UserLimits* uBStepLimits = new G4UserLimits(pBStepLimit);

	G4LogicalVolume* LDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"LDipoleLogical",0,0,LarmStepLimits);
//	LDipoleLogical->SetUserLimits(uBStepLimits);

	G4LogicalVolume* RDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"RDipoleLogical",0,0,LarmStepLimits);

        LDipoleLogical->SetVisAttributes(OrangeVisAtt); 
        RDipoleLogical->SetVisAttributes(OrangeVisAtt); 
//        LDipoleLogical->SetVisAttributes(G4VisAttributes::Invisible); 
//        RDipoleLogical->SetVisAttributes(G4VisAttributes::Invisible); 

	G4RotationMatrix *pRotDipoleInContainer=new G4RotationMatrix();
	G4RotationMatrix *pRotDipoleFringe1InContainer=new G4RotationMatrix();
	G4RotationMatrix *pRotDipoleFringe2InContainer=new G4RotationMatrix();
	
	pRotDipoleInContainer->rotateY(90*deg);
	//pRotDipoleFringeInContainer->rotateY(-90*deg);
	pRotDipoleFringe1InContainer->rotateY(180*deg);
	pRotDipoleFringe1InContainer->rotateZ(90*deg);
	pRotDipoleFringe1InContainer->rotateX(90*deg); 
	pRotDipoleFringe2InContainer->rotateX(90*deg);
	pRotDipoleFringe2InContainer->rotateY(90*deg + 45.*deg); 

	G4RotationMatrix *pRot_den=new G4RotationMatrix();
	G4RotationMatrix *pRot_dex=new G4RotationMatrix();	
	pRot_den->rotateY(90*deg);
	pRot_den->rotateX(-60*deg);
	pRot_dex->rotateY(90*deg);
	pRot_dex->rotateX(-90*deg + 105 * deg);

	//pRotDipoleFringe2InContainer->rotateX(90*deg); 
	//if(0)
	if(mSetupLHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    LDipoleLogical,"LDipolePhys",LHRSContainerLogical,0,0,0);
	}
	//if(0)
	if(mSetupRHRS>=4)
	{
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    RDipoleLogical,"RDipolePhys",RHRSContainerLogical,0,0,0);
	}

	


	/////////////////////////
	// HRS Q3 
	/////////////////////////
	double pQ3Rin=30.0*cm;
	double pQ3Rout=75.0*cm;
	double pQ3Length;//=180*cm;//NIM
	//pQ3Length = ( IS_NIM == 1 ) ? 180.0 * cm : 182.68 * cm; 
	pQ3Length = 182.68 * cm; 
	double pQ3CenterY;//=pDipoleR*(1-cos(pDipoleBendAngle))+2.4*m*sin(pDipoleBendAngle);//NIM
	double pQ3CenterZ;//=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.4*m*cos(pDipoleBendAngle);//NIM
	//pQ3CenterY = ( IS_NIM == 1 ) ? pDipoleR * ( 1 - cos( pDipoleBendAngle ) ) + 2.4 * m * sin( pDipoleBendAngle ) : 3.5853101 * m  + pQ3Length / sqrt(2.) / 2.;
	//pQ3CenterZ = ( IS_NIM == 1 ) ? 9.96 * m + pDipoleR * sin( pDipoleBendAngle ) + 2.4 * m * cos( pDipoleBendAngle ) : 17.0257042 * m  + pQ3Length / sqrt(2) / 2;
	pQ3CenterY =  3.58637 * m  + pQ3Length / sqrt(2.) / 2.;
	pQ3CenterZ = 17.0267042 * m  + pQ3Length / sqrt(2.) / 2;

	double pQ3LengthMag=pQ3Length * fringe_extension;
	//double pQ3LengthMag=182.68*cm;//SNAKE

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q3Solid = new G4Tubs("Q3Tub",pQ3Rin,pQ3Rout,pQ3Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"LQ3Logical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"RQ3Logical",0,0,LarmStepLimits);


//	LQ3Logical->SetVisAttributes(G4VisAttributes::Invisible); 
	LQ3Logical->SetVisAttributes(IronVisAtt); 
	RQ3Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3InContainer=new G4RotationMatrix();
	pRotQ3InContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			LQ3Logical,"LQ3Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			RQ3Logical,"RQ3Phys",RHRSContainerLogical,0,0,0);
	}


	//vacuum pipe
      double pQ3vacRin  = 30.0*cm;
//        double pQ3vacRin  = 35.3*cm;  //Vardan
	double pQ3vacRout = 45.0*cm;
        double pQ3vacLength = 566. * mm;   // Vardan
//old   double pQ3vacLength = 575. * mm; 
	double pQ3vacCenterY = pQ3CenterY + pQ3Length / 2. / sqrt( 2 ) + pQ3vacLength / 2. / sqrt( 2 );
	double pQ3vacCenterZ = pQ3CenterZ + pQ3Length / 2. / sqrt( 2 ) + pQ3vacLength / 2. / sqrt( 2 );
	G4VSolid* Q3vacSolid = new G4Tubs("Q3vacTub",pQ3vacRin,pQ3vacRout,pQ3vacLength/2.0,0.0,360.0*deg);
//??delete
	G4LogicalVolume* LQ3vacLogical = new G4LogicalVolume(Q3vacSolid,
		mMaterialManager->siliconsteel,"LQ3vacLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ3vacLogical = new G4LogicalVolume(Q3vacSolid,
		mMaterialManager->siliconsteel,"RQ3vacLogical",0,0,LarmStepLimits);

//	LQ3vacLogical->SetVisAttributes(IronVisAtt); 
//	RQ3vacLogical->SetVisAttributes(IronVisAtt); 

	LQ3vacLogical->SetVisAttributes(G4Colour::Green()); 
	RQ3vacLogical->SetVisAttributes(G4Colour::Green()); 

	G4RotationMatrix *pRotQ3vacInContainer=new G4RotationMatrix();
	pRotQ3vacInContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3vacInContainer,G4ThreeVector(0,-pQ3vacCenterZ,pQ3vacCenterY),
			LQ3vacLogical,"LQ3vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3vacInContainer,G4ThreeVector(0,-pQ3vacCenterZ,pQ3vacCenterY),
			RQ3vacLogical,"RQ3vacPhys",RHRSContainerLogical,0,0,0);
	}
	//vacuum pipe
        double pQ3vacRin1  = 30.3276*cm;  //Vardan
	double pQ3vacRout1 = 45.0*cm;
        double pQ3vacLength1 = 566. * mm;   // Vardan
	double pQ3vac2CenterY = pQ3CenterY - pQ3Length / 2. / sqrt( 2 ) - pQ3vacLength1 / 2. / sqrt( 2 );
	double pQ3vac2CenterZ = pQ3CenterZ - pQ3Length / 2. / sqrt( 2 ) - pQ3vacLength1 / 2. / sqrt( 2 );
	G4VSolid* Q3vac2Solid = new G4Tubs("Q3vac2Tub",pQ3vacRin1,pQ3vacRout1,pQ3vacLength1/2.0,0.0,360.0*deg);
//??delete
	G4LogicalVolume* LQ3vac2Logical = new G4LogicalVolume(Q3vac2Solid,
		mMaterialManager->siliconsteel,"LQ3vac2Logical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ3vac2Logical = new G4LogicalVolume(Q3vac2Solid,
		mMaterialManager->siliconsteel,"RQ3vac2Logical",0,0,LarmStepLimits);

	LQ3vac2Logical->SetVisAttributes(IronVisAtt); 
	RQ3vac2Logical->SetVisAttributes(IronVisAtt); 
	LQ3vac2Logical->SetVisAttributes(G4Colour::Red()); 
	RQ3vac2Logical->SetVisAttributes(G4Colour::Green()); 

	G4RotationMatrix *pRotQ3vac2InContainer=new G4RotationMatrix();
	pRotQ3vac2InContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3vac2InContainer,G4ThreeVector(0,-pQ3vac2CenterZ,pQ3vac2CenterY),
			LQ3vac2Logical,"LQ3vac2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3vac2InContainer,G4ThreeVector(0,-pQ3vac2CenterZ,pQ3vac2CenterY),
			RQ3vac2Logical,"RQ3vac2Phys",RHRSContainerLogical,0,0,0);
	}



	//////////////////////////////////////////////////////////
	//#ifdef G4DEBUG_GEOMETRY
	//plot the tunnel to view
	//////////////////////////
	//Q1Q2 tunnel
	//////////////////////////

	const int kNPlane_Q1Q2Tunnel=4;
	double rInner_Q1Q2Tunnel[] = {0,0,0,0};
	double rOuter_Q1Q2Tunnel[] = {pQ1Rin-0.1*mm,pQ1Rin-0.1*mm,pQ2Rin-0.1*mm,pQ2Rin-0.1*mm};
	double zPlane_Q1Q2Tunnel[] = {pHRSContainerRin+6*mm,pHallCenter2RQ1Face+pQ1PosAct+10.0*cm,
		pHallCenter2RQ1Face+pQ1PosAct+30*cm,9.7*m};
	G4Polycone* Q1Q2TunnelSolid = new G4Polycone("Q1Q2TunnelPolycone",0.0,360.0*deg,
		kNPlane_Q1Q2Tunnel,zPlane_Q1Q2Tunnel,rInner_Q1Q2Tunnel,rOuter_Q1Q2Tunnel);

	G4LogicalVolume* LQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"LQ1Q2TunnelLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"RQ1Q2TunnelLogical",0,0,LarmStepLimits);


//	LQ1Q2TunnelLogical->SetVisAttributes(LightYellowVisAtt); 
	LQ1Q2TunnelLogical->SetVisAttributes(G4VisAttributes::Invisible); 
	RQ1Q2TunnelLogical->SetVisAttributes(LightYellowVisAtt); 
	/*
	if(mSetupLHRS>=3){//This is the old Jixie way
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			    LQ1Q2TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3){//This is the old Jixie way
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			    RQ1Q2TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}
	*/
 	//This is a copy of Q1, but of smaller radius.//This is the new Nickie way
	//G4VSolid* Q1SolidMag = new G4Tubs("Q1TubMag",0,pQ1Rin,pQ1Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q1SolidMag = new G4Tubs("Q1TubMag",0,pQ1Rin,pQ1LengthMag/2.0,0.0,360.0*deg);//for fringe fields, note longer length

	//build 2 copy since there are different fields involved in it
	G4LogicalVolume* LQ1MagLogical = new G4LogicalVolume(Q1SolidMag,
		mMaterialManager->vacuum,"LQ1MagLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ1MagLogical = new G4LogicalVolume(Q1SolidMag,
		mMaterialManager->vacuum,"RQ1MagLogical",0,0,LarmStepLimits);

	//Nickie's Q1 field
	G4FieldManager* LQ1FieldManager = mEMFieldSetup->GetFieldManagerFZBL1();
	G4FieldManager* RQ1FieldManager = mEMFieldSetup->GetFieldManagerFZBR1();
	G4bool allLocal = true;
	//G4bool allLocal = false;
	LQ1MagLogical->SetFieldManager(LQ1FieldManager,allLocal);
	RQ1MagLogical->SetFieldManager(RQ1FieldManager,allLocal);

	LQ1MagLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ1MagLogical->SetVisAttributes(LightYellowVisAtt); 


	//Nickie's old way
	if(mSetupLHRS>=2){
	  //put it in the container, which also center at the hall center
	  //therefore only the z_at_lab position need to be considered
	  double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1PosAct/2.0 + q1shift);//NIM
	  //double pLQ1Pos_Z=(pQ1PosAct + q1shift);//NIM
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			    LQ1MagLogical,"LQ1MagPhys",LHRSContainerLogical,0,0,0);
	}if(mSetupRHRS>=2){
	  double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1PosAct/2.0 + q1shift);//NIM
	  //double pRQ1Pos_Z=(pQ1PosAct + q1shift);//NIM
	  //double pRQ1Pos_Z=(pHallCenter2RQ1FaceMag+pQ1LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			    RQ1MagLogical,"RQ1MagPhys",RHRSContainerLogical,0,0,0);
	}
	
	//G4VSolid* Q2SolidMag = new G4Tubs("Q2TubMag",0,pQ2Rin,pQ2Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q2SolidMag = new G4Tubs("Q2TubMag",0,pQ2Rin,pQ2LengthMag/2.0,0.0,360.0*deg);//SNAKE

	G4LogicalVolume* LQ2MagLogical = new G4LogicalVolume(Q2SolidMag,
		mMaterialManager->vacuum,"LQ2MagLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ2MagLogical = new G4LogicalVolume(Q2SolidMag,
		mMaterialManager->vacuum,"RQ2MagLogical",0,0,LarmStepLimits);

	//Nickie's Q2 field
	G4FieldManager* LQ2FieldManager = mEMFieldSetup->GetFieldManagerFZBL2();
	G4FieldManager* RQ2FieldManager = mEMFieldSetup->GetFieldManagerFZBR2();
	
	LQ2MagLogical->SetFieldManager(LQ2FieldManager,allLocal);
	RQ2MagLogical->SetFieldManager(RQ2FieldManager,allLocal);

	LQ2MagLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ2MagLogical->SetVisAttributes(LightYellowVisAtt); 

	 //Nickie's old way
	if(mSetupLHRS>=3){
	  //put it in the container, which also center at the hall center
	  double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);//NIM
	  //double pLQ2Pos_Z=(pHallCenter2LQ2FaceMag+pQ2LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			    LQ2MagLogical,"LQ2MagPhys",LHRSContainerLogical,0,0,0);
	}if(mSetupRHRS>=3){
	  double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);//NIM
	  //double pRQ2Pos_Z=(pHallCenter2RQ2FaceMag+pQ2LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			    RQ2MagLogical,"RQ2MagPhys",RHRSContainerLogical,0,0,0);
	}
	

	//////////////////////////
	//Dipole tunnel
	//////////////////////////
	//G4VSolid* DipoleTunnelSolid = new G4Tubs("DipoleTunnelSolid",pDipoleR-0.4*m+0.1*mm,
	//	pDipoleR+0.4*m-0.1*mm,0.125*m-0.1*mm,180*deg,pDipoleBendAngle);
	double pLQ1Pos_Z_en2=(pHallCenter2LQ1Face);//NIM
	double pLQ1Pos_Z_ex2=(pHallCenter2LQ1Face + pQ1PosAct);//NIM
	double pLQ2Pos_Z_en2=(pHallCenter2LQ2Face);//NIM
	double pLQ2Pos_Z_ex2=(pHallCenter2LQ2Face + pQ2Length);//NIM
	double pRQ1Pos_Z_en2=(pHallCenter2RQ1Face);//NIM
	double pRQ1Pos_Z_ex2=(pHallCenter2RQ1Face + pQ1PosAct);//NIM
	double pRQ2Pos_Z_en2=(pHallCenter2RQ2Face);//NIM
	double pRQ2Pos_Z_ex2=(pHallCenter2RQ2Face + pQ2Length);//NIM
	double pLDPos_Z_en2 = pLQ2Pos_Z_ex2 + 4.42 * m;
	double pRDPos_Z_en2 = pRQ2Pos_Z_ex2 + 4.42 * m;

	double rInner_DipoleTunnel2[] = {0.,0.,0.,0.};
	double rOuter_DipoleTunnel2[] = {pDipoleR-0.433*m,pDipoleR+0.4763*m,pDipoleR+0.4763*m,pDipoleR-0.433*m};

	G4Polycone* DipoleTunnelCone1 = new G4Polycone("DipoleTunnelCone1",
	            150.0*deg,105.0*deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel2,rOuter_DipoleTunnel2);
	G4Polycone* DipoleTunnelCone2 = new G4Polycone("DipoleTunnelCone2",
	            175.0*deg,pDipoleBendAngle + 10 * deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4Polycone* DipoleTunnelSolid = new G4Polycone("DipoleTunnelSolid",
		180.0*deg,pDipoleBendAngle,
		kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4double psiadj = 180. * deg - 22.5 * deg - 30. * deg;
	G4double badj = 8.4 * m * sin( 30.0 * deg ) / sin( psiadj );
	G4double xadj = badj * sin( 22.5 * deg );
	G4double yadj = badj * cos( 22.5 * deg );
	G4IntersectionSolid* DipoleTunnelCone3 = new G4IntersectionSolid("DipoleTunnelCone3",
			 DipoleTunnelCone2, DipoleTunnelCone1, 0, G4ThreeVector(-yadj, -xadj, 0));
	

	//DipoleTunnelSolid is Jixie with flat entrance
	//DipoleTunnelCone3 is Nickie with 30 deg entrance

	//But that is not good enough. We also have to add extensions to have fringe fields.
	//This is just going to be a trapezoidal extensionpRotDipoleInContaine
	//G4Trd* DipoleFringeSolid = new G4Trd("DipoleFringeSolid", 0.125*m+dy, 0.125*m-dy, 1.4 * m + 0.4 * m * 2. * tan( pi / 6 ), 1.4 * m, 0.4 * m);
	//G4cout << 0.125*m+dy << " and " << 0.125*m-dy << " better be greater than zero!" << G4endl;
	G4Trap* DipoleFringeSolid1 = new G4Trap("DipoleFringeSolid1",
					       0.4763 * m , -pi / 12. ,
					       0.0     , 0.125*m-dy_dn,
					       2.05 * m + 0.4763 * m * 2. * tan( pi / 12. ), 2.05 * m + 0.4763 * m * 2. * tan( pi / 12. ),
					       0.0     , 0.125*m+dy_dn,
					       2.05 * m, 2.05 * m  ,
					       0.0   );
	G4Trap* DipoleFringeSolid2 = new G4Trap("DipoleFringeSolid2",
					       0.4763 * m , -pi / 12. ,
					       0.0     , 0.125*m-dy_dn,
					       0.65 * m + 0.4763 * m * 2. * tan( pi / 12. ), 0.65 * m + 0.4763 * m * 2. * tan( pi / 12. ),
					       0.0     , 0.125*m+dy_dn,
					       0.65 * m, 0.65 * m  ,
					       0.0    );
	//pDz  , pTheta,
	//pPhi , pDy1,
	//pDx1 , pDx2,
	//pAlp1, pDy2,
	//pDx3 , pDx4,
	//pAlp2,

	//G4UnionSolid* DipoleTunnelCone4 = new G4UnionSolid("DipoleTunnelCone4",
	//DipoleTunnelCone3, DipoleFringeSolid, pRotDipoleInContainer, G4ThreeVector(-8.4 * m, 0., 0.));
	G4UnionSolid* DipoleTunnelCone4 = new G4UnionSolid("DipoleTunnelCone4",
		   DipoleTunnelCone3, DipoleFringeSolid1, pRotDipoleFringe1InContainer, G4ThreeVector(-8.4 * m, +2.05 * m, 0.));
	G4UnionSolid* DipoleTunnelCone5 = new G4UnionSolid("DipoleTunnelCone4",
		   DipoleTunnelCone4, DipoleFringeSolid2, pRotDipoleFringe2InContainer,
							   G4ThreeVector(-8.4 * m * cos( pi / 4. ) + 0.3 * m / sqrt( 2. ),
									 -8.4 * m * sin( pi / 4. ) - 0.3 * m / sqrt( 2. ), 0.));

	G4LogicalVolume* LDipoleFringe1Logical = new G4LogicalVolume(DipoleFringeSolid1,
	        mMaterialManager->vacuum,"LDipoleFringe1Logical",0,0,LarmStepLimits);
	G4LogicalVolume* RDipoleFringe1Logical = new G4LogicalVolume(DipoleFringeSolid1,
	        mMaterialManager->vacuum,"RDipoleFringe1Logical",0,0,LarmStepLimits);
	G4LogicalVolume* LDipoleFringe2Logical = new G4LogicalVolume(DipoleFringeSolid2,
	        mMaterialManager->vacuum,"LDipoleFringe2Logical",0,0,LarmStepLimits);
	G4LogicalVolume* RDipoleFringe2Logical = new G4LogicalVolume(DipoleFringeSolid2,
	        mMaterialManager->vacuum,"RDipoleFringe2Logical",0,0,LarmStepLimits);

	G4LogicalVolume* LDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone3,
		mMaterialManager->vacuum,"LDipoleTunnelLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone3,
		mMaterialManager->vacuum,"RDipoleTunnelLogical",0,0,LarmStepLimits);


        LDipoleTunnelLogical->SetVisAttributes(G4VisAttributes::Invisible);//G4Colour::Green()); 
        RDipoleTunnelLogical->SetVisAttributes(G4VisAttributes::Invisible);//YellowVisAtt); 
        LDipoleFringe1Logical->SetVisAttributes(G4Colour::Red()); // was invisible
        RDipoleFringe1Logical->SetVisAttributes(YellowVisAtt); 
        LDipoleFringe2Logical->SetVisAttributes(G4Colour::Red()); // was invisible
        RDipoleFringe2Logical->SetVisAttributes(YellowVisAtt); 

//        LDipoleTunnelLogical->SetVisAttributes(G4VisAttributes::Invisible); 
//        RDipoleTunnelLogical->SetVisAttributes(G4VisAttributes::Invisible); 
//        LDipoleFringe1Logical->SetVisAttributes(G4VisAttributes::Invisible); // was invisible
//        RDipoleFringe1Logical->SetVisAttributes(G4VisAttributes::Invisible); 
//        LDipoleFringe2Logical->SetVisAttributes(G4VisAttributes::Invisible); // was invisible
//        RDipoleFringe2Logical->SetVisAttributes(G4VisAttributes::Invisible); 

	//Nickie's dipole field
	G4FieldManager* LdipoleFieldManager = mEMFieldSetup->GetFieldManagerFZBL3();
	G4FieldManager* RdipoleFieldManager = mEMFieldSetup->GetFieldManagerFZBR3();
	//G4double minEps= 1.0e-9;  //   Minimum & value for smallest steps
	//G4double maxEps= 1.0e-8;  //   Maximum & value for largest steps

	//dipoleFieldManager->SetMinimumEpsilonStep( minEps );
	//dipoleFieldManager->SetMaximumEpsilonStep( maxEps );
	LdipoleFieldManager->SetDeltaOneStep( 0.001 );
	LdipoleFieldManager->GetChordFinder()->SetDeltaChord( 0.01 * mm );
	LDipoleTunnelLogical->SetFieldManager(LdipoleFieldManager,allLocal);
	RDipoleTunnelLogical->SetFieldManager(RdipoleFieldManager,allLocal);

        //double pFringeX_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) + 0.75 * m / sqrt( 2 ) :   2.4603032 * m + 0.75 * m / sqrt( 2 );
        //double pFringeZ_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) - 0.75 * m / sqrt( 2 ) : -15.9006973 * m - 0.75 * m / sqrt( 2 );

	if(mSetupLHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    LDipoleTunnelLogical,"LDipoleTunnelPhys",LHRSContainerLogical,0,0,0);
	  /*
	  new G4PVPlacement(pRotDipoleFringe1InContainer,
			    G4ThreeVector(0., -7.75 * m, 0.),
			    LDipoleFringe1Logical,"LDipoleFringe1Phys",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotDipoleFringe2InContainer,
			    G4ThreeVector(0., pFringeZ_ex, pFringeX_ex),
			    LDipoleFringe2Logical,"LDipoleFringe2Phys",LHRSContainerLogical,0,0,0);
	  */
	}
	if(mSetupRHRS>=4)
	{
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    RDipoleTunnelLogical,"RDipoleTunnelPhys",RHRSContainerLogical,0,0,0);
	  /*
	  new G4PVPlacement(pRotDipoleFringe1InContainer,
			    G4ThreeVector(0., -7.75 * m, 0.),
			    RDipoleFringe1Logical,"RDipoleFringe1Phys",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotDipoleFringe2InContainer,
			    G4ThreeVector(0., pFringeZ_ex, pFringeX_ex),
			    RDipoleFringe2Logical,"RDipoleFringe2Phys",RHRSContainerLogical,0,0,0);

	  */
	}


	//////////////////////////
	//Q3 tunnel
	//////////////////////////

	G4VSolid* Q3TunnelSolid = new G4Tubs("Q3TunnelTub",0,pQ3Rin-0.1*mm,
		pQ3Length,0.0,360.0*deg);
	G4LogicalVolume* LQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"LQ3TunnelLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"RQ3TunnelLogical",0,0,LarmStepLimits);

//	LQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 
	LQ3TunnelLogical->SetVisAttributes(G4VisAttributes::Invisible); 
	RQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 

	//double pQ3TunnelPos_Y=pDipoleR*(1-cos(pDipoleBendAngle))+2.0*m*sin(pDipoleBendAngle);
	//double pQ3TunnelPos_Z=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	//double pQ3TunnelPos_Z=pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	//pQ3TunnelPos_Z += ( IS_NIM == 1 ) ? 9.96*m : pDipoleRCenterZ; 
	//if(mSetupLHRS>=4)
	//{
		//put it in the container, which also center at the hall center
		//new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
				  //LQ3TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	//}
	//if(mSetupRHRS>=4)
	//{
	  //new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
	  //RQ3TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	//}

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	//G4VSolid* Q3MagSolid = new G4Tubs("Q3MagTub",0,pQ3Rin,pQ3Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q3MagSolid = new G4Tubs("Q3MagTub",0,pQ3Rin,pQ3LengthMag/2.0,0.0,360.0*deg);//SNAKE

	G4LogicalVolume* LQ3MagLogical = new G4LogicalVolume(Q3MagSolid,
		mMaterialManager->vacuum,"LQ3MagLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RQ3MagLogical = new G4LogicalVolume(Q3MagSolid,
		mMaterialManager->vacuum,"RQ3MagLogical",0,0,LarmStepLimits);

	//Nickie's Q3 field
	G4FieldManager* LQ3FieldManager = mEMFieldSetup->GetFieldManagerFZBL4();
	G4FieldManager* RQ3FieldManager = mEMFieldSetup->GetFieldManagerFZBR4();
	LQ3MagLogical->SetFieldManager(LQ3FieldManager,allLocal);
	RQ3MagLogical->SetFieldManager(RQ3FieldManager,allLocal);

	LQ3MagLogical->SetVisAttributes(MagFieldVisAtt); 
	RQ3MagLogical->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotQ3MagInContainer=new G4RotationMatrix();
	pRotQ3MagInContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  //new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZMag,pQ3CenterYMag),//SNAKE
	  //LQ3MagLogical,"LQ3MagPhys",LHRSContainerLogical,0,0,0);//SNAKE
	  new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),//NIM
			    LQ3MagLogical,"LQ3MagPhys",LHRSContainerLogical,0,0,0);//NIM
	}
	if(mSetupRHRS>=4)
	{
	  //new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZMag,pQ3CenterYMag),//SNAKE
	  //RQ3MagLogical,"RQ3MagPhys",RHRSContainerLogical,0,0,0);//SNAKE
	  new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),//NIM
	    		    RQ3MagLogical,"RQ3MagPhys",RHRSContainerLogical,0,0,0);//NIM
	}


	//Nickie's sensitive detector, at the focal plane, virtual boundary, vb, fp
	//double pFPR      = pDipoleR; //radius of curvature of dipole
	//double pFPA      = 9.96 * m; //distance from pivot to dipole entrance
	//double pFPB      = 2.4  * m; //distance from center of dipole to center of Q3
	//
	//double pFPB      = 2.4  * m + pQ3Length / 2.0 + 3.57 * m + 1.43 * m;
	//double pFPB      = 2.4  * m + pQ3Length / 2.0 + 3.57 * m + 0.7 * m;
	//double pFPCenterY=       pFPR * ( 1 - cos( pDipoleBendAngle ) ) + pFPB * sin( pDipoleBendAngle );
	//double pFPCenterZ=pFPA + pFPR *       sin( pDipoleBendAngle )   + pFPB * cos( pDipoleBendAngle );
        //double pFPA      = 9.96 * m;//distance from pivot to entrance of dipole
        //double pFPCenterX=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) ) *-sin( mRHRSAngle );
        //double pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) ) * cos( mRHRSAngle );

	double pFPCenterZ , pFPCenterY ;
	double vdc1Z , vdc1Y ;
	double pVDCCenterZ, pVDCCenterY;
	double pQZ1CenterZ, pQZ1CenterY;
	double pQZ2CenterZ, pQZ2CenterY;
	//if( IS_NIM == 1 ){
	//double pFPR      = 8.4  * m;//radius of curvature of dipole
	//double pFPA      = ( IS_NIM == 1 ) ? 9.96 * m : pDipoleRCenterZ;//distance from pivot to entrance of dipole
	//double pFPH      = pFPR * tan ( 45. / 2. * deg );
	//pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) );
	//pFPCenterY=(                 pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2);
	//}else{
//	vdc1Z=pQ3CenterZ + ( pQ3Length / 2. + 3.57 * m ) / sqrt(2);
//	vdc1Y=pQ3CenterY + ( pQ3Length / 2. + 3.57 * m ) / sqrt(2);
	vdc1Z=pQ3CenterZ + ( pQ3Length / 2. + 3.4523 * m ) / sqrt(2);
	vdc1Y=pQ3CenterY + ( pQ3Length / 2. + 3.4523 * m ) / sqrt(2);
	cout<<"vdc1Z="<<vdc1Z<<"    vdc1Y="<<vdc1Y<<endl;

//SIMC
//            Q1 center position	= 207.064 cm
//            Q2 center position	= 462.064 cm
//            Q3 center position	= 421.038 cm before vacuum window.

//            Q1 mechanical length 	= 165.40  cm
//            Q2/Q3 mechanical length	= 295.20  cm

//            Q1 optical length	=  94.13  cm
//            Q2/Q3 optical length	= 182.66  cm	!COSY Q3 is actually 182.68.
//
//            Focal plane		= 436.56  cm after Q3 optical exit.
//            Vacuum window		= 273.44  cm after Q3 mechanical exit.
//
//            So taking focal plane = 2342.99446 from COSY, we work backwards:
//
//            Q3 center	= 2342.99446 - 436.56	= 1906.43446
//            Q3 optical exit	= 1906.43446 + 182.66/2	= 1997.76446
//            Q3 mech. exit	= 1906.43446 + 295.20/2	= 2054.03446
//            Vacuum window	= 1906.43446 + 421.038	= 2327.47246 = 15.522 cm before FP.
//            Rectangle before
//            vacuum window	= 2327.47246 - 247.085	= 2080.38746


	pVDCCenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pVDCCenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pQZ1CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ1CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ2CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pQZ2CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pFPCenterZ  = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	pFPCenterY  = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	//}
	
	G4double vb_thickness = .01 * mm;
	G4VSolid* FPSolid = new G4Tubs("FPTub",0,pQ3Rout * 2,vb_thickness,0.0,360.0*deg);
	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,pQ1Rout,vb_thickness,0.0,360.0*deg); //circles
	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,pQ2Rout,vb_thickness,0.0,360.0*deg); //circles

	G4VSolid* vdc1PlaneSolid = new G4Box("vdc1Plane",3*288*mm,2118*mm,0.01*mm); //box

	double X_S0 = 29.5 * cm;
	double Y_S0 = 35.5 * cm;
	double Z_S0 = vb_thickness;

	double X_S1 = 29.5 * cm;
	double Y_S1 = 35.5 * cm;
	double Z_S1 = vb_thickness;

	double X_Q1 = 29.5 * cm;
	double Y_Q1 = 35.5 * cm;
	double Z_Q1 = vb_thickness;

	double X_Q2 = 29.5 * cm;
	double Y_Q2 = 35.5 * cm;
	double Z_Q2 = vb_thickness;

	G4VSolid* Scint0  = new G4Box("Scint0" , X_S0 / 2.0, Y_S0 / 2.0, Z_S0 / 2.0);
	G4VSolid* Scint1  = new G4Box("Scint1" , X_S1 / 2.0, Y_S1 / 2.0, Z_S1 / 2.0);
	G4VSolid* Quartz1 = new G4Box("Quartz1", X_Q1 / 2.0, Y_Q1 / 2.0, Z_Q1 / 2.0);
	G4VSolid* Quartz2 = new G4Box("Quartz2", X_Q2 / 2.0, Y_Q2 / 2.0, Z_Q2 / 2.0);

	G4VSolid* magneticSolid = new G4Box("magneticBox",mFieldX/2.0,mFieldY/2.0,mFieldZ/2.0);

	G4LogicalVolume* vdc1LLogical = new G4LogicalVolume(vdc1PlaneSolid,
		mMaterialManager->vacuum,"vdc1Logical",0,0,LarmStepLimits);
	G4LogicalVolume* vdc1RLogical = new G4LogicalVolume(vdc1PlaneSolid,
		mMaterialManager->vacuum,"vdc1Logical",0,0,LarmStepLimits);


	G4LogicalVolume* FocPlLLogical = new G4LogicalVolume(vdc1PlaneSolid,
		mMaterialManager->vacuum,"FocPlLLogical",0,0,LarmStepLimits);


	G4LogicalVolume* LFPLogical = new G4LogicalVolume(FPSolid,
		mMaterialManager->vacuum,"LFPLogical",0,0,LarmStepLimits);
	G4LogicalVolume* RFPLogical = new G4LogicalVolume(FPSolid,
		mMaterialManager->vacuum,"RFPLogical",0,0,LarmStepLimits);

	G4LogicalVolume* LPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		mMaterialManager->vacuum,"LPlaneLogical1",0,0,LarmStepLimits);
	G4LogicalVolume* RPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		mMaterialManager->vacuum,"RPlaneLogical1",0,0,LarmStepLimits);
	G4LogicalVolume* LPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		mMaterialManager->vacuum,"LPlaneLogical1",0,0,LarmStepLimits);
	G4LogicalVolume* RPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		mMaterialManager->vacuum,"RPlaneLogical1",0,0,LarmStepLimits);

	vdc1LLogical->SetVisAttributes(IronVisAtt); 
	vdc1RLogical->SetVisAttributes(IronVisAtt); 
//	FocPlLLogical->SetVisAttributes(G4Colour::Red()); 
	FocPlLLogical->SetVisAttributes(G4VisAttributes::Invisible); 
	LFPLogical->SetVisAttributes(G4VisAttributes::Invisible);
	RFPLogical->SetVisAttributes(G4VisAttributes::Invisible); 
	LPlaneLogical1->SetVisAttributes(G4VisAttributes::Invisible); 
	RPlaneLogical1->SetVisAttributes(G4VisAttributes::Invisible); 
	LPlaneLogical2->SetVisAttributes(G4VisAttributes::Invisible); 
	RPlaneLogical2->SetVisAttributes(G4VisAttributes::Invisible); 

/*
	RFPLogical->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 
*/
	G4RotationMatrix *pRotVDCInContainer=new G4RotationMatrix();
	pRotVDCInContainer->rotateX(0.*deg); 
	G4RotationMatrix *pRotFPInContainer=new G4RotationMatrix();
	pRotFPInContainer->rotateX(-45*deg); 

//	G4RotationMatrix *pRotvdc1=new G4RotationMatrix();
//	pRotbvdc1->rotateX(-45*deg); 

	G4RotationMatrix *pRotFocPlL=new G4RotationMatrix();
	pRotFocPlL->rotateX(fp_angle);
	cout<<"focal plane angle is "<<fp_angle/deg<<" deg"<<endl;
/*
	pRotFocPlL->rotateX(-50.143);   //p=0.72
	pRotFocPlL->rotateX(-49.2431);  //p=0.725
	pRotFocPlL->rotateX(-48.3605);  //p=0.73
	pRotFocPlL->rotateX(-47.4966);  //p=0.735
	pRotFocPlL->rotateX(-46.6518);  //p=0.74
	pRotFocPlL->rotateX(-45.8258);  //p=0.745
	pRotFocPlL->rotateX(-45.0201);  //p=0.75
	pRotFocPlL->rotateX(-44.234);   //p=0.755
	pRotFocPlL->rotateX(-43.4677);  //p=0.76
	pRotFocPlL->rotateX(-42.7212);  //p=0.765
	pRotFocPlL->rotateX(-41.9949);  //p=0.77
	pRotFocPlL->rotateX(-41.2885);  //p=0.775
	pRotFocPlL->rotateX(-40.6025);  //p=0.78
*/

	//if(mSetupLHRS>=4){
	if(mSnakeModel == 49 || mSnakeModel == 48 || mSnakeModel > 51 ){
	  double pSeptumX      = 140.0  * cm;
	  double pSeptumY      = 84.4   * cm;
	  double pSeptumZ      = 74.0   * cm;
	  //double pSeptumPlaceZ = 70.414 * cm;
	  double pSeptumPlaceZ = 69.99937 * cm;
	  new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ - 0.5 * pSeptumZ + 2 * vb_thickness),
			    LPlaneLogical2,"virtualBoundaryPhys_sen",motherLogical,0,0,0);//sen
	  new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ),
			    LPlaneLogical2,"virtualBoundaryPhys_sm",motherLogical,0,0,0);//sm
	  new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ + 0.5 * pSeptumZ),
			    LPlaneLogical2,"virtualBoundaryPhys_sex",motherLogical,0,0,0);//sex
	  new G4PVPlacement(0,G4ThreeVector(0,0,36. * cm),
			    LPlaneLogical2,"virtualBoundaryPhys_coil",motherLogical,0,0,0);//coil
	  new G4PVPlacement(0,G4ThreeVector(0,0,-50. * cm),
			    LPlaneLogical2,"virtualBoundaryPhys_mid",motherLogical,0,0,0);//mid
	  
	}
	//double pLQ1Pos_Z_en=(pHallCenter2LQ1Face);//NIM
	//double pLQ1Pos_Z_ex=(pHallCenter2LQ1Face + pQ1Length);//NIM
	//double pLQ2Pos_Z_en=(pHallCenter2LQ2Face);//NIM
	//double pLQ2Pos_Z_ex=(pHallCenter2LQ2Face + pQ2Length);//NIM
	//double pRQ1Pos_Z_en=(pHallCenter2RQ1Face);//NIM
	//double pRQ1Pos_Z_ex=(pHallCenter2RQ1Face + pQ1Length);//NIM
	//double pRQ2Pos_Z_en=(pHallCenter2RQ2Face);//NIM
	//double pRQ2Pos_Z_ex=(pHallCenter2RQ2Face + pQ2Length);//NIM
	//double pLDPos_Z_en = pLQ2Pos_Z_ex + 4.42 * m;
	//double pRDPos_Z_en = pRQ2Pos_Z_ex + 4.42 * m;
	//double pLDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
	//double pRDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
	//double pLDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;
	//double pRDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;

	double pHallCenter2Col = 1.38 * m;
	double pPaulColT        = 0.01 * m;
	double pPaulX = ( - pHallCenter2Col - pPaulColT * 2. ) * cos(mLHRSAngle) ;
	double pPaulY = ( - pHallCenter2Col - pPaulColT * 2. ) * sin(mLHRSAngle);
	if(mSnakeModel == 49 || mSnakeModel == 48 || mSnakeModel > 50 ){
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE
	  new G4PVPlacement(pRotXLHRSdeg,G4ThreeVector(pPaulY,0,-pPaulX),
			    LPlaneLogical1,"virtualBoundaryPhys_col_LHRS",motherLogical,0,0,0);//q1en
	  /*
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_en,0),
			    LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LHRSContainerLogical,0,0,0);//q1en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_ex,0),
			    LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LHRSContainerLogical,0,0,0);//q1ex
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_en,0),
			    LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LHRSContainerLogical,0,0,0);//q2en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_ex,0),
			    LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LHRSContainerLogical,0,0,0);//q2ex
	  */
	  new G4PVPlacement(0,G4ThreeVector(0, 0,-pQ1Length/2.),
			    LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LQ1MagLogical,0,0,0);//q1en
	  new G4PVPlacement(0,G4ThreeVector(0, 0, pQ1Length/2.),
			    LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LQ1MagLogical,0,0,0);//q1ex
	  new G4PVPlacement(0,G4ThreeVector(0, 0,-pQ2Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LQ2MagLogical,0,0,0);//q2en
	  new G4PVPlacement(0,G4ThreeVector(0, 0, pQ2Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LQ2MagLogical,0,0,0);//q2ex
	  //new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pLDPos_Z_en,0),
	  //LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LHRSContainerLogical,0,0,0);//den
	  //new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pLDPos_Z_ex,pLDPos_X_ex),
	  //LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LHRSContainerLogical,0,0,0);//dex
	  
	  new G4PVPlacement(pRot_den,G4ThreeVector(-8.4*m,0,0),
			    LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LDipoleTunnelLogical,0,0,0);//den
	  new G4PVPlacement(pRot_dex,G4ThreeVector(-8.4 * m * cos( pi / 4. ), -8.4 * m * sin( pi / 4. ),0),
			    LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LDipoleTunnelLogical,0,0,0);//dex
	  
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ + pQ3Length / sqrt(2) / 2,pQ3CenterY - pQ3Length / sqrt(2) / 2),
	  //LPlaneLogical2,"virtualBoundaryPhys_q3en_LHRS",LHRSContainerLogical,0,0,0);//q3en
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ - pQ3Length / sqrt(2) / 2,pQ3CenterY + pQ3Length / sqrt(2) / 2),
	  //LPlaneLogical2,"virtualBoundaryPhys_q3ex_LHRS",LHRSContainerLogical,0,0,0);//q3ex
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ3Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q3en_LHRS",LQ3MagLogical,0,0,0);//q3en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ3Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q3ex_LHRS",LQ3MagLogical,0,0,0);//q3ex
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pVDCCenterZ,pVDCCenterY),
			    LFPLogical,"virtualBoundaryPhys_vdc_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pFPCenterZ,pFPCenterY),
			    LFPLogical,"virtualBoundaryPhys_fp_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pQZ1CenterZ,pQZ1CenterY),
			    LFPLogical,"virtualBoundaryPhys_qz1_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pQZ2CenterZ,pQZ2CenterY),
			    LFPLogical,"virtualBoundaryPhys_qz2_LHRS",LHRSContainerLogical,0,0,0);

	}
	//if(mSetupRHRS>=4){
	if(mSnakeModel == 49 || mSnakeModel == 48 || mSnakeModel > 50 ){
	  new G4PVPlacement(pRotXRHRSdeg,G4ThreeVector(-pPaulY, 0, -pPaulX),
			    LPlaneLogical1,"virtualBoundaryPhys_col_RHRS",motherLogical,0,0,0);//q1en

	  //double pRQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE
	  /*
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z_en,0),
			    RPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RHRSContainerLogical,0,0,0);//q1en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z_ex,0),
			    RPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RHRSContainerLogical,0,0,0);//q1ex
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_en,0),
			    RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RHRSContainerLogical,0,0,0);//q2en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_ex,0),
			    RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RHRSContainerLogical,0,0,0);//q2ex
	  */
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ1Length/2.),
			    RPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RQ1MagLogical,0,0,0);//q1en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ1Length/2.),
			    RPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RQ1MagLogical,0,0,0);//q1ex
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ2Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RQ2MagLogical,0,0,0);//q2en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ2Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RQ2MagLogical,0,0,0);//q2ex
	  //new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pRDPos_Z_en,0),
	  //RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RHRSContainerLogical,0,0,0);//den
	  //new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pRDPos_Z_ex,pRDPos_X_ex),
	  //RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RHRSContainerLogical,0,0,0);//dex
	  
	  new G4PVPlacement(pRot_den,G4ThreeVector(-8.4*m,0,0),
			    RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RDipoleTunnelLogical,0,0,0);//den
	  new G4PVPlacement(pRot_dex,G4ThreeVector(-8.4 * m * cos( pi / 4. ), -8.4 * m * sin( pi / 4. ),0),
			    RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RDipoleTunnelLogical,0,0,0);//dex
	  
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ + pQ3Length / sqrt(2) / 2,pQ3CenterY - pQ3Length / sqrt(2) / 2),
	  //RPlaneLogical2,"virtualBoundaryPhys_q3en_RHRS",RHRSContainerLogical,0,0,0);//q3en
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ - pQ3Length / sqrt(2) / 2,pQ3CenterY + pQ3Length / sqrt(2) / 2),
	  //RPlaneLogical2,"virtualBoundaryPhys_q3ex_RHRS",RHRSContainerLogical,0,0,0);//q3ex
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ3Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q3en_RHRS",RQ3MagLogical,0,0,0);//q3en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ3Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q3ex_RHRS",RQ3MagLogical,0,0,0);//q3ex


//	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-vdc1Z,vdc1Y),
//			    vdc1LLogical,"vdc1_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-vdc1Z,vdc1Y),
			    vdc1LLogical,"vdc1_LHRS",LHRSContainerLogical,0,0,0);

/*
          new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-(vdc1Z-150.*cm), 505.*cm),
                            vdc1LLogical,"Q3_L1",LHRSContainerLogical,0,0,0);

	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-(vdc1Z-150.*cm), 515.*cm),
                            vdc1LLogical,"Q3_L2",LHRSContainerLogical,0,0,0);

	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-(vdc1Z-150.*cm), 495.*cm),
                            vdc1LLogical,"Q3_L3",LHRSContainerLogical,0,0,0);
*/
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-vdc1Z,vdc1Y),
			    vdc1RLogical,"vdc1_RHRS",RHRSContainerLogical,0,0,0);

	  new G4PVPlacement(pRotFocPlL,G4ThreeVector(0,-vdc1Z,vdc1Y),
			    FocPlLLogical,"FocPl_LHRS",LHRSContainerLogical,0,0,0);


	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pVDCCenterZ,pVDCCenterY),
			    RFPLogical,"virtualBoundaryPhys_vdc_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pFPCenterZ,pFPCenterY),
			    RFPLogical,"virtualBoundaryPhys_fp_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pQZ1CenterZ,pQZ1CenterY),
			    RFPLogical,"virtualBoundaryPhys_qz1_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pQZ2CenterZ,pQZ2CenterY),
			    RFPLogical,"virtualBoundaryPhys_qz2_RHRS",RHRSContainerLogical,0,0,0);

	}

	//#endif
	//////////////////////////////////////////////////////////

	return theHRSPhys;
}


/////////////////////////////////////////////////////////////////////


G4VPhysicalVolume* HRSDetectorConstruction::ConstructRadiator(G4LogicalVolume* motherLogical)
{
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4VSensitiveDetector* virtualDetectorSD=new HRSStdSD(SDname="virtualDetector_Rad");

	// sensitive detectors   
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4VPhysicalVolume* theRadiatorPhys = 0;

	////////////////////////
	G4Material *theRadMaterial = mMaterialManager->copper;
	if(mRadiatorType == 1) theRadMaterial = mMaterialManager->iron;
	else if(mRadiatorType == 2) theRadMaterial = mMaterialManager->copper;
	else if(mRadiatorType == 3) theRadMaterial = mMaterialManager->tantalum;
	else if(mRadiatorType == 4) theRadMaterial = mMaterialManager->tungsten;


	/////////////////////////
	//the radiator
	/////////////////////////
	double pRadX=4*cm, pRadY=4*cm;
	double pRadZ=mRaditorThickness;

	G4VSolid* radiatorSolid = new G4Box("radiatorBox",pRadX/2.0,pRadY/2.0,pRadZ/2.0);
	G4LogicalVolume* radiatorLogical = new G4LogicalVolume(radiatorSolid,
		theRadMaterial,"radiatorLogical",0,0,LarmStepLimits);
	radiatorLogical->SetVisAttributes(OrangeVisAtt);


	double pRadPos_X=mPivotXOffset, pRadPos_Y=mPivotYOffset;
	double pRadPos_Z=mPivotYOffset + mRadiator2Pivot;
	theRadiatorPhys = new G4PVPlacement(0,
		G4ThreeVector(pRadPos_X,pRadPos_Y,pRadPos_Z),
		radiatorLogical,"radiatorPhys",motherLogical,0,0);


	/////////////////////////
	//the radiator virtual ditector
	/////////////////////////
	if(mSetupRadiatorVD)
	if (0)
	{
		double pRadVDX=20*cm, pRadVDY=20*cm, pRadVDZ=2*mm;

		G4VSolid* radiatorVDSolid = new G4Box("radiatorVDBox",
			pRadVDX/2.0,pRadVDY/2.0,pRadVDZ/2.0);
		G4LogicalVolume* radiatorVDLogical = new G4LogicalVolume(radiatorVDSolid,
			mMaterialManager->air,"radiatorVDLogical",0,0,LarmStepLimits);

		SDman->AddNewDetector(virtualDetectorSD);
		radiatorVDLogical->SetSensitiveDetector(virtualDetectorSD);
		radiatorVDLogical->SetVisAttributes(LightYellowVisAtt);


		//place the VD 2 cm downstream from the downstream end plane of the radiator
		double pRadVDPos_X=pRadPos_X, pRadVDPos_Y=pRadPos_Y;
		double pRadVDPos_Z=pRadPos_Z + mRaditorThickness/2 + 2*cm;
		G4String VDPhysName = (mSetupRadiatorVD==2) ? "virtualBoundaryPhys_Rad" : "radiatorVDPhys";
		new G4PVPlacement(0,G4ThreeVector(pRadVDPos_X,pRadVDPos_Y,pRadVDPos_Z),
			radiatorVDLogical,VDPhysName,motherLogical,0,0);
	}

	return theRadiatorPhys;
}

/////////////////////////////////////////////////////////////////////

