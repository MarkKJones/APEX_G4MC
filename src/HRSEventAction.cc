// ********************************************************************
//
// $Id: HRSEventAction.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#include "math.h"
#include "HRSEventAction.hh"
#include "HRSEventActionMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "HRSAnalysisManager.hh"
#endif // G4ANALYSIS_USE

#include "HRSTrajectory.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4PrimaryParticle.hh"
#include "G4ThreeVector.hh"
#include "HRSStdHit.hh"
#include "HRSDCHit.hh"
#include "HRSCalorimeterHit.hh" 
//#include "HRSRootTree.hh"
#include "UsageManager.hh"

//I am not using this enum list any longer
//remember to add your SD name into this emum list
enum SensitiveDetectorList {  
	FZB1VD=1,FZB2VD=2,
	sieveSlit=11,septumWindow=12,
	thirdArmShielding=30, thirdArmSC1=31, thirdArmSC2=32, 
	virtualDetector=40,
	BBSC1=51,BBSC2=52,	
	VETO=60,NDSC1=61,NDSC2=62,NDSC3=63,NDSC4=64,	
	CREXUpBlock=70,CREXTarget=71,CREXDownBlock=72,
	virtualBoundary=99,
	RTPC=200
};

//extern  HRSRootTree* gHRSTree;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//verbose level: 
//0:  nothing will be printed
//1:  HRSStdHit::Print() will be called
//2:  Event header will be printed, as well as all primary particle info (print togather). 
//3:  HRSTrajectory::ShowTrajectory() will be called;
//4:  HRSEventAction::PrintPrimary() will be call;

HRSEventAction::HRSEventAction()
{
	verboseLevel = 0;
	messenger = new HRSEventActionMessenger(this);

	//G4cout<<"HRSEventAction() construction done!"<<G4endl;
}

HRSEventAction::~HRSEventAction()
{
	delete messenger;
	//G4cout<<"delete HRSEventAction ... done!"<<G4endl;
}

void HRSEventAction::BeginOfEventAction(const G4Event* evt)
{
  cout<<"evtNo:"<<i_ev<<endl;
  i_ev++;
}

void HRSEventAction::ProcessDCSDHits(G4HCofThisEvent* HCE)
{
}


void HRSEventAction::ProcessStdSDHits(G4HCofThisEvent* HCE)
{
}


void HRSEventAction::ProcessCalorimeterSDHits(G4HCofThisEvent* HCE)
{
}


void HRSEventAction::ProcessSDHits(G4HCofThisEvent* HCE)
{
}


void HRSEventAction::EndOfEventAction(const G4Event* evt)
{
}


void HRSEventAction::PrintPrimary(G4PrimaryParticle* pp,G4int ind)
{
}

