// ********************************************************************
//
// $Id: HRSSteppingAction.cc,v3.1 2008/3/16 HRS Exp $
//
//..............................................................................
#include <iomanip>
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4String.hh"
#include "HRSSteppingAction.hh"
#include "G4FieldManager.hh"
#include "G4SteppingManager.hh"
#include "HRSSteppingActionMessenger.hh"
#include "HRSTrackInformation.hh"
#include "UsageManager.hh"
#include "HRSTransform_TCSNHCS.hh"
#include "HRSEMField.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "HRSPrimaryGeneratorAction.hh"

//..............................................................................

//#define STEPPING_DEBUG 2
extern UsageManager* gConfig;
//extern bool CreateFileName(char *instr,char *outstr,int type);
HRSSteppingAction::HRSSteppingAction()
{
	PrintHeadFlag=true;
	//verboseLevel=2;
	verboseLevel=0;
	vPrintList.push_back("all"); 

	int pBookTrees=0, pBookHistos=0, pBookTxt=0;
	gConfig->GetParameter("BookTrees",pBookTrees);
        gConfig->GetParameter("SeptumOn",mSeptumOn);
        gConfig->GetParameter("LHRSAngle",mLHRSAngle);
        mLHRSAngle*=deg;
        //    mLHRSAngle=45.*deg;
//        ang_fp=45.032*deg;
	gConfig->GetParameter("LFocalPlaneAngle",ang_fp);
	ang_fp*=deg;
	gConfig->GetParameter("BookHistos",pBookHistos);
	gConfig->GetParameter("BookTxt",pBookTxt);

	CreateTxt=(pBookTxt>0);
	CreateRootNt=(pBookTrees>0 || pBookHistos>0);

	messenger = new HRSSteppingActionMessenger(this);

	//G4cout<<"HRSSteppingAction() construction done!"<<G4endl;
        prevevt_id=-1;

	//G4cout<<"HRSSteppingAction() construction done!"<<G4endl;

//        output_vdc.open("vdc_coord.C", ios::out | ios::app );
//        output_vdc <<" {"<<endl;
//        output_vdc <<" TH2D * vdc = new TH2D(\"vdc xy\", \"vdc xy\", 400, -40., 40., 400, -40., 40.);"<<endl;

        output.open("coord.C", ios::out | ios::app );
        output <<"{"<<endl;
        output <<"double x[5000];"<<endl;
        output <<"double y[5000];"<<endl;

        output1.open("coord_xz.C", ios::out | ios::app );
        output1 <<"{"<<endl;
        output1 <<"double x[5000];"<<endl;
        output1 <<"double y[5000];"<<endl;

        output_focpl.open("FocalPlanecoords.txt", ios::out | ios::app );

        i_st=0;

//	cout<<"stepping test 1"<<endl;
//        G4FieldManager *theFieldManager;// = theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetFieldManager();
//	cout<<"stepping test 1"<<endl;
//	if(!theFieldManager) theFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
//	if(G4TransportationManager::GetTransportationManager()->GetFieldManager()->DoesFieldExist());
//	double Point[4]={700.*cm,0.*cm,700.*cm,0.};
//	double fieldArr[6]={0,0,0,0,0,0};

	P1x=0.;
	P1y=0.;
	P1z=0.;
	P2x=0.;
	P2y=0.;
	P2z=0.;
        if (1==2)
        {
          HRSEMField * ffield;
          ffield = new HRSEMField();
          output1.open("ffield.C", ios::out | ios::app );
          output1 <<"{"<<endl;
          output1 <<"TH2D * Bx = new TH2D(\"Bx\", \"Bx\", 2300, 0, 2300, 1000, -200, 800);"<<endl;
	  for (double xz=1.*cm; xz<=2300.*cm; xz+=1.*cm)
	  for (double yy=-200.*cm; yy<=800.*cm; yy+=1.*cm)
	  {

                  G4double coord_d[4]={xz*sin(mLHRSAngle),yy,xz*cos(mLHRSAngle),0};
                  G4double bdasht[3]={0,0,0};
                  ffield->GetFieldValue(coord_d, bdasht) ;
                  if ((bdasht[0]!=0.) && (bdasht[1]!=0) && (bdasht[2]!=0.) && (fabs(bdasht[0])<=1000.*tesla) && (fabs(bdasht[1])<=1000.*tesla) && (fabs(bdasht[2])<=1000.*tesla) )
                  {
                    double angg=acos(bdasht[2]/sqrt(bdasht[2]*bdasht[2]+bdasht[0]*bdasht[0]))*fabs(bdasht[0])/bdasht[0]-mLHRSAngle;
                    double b_new=sqrt(bdasht[2]*bdasht[2]+bdasht[0]*bdasht[0])*sin(angg);
//                  cout<<"dip   x="<<bdasht[0]<<"  y="<<bdasht[1]<<"  z="<<bdasht[2]<<endl;
                    output1 << "  Bx->Fill("<<xz/cm<<", "<< yy/cm <<", "<<b_new/tesla<< ");"<<endl;
                  }
                  else 
                  {
                    output1 << "  Bx->Fill("<<xz/cm<<", "<< yy/cm <<", "<<-5999.<< ");"<<endl;
                  }
          }
          output1 <<"  Bx->Draw(\"colz\");"<<endl;
          output1 <<"}"<<endl;
          output1 <<endl;
          output1.close();
        }


/*
        G4FieldManager *globalFieldManager;// = theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetFieldManager();
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

//..............................................................................

HRSSteppingAction::~HRSSteppingAction()
{
        output <<"int npoints="<<i_st<<";"<<endl;
        output <<"TGraph *gr1 = new TGraph(npoints,x,y);"<<endl;
        output <<"gr1->SetMarkerStyle(8);"<<endl;
        output <<"gr1->Draw(\"APL\");"<<endl;
        output <<"}"<<endl;
        output <<endl;
        output.close();

        output1 <<"int npoints="<<i_st<<";"<<endl;
        output1 <<"TGraph *gr1 = new TGraph(npoints,x,y);"<<endl;
        output1 <<"gr1->SetMarkerStyle(8);"<<endl;
        output1 <<"gr1->Draw(\"APL\");"<<endl;
        output1 <<"}"<<endl;
        output1 <<endl;
        output1.close();

//        output_vdc <<" vdc->Draw(\"colz\");"<<endl;
//        output_vdc <<" }"<<endl;
//        output_vdc <<endl;
//        output_vdc.close();


	if (CreateTxt)
	{
		OutTxt.close();
	}
	vPrintList.clear();
	delete messenger;
	//G4cout<<"delete HRSSteppingAction ... done!"<<G4endl;
}
//..............................................................................

void HRSSteppingAction::InitOutTxt()
{
	if (!CreateTxt) return;

	char strFileName[200], strRawName[200];
	
	std::string pOutFileName=gConfig->GetArgument("OutFileName");
/*
	if(gHRSTree)
	{
		pOutFileName=gHRSTree->GetRootFileName();
	}
	if(pOutFileName.length()>3) 
	{
		string tmpStr(pOutFileName);
		size_t pos=tmpStr.rfind(".root");
		if(pos!=string::npos) tmpStr.replace(pos,5,".txt"); 
		strcpy(strFileName,tmpStr.c_str());
	}
	else
	{
		//sprintf(strRawName,"G4Sim_txt_%02d.txt",gHRSTree->iRunNumber);
		//CreateFileName(strRawName,strFileName);
		//in order to make the txt output has the same name with the root output, I check 
		//the existance of the ntuple root file instead of txt file
		sprintf(strRawName,"G4Sim_nt_%02d.root",gHRSTree->iRunNumber);
		CreateFileName(strRawName,strFileName,2);
	}
	OutTxt.open (strFileName, fstream::out | fstream::app);
*/
}

void HRSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
//  cout << "WE ARE IN THE USER STEPPING ACTION DOOOOOOOOOOOOOOOOOOOOP!" << endl;
        G4VPhysicalVolume* curPV  = theStep->GetPreStepPoint()->GetPhysicalVolume();
        G4String name = curPV->GetName();
        name.assign(name,0,70);
        G4String name1 = curPV->GetName();
        name1.assign(name1,0,40);
        G4String name2 = curPV->GetName();
        name2.assign(name2,0,7);

        double edeposit   =  theStep->GetTotalEnergyDeposit();
        //e_sum+=edeposit/1000.;
//        G4int eID = 0;
//        const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
//        if(evt) eID = evt->GetEventID();

        G4int evtNb = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
        evt_id=evtNb;
        if (prevevt_id != evt_id)
        {
          track_id_rad_tail_foil=-1*evt_id;
          px_rad_tail_foil=-999.1;
          py_rad_tail_foil=-999.2;
          pz_rad_tail_foil=-999.3 ;
          track_id_sieve_back==-1*evt_id;
          px_sieve_back=-999.1;
          py_sieve_back=-999.2;
          pz_sieve_back=-999.3;
          x_sieve_back=-1999.1;
          y_sieve_back=-1999.2;
          z_sieve_back=-1999.3;

          write_trig=true;
          prevevt_id =evt_id;
//          cout<<"mass    "<<E1v<<setw(18)<<P1x<<setw(18)<<P1y<<setw(18)<<P1z<<setw(18)<<E2v<<setw(18)<<P2x<<setw(18)<<P2y<<setw(18)<<P2z<<endl;
          if (pz_tmp_mass != 0)
          {
/*
            if (P1z == 0)
            {
              P1x=px_tmp_mass;
	      P1y=py_tmp_mass;
	      P1z=pz_tmp_mass;
	      E1v=e_tmp_mass;
	    }
*/
//	    else if ((P2z == 0) && (P1z != pz_tmp_mass) )
	    {
//              P2x=px_tmp_mass;
//	      P2y=py_tmp_mass;
//	      P2z=pz_tmp_mass;
//	      E2v=e_tmp_mass;
	      if (el_n_mass == pos_n_mass)
	      {
                std::ofstream output_mass;
                output_mass.open( "mass.txt", ios::out | ios::app );
                output_mass<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()
                <<setw(15)<<E1v<<setw(18)<<P1x<<setw(18)<<P1y<<setw(18)<<P1z
                <<setw(18)<<E2v<<setw(18)<<P2x<<setw(18)<<P2y<<setw(18)<<P2z<<endl;
//              cout<<"mass    "<<E1v<<setw(18)<<P1x<<setw(18)<<P1y<<setw(18)<<P1z<<setw(18)<<E2v<<setw(18)<<P2x<<setw(18)<<P2y<<setw(18)<<P2z<<endl;
                output_mass.close();
                el_n_mass =-1;
                pos_n_mass=-2;
              }
            }
	  }
//          P1x=0.;
//	  P1y=0.;
//	  P1z=0.;
//	  E1v=0.;
//          P2x=0.;
//	  P2y=0.;
//	  P2z=0.;
//	  E2v=0.;
//          e_tmp_mass=0.;
//          px_tmp_mass=0.;
//	  py_tmp_mass=0.;
//	  pz_tmp_mass=0.;
        }


        G4double xxx=theStep->GetTrack()->GetPosition().x();
        G4double yyy=theStep->GetTrack()->GetPosition().y();
        G4double zzz=theStep->GetTrack()->GetPosition().z();
        G4double z=theStep->GetTrack()->GetPosition().z()/cm;
        G4double x=theStep->GetTrack()->GetPosition().x()/cm;
        G4double y=theStep->GetTrack()->GetPosition().y()/cm;


        double y_TCS=y;
        double thet_tcs=acos(z/sqrt(x*x+z*z))-mLHRSAngle;
    //    cout<<"thet_tcs="<<thet_tcs/deg<<endl;
        double x_TCS=sqrt(x*x+z*z)*sin(thet_tcs);
        double z_TCS=sqrt(x*x+z*z)*cos(thet_tcs);
//        cout<<"Target CS:"<<"  angle="<<thet_tcs<<"     x="<<x_TCS<<"   y="<<y_TCS<<"   z="<<z_TCS<<endl;

        G4double x_vert=theStep->GetTrack()->GetVertexPosition().x()/cm;
        G4double y_vert=theStep->GetTrack()->GetVertexPosition().y()/cm;
        G4double z_vert=theStep->GetTrack()->GetVertexPosition().z()/cm;

        float r=sqrt(x*x + y*y);
        float ang_phi=acos(x/r) * y/fabs(y);
        float ang_theta = atan(r/z);
        float eta_tmp = -1*log(fabs(tan(ang_theta/2)))*z/fabs(z);
        double xcoord=theStep->GetPreStepPoint()->GetPosition().x();
        double ycoord=theStep->GetPreStepPoint()->GetPosition().y();
        double zcoord=theStep->GetPreStepPoint()->GetPosition().z();
        double px_step=theStep->GetPreStepPoint()->GetMomentum().x()/GeV;
        double py_step=theStep->GetPreStepPoint()->GetMomentum().y()/GeV;
        double pz_step=theStep->GetPreStepPoint()->GetMomentum().z()/GeV;
        double pt=sqrt(theStep->GetTrack()->GetMomentum().x()*theStep->GetTrack()->GetMomentum().x() + theStep->GetTrack()->GetMomentum().y()*theStep->GetTrack()->GetMomentum().y());

        double p_HCS=sqrt(px_step*px_step + py_step*py_step + pz_step*pz_step);

        double p_thet=asin(px_step/sqrt(px_step*px_step + pz_step*pz_step))-mLHRSAngle;
        double p_x_TCS=sqrt(px_step*px_step + pz_step*pz_step)*sin(p_thet);
        double p_y_TCS=py_step;
        double p_z_TCS=sqrt(px_step*px_step + pz_step*pz_step)*cos(p_thet);
        double p_phi =asin(py_step/sqrt(p_HCS*p_HCS-p_x_TCS*p_x_TCS));
        
        double px_dir_vert=theStep->GetTrack()->GetVertexMomentumDirection().x();
        double py_dir_vert=theStep->GetTrack()->GetVertexMomentumDirection().y();
        double pz_dir_vert=theStep->GetTrack()->GetVertexMomentumDirection().z();
        double thet_vert = asin(px_dir_vert/sqrt(px_dir_vert*px_dir_vert + pz_dir_vert*pz_dir_vert))-mLHRSAngle;
        double p_z_vert_tcs=sqrt(px_dir_vert*px_dir_vert + pz_dir_vert*pz_dir_vert)*cos(thet_vert);
        double phi_vert = asin(py_dir_vert/p_z_vert_tcs);



	double t_targ=sqrt(x*x+z*z)*sin(acos(z/sqrt(x*x+z*z))-mLHRSAngle);
//	output1 <<"x["<<i_st<<"]="<<sqrt(x*x+z*z)<<";  y["<<i_st<<"]="<<t_targ<<";"<<endl;
//	output <<"x["<<i_st<<"]="<<sqrt(x*x+z*z)<<";  y["<<i_st<<"]="<<y<<";"<<endl;
	i_st++;
	double tmpField[6]={0,0,0,0,0,0};
	//By Jixie: In case there is a local field, this is the correct way
	G4FieldManager *theFieldManager = theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetFieldManager();   
	if(!theFieldManager) theFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();

	if(theFieldManager)
	{
		double tmpPoint[4]={xxx*mm,yyy*mm,zzz*mm,0};
		theFieldManager->GetDetectorField()->GetFieldValue(tmpPoint,tmpField);
//        	cout<<tmpField[0]/tesla<<endl;
//        	cout<<tmpField[1]/tesla<<endl;
//        	cout<<tmpField[2]/tesla<<endl;
//        	cout<<tmpField[3]<<endl;
//        	cout<<tmpField[4]<<endl;
//        	cout<<tmpField[5]<<endl;
	}




//    	G4FieldManager bfld=GetDetectorField()->GetFieldValue();
//        if (1==2) cout<<evtNb<<"x="<<xcoord<<"    px_step="<<px_step<<"    y="<<ycoord<<"    py_step="<<py_step<<"    z="<<zcoord<<"    pz_step="<<pz_step<<endl;
//        if (1==2) cout<<"   material name is             "<< theStep->GetPreStepPoint()->GetMaterial()->GetName()<<"       z = "<< theStep->GetPreStepPoint()->GetPosition().z()/10<<"      (x, y) = ("<<x<<","<<y<<")"<<endl<<endl;
//        std::cout << "Event:"<<evtNb<<"   name: "<<name<<"     "<<"     Track:"<< theStep->GetTrack()->GetTrackID() <<"    "<<theStep->GetTrack()->GetDefinition()->GetParticleName()<< "   x=" <<theStep->GetTrack()->GetPosition().x()/10 << "   y=" <<theStep->GetTrack()->GetPosition().y()/10 << "   z=" <<theStep->GetTrack()->GetPosition().z()/10<<"    energy = ?? add later" <<"    e_dep = " << edeposit/1000. <<std::endl;
	  double z00=9.961*m;
	  double y00=8.4*m;
	  double R_new= sqrt( (y00-yyy)*(y00-yyy) + (zzz - z00) * (zzz - z00) );

//          if (z < -0.)
//          if (z > -4.9)
          if (0)
          std::cout << "Event:"<<evtNb<<"  name: "<<setw(30)<<name<<"   Track:"<< theStep->GetTrack()->GetTrackID() 
                    <<"  "<<theStep->GetTrack()->GetDefinition()->GetParticleName()<<"  P="<<p_HCS<<"  E0="<<0.00051+theStep->GetTrack()->GetVertexKineticEnergy()/GeV<<"  x=" <<setw(10)<<x << ", y="<<setw(10)<<y << ", z=" <<setw(10)<<z
                    <<",     theta="<<setw(10)<<atan(px_step/pz_step)/deg<<"      phi="<<atan(py_step/pz_step)/deg
                    <<"            tar:"<<setw(10)<<x_TCS<<", "<<setw(10)<<y_TCS<<", "<<setw(10)<<z_TCS
                    <<" p="<<setw(6)<<p_HCS<<"    phi_tar="<<setw(10)<<p_phi/deg<<"   p_thet_tar="<<setw(10)<<p_thet
                    <<",   px=" <<setw(10)<<px_step<<", py=" <<setw(10)<<py_step<<", pz=" <<setw(10)<< pz_step 
                    <<",  B="<<setw(10)<<tmpField[0]/tesla<<", "<<setw(10)<<tmpField[1]/tesla<<", "<<setw(10)<<tmpField[2]/tesla<<");"<<"  e_dep:" <<setw(10)<< edeposit/GeV 
                    <<"   material: "<< theStep->GetPreStepPoint()->GetMaterial()->GetName()<<"   "<<theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<endl;
          if (name=="targtrans")
          {
            std::ofstream output_e_loss;
            output_e_loss.open( "e_loss_vs_target.txt", ios::out | ios::app );
            output_e_loss<<evtNb<<setw(15)<<z<<setw(20)<<px_step<<setw(20)<<py_step<<setw(20)<<pz_step<<endl;
            output_e_loss.close();
          }
//          std::cout << "Event:"<<evtNb<<"  name: "<<name<<"   Track:"<< theStep->GetTrack()->GetTrackID() <<"  "<<theStep->GetTrack()->GetDefinition()->GetParticleName()<< "  x=" <<x << ", y=" <<y << ", z=" <<z<<", px=" <<px_step<<", py=" <<py_step<<", pz=" << pz_step <<",    bx="<<tmpField[0]/tesla<<", by="<<tmpField[1]/tesla<<", bz="<<tmpField[2]/tesla<<"     r_new="<<(R_new-y00)/cm<<"     Vertex[cm]:("<<x_vert<<","<<y_vert<<","<<z_vert<<")  direction:"<<theStep->GetTrack()->GetVertexMomentumDirection()<<"   material name is             "<< theStep->GetPreStepPoint()->GetMaterial()->GetName()<<endl;
//          cout<<
//          std::cout << "Event:"<<evtNb<<"  name: "<<name<<"   Track:"<< theStep->GetTrack()->GetTrackID() <<"  "<<theStep->GetTrack()->GetDefinition()->GetParticleName()<< "  x=" <<x << ", y=" <<y << ", z=" <<z<<", px=" <<px_step<<", py=" <<py_step<<", pz=" << pz_step <<",    bx="<<tmpField[0]/tesla<<", by="<<tmpField[1]/tesla<<", bz="<<tmpField[2]/tesla<<"     r_new="<<(R_new-y00)/cm<<std::endl;
//   -200., -100.,-200.      ;minbox x,y,z
//   +200., +1600.,+200.      ;maxbox

        if (z_TCS>300.)
	if(0)
        {
          cout<<"the track is killed after Q1. Go to SIMC simulations"<<endl;
          theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }

        if ( (z_TCS>620.) && (z_TCS<630.) )
        if ( (sqrt(x_TCS*x_TCS+y_TCS*y_TCS)>31.) )
        {
            cout<<"the track is killed Q2 ex"<<endl;
            theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }
/*
        if (theStep->GetTrack()->GetTrackID()!=1)
        {
            cout<<"the track is killed - secondary particle"<<endl;
            theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }
*/
        double pQ3Cen_y_tmp =  358.637  + 182./sqrt(2.)/2.;
        double pQ3Cen_z_tmp = 1702.67042  + 182./sqrt(2.)/2.;

        if ( (z_TCS>pQ3Cen_z_tmp) && (z_TCS<pQ3Cen_z_tmp+10.) )
        if (sqrt( x_TCS*x_TCS+(y_TCS-pQ3Cen_y_tmp)*(y_TCS-pQ3Cen_y_tmp) ) > 45. )
        {
          cout<<"the track is killed Q3 cen"<<endl;
          theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }

        if ( (z_TCS>-10.) && (z_TCS<160.) )
        if ( (fabs(x_TCS)>20.) || (fabs(y_TCS)>20.) )
        {
            cout<<"the track is killed before Q1"<<endl;
            theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }
        if (fabs(y)>820.)
        {
            cout<<"the track is killed y > 820"<<endl;
            theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }


          //septum vacuuum box
          if (mSeptumOn)
          if ((z>-109.) && z<(160.))
          {
            double inch           = 2.54;
            double zlength        = 173.939;
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
            double length_max_1   = 50.19;//* cm; // 19.76*inch
            double length_min_1   = 52.1 ;//* cm; // 19.76*inch
            double length_min_2   = 59.82;//* cm; // 23.55*inch
            double length_max_2   = 60.3 ;//* cm; // 23.74*inch

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

//            cout<<"kuku  "<<z<<"     "<<x<<"        "<<xsep_tmp<<endl;

/*
            if ((z>z_sept_en_min1) && (z<z_sept_ex_min1))
            {
              double xsep_tmp = xmin_sep_en1 + (xmin_sep_ex1-xmin_sep_en1)*(z-z_sept_en_min1)/(z_sept_ex_min1-z_sept_en_min1);
              if (xsep_tmp > x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
            }

            if ((z>z_sept_en_min2) && (z<z_sept_ex_min2))
            {
              double xsep_tmp = xmin_sep_en2 + (xmin_sep_ex2-xmin_sep_en2)*(z-z_sept_en_min2)/(z_sept_ex_min2-z_sept_en_min2);
//              cout<<"kuku  "<<z<<"     "<<x<<"        "<<xsep_tmp<<endl;
              if (xsep_tmp > x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
            }

            if ((z>z_sept_en_max1) && (z<z_sept_ex_max1))
            {
              double xsep_tmp = xmax_sep_en1 + (xmax_sep_ex1-xmax_sep_en1)*(z-z_sept_en_max1)/(z_sept_ex_max1-z_sept_en_max1);
              if (xsep_tmp < x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
            }
            if ((z>z_sept_en_max2) && (z<z_sept_ex_max2))
            {
              double xsep_tmp = xmax_sep_en2 + (xmax_sep_ex2-xmax_sep_en2)*(z-z_sept_en_max2)/(z_sept_ex_max2-z_sept_en_max2);
              if (xsep_tmp < x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
            }
*/

            if ( (z > z_sept_en_max1) && (z < z_sept_ex_min2) )
            {
              double y_en = 2.44*inch;
              double y_ex = 4.7 *inch;
              double ysep_tmp = y_en + (z-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
//              cout<<z<<"     "<<ysep_tmp<<endl;
              if (fabs(y) > ysep_tmp )
              {
                cout<<"the track is killed (y > ysep_tmp)"<<endl;
                theStep->GetTrack()->SetTrackStatus(fStopAndKill);
              }
            }
          }




        if (name == "virtual")
        {
          if ( fabs(y_TCS-246.) < 10.)
          if ( fabs(z_TCS-1590.) < 40.)
          if ( fabs(x_TCS) < 15.)
          {
            cout<<"   look here ------------------> "<<" tar:"<<x_TCS<<", "<<y_TCS<<", "<<z_TCS<<"   p="<<p_HCS<<"   phi="<<p_phi/deg<<"   p_thet="<<p_thet<<endl;
            cout<<"   look here ------------------> the distance is "<<sqrt((z_TCS-1590.)*(z_TCS-1590.) + (y_TCS-246.)*(y_TCS-246.))<<endl;
          }
        }

//        if (theStep->GetTrack()->GetTrackID()==1)
        if (p_HCS<1.27)
        {
          cout<<"the track is killed (low momentum track, p="<<p_HCS<<" GeV)"<<endl;
          theStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }


        if (p_HCS>0.5)
        if ((name2 == "RSvSlBa") || (name2 == "LSvSlBa"))
        {
          double e_kin = theStep->GetTrack()->GetVertexKineticEnergy()+0.511;
          double px_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().x();
          double py_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().y();
          double pz_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().z();
          std::ofstream output_q1;
          output_q1.open( "Sieve_back.dat", ios::out | ios::app );
          output_q1<<evtNb<<setw(18)<<theStep->GetTrack()->GetTrackID()<<setw(18)<<p_HCS<<setw(18)<<x<<setw(18)<<y<<setw(18)<<z<<setw(18)<<px_step<<setw(18)<<py_step<<setw(18)<<pz_step<<setw(18)<<px_vtx_tmp<<setw(18)<<py_vtx_tmp<<setw(18)<<pz_vtx_tmp<<setw(18)<<theStep->GetTrack()->GetDefinition()->GetParticleName()<<setw(18)<<track_id_rad_tail_foil<<setw(18)<<px_rad_tail_foil<<setw(18)<<py_rad_tail_foil<<setw(18)<<pz_rad_tail_foil<<endl;
//          output_q1<<atan(px_vtx_tmp/pz_vtx_tmp)*180./3.1415927<<endl;
          if ( pz_step>pz_sieve_back )
          {
            track_id_sieve_back = theStep->GetTrack()->GetTrackID();
            x_sieve_back = x;
            y_sieve_back = y;
            z_sieve_back = z;
            px_sieve_back = px_step;
            py_sieve_back = py_step;
            pz_sieve_back = pz_step;
          }
          output_q1.close();
        }

//        if (theStep->GetTrack()->GetTrackID()==1)
        if (p_HCS>0.5)
        if (name2 == "LQ1Fron")
        if (theStep->GetTrack()->GetDefinition()->GetParticleName()=="e-")
        {
          double e_kin = theStep->GetTrack()->GetVertexKineticEnergy()+0.511;
          double px_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().x();
          double py_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().y();
          double pz_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().z();
          std::ofstream output_q1;
          output_q1.open( "Q1_front.dat", ios::out | ios::app );
          output_q1<<evtNb<<setw(18)<<theStep->GetTrack()->GetTrackID()<<setw(18)<<z_TCS<<setw(18)<<p_HCS<<setw(18)<<p_x_TCS/p_z_TCS<<setw(18)<<p_y_TCS/p_z_TCS<<setw(18)<<x_TCS<<setw(18)<<y_TCS<<setw(18)<<px_vtx_tmp<<setw(18)<<py_vtx_tmp<<setw(18)<<pz_vtx_tmp<<setw(18)<<theStep->GetTrack()->GetDefinition()->GetParticleName()<<setw(18)<<track_id_rad_tail_foil<<setw(18)<<px_rad_tail_foil<<setw(18)<<py_rad_tail_foil<<setw(18)<<pz_rad_tail_foil<<setw(18)<<track_id_sieve_back<<setw(18)<<x_sieve_back<<setw(18)<<y_sieve_back<<setw(18)<<z_sieve_back<<setw(18)<<px_sieve_back<<setw(18)<<py_sieve_back<<setw(18)<<pz_sieve_back<<endl;
          output_q1.close();
        }
        if (name=="RadTailDet")
        {
          double e_kin = theStep->GetTrack()->GetVertexKineticEnergy()+0.511;
          double px_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().x();
          double py_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().y();
          double pz_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().z();


          if ( pz_step>pz_rad_tail_foil )
          {
            track_id_rad_tail_foil = theStep->GetTrack()->GetTrackID();
            px_rad_tail_foil = px_step;
            py_rad_tail_foil = py_step;
            pz_rad_tail_foil = pz_step;
          }
          std::ofstream rad_tail;
          rad_tail.open( "after_W_foil.dat", ios::out | ios::app );
          rad_tail<<evtNb<<setw(18)<<theStep->GetTrack()->GetTrackID()<<setw(18)<<p_HCS<<setw(18)<<x<<setw(18)<<y<<setw(18)<<z<<setw(18)<<px_step<<setw(18)<<py_step<<setw(18)<<pz_step<<setw(18)<<px_vtx_tmp<<setw(18)<<py_vtx_tmp<<setw(18)<<pz_vtx_tmp<<setw(18)<<theStep->GetTrack()->GetDefinition()->GetParticleName()<<endl;
          rad_tail.close();
          if( p_HCS > 1.8 )
          {
            cout<<"the track has too high energy after rad_tail_foil. Killed!"<<endl;
            theStep->GetTrack()->SetTrackStatus(fStopAndKill);
          }
        }

        if (name2 == "FocPl_L")
        {
          double vdc1Z=2075.965;
          double vdc1Y=731.9258;
          double p_y_fp = p_x_TCS;
          cout<<"Focal plane TCS:"<<" px="<<p_x_TCS<<"   py="<<p_y_TCS<<"  pz="<<p_z_TCS<<"     x="<< x_TCS<<"   y="<<y_TCS<<"   z="<<z_TCS<<"   vdc_y="<<vdc1Y<<"   vdc_z="<<vdc1Z<<endl;
          double p_thet_fp=3.1415927/2-(acos(p_z_TCS/sqrt(p_y_TCS*p_y_TCS+p_z_TCS*p_z_TCS))+fabs(ang_fp));
          double p_z_fp = sqrt(p_y_TCS*p_y_TCS+p_z_TCS*p_z_TCS) * cos(p_thet_fp);
          double p_x_fp = sqrt(p_y_TCS*p_y_TCS+p_z_TCS*p_z_TCS) * sin(p_thet_fp);
          double p_phi_fp = asin(p_y_fp/sqrt(p_y_fp*p_y_fp+p_z_fp*p_z_fp));
          double y_fp    = x_TCS;
          double z_fp  = 0;
          double dist_fp = sqrt( y_fp*y_fp + (y_TCS-vdc1Y)*(y_TCS-vdc1Y) + (z_TCS-vdc1Z)*(z_TCS-vdc1Z) );  //y_TCS - is the x_TCS
          double x_fp    = sqrt( dist_fp * dist_fp - y_fp*y_fp) * (-1.)*fabs(y-vdc1Y)/(y-vdc1Y);
          cout<<"Focal plane :"<<" px="<<p_x_fp<<"   py="<<p_y_fp<<"  pz="<<p_z_fp<<"     x="<< x_fp <<"   y="<<y_fp<<"       theta="<<p_thet_fp<<"     phi="<<p_phi_fp<<endl;
//          double angg=acos(z_transport/sqrt(z_transport*z_transport+x_transport*x_transport))*fabs(x_transport)/x_transport-45.*deg;
//          output_vdc<<"  vdc->Fill("<<x_TrCS<<","<<y_TrCS<<");"<<"     //"<<x_vert<<"  "<<y_vert<<"   "<<z_vert<<"     "<<theStep->GetTrack()->GetVertexMomentumDirection().x()<<"   "<<theStep->GetTrack()->GetVertexMomentumDirection().y()<<"   "<<theStep->GetTrack()->GetVertexMomentumDirection().z()<<endl;
//          output_vdc<<x_TrCS<<"    "<<y_TrCS<<"                     "<<x_vert<<"  "<<y_vert<<"   "<<z_vert<<"                  "<<theStep->GetTrack()->GetVertexMomentumDirection().x()<<"   "<<theStep->GetTrack()->GetVertexMomentumDirection().y()<<"   "<<theStep->GetTrack()->GetVertexMomentumDirection().z()<<"               "<<px_step<<"   " <<py_step<<"   " << pz_step <<endl;
          std::ofstream output_fp;
          output_fp.open( "fp_coord.dat", ios::out | ios::app );
          output_fp<<x_fp<<setw(18)<<y_fp<<setw(18)<<p_thet_fp<<setw(18)<<p_phi_fp<<setw(18)<< p_x_fp <<setw(18)<< p_y_fp <<setw(18)<< p_z_fp << setw(18) << px_step << setw(18) << py_step << setw(18) << pz_step <<setw(18) << p_x_TCS << setw(18) << p_y_TCS << setw(18) << p_z_TCS <<setw(18)<<p_HCS<<"       0.0        0.0"<<setw(18)<<-1.*phi_vert<<setw(18)<<thet_vert<<endl;
          output_fp.close();
        }




        if (name2 == "vdc1_LH")
        {
//          HRSPrimaryGeneratorAction * root_ent_no;
//          root_ent_no = new HRSPrimaryGeneratorAction;
//          int root_evt_no = root_ent_no->root_nr(evtNb);
          int root_evt_no = 0;//root_ent_no->root_nr(evtNb);

          cout<<"stepping "<<root_evt_no<<endl;
//          output_vdc.open("vdc_coord.C", ios::out | ios::app );
          output_vdc.open("vdc_pairs_coord.txt", ios::out | ios::app );

//          double vdc1Z=2084.18;
//          double vdc1Y=740.142;
          double vdc1Z=2075.965;
          double vdc1Y=731.9258;
          double x0_vdc=vdc1Z*sin(mLHRSAngle);
          double z0_vdc=vdc1Z*cos(mLHRSAngle);
          double x_transport=x-x0_vdc;
          double z_transport=z-z0_vdc;
          double thet_dcs=acos(z/sqrt(x*x+z*z))-mLHRSAngle;
          
//          double Y_TCS=y;
//          double thet_tcs=acos(z/sqrt(x*x+z*z))-45.*deg;
//          cout<<"thet_tcs="<<thet_tcs<<endl;
//          double x_TCS=sqrt(x*x+z*z)*sin(thet_tcs)
          double e_kin = theStep->GetTrack()->GetVertexKineticEnergy()+0.511;
          double px_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().x();
          double py_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().y();
          double pz_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().z();


          double th_vtx_tmp=180./3.1415926*atan(px_vtx_tmp/pz_vtx_tmp);
          double ph_vtx_tmp=180./3.1415926*atan(py_vtx_tmp/sqrt(pz_vtx_tmp*pz_vtx_tmp+px_vtx_tmp*px_vtx_tmp));
          output_vdc<<"    "<<root_evt_no<<"  "<<evtNb<<"  11  "<<x_TCS<<"    "<<z_TCS-vdc1Z
                                 <<"        "<<x_vert<<"  "<<y_vert<<"   "<<z_vert<<"   "<<th_vtx_tmp<<"    "<<ph_vtx_tmp
                                 <<"        "<<px_vtx_tmp<<"   "<<py_vtx_tmp<<"   "<<pz_vtx_tmp
                                 <<"        "<<px_step<<"    " <<py_step<<"   " << pz_step << "   "<<e_kin <<endl;
          double x_TrCS=x_TCS;
          double z_TrCS=0;
          double dist_TrCS=sqrt( x_TrCS*x_TrCS + (y-vdc1Y)*(y-vdc1Y) + (z_TCS-vdc1Z)*(z_TCS-vdc1Z) );
          double y_TrCS=sqrt( dist_TrCS * dist_TrCS - x_TrCS*x_TrCS) * fabs(y-vdc1Y)/(y-vdc1Y);
          cout<<"left  transport_CS:"<<" p="<<p_HCS<<"   phi="<<p_phi/deg<<"   p_thet="<<p_thet<<"     x_trans="<<x_TrCS<<"      y_trans="<<y_TrCS<<"     z_trans="<<z_TrCS<<"     Vertex[cm]:("<<x_vert<<","<<y_vert<<","<<z_vert<<")  direction:"<<theStep->GetTrack()->GetVertexMomentumDirection()<<endl<<endl;
          double angg=acos(z_transport/sqrt(z_transport*z_transport+x_transport*x_transport))*fabs(x_transport)/x_transport-mLHRSAngle;
//          output_vdc<<"  vdc->Fill("<<x_TrCS<<","<<y_TrCS<<");"<<"     //"<<x_vert<<"  "<<y_vert<<"   "<<z_vert<<"     "<<px_vtx_tmp<<"   "<<py_vtx_tmp<<"   "<<pz_vtx_tmp<<endl;
          if (0)
          output_vdc<<evtNb<<"  11  "<<x_TrCS<<"    "<<y_TrCS<<"                     "<<x_vert<<"  "<<y_vert<<"   "<<z_vert
                                     <<"            "<<px_vtx_tmp<<"   "<<py_vtx_tmp<<"   "<<pz_vtx_tmp
                                     <<"            "<<px_step<<"   " <<py_step<<"   " << pz_step <<endl;
//          output_vdc <<endl;

          if ( ( fabs(x_TrCS)<14.4 ) && ( fabs(y_TrCS) < 0.5*211.8) )
          if ( p_HCS > 0.97 * ( 0.00051 + theStep->GetTrack()->GetVertexKineticEnergy()/GeV))
          {
            px_tmp_mass = px_vtx_tmp;
            py_tmp_mass = py_vtx_tmp;
            pz_tmp_mass = pz_vtx_tmp;
            e_tmp_mass  =  0.51 + theStep->GetTrack()->GetVertexKineticEnergy();

            P1x=px_tmp_mass;
  	    P1y=py_tmp_mass;
 	    P1z=pz_tmp_mass;
	    E1v=e_tmp_mass;
            
            
            cout<<px_tmp_mass<<"  "<<py_tmp_mass<<"  "<<pz_tmp_mass<<"  "<<e_tmp_mass<<endl;
            el_n_mass=root_evt_no;

          }
          output_vdc.close();
        }

        if (name2 == "vdc1_RH")
        {
//          HRSPrimaryGeneratorAction * root_ent_no;
//          root_ent_no = new HRSPrimaryGeneratorAction;
//          int root_evt_no = root_ent_no->root_nr(evtNb);
          int root_evt_no = 0;//root_ent_no->root_nr(evtNb);
          output_vdc.open("vdc_pairs_coord.txt", ios::out | ios::app );

//          double vdc1Z=2084.18;
//          double vdc1Y=740.142;
          double vdc1Z=2075.965;
          double vdc1Y=731.9258;
          double x0_vdc=vdc1Z*sin(mLHRSAngle);
          double z0_vdc=vdc1Z*cos(mLHRSAngle);
          double x_transport=x-x0_vdc;
          double z_transport=z-z0_vdc;
          double thet_dcs=acos(z/sqrt(x*x+z*z))-mLHRSAngle;
          
//          double Y_TCS=y;
//          double thet_tcs=acos(z/sqrt(x*x+z*z))-45.*deg;
//          cout<<"thet_tcs="<<thet_tcs<<endl;
//          double x_TCS=sqrt(x*x+z*z)*sin(thet_tcs)
          double e_kin = theStep->GetTrack()->GetVertexKineticEnergy()+0.511;
          double px_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().x();
          double py_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().y();
          double pz_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().z();


          double th_vtx_tmp=180./3.1415926*atan(px_vtx_tmp/pz_vtx_tmp);
          double ph_vtx_tmp=180./3.1415926*atan(py_vtx_tmp/sqrt(pz_vtx_tmp*pz_vtx_tmp+px_vtx_tmp*px_vtx_tmp));
          output_vdc<<"    "<<root_evt_no<<"  "<<evtNb<<"  -11  "<<x_TCS<<"    "<<z_TCS-vdc1Z
                                <<"         "<<x_vert<<"  "<<y_vert<<"   "<<z_vert<<"   "<<th_vtx_tmp<<"    "<<ph_vtx_tmp
                                <<"         "<<px_vtx_tmp<<"   "<<py_vtx_tmp<<"   "<<pz_vtx_tmp
                                <<"         "<<px_step<<"    " <<py_step<<"   " << pz_step << "   "<<e_kin << endl;
          double x_TrCS=x_TCS;
          double z_TrCS=0;
          double dist_TrCS=sqrt( x_TrCS*x_TrCS + (y-vdc1Y)*(y-vdc1Y) + (z_TCS-vdc1Z)*(z_TCS-vdc1Z) );
          double y_TrCS=sqrt( dist_TrCS * dist_TrCS - x_TrCS*x_TrCS) * fabs(y-vdc1Y)/(y-vdc1Y);
          cout<<"right transport_CS:"<<" p="<<p_HCS<<"   phi="<<p_phi/deg<<"   p_thet="<<p_thet<<"     x_trans="<<x_TrCS<<"      y_trans="<<y_TrCS<<"     z_trans="<<z_TrCS<<"     Vertex[cm]:("<<x_vert<<","<<y_vert<<","<<z_vert<<")  direction:"<<theStep->GetTrack()->GetVertexMomentumDirection()<<endl<<endl;
          double angg=acos(z_transport/sqrt(z_transport*z_transport+x_transport*x_transport))*fabs(x_transport)/x_transport-mLHRSAngle;
//          output_vdc<<"  vdc->Fill("<<x_TrCS<<","<<y_TrCS<<");"<<"     //"<<x_vert<<"  "<<y_vert<<"   "<<z_vert<<"     "<<px_vtx_tmp<<"   "<<py_vtx_tmp<<"   "<<pz_vtx_tmp<<endl;
          if (0)
          output_vdc<<evtNb<<"  -11  "<<x_TrCS<<"    "<<y_TrCS<<"                     "<<x_vert<<"  "<<y_vert<<"   "<<z_vert
                                      <<"            "<<px_vtx_tmp<<"   "<<py_vtx_tmp<<"   "<<pz_vtx_tmp
                                      <<"            "<<px_step<<"   " <<py_step<<"   " << pz_step <<endl;
//          output_vdc <<endl;

          if ( ( fabs(x_TrCS)<14.4 ) && ( fabs(y_TrCS) < 0.5*211.8) )
          if ( p_HCS > 0.97 * ( 0.00051 + theStep->GetTrack()->GetVertexKineticEnergy()/GeV))
          {
            px_tmp_mass = px_vtx_tmp;
            py_tmp_mass = py_vtx_tmp;
            pz_tmp_mass = pz_vtx_tmp;
            e_tmp_mass  = 0.51 + theStep->GetTrack()->GetVertexKineticEnergy();
            P2x=px_tmp_mass;
  	    P2y=py_tmp_mass;
 	    P2z=pz_tmp_mass;
	    E2v=e_tmp_mass;
            cout<<px_tmp_mass<<"  "<<py_tmp_mass<<"  "<<pz_tmp_mass<<"  "<<e_tmp_mass<<endl;
            pos_n_mass=root_evt_no;
          }
          output_vdc.close();
        }

        //check if evt is good
        {
          double y_dip_end=840.*(1.-sqrt(2.)/2.);
          if (y>y_dip_end)
          {
            if (1==2) cout<<"y is ok"<<endl;
            double z_dip_end=960.+840.*sqrt(2.)/2.;
            double z_tr_tmp = sqrt(x*x + z*z);
            if (z_tr_tmp>z_dip_end)
            {
              if (1==2) cout<<"z is also ok"<<endl;
              double p_fp=sqrt(px_step*px_step + pz_step*pz_step);
              if ((0.9*p_fp < py_step) && (0.9*p_fp < py_step))
              {
                if (1==2) cout<<"momentum check is also ok"<<endl;
                if (1==2) cout<<"trig is "<<write_trig<<endl;
                if (write_trig)
                {
                  write_trig=false;
                  std::ofstream output;
                  output.exceptions(std::ios_base::badbit | std::ios_base::failbit);
                  output.open("good_events", ios::out | ios::app );
                  output << "Event"<<evtNb<<"  name: "<<name<<"   Track:"<< theStep->GetTrack()->GetTrackID() <<"  "<<theStep->GetTrack()->GetDefinition()->GetParticleName()<< "  x=" <<x << ", y=" <<y << ", z=" <<z<<", px=" <<px_step<<", py=" <<py_step<<", pz=" << pz_step <<",    bx="<<tmpField[0]/tesla<<", by="<<tmpField[1]/tesla<<", bz="<<tmpField[2]/tesla<<"  Vertex[cm]:("<<x_vert<<","<<y_vert<<","<<z_vert<<")  direction:"<<theStep->GetTrack()->GetVertexMomentumDirection()<<std::endl;
                  output.close();
                }
              }
            }
          }
        }

	if (1==2)
	{
	  double sum_max=0;
	  cout<<"mtnum es?"<<endl;

  	  for (G4double xit=100.*cm; xit<=700.*cm; xit+=3.*cm)
	  for (G4double yit=0.*cm; yit<=500.*cm; yit+=3.*cm)
	  for (G4double zit=xit-10.*cm; zit<=xit+10.*cm; zit+=3.*cm)
	  {
	    double tmpPoint1[4]={xit,yit,zit,0};
	    theFieldManager->GetDetectorField()->GetFieldValue(tmpPoint1,tmpField);
//	    cout<<setw(15)<<xit<<setw(15)<<yit<<setw(15)<<zit<<setw(15)<<tmpField[0]<<setw(15)<<tmpField[1]<<setw(15)<<tmpField[2]<<setw(15)<<tmpField[3]<<setw(15)<<tmpField[4]<<setw(15)<<tmpField[5]<<"  yohoo"<<endl;
	    double sum=0;
	  
	    sum=sqrt(tmpField[0]*tmpField[0] + tmpField[1]*tmpField[1] + tmpField[2]*tmpField[2] + tmpField[3]*tmpField[3] + tmpField[4]*tmpField[4] + tmpField[5]*tmpField[5]);
	    if (sum>0.00001)
	    {
	      sum_max=sum;
//	    if (sum>=0.00001)
	      cout<<setw(15)<<xit<<setw(15)<<yit<<setw(15)<<zit<<setw(15)<<tmpField[0]/tesla<<setw(15)<<tmpField[1]/tesla<<setw(15)<<tmpField[2]/tesla<<setw(15)<<tmpField[3]<<setw(15)<<tmpField[4]<<setw(15)<<tmpField[5]<<"  yohoo"<<endl;
	    }
	  }
	}
}

//..............................................................................

void HRSSteppingAction::DoPrint(const G4Step* theStep)
{
}


//..............................................................................

void HRSSteppingAction::PrintHead(const G4Step* theStep,ostream& pOut)
{
}
//..............................................................................

void HRSSteppingAction::PrintStep(const G4Step* theStep,ostream& pOut)
{	
}
//..............................................................................

void HRSSteppingAction::FillRootArray(const G4Step* theStep)
{
}

