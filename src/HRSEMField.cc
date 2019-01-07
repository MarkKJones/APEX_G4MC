// ********************************************************************
//
// $Id: HRSEMField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//   User Field class Setup implementation.
//
//  
#include "HRSEMField.hh"
#include "HRSEMFieldMessenger.hh"
#include "UsageManager.hh"

#include "BField_Dipole_Fringe.hh"
#include "BField_Fringe_Q.hh"
#include "BField_Quad_Snake.hh"

extern UsageManager* gConfig;

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

HRSEMField::HRSEMField():mBField_Helm(0), mBField_Septum(0)
{

        gConfig->GetArgument("LHRSMomentum",mLHRSMomentum);
        gConfig->GetArgument("RHRSMomentum",mRHRSMomentum);
        gConfig->GetParameter("FringeField",mFringeField);
        gConfig->GetParameter("LHRSAngle",mLHRSAngle);
        mLHRSAngle*=deg;
        gConfig->GetParameter("RHRSAngle",mRHRSAngle);
        mRHRSAngle*=deg;
        gConfig->GetParameter("SeptumOn",Septum_On);
        gConfig->GetParameter("SeptumNew",SeptumNew);
        gConfig->GetParameter("SeptumFieldScale",SeptumFieldScale);
        cout<<"SeptumFieldScale="<<SeptumFieldScale<<endl;
        gConfig->GetParameter("Q1Sos",nQ1Sos);
        KAPPA1=0.243838 * tesla * (mLHRSMomentum/GeV)/0.83756;// / snakemagnumber / .1492 / m; //STD tune
        KAPPA2=-0.1934 * tesla * (mLHRSMomentum/GeV)/0.83756;// / snakemagnumber / .300  / m; //STD tune
        KAPPA3=0.17892 * tesla * (mLHRSMomentum/GeV)/0.83756;// / snakemagnumber / .300  / m; //STD tune
        DipField= -0.4192 * tesla * mLHRSMomentum / (1.063 * GeV);
        if (nQ1Sos)
        {
          double a_sup = 0.1492;                    //Bore radius, in meters
          double l_sup = 0.9413;                    //Length of quad, in meters
          double a_sos = 0.12827;
          double l_sos = 0.70;
          KAPPA1  *= l_sup / l_sos * a_sos / a_sup;
        }
        cout<<"     KAPPA1="<<KAPPA1*0.83756/mLHRSMomentum<<"     KAPPA2="<<KAPPA2*0.83756/mLHRSMomentum<<"     KAPPA3="<<KAPPA3*0.83756/mLHRSMomentum<<"     DipField="<<DipField<<endl;
/*
	messenger = new HRSEMFieldMessenger(this); 

	UsageManager* gConfig=UsageManager::GetUsageManager();

	//currently target field is always on
	//determine whether helm field is needed
	string pTargetFieldIni=gConfig->GetArgument("TargetFieldIni");
	string pTargetFieldMap=gConfig->GetArgument("TargetFieldMap");
	mBField_Helm = new BField_Helm(pTargetFieldIni.c_str() ,pTargetFieldMap.c_str());

	//determine whether septum field is needed
	//if HRS is built, septum will be loaded
	//There is a potential problem here: if one set SetupLHRS=SetupRHRS=0 in Detector.ini
	//and later on turn them on in G2P or CREX, this class will only response to Detector.ini
	//but not other config files, therefore the septum field will not setup properly
	int pSetupLHRS=0,pSetupRHRS=0;
	gConfig->GetParameter("SetupLHRS",pSetupLHRS); 
	gConfig->GetParameter("SetupRHRS",pSetupRHRS);
	if(pSetupLHRS || pSetupRHRS)
	{
		double pLHRSMomentum,pRHRSMomentum;
		gConfig->GetArgument("LHRSMomentum",pLHRSMomentum); 
		gConfig->GetArgument("RHRSMomentum",pRHRSMomentum);
		//cout << "LHRSMomentum and RHRSMomentum are: " << pLHRSMomentum << " and " << pRHRSMomentum << endl;
		string pSeptumFieldIni=gConfig->GetArgument("SeptumFieldIni");
		string pSeptumFieldMap=gConfig->GetArgument("SeptumFieldMap");	
		mBField_Septum = new BField_Septum(pLHRSMomentum,pRHRSMomentum,
						   pSeptumFieldIni.c_str(),pSeptumFieldMap.c_str());
		cout << "LHRSMomentum and RHRSMomentum are: " << pLHRSMomentum << " and " << pRHRSMomentum << " when making the septum" << endl;
	}

	ErDC = 0 *kilovolt/cm; 
	ErInner = 0 *kilovolt/cm; 

	bUseUniformEField=true;
	bUseUniformBField=false;
	EField3V.set(0,0,0);
	BField3V.set(0,0,0);
*/
        if (Septum_On)
        {
          if (SeptumNew==0)
            b_septum = new BField_Septum( mLHRSMomentum, mLHRSMomentum, "S80-map.table" );
          if (SeptumNew!=0)
            b_septum_new = new BField_Septum_New( mLHRSMomentum, mLHRSMomentum, "Septa-JB_map.table" );
        }
//        cout<<"HRSEMField construct"<<endl;

}

//////////////////////////////////////////////////////////////////////////
//
//  Deconstructors:
HRSEMField::~HRSEMField()
{
//	delete messenger;
        if (Septum_On)
        {
          delete b_septum;
          delete b_septum_new;
	}
	if(mBField_Helm)   delete mBField_Helm;
	if(mBField_Septum) delete mBField_Septum;
}


////////////////////////////////////////////////////////////////////////////
//input Point[4] (x,y,z,t) 
//
inline void HRSEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
//        cout<<"HRSEMField body"<<endl;
	//////////////////////////////////////////////////////////
	//get BField
	if(this->bUseUniformBField) 
	{
		Bfield[0]=BField3V.x();
		Bfield[1]=BField3V.y();
		Bfield[2]=BField3V.z();
	}
	else
	{
		double pB[3],pPos[3]={Point[0]/cm,Point[1]/cm,Point[2]/cm};  //turn into cm
		for(int i=0;i<3;i++) Bfield[i]=0.0;  //reset

		//target field read from map
		if(mBField_Helm)
		{
		if (1==2)
			for(int i=0;i<3;i++) pB[i]=0.0;  //reset
			if (! mBField_Helm->IsUniformField() )  mBField_Helm->GetBField(pPos,pB); 
			else  mBField_Helm->GetUniformField(pB); 
			for(int i=0;i<3;i++) Bfield[i]+=pB[i]*tesla;
		}


	}

	//////////////////////////////////////////////////////////
	//get EFiled,
	//g2p has no Efiled, 

	double  *Efield=&Bfield[3];
	if(this->bUseUniformEField) 
	{
//		Efield[0]=EField3V.x();
//		Efield[1]=EField3V.y();
//		Efield[2]=EField3V.z();
		for(int i=0;i<3;i++) Efield[i]=0.;
	}
	else
	{
		for(int i=0;i<3;i++) Efield[i]=0.;
	}

	{
		Bfield[0]=0;
		Bfield[1]=0;
		Bfield[2]=0;
	}

          if (Septum_On)
          {
            if (Point[2]>-40.*cm)
            if (Point[2]<170*cm)
            if (SeptumNew==0)
            {
              G4double pos_sept[3]={Point[0], Point[1], Point[2]};
              G4double B_sept[3]={0,0,0};
              b_septum->GetBField(pos_sept, B_sept);
              B_sept[0]*=2200./2178.432;
              B_sept[1]*=2200./2178.432;
              B_sept[2]*=2200./2178.432;
              Bfield[0]+=B_sept[0]*(mLHRSMomentum/GeV)/2.2;
              Bfield[1]+=B_sept[1]*(mLHRSMomentum/GeV)/2.2;
              Bfield[2]+=B_sept[2]*(mLHRSMomentum/GeV)/2.2;
              //VERTEX_Z=-110.0 cm;
            }

            if (Point[2]>-100.*cm)
            if (Point[2]<250*cm)
            if (SeptumNew==1)
            {

              G4double pos_sept[3]={Point[0], Point[1], Point[2]};
              G4double B_sept[3]={0,0,0};
              b_septum_new->GetBField(pos_sept, B_sept);
              B_sept[0]*=2.2/2.140045;
              B_sept[1]*=2.2/2.140045;
              B_sept[2]*=2.2/2.140045;
//              Bfield[0]+=B_sept[0]*(mLHRSMomentum/GeV)/2.2;
//              Bfield[1]+=B_sept[1]*(mLHRSMomentum/GeV)/2.2;
//              Bfield[2]+=B_sept[2]*(mLHRSMomentum/GeV)/2.2;
              Bfield[0]+=SeptumFieldScale*B_sept[0]*(mLHRSMomentum/GeV)/2.2;
              Bfield[1]+=SeptumFieldScale*B_sept[1]*(mLHRSMomentum/GeV)/2.2;
              Bfield[2]+=SeptumFieldScale*B_sept[2]*(mLHRSMomentum/GeV)/2.2;
//              cout<<"septum field: "<<Point[0]<<"      "<<Point[1]<<"       "<<Point[2]<<"       "<<SeptumFieldScale*B_sept[0]*(mLHRSMomentum/GeV)/2.2<<"     "<<SeptumFieldScale*B_sept[1]*(mLHRSMomentum/GeV)/2.2<<"      "<<SeptumFieldScale*B_sept[2]*(mLHRSMomentum/GeV)/2.2<<endl;
              //VERTEX_Z=-108.1365;
            }
          }





	  if (mFringeField==1)
//	  if (0)
	  {
            bool debug = true;
            debug=false;
            double pTarget   =             0.0  * cm;
            double pQ1en     = pTarget + 160.0  * cm;
            double pQ1Length = 94.13*cm;//SNAKE
            double pQ1ex     = pQ1en   +  pQ1Length;
            double pQ2en     = pQ1ex   + 115.58 * cm;
            double pQ1Radius = 0.1492*m;//12.827 * cm;
            double pQ2Length = 182.66*cm;//SNAKE
            double pQ2Radius =  0.3 * m;//SNAKE
            double pQ2ex     = pQ2en   + pQ2Length;
            double pQ3Radius = 0.3 * m;//SNAKE
            double pQ3Length = 182.68*cm;//SNAKE

            if (nQ1Sos)
            {
              pQ1en   = pTarget + 171.1*cm;
              pQ1Radius = 12.827 * cm;
              pQ1ex   = pQ1en   +  70.0 * cm;
            }


            double pLQ3_Fr_EnX = (17.0267042 * m  ) *  sin( mLHRSAngle );//SNAKE
            double pLQ3_Fr_EnZ = (17.0267042 * m  ) *  cos( mLHRSAngle );//SNAKE
            double pLQ3_Fr_EnY = ( 3.58637   * m  );//SNAKE	
            double pLQ3_Fr_ExX = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  sin( mLHRSAngle );//SNAKE
            double pLQ3_Fr_ExZ = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  cos( mLHRSAngle );//SNAKE
            double pLQ3_Fr_ExY = ( 3.58637   * m  + pQ3Length / sqrt(2.) );//SNAKE	

            double pRQ3_Fr_EnX = (17.0267042 * m  ) *  sin( mRHRSAngle );//SNAKE
            double pRQ3_Fr_EnZ = (17.0267042 * m  ) *  cos( mRHRSAngle );//SNAKE
            double pRQ3_Fr_EnY = ( 3.58637   * m  );//SNAKE	
            double pRQ3_Fr_ExX = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  sin( mRHRSAngle );//SNAKE
            double pRQ3_Fr_ExZ = (17.0267042 * m  + pQ3Length / sqrt(2.) ) *  cos( mRHRSAngle );//SNAKE
            double pRQ3_Fr_ExY = ( 3.58637   * m  + pQ3Length / sqrt(2.) );//SNAKE	

//            cout<<"noric tes "<<pLQ3_Fr_EnX<<"    "<<pLQ3_Fr_EnY<<"     "<<pLQ3_Fr_EnZ<<endl;
//            cout<<"noric tes "<<0.5*(pLQ3_Fr_EnX+pLQ3_Fr_ExX)<<"    "<<0.5*(pLQ3_Fr_EnY+pLQ3_Fr_ExY)<<"     "<<0.5*(pLQ3_Fr_EnZ+pLQ3_Fr_ExZ)<<endl;
//            cout<<"noric tes "<<pLQ3_Fr_ExX<<"    "<<pLQ3_Fr_ExY<<"     "<<pLQ3_Fr_ExZ<<endl<<endl;
//            q3 cent coord = 3825.05       4232.24      17253.7
//            G4ThreeVector     LORIGINQ3(0.5*(pLQ3_Fr_ExX+pLQ3_Fr_EnX), 0.5*(pLQ3_Fr_ExY+pLQ3_Fr_EnY), 0.5*(pLQ3_Fr_ExZ+pLQ3_Fr_EnZ));

//            double KAPPA2    = -0.1939*tesla*mLHRSMomentum/GeV;
//            double KAPPA1    = 0.24450*tesla*mLHRSMomentum/GeV;
//            double KAPPA3    = -0.1794 * tesla * mLHRSMomentum/GeV;

	    if (Point[0] > 0)
	    //if (0)
	    {
      	      // q1fringe enter
	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ1en+0*cm))
	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ1en-100*cm))  //60 cm is the fringe field length: 50cm+10cm 
	      {

                G4ThreeVector     LORIGINQ1(pQ1en * sin(mLHRSAngle), 0., pQ1en * cos(mLHRSAngle));
                G4RotationMatrix* LROTATEQ1 = new G4RotationMatrix;
                LROTATEQ1->rotateY( mLHRSAngle);

                BField_Fringe_Q  * fMagFieldFringeQ1;
                fMagFieldFringeQ1 = new BField_Fringe_Q(1, KAPPA1, pQ1Radius, LORIGINQ1, LROTATEQ1, 1);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad en "<<endl;
                fMagFieldFringeQ1->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
     	          Bfield[0]+=yBBB[0];
  	          Bfield[1]+=yBBB[1];
  	          Bfield[2]+=yBBB[2];
                  if (debug)
  	            cout<<"Q1 entrance   "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;

  	        }
  	        //delete LROTATEQ1;
  	        delete fMagFieldFringeQ1;
	      }







              // q1 inside  exit
              if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ1en+0*cm))
              if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ1ex+0*cm))  //60 cm is the fringe field length: 50cm+10cm 
              {
                double q1_cent=(pQ1en+pQ1ex)/2.;

                G4ThreeVector     LORIGINQ1(q1_cent * sin(mLHRSAngle), 0., q1_cent * cos(mLHRSAngle));
                G4RotationMatrix* LROTATEQ1 = new G4RotationMatrix;
                LROTATEQ1->rotateY( mLHRSAngle);

                BField_Quad_Snake  * fMagFieldQ1_sn;
                fMagFieldQ1_sn = new BField_Quad_Snake(pQ1ex-pQ1en, KAPPA1, pQ1Radius, LORIGINQ1, LROTATEQ1, 1);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad en "<<endl;
                fMagFieldQ1_sn->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
                  Bfield[0]+=yBBB[0];
                  Bfield[1]+=yBBB[1];
                  Bfield[2]+=yBBB[2];
                  if (debug)
                    cout<<"Q1 inside   "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
                }
                //delete LROTATEQ1;
                delete fMagFieldQ1_sn;
              }








  	      // q1fringe exit

  	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ1ex-0*cm) )
  	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ1ex+100*cm) )  //60 cm is the fringe field length: 50cm+10cm 
  	      //if (0)
  	      {
                G4ThreeVector     LORIGINQ1(pQ1ex * sin(mLHRSAngle), 0., pQ1ex * cos(mLHRSAngle));
                G4RotationMatrix* LROTATEQ1 = new G4RotationMatrix;
                LROTATEQ1->rotateY( mLHRSAngle);

                BField_Fringe_Q  * fMagFieldFringeQ1;
                fMagFieldFringeQ1 = new BField_Fringe_Q(3, KAPPA1, pQ1Radius, LORIGINQ1, LROTATEQ1, 1);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad q1 ex "<<endl;
                fMagFieldFringeQ1->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe q1 ex  "<<endl<<endl<<endl;
//                    cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
     	          Bfield[0]+=yBBB[0];
    	          Bfield[1]+=yBBB[1];
  	          Bfield[2]+=yBBB[2];
                  if (debug)
    	            cout<<"Q1 exit, Distance "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
  	        }
  	        //delete LROTATEQ1;
  	        delete fMagFieldFringeQ1;
	      }
	    }


	    if (Point[0] < 0.)
	    //if (0)
	    {
      	      // q1fringe enter
	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ1en+0*cm))
	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ1en-100*cm))  //60 cm is the fringe field length: 50cm+10cm 
	      {

                G4ThreeVector     RORIGINQ1(pQ1en * sin(mRHRSAngle), 0., pQ1en * cos(mRHRSAngle));
                G4RotationMatrix* RROTATEQ1 = new G4RotationMatrix;
                RROTATEQ1->rotateY( mRHRSAngle);

                BField_Fringe_Q  * fMagFieldFringeQ1;
                fMagFieldFringeQ1 = new BField_Fringe_Q(1, KAPPA1, pQ1Radius, RORIGINQ1, RROTATEQ1, 1);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad en "<<endl;
                fMagFieldFringeQ1->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
     	          Bfield[0]+=-1.*yBBB[0];
  	          Bfield[1]+=-1.*yBBB[1];
  	          Bfield[2]+=-1.*yBBB[2];
  	        }
  	        //delete RROTATEQ1;
  	        delete fMagFieldFringeQ1;
	      }

  	      // q1fringe exit

  	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ1ex-0*cm) )
  	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ1ex+100*cm) )  //60 cm is the fringe field length: 50cm+10cm 
  	      //if (0)
  	      {
                G4ThreeVector     RORIGINQ1(pQ1ex * sin(mRHRSAngle), 0., pQ1ex * cos(mRHRSAngle));
                G4RotationMatrix* RROTATEQ1 = new G4RotationMatrix;
                RROTATEQ1->rotateY( mRHRSAngle);

                BField_Fringe_Q  * fMagFieldFringeQ1;
                fMagFieldFringeQ1 = new BField_Fringe_Q(3, KAPPA1, pQ1Radius, RORIGINQ1, RROTATEQ1, 1);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad q1 ex "<<endl;
                fMagFieldFringeQ1->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe q1 ex  "<<endl<<endl<<endl;
//                    cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
     	          Bfield[0]+=-1.*yBBB[0];
    	          Bfield[1]+=-1.*yBBB[1];
  	          Bfield[2]+=-1.*yBBB[2];
  	        }
  	        //delete RROTATEQ1;
  	        delete fMagFieldFringeQ1;
	      }
	    }





	    // q2fringe enter

	    if (Point[0] > 0)
	    //if (0)
	    {
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ2en+0*cm) )
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ2en-100*cm) )  //60 cm is the fringe field length: 50cm+10cm 
	     //if (0)
	     {

              G4ThreeVector     LORIGINQ2(pQ2en * sin(mLHRSAngle), 0., pQ2en * cos(mLHRSAngle));
              G4RotationMatrix* LROTATEQ2 = new G4RotationMatrix;
              LROTATEQ2->rotateY( mLHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ2;
              fMagFieldFringeQ2 = new BField_Fringe_Q(1, KAPPA2, pQ2Radius, LORIGINQ2, LROTATEQ2, 2);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad2 en "<<endl;
              fMagFieldFringeQ2->GetFieldValue(yyyy, yBBB) ;
//              cout<<"   end fringe2 en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
     	        Bfield[0]+=yBBB[0];
  	        Bfield[1]+=yBBB[1];
  	        Bfield[2]+=yBBB[2];
                if (debug)
                  cout<<"Q2 entrance, Distance "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
  	      }
    	      //delete LROTATEQ2;
  	      delete fMagFieldFringeQ2;
	     }





              // q2 inside  exit
              if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ2en+0*cm))
              if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ2ex+0*cm))  //60 cm is the fringe field length: 50cm+10cm 
              {
                double q2_cent=(pQ2en+pQ2ex)/2.;

                G4ThreeVector     LORIGINQ2(q2_cent * sin(mLHRSAngle), 0., q2_cent * cos(mLHRSAngle));
                G4RotationMatrix* LROTATEQ2 = new G4RotationMatrix;
                LROTATEQ2->rotateY( mLHRSAngle);

                BField_Quad_Snake  * fMagFieldQ2_sn;
                fMagFieldQ2_sn = new BField_Quad_Snake(pQ2ex-pQ2en, KAPPA2, pQ2Radius, LORIGINQ2, LROTATEQ2, 2);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad en "<<endl;
                fMagFieldQ2_sn->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
                  Bfield[0]+=yBBB[0];
                  Bfield[1]+=yBBB[1];
                  Bfield[2]+=yBBB[2];
                  if (debug)
                    cout<<"Q2 inside   "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
                }
                //delete LROTATEQ1;
                delete fMagFieldQ2_sn;
              }





	     // q2fringe exit

	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ2ex-0*cm) )
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ2ex+100*cm) )  //60 cm is the fringe field length: 50cm+10cm 
//	     if (0)
	     {
              G4ThreeVector     LORIGINQ2(pQ2ex * sin(mLHRSAngle), 0., pQ2ex * cos(mLHRSAngle));
              G4RotationMatrix* LROTATEQ2 = new G4RotationMatrix;
              LROTATEQ2->rotateY( mLHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ2;
              fMagFieldFringeQ2 = new BField_Fringe_Q(3, KAPPA2, pQ2Radius, LORIGINQ2, LROTATEQ2, 2);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad2 ex "<<endl;
              fMagFieldFringeQ2->GetFieldValue(yyyy, yBBB) ;
//              cout<<"   end fringe ex  "<<endl<<endl<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
//                cout<<"   xyz=("<<Point[0]/cm<<","<<Point[1]/cm<<","<<Point[2]/cm<<")"<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                Bfield[0]+=yBBB[0];
   	        Bfield[1]+=yBBB[1];
   	        Bfield[2]+=yBBB[2];
                if (debug)
                  cout<<"Q2 exit  "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
   	      }
    	      //delete LROTATEQ2;
  	      delete fMagFieldFringeQ2;
  	     }
	    }


	    if (Point[0] < 0)
	    //if (0)
	    {
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ2en+0*cm) )
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ2en-100*cm) )  //60 cm is the fringe field length: 50cm+10cm 
	     {

              G4ThreeVector     RORIGINQ2(pQ2en * sin(mRHRSAngle), 0., pQ2en * cos(mRHRSAngle));
              G4RotationMatrix* RROTATEQ2 = new G4RotationMatrix;
              RROTATEQ2->rotateY( mRHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ2;
              fMagFieldFringeQ2 = new BField_Fringe_Q(1, KAPPA2, pQ2Radius, RORIGINQ2, RROTATEQ2, 2);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad2 en "<<endl;
              fMagFieldFringeQ2->GetFieldValue(yyyy, yBBB) ;
//              cout<<"   end fringe2 en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
     	        Bfield[0]+=-1.*yBBB[0];
  	        Bfield[1]+=-1.*yBBB[1];
  	        Bfield[2]+=-1.*yBBB[2];
  	      }
    	      //delete RROTATEQ2;
  	      delete fMagFieldFringeQ2;
	     }

	     // q2fringe exit

	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > (pQ2ex-0*cm) )
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) < (pQ2ex+100*cm) )  //60 cm is the fringe field length: 50cm+10cm 
	     {
              G4ThreeVector     RORIGINQ2(pQ2ex * sin(mRHRSAngle), 0., pQ2ex * cos(mRHRSAngle));
              G4RotationMatrix* RROTATEQ2 = new G4RotationMatrix;
              RROTATEQ2->rotateY( mRHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ2;
              fMagFieldFringeQ2 = new BField_Fringe_Q(3, KAPPA2, pQ2Radius, RORIGINQ2, RROTATEQ2, 2);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad2 ex "<<endl;
              fMagFieldFringeQ2->GetFieldValue(yyyy, yBBB) ;
//              cout<<"   end fringe ex  "<<endl<<endl<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
//                cout<<"   xyz=("<<Point[0]/cm<<","<<Point[1]/cm<<","<<Point[2]/cm<<")"<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                Bfield[0]+=-1.*yBBB[0];
   	        Bfield[1]+=-1.*yBBB[1];
   	        Bfield[2]+=-1.*yBBB[2];
   	      }
    	      //delete RROTATEQ2;
  	      delete fMagFieldFringeQ2;
  	     }
	    }








	    if (Point[0] > 0)
	    {
              //dipole fringe
              double rn = sqrt(Point[0]*Point[0] + Point[2]*Point[2]);
              if (rn > 8.*m)
              {
                  G4double dipoleField = DipField;
                  G4RotationMatrix* LROTATED = new G4RotationMatrix;
                  LROTATED->rotateY( mLHRSAngle );
                  double pDipoleRCenterY=8.4  * m;
                  double pDipoleRCenterZ=9.961 * m;//SNAKE
/*
{
                              G4ThreeVector     LORIGIND(pDipoleRCenterZ * sin( mLHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mLHRSAngle ));
                              BField_Dipole_Fringe *    fMagFieldFZBL3_fringe ;
                              fMagFieldFZBL3_fringe = new BField_Dipole_Fringe( dipoleField, LORIGIND, LROTATED );
                              G4double coord_fr[7]={Point[0],Point[1],Point[2],0,0,0,0};
                              G4double bdasht_fr[3]={0,0,0};
                              fMagFieldFZBL3_fringe->GetFieldValue(coord_fr, bdasht_fr) ;
}
*/

                  {
	            //double rn = sqrt(Point[0]*Point[0] + Point[2]*Point[2]);
	            double thet=asin(Point[0]/rn);
	            double thet_new=thet-mLHRSAngle;
	            double xn=sin(thet_new)*rn;
	            double yn=Point[1];
	            double zn=cos(thet_new)*rn;
	            {
	              double x_width=0.30*m;
	              double z0=9.961*m;
	              double y0=8.4*m;
	              double y1=y0*(1-sin(45.*deg));     //dipole end plane, 45deg is the bending angle 
	              double z1=z0+y0*sin(45.*deg);    //dipole end plane, 45deg is the bending angle 

	              if (fabs(xn) < 0.5*x_width)    //dipole x is ok
	              {
	                if ( (zn) < (z0 + (yn)*tan(30*deg)) )   //front plane
	                {
//                          G4double dipoleField = DipField;
	                  double R_new= sqrt( (y0-yn)*(y0-yn) + (zn - z0) * (zn - z0) );
//  	                  if (fabs(R_new-y0)/m<(0.4))    //tunnel radius
  	                  if (fabs(Point[1])/m < 0.477 )    //tunnel radius
  	                  {
  	                    {
                              G4ThreeVector     LORIGIND(pDipoleRCenterZ * sin( mLHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mLHRSAngle ));
                              BField_Dipole_Fringe *    fMagFieldFZBL3_fringe ;
                              fMagFieldFZBL3_fringe = new BField_Dipole_Fringe( dipoleField, LORIGIND, LROTATED);
                              G4double coord_fr[7]={Point[0],Point[1],Point[2],0,0,0,0};
                              G4double bdasht_fr[3]={0,0,0};
                              fMagFieldFZBL3_fringe->GetFieldValue(coord_fr, bdasht_fr) ;
//                              cout<<"Mtnum es che stex?   point:"<<Point[0]<<", "<< Point[1] <<", "<< Point[2] <<endl;
                              if (bdasht_fr[0]>-1000.)
                              if (bdasht_fr[1]>-1000.)
                              if (bdasht_fr[2]>-1000.)
                              if (bdasht_fr[0]<1000.)
                              if (bdasht_fr[1]<1000.)
                              if (bdasht_fr[2]<1000.)
                              {
                                Bfield[0]+=bdasht_fr[0];
     	                        Bfield[1]+=bdasht_fr[1];
                       	        Bfield[2]+=bdasht_fr[2];
  	                      }
                              delete fMagFieldFZBL3_fringe;
  	                    }
                          }
	                }
                        double d1=sqrt((zn-z1)*(zn-z1) + (yn-y1)*(yn-y1));
                        if ( (d1 < 1.*m) && fabs(Point[1]-(y1+(zn-z1)))<47.7*cm*sqrt(2.) )
  	                if ( yn > ( (z1 - zn ) * tan (15.*deg) + y1) )
  	                {
                            G4ThreeVector     LORIGIND(pDipoleRCenterZ * sin( mLHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mLHRSAngle ));
                            BField_Dipole_Fringe *    fMagFieldFZBL3_fringe ;
                            fMagFieldFZBL3_fringe = new BField_Dipole_Fringe( dipoleField, LORIGIND, LROTATED);
                            G4double coord_fr[7]={Point[0],Point[1],Point[2],0,0,0,0};
                            G4double bdasht_fr[3]={0,0,0};
                            fMagFieldFZBL3_fringe->GetFieldValue(coord_fr, bdasht_fr) ;

//                            cout<<"-Mtnum es che stex?   point:"<<Point[0]<<", "<< Point[1] <<", "<< Point[2] <<endl;
                            if (bdasht_fr[0]>-1000.)
                            if (bdasht_fr[1]>-1000.)
                            if (bdasht_fr[2]>-1000.)
                            if (bdasht_fr[0]<1000.)
                            if (bdasht_fr[1]<1000.)
                            if (bdasht_fr[2]<1000.)
                            {
                              Bfield[0]=+bdasht_fr[0];
   	                      Bfield[1]=+bdasht_fr[1];
                     	      Bfield[2]=+bdasht_fr[2];
                     	    }
                     	    delete fMagFieldFZBL3_fringe;
                        }
                      }
	            }
                  }
                  //delete LROTATED;
              }
            }



	    if (Point[0] < 0)
	    {
              //dipole fringe
              double rn = sqrt(Point[0]*Point[0] + Point[2]*Point[2]);
              if (rn > 8.*m)
              {
                  G4double dipoleField = DipField;
                  G4RotationMatrix* RROTATED = new G4RotationMatrix;
                  RROTATED->rotateY( mRHRSAngle );
                  double pDipoleRCenterY=8.4  * m;
                  double pDipoleRCenterZ=9.961 * m;//SNAKE

                  {
	            //double rn = sqrt(Point[0]*Point[0] + Point[2]*Point[2]);
	            double thet=asin(Point[0]/rn);
	            double thet_new=thet-mRHRSAngle;
	            double xn=sin(thet_new)*rn;
	            double yn=Point[1];
	            double zn=cos(thet_new)*rn;
	            {
	              double x_width=0.30*m;
	              double z0=9.961*m;
	              double y0=8.4*m;
	              double y1=y0*(1-sin(45.*deg));     //dipole end plane, 45deg is the bending angle 
	              double z1=z0+y0*sin(45.*deg);    //dipole end plane, 45deg is the bending angle 

	              if (fabs(xn) < 0.5*x_width)    //dipole x is ok
	              {
	                if ( (zn) < (z0 + (yn)*tan(30*deg)) )   //front plane
	                {
//                          G4double dipoleField = DipField;
	                  double R_new= sqrt( (y0-yn)*(y0-yn) + (zn - z0) * (zn - z0) );
//  	                  if (fabs(R_new-y0)/m<(0.4))    //tunnel radius
  	                  if (fabs(Point[1])/m < 0.477 )    //tunnel radius
  	                  {
  	                    {
                              G4ThreeVector     RORIGIND(pDipoleRCenterZ * sin( mRHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mRHRSAngle ));
                              BField_Dipole_Fringe *    fMagFieldFZBR3_fringe ;
                              fMagFieldFZBR3_fringe = new BField_Dipole_Fringe( dipoleField, RORIGIND, RROTATED);
                              G4double coord_fr[7]={Point[0],Point[1],Point[2],0,0,0,0};
                              G4double bdasht_fr[3]={0,0,0};
                              fMagFieldFZBR3_fringe->GetFieldValue(coord_fr, bdasht_fr) ;
//                              cout<<"Mtnum es che stex?   point:"<<Point[0]<<", "<< Point[1] <<", "<< Point[2] <<endl;
                              if (bdasht_fr[0]>-1000.)
                              if (bdasht_fr[1]>-1000.)
                              if (bdasht_fr[2]>-1000.)
                              if (bdasht_fr[0]<1000.)
                              if (bdasht_fr[1]<1000.)
                              if (bdasht_fr[2]<1000.)
                              {
                                Bfield[0]=-1.*bdasht_fr[0];
     	                        Bfield[1]=-1.*bdasht_fr[1];
                       	        Bfield[2]=-1.*bdasht_fr[2];
  	                      }
  	                      delete fMagFieldFZBR3_fringe;
  	                    }
                          }
	                }
                        double d1=sqrt((zn-z1)*(zn-z1) + (yn-y1)*(yn-y1));
                        if ( (d1 < 1.*m) && fabs(Point[1]-(y1+(zn-z1)))<47.7*cm*sqrt(2.) )
  	                if ( yn > ( (z1 - zn ) * tan (15.*deg) + y1) )
  	                {
                            G4ThreeVector     RORIGIND(pDipoleRCenterZ * sin( mRHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mRHRSAngle ));
                            BField_Dipole_Fringe *    fMagFieldFZBR3_fringe ;
                            fMagFieldFZBR3_fringe = new BField_Dipole_Fringe( dipoleField, RORIGIND, RROTATED);
                            G4double coord_fr[7]={Point[0],Point[1],Point[2],0,0,0,0};
                            G4double bdasht_fr[3]={0,0,0};
                            fMagFieldFZBR3_fringe->GetFieldValue(coord_fr, bdasht_fr) ;

//                            cout<<"+Mtnum es che stex?   point:"<<Point[0]<<", "<< Point[1] <<", "<< Point[2] <<endl;
                            if (bdasht_fr[0]>-1000.)
                            if (bdasht_fr[1]>-1000.)
                            if (bdasht_fr[2]>-1000.)
                            if (bdasht_fr[0]<1000.)
                            if (bdasht_fr[1]<1000.)
                            if (bdasht_fr[2]<1000.)
                            {
                              Bfield[0]+=-1.*bdasht_fr[0];
   	                      Bfield[1]+=-1.*bdasht_fr[1];
                     	      Bfield[2]+=-1.*bdasht_fr[2];
                     	    }
                     	    delete fMagFieldFZBR3_fringe;
                        }
                      }
	            }
                  }
                  //delete RROTATED;
              }
            }








	    if (Point[0] > 0)
	    //if (0)
	    if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > 15.*m )
	    {
              double z0_en_cent=17.0267042 * m;
              double y0_en_cent=3.58637 * m;
              double z0_ex_cent=z0_en_cent + pQ3Length / sqrt(2.);
              double y0_ex_cent=y0_en_cent + pQ3Length / sqrt(2.);
              double z_cur_coord=sqrt(Point[0]*Point[0]+Point[2]*Point[2]);
              double y_cur_coord=fabs(Point[1]);
//              if (fabs(y_cur_coord) > (-1.*(z_cur_coord-z0_ex_cent)+y0_ex_cent))    cout<<"q3 exit"<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<endl;

	     // Q3fringe enter
             if (fabs(y_cur_coord) < (-1.*(z_cur_coord-z0_en_cent)+y0_en_cent)+5.*cm)
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > 15.*m )
	     {
              G4ThreeVector     LORIGINQ3_Fr_En(pLQ3_Fr_EnX, pLQ3_Fr_EnY, pLQ3_Fr_EnZ);
              G4RotationMatrix* LROTATEQ3_Fr_En = new G4RotationMatrix;
              LROTATEQ3_Fr_En->rotateX(-45.0 * deg);
              LROTATEQ3_Fr_En->rotateY( mLHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ3;
              fMagFieldFringeQ3 = new BField_Fringe_Q(1, KAPPA3, pQ3Radius, LORIGINQ3_Fr_En, LROTATEQ3_Fr_En, 3);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad en "<<endl;
              fMagFieldFringeQ3->GetFieldValue(yyyy, yBBB) ;
//              cout<<"   end fringe en "<<endl<<endl<<endl;
//              cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
     	        Bfield[0]+=yBBB[0];
  	        Bfield[1]+=yBBB[1];
  	        Bfield[2]+=yBBB[2];
                if (debug)
                  cout<<"Q3 entrance   "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
	      }
    	      //delete LROTATEQ3_Fr_En;
  	      delete fMagFieldFringeQ3;
	    }





              // q3 inside
              if (fabs(y_cur_coord) < (-1.*(z_cur_coord-z0_ex_cent)+y0_ex_cent)+5.*cm)
              if (fabs(y_cur_coord) > (-1.*(z_cur_coord-z0_en_cent)+y0_en_cent)-5.*cm)
 	      if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > 16.*m )
              {
                G4ThreeVector     LORIGINQ3(0.5*(pLQ3_Fr_ExX+pLQ3_Fr_EnX), 0.5*(pLQ3_Fr_ExY+pLQ3_Fr_EnY), 0.5*(pLQ3_Fr_ExZ+pLQ3_Fr_EnZ));
//                cout<<"     q3 cent coord = "<<0.5*(pLQ3_Fr_ExX+pLQ3_Fr_EnX)<<"       "<<0.5*(pLQ3_Fr_ExY+pLQ3_Fr_EnY)<<"      "<<0.5*(pLQ3_Fr_ExZ+pLQ3_Fr_EnZ)<<endl;
                G4RotationMatrix* LROTATEQ3 = new G4RotationMatrix;
                LROTATEQ3->rotateX(-45.0 * deg);
                LROTATEQ3->rotateY( mLHRSAngle);

                BField_Quad_Snake  * fMagFieldQ3_sn;
                fMagFieldQ3_sn = new BField_Quad_Snake(pQ3Length , KAPPA3, pQ3Radius, LORIGINQ3, LROTATEQ3, 3);

                G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
                G4double yBBB[3]={0,0,0};
//                cout<<"   finge  quad en "<<endl;
                fMagFieldQ3_sn->GetFieldValue(yyyy, yBBB) ;
//                cout<<"   end fringe en "<<endl<<endl<<endl;
//                  cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
                if (yBBB[0]>-1000.)
                if (yBBB[1]>-1000.)
                if (yBBB[2]>-1000.)
                if (yBBB[0]<1000.)
                if (yBBB[1]<1000.)
                if (yBBB[2]<1000.)
                {
                  Bfield[0]+=yBBB[0];
                  Bfield[1]+=yBBB[1];
                  Bfield[2]+=yBBB[2];
                  if (debug)
                    cout<<"Q3 inside   "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
                }
                //delete LROTATEQ1;
                delete fMagFieldQ3_sn;
              }





	    // Q3fringe exit
            if (fabs(y_cur_coord) > (-1.*(z_cur_coord-z0_ex_cent)+y0_ex_cent)-5.*cm)
	    if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > 16.*m )
	    {

              G4ThreeVector     LORIGINQ3_Fr_Ex(pLQ3_Fr_ExX, pLQ3_Fr_ExY, pLQ3_Fr_ExZ);
              G4RotationMatrix* LROTATEQ3_Fr_Ex = new G4RotationMatrix;
              LROTATEQ3_Fr_Ex->rotateX(-45.0 * deg);
              LROTATEQ3_Fr_Ex->rotateY( mLHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ3;
              fMagFieldFringeQ3 = new BField_Fringe_Q(3, KAPPA3, pQ3Radius, LORIGINQ3_Fr_Ex, LROTATEQ3_Fr_Ex, 3);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad ex "<<endl;
              fMagFieldFringeQ3->GetFieldValue(yyyy, yBBB);
//              cout<<"   end fringe ex  "<<endl<<endl<<endl;
//              cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
     	        Bfield[0]+=yBBB[0];
  	        Bfield[1]+=yBBB[1];
  	        Bfield[2]+=yBBB[2];
                if (debug)
                  cout<<"Q3 exit   "<<Point[0]/cm<<"      "<<Point[1]/cm<<"       "<<Point[2]/cm<<"       "<<yBBB[0]/tesla<<"       "<<yBBB[0]<<"     "<<yBBB[1]<<"      "<<yBBB[2]<<endl;
  	      }
    	      //delete LROTATEQ3_Fr_Ex;
  	      delete fMagFieldFringeQ3;
	     }
	    }






	    if (Point[0] < 0)
	    //if (0)
	    {
	     // Q3fringe enter
	     if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > 16.*m )
	     {
              G4ThreeVector     RORIGINQ3_Fr_En(pRQ3_Fr_EnX, pRQ3_Fr_EnY, pRQ3_Fr_EnZ);
              G4RotationMatrix* RROTATEQ3_Fr_En = new G4RotationMatrix;
              RROTATEQ3_Fr_En->rotateX(-45.0 * deg);
              RROTATEQ3_Fr_En->rotateY( mRHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ3;
              fMagFieldFringeQ3 = new BField_Fringe_Q(1, KAPPA3, pQ3Radius, RORIGINQ3_Fr_En, RROTATEQ3_Fr_En, 3);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad en "<<endl;
              fMagFieldFringeQ3->GetFieldValue(yyyy, yBBB) ;
//              cout<<"   end fringe en "<<endl<<endl<<endl;
//              cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
     	        Bfield[0]+=-1.*yBBB[0];
  	        Bfield[1]+=-1.*yBBB[1];
  	        Bfield[2]+=-1.*yBBB[2];
	      }
    	      //delete RROTATEQ3_Fr_En;
  	      delete fMagFieldFringeQ3;
	    }

	    // Q3fringe exit
	    if ( sqrt(Point[0]*Point[0]+Point[2]*Point[2]) > 16.*m )
	    {

              G4ThreeVector     RORIGINQ3_Fr_Ex(pRQ3_Fr_ExX, pRQ3_Fr_ExY, pRQ3_Fr_ExZ);
              G4RotationMatrix* RROTATEQ3_Fr_Ex = new G4RotationMatrix;
              RROTATEQ3_Fr_Ex->rotateX(-45.0 * deg);
              RROTATEQ3_Fr_Ex->rotateY( mRHRSAngle);

              BField_Fringe_Q  * fMagFieldFringeQ3;
              fMagFieldFringeQ3 = new BField_Fringe_Q(3, KAPPA3, pQ3Radius, RORIGINQ3_Fr_Ex, RROTATEQ3_Fr_Ex, 3);

              G4double yyyy[7]={Point[0],Point[1],Point[2],0,0,0,0};
              G4double yBBB[3]={0,0,0};
//              cout<<"   finge  quad ex "<<endl;
              fMagFieldFringeQ3->GetFieldValue(yyyy, yBBB);
//              cout<<"   end fringe ex  "<<endl<<endl<<endl;
//              cout<<"   thet_new"<<thet_new<<" rad or "<<thet_new/3.1415926*180.<<"     B=("<<yBBB[0]/tesla<<","<<yBBB[1]/tesla<<","<<yBBB[2]/tesla<<")"<<endl;
              if (yBBB[0]>-1000.)
              if (yBBB[1]>-1000.)
              if (yBBB[2]>-1000.)
              if (yBBB[0]<1000.)
              if (yBBB[1]<1000.)
              if (yBBB[2]<1000.)
              {
     	        Bfield[0]+=-1.*yBBB[0];
  	        Bfield[1]+=-1.*yBBB[1];
  	        Bfield[2]+=-1.*yBBB[2];
  	      }
    	      //delete RROTATEQ3_Fr_Ex;
  	      delete fMagFieldFringeQ3;
	     }
	    }
	  }

}
