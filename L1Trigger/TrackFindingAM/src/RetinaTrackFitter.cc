#include "../interface/RetinaTrackFitter.h"

RetinaTrackFitter::RetinaTrackFitter():TrackFitter(0){
  initialize();
}

RetinaTrackFitter::RetinaTrackFitter(int nb):TrackFitter(nb){
  initialize();
}

RetinaTrackFitter::~RetinaTrackFitter(){
}


void RetinaTrackFitter::fit(vector<Hit*> hits_){
  if(hits_.size()>1024){
    cout<<"ERROR : too many stubs for fitting!"<<endl;
    return;
  }


  // --- Determine in which eta range the trigger tower is:
  const int eta_range = (int) sector_id / 8;


  // --- Constants used in X+-X- transformation:
  double y0 =  0.0217391304347826081; // 0.5/23.
  double y1 =  0.0046296296296296294; // 0.5/108.


  // --- Fill the stubs vector:
  vector <Hit_t*> hits;
  vector <Hit_t*> hits_RZ;

  for(unsigned int ihit=0; ihit<hits_.size(); ihit++){
  
    Hit_t* hit = new Hit_t();
    hit->x     = hits_[ihit]->getX();
    hit->y     = hits_[ihit]->getY();
    hit->z     = hits_[ihit]->getZ();
    hit->rho   = sqrt( hits_[ihit]->getX()*hits_[ihit]->getX() +
		       hits_[ihit]->getY()*hits_[ihit]->getY() );
    hit->layer = hits_[ihit]->getLayer();
    hit->id    = hits_[ihit]->getID();
    
    hits.push_back(hit);

    //cout << " x = " << hits_[ihit]->getX() << "   "
    //	 << " y = " << hits_[ihit]->getY() << "   "
    //	 << " z = " << hits_[ihit]->getZ() << " ---> " 
    //	 << " R = " << sqrt(hits_[ihit]->getX()*hits_[ihit]->getX()+
    //			    hits_[ihit]->getY()*hits_[ihit]->getY())
    //	 << "  -  Layer = " <<  hits_[ihit]->getLayer()
    //	 << endl;
   
  }



  // ===========================================================================
  //  XY fit
  // ===========================================================================

  // --- Phi sector rotation:
  rotateHits(hits, rot_angle[sector_id%8]);


  // --- Conformal transformation: 
  confTrans(hits);


  //
  // --- First step ------------------------------------------------------------
  //

  // --- Setup the retina:
  double pbins_step1 = config[eta_range]["xy_pbins_step1"];
  double qbins_step1 = config[eta_range]["xy_qbins_step1"];
  double pmin_step1  = config[eta_range]["xy_pmin_step1"];
  double pmax_step1  = config[eta_range]["xy_pmax_step1"];
  double qmin_step1  = config[eta_range]["xy_qmin_step1"];
  double qmax_step1  = config[eta_range]["xy_qmax_step1"];
  
  double minWeight_step1 = config[eta_range]["xy_threshold_step1"];

  double pstep_step1 = (pmax_step1-pmin_step1)/pbins_step1;
  double qstep_step1 = (qmax_step1-qmin_step1)/qbins_step1;

  vector <double> sigma_step1(8,0.25*sqrt(pstep_step1*pstep_step1+qstep_step1*qstep_step1));
  if ( config[eta_range]["xy_sigma1_step1"] != 0. ) 
    sigma_step1[0] = config[eta_range]["xy_sigma1_step1"];
  if ( config[eta_range]["xy_sigma2_step1"] != 0. ) 
    sigma_step1[1] = config[eta_range]["xy_sigma2_step1"];
  if ( config[eta_range]["xy_sigma3_step1"] != 0. ) 
    sigma_step1[2] = config[eta_range]["xy_sigma3_step1"];
  if ( config[eta_range]["xy_sigma4_step1"] != 0. ) 
    sigma_step1[3] = config[eta_range]["xy_sigma4_step1"];
  if ( config[eta_range]["xy_sigma5_step1"] != 0. ) 
    sigma_step1[4] = config[eta_range]["xy_sigma5_step1"];
  if ( config[eta_range]["xy_sigma6_step1"] != 0. ) 
    sigma_step1[5] = config[eta_range]["xy_sigma6_step1"];
  if ( config[eta_range]["xy_sigma7_step1"] != 0. ) 
    sigma_step1[6] = config[eta_range]["xy_sigma7_step1"];
  if ( config[eta_range]["xy_sigma8_step1"] != 0. ) 
    sigma_step1[7] = config[eta_range]["xy_sigma8_step1"];


  Retina retinaXY_step1(hits, pbins_step1+2, qbins_step1+2, 
			pmin_step1-pstep_step1, pmax_step1+pstep_step1, 
			qmin_step1-qstep_step1, qmax_step1+qstep_step1, 
			sigma_step1, minWeight_step1, 1, XY);

  // --- Fill the retina and find maxima:
  retinaXY_step1.fillGrid();
  retinaXY_step1.findMaxima();
  //retinaXY_step1.dumpGrid(eventID.event(),1);
  //retinaXY_step1.printMaxima();


  // --- Get first step maxima:
  vector <pqPoint> maximaXY_step1 = retinaXY_step1.getMaxima();


  //
  // --- Second step -----------------------------------------------------------
  //

  double pbins_step2 = config[eta_range]["xy_pbins_step2"];
  double qbins_step2 = config[eta_range]["xy_qbins_step2"];


  // --- Zoom around first step maxima:
  for (unsigned int imax=0; imax<maximaXY_step1.size(); ++imax){
    if (maximaXY_step1[imax].w == -1.) continue;

    // --- Retina setup:
    double pmin_step2 = maximaXY_step1[imax].p - config[eta_range]["xy_zoom_step2"]*pstep_step1;
    double pmax_step2 = maximaXY_step1[imax].p + config[eta_range]["xy_zoom_step2"]*pstep_step1;
    double qmin_step2 = maximaXY_step1[imax].q - config[eta_range]["xy_zoom_step2"]*qstep_step1;
    double qmax_step2 = maximaXY_step1[imax].q + config[eta_range]["xy_zoom_step2"]*qstep_step1;
   
    double pstep_step2 = (pmax_step2-pmin_step2)/pbins_step2;
    double qstep_step2 = (qmax_step2-qmin_step2)/qbins_step2;
    
    double minWeight_step2 = config[eta_range]["xy_threshold_step2"];

    vector <double> sigma_step2(8,0.25*sqrt(pstep_step2*pstep_step2+qstep_step2*qstep_step2));
    if ( config[eta_range]["xy_sigma1_step2"] != 0. ) 
      sigma_step2[0] = config[eta_range]["xy_sigma1_step2"];
    if ( config[eta_range]["xy_sigma2_step2"] != 0. ) 
      sigma_step2[1] = config[eta_range]["xy_sigma2_step2"];
    if ( config[eta_range]["xy_sigma3_step2"] != 0. ) 
      sigma_step2[2] = config[eta_range]["xy_sigma3_step2"];
    if ( config[eta_range]["xy_sigma4_step2"] != 0. ) 
      sigma_step2[3] = config[eta_range]["xy_sigma4_step2"];
    if ( config[eta_range]["xy_sigma5_step2"] != 0. ) 
      sigma_step2[4] = config[eta_range]["xy_sigma5_step2"];
    if ( config[eta_range]["xy_sigma6_step2"] != 0. ) 
      sigma_step2[5] = config[eta_range]["xy_sigma6_step2"];
    if ( config[eta_range]["xy_sigma7_step2"] != 0. ) 
      sigma_step2[6] = config[eta_range]["xy_sigma7_step2"];
    if ( config[eta_range]["xy_sigma8_step2"] != 0. ) 
      sigma_step2[7] = config[eta_range]["xy_sigma8_step2"];


    Retina retinaXY_step2(hits, pbins_step2+2, qbins_step2+2, 
			  pmin_step2-pstep_step2, pmax_step2+pstep_step2, 
			  qmin_step2-qstep_step2, qmax_step2+qstep_step2, 
			  sigma_step2, minWeight_step2, 1, XY);

    // --- Fill the retina and find maxima:
    retinaXY_step2.fillGrid();
    retinaXY_step2.findMaxima();
    //retinaXY_step2.dumpGrid(eventID.event(),2,imax);
    //retinaXY_step2.printMaxima();

    pqPoint bestpqXY_step2 = retinaXY_step2.getBestPQ();
    if ( bestpqXY_step2.w == -1. ) continue;


    // --- Invert the X+-X- transformation:
    double p = 0.5*(y1 - y0)/bestpqXY_step2.q;
    double q = y0 - p*(bestpqXY_step2.p-bestpqXY_step2.q);


    // --- Associate stubs to this maxumum:
    hits_RZ.clear();
    for (unsigned int ihit=0; ihit<hits.size(); ++ihit){
      
      double dist   = fabs(hits[ihit]->y-p*hits[ihit]->x-q)/sqrt(1.+p*p);
      double weight = exp(-0.5*dist*dist/(sigma_step2[0]*sigma_step2[0]));

      if ( weight>0.5 ){
    	hits_RZ.push_back(hits[ihit]);
      }
      //cout << ihit << " - " << dist << "  " << weight << endl;

    }


    // --- Rotate back the original phi sector:
    q = q/(cos(rot_angle[sector_id%8])+p*sin(rot_angle[sector_id%8]));
    p = (p*cos(rot_angle[sector_id%8])-sin(rot_angle[sector_id%8]))/
        (cos(rot_angle[sector_id%8])+p*sin(rot_angle[sector_id%8]));


    // --- Invert the conformal transformation and get the track parameters:
    double a = -0.5*p/q;
    double b =  0.5/q;
    
    double c   = 1./sqrt(a*a+b*b);
    double phi = atan(p);
    //if (phi<0.)
    //  phi += TMath::Pi();


    // =========================================================================
    //  RZ fit
    // =========================================================================

    y0 = 0.5/y0;
    y1 = 0.5/y1;

    double eta = -9999.;
    double z0  = -9999.;
    
    //
    // --- First step ----------------------------------------------------------
    //

    pbins_step1 = config[eta_range]["rz_pbins_step1"];
    qbins_step1 = config[eta_range]["rz_qbins_step1"];
    pmin_step1  = config[eta_range]["rz_pmin_step1"];
    pmax_step1  = config[eta_range]["rz_pmax_step1"];
    qmin_step1  = config[eta_range]["rz_qmin_step1"];
    qmax_step1  = config[eta_range]["rz_qmax_step1"];

    minWeight_step1 = config[eta_range]["rz_threshold_step1"];

    pstep_step1 = (pmax_step1-pmin_step1)/pbins_step1;
    qstep_step1 = (qmax_step1-qmin_step1)/qbins_step1;

    for (unsigned int ilayer=0; ilayer<8; ++ilayer)
      sigma_step1[ilayer] = 0.5*sqrt(pstep_step1*pstep_step1+qstep_step1*qstep_step1);

    if ( config[eta_range]["rz_sigma1_step1"] != 0. ) 
      sigma_step1[0] = config[eta_range]["rz_sigma1_step1"];
    if ( config[eta_range]["rz_sigma2_step1"] != 0. ) 
      sigma_step1[1] = config[eta_range]["rz_sigma2_step1"];
    if ( config[eta_range]["rz_sigma3_step1"] != 0. ) 
      sigma_step1[2] = config[eta_range]["rz_sigma3_step1"];
    if ( config[eta_range]["rz_sigma4_step1"] != 0. ) 
      sigma_step1[3] = config[eta_range]["rz_sigma4_step1"];
    if ( config[eta_range]["rz_sigma5_step1"] != 0. ) 
      sigma_step1[4] = config[eta_range]["rz_sigma5_step1"];
    if ( config[eta_range]["rz_sigma6_step1"] != 0. ) 
      sigma_step1[5] = config[eta_range]["rz_sigma6_step1"];
    if ( config[eta_range]["rz_sigma7_step1"] != 0. ) 
      sigma_step1[6] = config[eta_range]["rz_sigma7_step1"];
    if ( config[eta_range]["rz_sigma8_step1"] != 0. ) 
      sigma_step1[7] = config[eta_range]["rz_sigma8_step1"];


    Retina retinaRZ_step1(hits_RZ, pbins_step1+2, qbins_step1+2, 
			  pmin_step1-pstep_step1, pmax_step1+pstep_step1, 
			  qmin_step1-qstep_step1, qmax_step1+qstep_step1, 
			  sigma_step1, minWeight_step1, 1, RZ);

    retinaRZ_step1.fillGrid();
    retinaRZ_step1.findMaxima();
    //retinaRZ_step1.dumpGrid(eventID.event(),imax);
    //retinaRZ_step1.printMaxima();


    // --- Get first step maximum:
    vector <pqPoint> maximaRZ_step1 = retinaRZ_step1.getMaxima();
    //pqPoint bestpqRZ_step1 = retinaRZ_step1.getBestPQ();

    
    //
    // --- Second step ---------------------------------------------------------
    //

    pbins_step2 = config[eta_range]["rz_pbins_step2"];
    qbins_step2 = config[eta_range]["rz_qbins_step2"];

    // Zoom around first step maxima
    for (unsigned int imax_RZ=0; imax_RZ<maximaRZ_step1.size(); ++imax_RZ){
      if (maximaRZ_step1[imax].w == -1.) continue;
 
      double pmin_step2 = maximaRZ_step1[imax_RZ].p - config[eta_range]["rz_zoom_step2"]*pstep_step1;
      double pmax_step2 = maximaRZ_step1[imax_RZ].p + config[eta_range]["rz_zoom_step2"]*pstep_step1;
      double qmin_step2 = maximaRZ_step1[imax_RZ].q - config[eta_range]["rz_zoom_step2"]*qstep_step1;
      double qmax_step2 = maximaRZ_step1[imax_RZ].q + config[eta_range]["rz_zoom_step2"]*qstep_step1;
   
      double pstep_step2 = (pmax_step2-pmin_step2)/pbins_step2;
      double qstep_step2 = (qmax_step2-qmin_step2)/qbins_step2;
    
      double minWeight_step2 = config[eta_range]["rz_threshold_step2"];


      vector <double> sigma_step2(8,0.5*sqrt(pstep_step2*pstep_step2+qstep_step2*qstep_step2));
      for (unsigned int ilayer=3; ilayer<6; ++ilayer)
	sigma_step2[ilayer] = 6.*sqrt(pstep_step2*pstep_step2+qstep_step2*qstep_step2);

      if ( config[eta_range]["rz_sigma1_step2"] != 0. ) 
	sigma_step2[0] = config[eta_range]["rz_sigma1_step2"];
      if ( config[eta_range]["rz_sigma2_step2"] != 0. ) 
	sigma_step2[1] = config[eta_range]["rz_sigma2_step2"];
      if ( config[eta_range]["rz_sigma3_step2"] != 0. ) 
	sigma_step2[2] = config[eta_range]["rz_sigma3_step2"];
      if ( config[eta_range]["rz_sigma4_step2"] != 0. ) 
	sigma_step2[3] = config[eta_range]["rz_sigma4_step2"];
      if ( config[eta_range]["rz_sigma5_step2"] != 0. ) 
	sigma_step2[4] = config[eta_range]["rz_sigma5_step2"];
      if ( config[eta_range]["rz_sigma6_step2"] != 0. ) 
	sigma_step2[5] = config[eta_range]["rz_sigma6_step2"];
      if ( config[eta_range]["rz_sigma7_step2"] != 0. ) 
	sigma_step2[6] = config[eta_range]["rz_sigma7_step2"];
      if ( config[eta_range]["rz_sigma8_step2"] != 0. ) 
	sigma_step2[7] = config[eta_range]["rz_sigma8_step2"];


      Retina retinaRZ_step2(hits_RZ, pbins_step2+2, qbins_step2+2, 
			    pmin_step2-pstep_step2, pmax_step2+pstep_step2, 
			    qmin_step2-qstep_step2, qmax_step2+qstep_step2, 
			    sigma_step2, minWeight_step2, 1, RZ);

      retinaRZ_step2.fillGrid();
      retinaRZ_step2.findMaxima();
      //retinaRZ_step2.dumpGrid(eventID.event(),2,10+imax_RZ);
      //retinaRZ_step2.printMaxima();

      pqPoint bestpqRZ_step2 = retinaRZ_step2.getBestPQ();
      if (bestpqRZ_step2.w == -1.) continue;


      // --- Invert the X+-X- transformation:
      double p = 0.5*(y1 - y0)/bestpqRZ_step2.q;
      double q = y0 - p*(bestpqRZ_step2.p-bestpqRZ_step2.q);


      // --- Get the track parameters:
      double theta = atan(p);
      eta = -log(tan(0.5*theta));
      z0  = -q/p;

      //cout << c << "  " << phi << "  -  " << eta << "  " << z0 << endl; 

      if ( !std::isnan(c) && !std::isnan(phi) && !std::isnan(eta) && !std::isnan(z0) &&
	   eta != -9999. && z0 != -9999. ){

	Track* trk = new Track(c, 0., phi, eta, z0);
      
	for(unsigned int ihit=0; ihit<hits_RZ.size(); ihit++)
	  trk->addStubIndex(hits_RZ[ihit]->id);
      
	tracks.push_back(trk);

      }
 

    } // imax_RZ loop


  } // imax loop


  // --- Clean-up pointers:
  for(vector<Hit_t*>::iterator it=hits.begin(); it!=hits.end(); ++it)
    delete *it;

  
}

void RetinaTrackFitter::fit(){

  vector<Hit*> activatedHits;

  //////// Get the list of unique stubs from all the patterns ///////////
  set<int> ids;
  int total=0;
  
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    total+=allHits.size();
    for(unsigned int j=0;j<allHits.size();j++){
      pair<set<int>::iterator,bool> result = ids.insert(allHits[j]->getID());
      if(result.second==true)
	activatedHits.push_back(allHits[j]);
    }
  }

  fit(activatedHits);
 
}


void RetinaTrackFitter::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}


void RetinaTrackFitter::mergeTracks(){
//  
//  if ( tracks.size()<2 ) return;   // There is nothing to merge
//
//  // --- Loop over tracks, count the stubs in common and discard the duplicates:
//  std::set<int> trk_index;
//
//  cout << ">>>>>>>>>>>>>>>>>>>>>>>> n track = " << tracks.size() << endl; 
//
//  for( unsigned int itrk=0; itrk<tracks.size(); ++itrk){
//    vector<int> stubs_i = tracks[itrk]->getStubs();
//    for( unsigned int jtrk=itrk; jtrk<tracks.size(); ++jtrk){
//      vector<int> stubs_j = tracks[jtrk]->getStubs();
//      
//      unsigned int nShared=0;
//      for ( unsigned int i=0; i<stubs_i.size() && nShared<2; i++){
//	for ( unsigned int j=0; j<stubs_j.size() && nShared<2; j++){
//	  
//	  cout << "               " << stubs_i[i] << " " << stubs_j[j] << endl;
//
//	  if ( stubs_i[i] == stubs_j[j] )
//	    nShared++;
//
//	}
//      }
//      
//      // --- If more than 1 stub in common flag the track as a duplicate
//      if ( nShared > 1 )
//	trk_index.insert(jtrk);
//      //stubs_i.size() > stubs_j.size() ? trk_index.insert(jtrk) : trk_index.insert(itrk);
//
//    } // loop over jtrk
//  } // loop over itrk
//
//  for (std::set<int>::iterator itrk=trk_index.begin(); itrk!=trk_index.end(); ++itrk) {
//    cout << " -- > " << *itrk << endl;
//  }
//
//
}


TrackFitter* RetinaTrackFitter::clone(){
  RetinaTrackFitter* fit = new RetinaTrackFitter(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

void RetinaTrackFitter::rotateHits(vector<Hit_t*> hits, double angle){
  
  for (unsigned int ihit=0; ihit<hits.size(); ihit++) {
    double x = hits[ihit]->x*cos(rot_angle[sector_id%8]) - hits[ihit]->y*sin(rot_angle[sector_id%8]);
    double y = hits[ihit]->x*sin(rot_angle[sector_id%8]) + hits[ihit]->y*cos(rot_angle[sector_id%8]);
    hits[ihit]->x = x;
    hits[ihit]->y = y;
  }

}

void RetinaTrackFitter::confTrans(vector<Hit_t*> hits){
  
  for (unsigned int ihit=0; ihit<hits.size(); ihit++) {
    double R2 = hits[ihit]->x*hits[ihit]->x + hits[ihit]->y*hits[ihit]->y;
    hits[ihit]->x /= R2;
    hits[ihit]->y /= R2;
  }

}

void RetinaTrackFitter::initialize(){

  // Enter all the retina parameters
  // (we refer to the detector geometry in
  //  http://sviret.web.cern.ch/sviret/Images/CMS/Upgrade/Eta6_Phi8.jpg)

  // --- eta range 1
  config[0]["xy_pbins_step1"]     = 40.;
  config[0]["xy_qbins_step1"]     = 40.;
  config[0]["xy_pmin_step1"]      = -0.05;
  config[0]["xy_pmax_step1"]      =  0.05;
  config[0]["xy_qmin_step1"]      = -0.05;
  config[0]["xy_qmax_step1"]      =  0.05;
  config[0]["xy_threshold_step1"] =  4.5;
  config[0]["xy_sigma1_step1"]    =  0.;
  config[0]["xy_sigma2_step1"]    =  0.;
  config[0]["xy_sigma3_step1"]    =  0.;
  config[0]["xy_sigma4_step1"]    =  0.;
  config[0]["xy_sigma2_step1"]    =  0.;
  config[0]["xy_sigma5_step1"]    =  0.;
  config[0]["xy_sigma6_step1"]    =  0.;
  config[0]["xy_sigma7_step1"]    =  0.;
  config[0]["xy_sigma8_step1"]    =  0.;
  config[0]["xy_pbins_step2"]     = 100.;
  config[0]["xy_qbins_step2"]     = 100.;
  config[0]["xy_zoom_step2"]      = 1.;
  config[0]["xy_threshold_step2"] =  4.5;
  config[0]["xy_sigma1_step2"]    =  0.;
  config[0]["xy_sigma2_step2"]    =  0.;
  
  config[0]["rz_pbins_step1"]     = 20.;
  config[0]["rz_qbins_step1"]     = 20.;
  config[0]["rz_pmin_step1"]      = -20.;
  config[0]["rz_pmax_step1"]      =  60.;
  config[0]["rz_qmin_step1"]      = -60.;
  config[0]["rz_qmax_step1"]      =  60.;
  config[0]["rz_threshold_step1"] =  4.5;
  config[0]["rz_sigma1_step1"]    =  0.;
  config[0]["rz_sigma2_step1"]    =  0.;
  config[0]["rz_sigma3_step1"]    =  0.;
  config[0]["rz_sigma4_step1"]    =  0.;
  config[0]["rz_sigma2_step1"]    =  0.;
  config[0]["rz_sigma5_step1"]    =  0.;
  config[0]["rz_sigma6_step1"]    =  0.;
  config[0]["rz_sigma7_step1"]    =  0.;
  config[0]["rz_sigma8_step1"]    =  0.;
  config[0]["rz_pbins_step2"]     = 80.;
  config[0]["rz_qbins_step2"]     = 80.;
  config[0]["rz_zoom_step2"]      = 1.5;
  config[0]["rz_threshold_step2"] =  4.5;
  config[0]["rz_sigma1_step2"]    =  0.;
  config[0]["rz_sigma2_step2"]    =  0.;
  
  // --- eta range 2
  config[1]["xy_pbins_step1"]     = 40.;
  config[1]["xy_qbins_step1"]     = 40.;
  config[1]["xy_pmin_step1"]      = -0.05;
  config[1]["xy_pmax_step1"]      =  0.05;
  config[1]["xy_qmin_step1"]      = -0.05;
  config[1]["xy_qmax_step1"]      =  0.05;
  config[1]["xy_threshold_step1"] =  4.5;
  config[1]["xy_sigma1_step1"]    =  0.;
  config[1]["xy_sigma2_step1"]    =  0.;
  config[1]["xy_sigma3_step1"]    =  0.;
  config[1]["xy_sigma4_step1"]    =  0.;
  config[1]["xy_sigma2_step1"]    =  0.;
  config[1]["xy_sigma5_step1"]    =  0.;
  config[1]["xy_sigma6_step1"]    =  0.;
  config[1]["xy_sigma7_step1"]    =  0.;
  config[1]["xy_sigma8_step1"]    =  0.;
  config[1]["xy_pbins_step2"]     = 100.;
  config[1]["xy_qbins_step2"]     = 100.;
  config[1]["xy_zoom_step2"]      = 1.;
  config[1]["xy_threshold_step2"] =  4.5;
  config[1]["xy_sigma1_step2"]    =  0.;
  config[1]["xy_sigma2_step2"]    =  0.;
  
  config[1]["rz_pbins_step1"]     = 20.;
  config[1]["rz_qbins_step1"]     = 20.;
  config[1]["rz_pmin_step1"]      = -20.;
  config[1]["rz_pmax_step1"]      =  60.;
  config[1]["rz_qmin_step1"]      = -60.;
  config[1]["rz_qmax_step1"]      =  60.;
  config[1]["rz_threshold_step1"] =  4.5;
  config[1]["rz_sigma1_step1"]    =  0.;
  config[1]["rz_sigma2_step1"]    =  0.;
  config[1]["rz_sigma3_step1"]    =  0.;
  config[1]["rz_sigma4_step1"]    =  0.;
  config[1]["rz_sigma2_step1"]    =  0.;
  config[1]["rz_sigma5_step1"]    =  0.;
  config[1]["rz_sigma6_step1"]    =  0.;
  config[1]["rz_sigma7_step1"]    =  0.;
  config[1]["rz_sigma8_step1"]    =  0.;
  config[1]["rz_pbins_step2"]     = 80.;
  config[1]["rz_qbins_step2"]     = 80.;
  config[1]["rz_zoom_step2"]      = 1.5;
  config[1]["rz_threshold_step2"] =  4.5;
  config[1]["rz_sigma1_step2"]    =  0.;
  config[1]["rz_sigma2_step2"]    =  0.;
  
  // --- eta range 3
  config[2]["xy_pbins_step1"]     = 40.;
  config[2]["xy_qbins_step1"]     = 40.;
  config[2]["xy_pmin_step1"]      = -0.05;
  config[2]["xy_pmax_step1"]      =  0.05;
  config[2]["xy_qmin_step1"]      = -0.05;
  config[2]["xy_qmax_step1"]      =  0.05;
  config[2]["xy_threshold_step1"] =  4.5;
  config[2]["xy_sigma1_step1"]    =  0.;
  config[2]["xy_sigma2_step1"]    =  0.;
  config[2]["xy_sigma3_step1"]    =  0.;
  config[2]["xy_sigma4_step1"]    =  0.;
  config[2]["xy_sigma2_step1"]    =  0.;
  config[2]["xy_sigma5_step1"]    =  0.;
  config[2]["xy_sigma6_step1"]    =  0.;
  config[2]["xy_sigma7_step1"]    =  0.;
  config[2]["xy_sigma8_step1"]    =  0.;
  config[2]["xy_pbins_step2"]     = 100.;
  config[2]["xy_qbins_step2"]     = 100.;
  config[2]["xy_zoom_step2"]      = 1.;
  config[2]["xy_threshold_step2"] =  4.5;
  config[2]["xy_sigma1_step2"]    =  0.;
  config[2]["xy_sigma2_step2"]    =  0.;

  config[2]["rz_pbins_step1"]     = 20.;
  config[2]["rz_qbins_step1"]     = 20.;
  config[2]["rz_pmin_step1"]      = -20.;
  config[2]["rz_pmax_step1"]      =  60.;
  config[2]["rz_qmin_step1"]      = -60.;
  config[2]["rz_qmax_step1"]      =  60.;
  config[2]["rz_threshold_step1"] =  4.5;
  config[2]["rz_sigma1_step1"]    =  0.;
  config[2]["rz_sigma2_step1"]    =  0.;
  config[2]["rz_sigma3_step1"]    =  0.;
  config[2]["rz_sigma4_step1"]    =  0.;
  config[2]["rz_sigma2_step1"]    =  0.;
  config[2]["rz_sigma5_step1"]    =  0.;
  config[2]["rz_sigma6_step1"]    =  0.;
  config[2]["rz_sigma7_step1"]    =  0.;
  config[2]["rz_sigma8_step1"]    =  0.;
  config[2]["rz_pbins_step2"]     = 80.;
  config[2]["rz_qbins_step2"]     = 80.;
  config[2]["rz_zoom_step2"]      = 1.5;
  config[2]["rz_threshold_step2"] =  4.5;
  config[2]["rz_sigma1_step2"]    =  0.;
  config[2]["rz_sigma2_step2"]    =  0.;
  
  // --- eta range 4
  config[3]["xy_pbins_step1"]     = 40.;
  config[3]["xy_qbins_step1"]     = 40.;
  config[3]["xy_pmin_step1"]      = -0.05;
  config[3]["xy_pmax_step1"]      =  0.05;
  config[3]["xy_qmin_step1"]      = -0.05;
  config[3]["xy_qmax_step1"]      =  0.05;
  config[3]["xy_threshold_step1"] =  4.;
  config[3]["xy_sigma1_step1"]    =  0.;
  config[3]["xy_sigma2_step1"]    =  0.;
  config[3]["xy_sigma3_step1"]    =  0.;
  config[3]["xy_sigma4_step1"]    =  0.;
  config[3]["xy_sigma2_step1"]    =  0.;
  config[3]["xy_sigma5_step1"]    =  0.;
  config[3]["xy_sigma6_step1"]    =  0.;
  config[3]["xy_sigma7_step1"]    =  0.;
  config[3]["xy_sigma8_step1"]    =  0.;
  config[3]["xy_pbins_step2"]     = 100.;
  config[3]["xy_qbins_step2"]     = 100.;
  config[3]["xy_zoom_step2"]      = 1.;
  config[3]["xy_threshold_step2"] =  4.5;
  config[3]["xy_sigma1_step2"]    =  0.;
  config[3]["xy_sigma2_step2"]    =  0.;

  config[3]["rz_pbins_step1"]     = 20.;
  config[3]["rz_qbins_step1"]     = 20.;
  config[3]["rz_pmin_step1"]      = -20.;
  config[3]["rz_pmax_step1"]      =  60.;
  config[3]["rz_qmin_step1"]      = -60.;
  config[3]["rz_qmax_step1"]      =  60.;
  config[3]["rz_threshold_step1"] =  4.;
  config[3]["rz_sigma1_step1"]    =  0.;
  config[3]["rz_sigma2_step1"]    =  0.;
  config[3]["rz_sigma3_step1"]    =  0.;
  config[3]["rz_sigma4_step1"]    =  0.;
  config[3]["rz_sigma2_step1"]    =  0.;
  config[3]["rz_sigma5_step1"]    =  0.;
  config[3]["rz_sigma6_step1"]    =  0.;
  config[3]["rz_sigma7_step1"]    =  0.;
  config[3]["rz_sigma8_step1"]    =  0.;
  config[3]["rz_pbins_step2"]     = 80.;
  config[3]["rz_qbins_step2"]     = 80.;
  config[3]["rz_zoom_step2"]      = 1.5;
  config[3]["rz_threshold_step2"] =  4.;
  config[3]["rz_sigma1_step2"]    =  0.;
  config[3]["rz_sigma2_step2"]    =  0.;
  
  // --- eta range 5
  config[4]["xy_pbins_step1"]     = 40.;
  config[4]["xy_qbins_step1"]     = 40.;
  config[4]["xy_pmin_step1"]      = -0.05;
  config[4]["xy_pmax_step1"]      =  0.05;
  config[4]["xy_qmin_step1"]      = -0.05;
  config[4]["xy_qmax_step1"]      =  0.05;
  config[4]["xy_threshold_step1"] =  4.5;
  config[4]["xy_sigma1_step1"]    =  0.;
  config[4]["xy_sigma2_step1"]    =  0.;
  config[4]["xy_sigma3_step1"]    =  0.;
  config[4]["xy_sigma4_step1"]    =  0.;
  config[4]["xy_sigma2_step1"]    =  0.;
  config[4]["xy_sigma5_step1"]    =  0.;
  config[4]["xy_sigma6_step1"]    =  0.;
  config[4]["xy_sigma7_step1"]    =  0.;
  config[4]["xy_sigma8_step1"]    =  0.;
  config[4]["xy_pbins_step2"]     = 100.;
  config[4]["xy_qbins_step2"]     = 100.;
  config[4]["xy_zoom_step2"]      = 1.;
  config[4]["xy_threshold_step2"] =  4.5;
  config[4]["xy_sigma1_step2"]    =  0.;
  config[4]["xy_sigma2_step2"]    =  0.;

  config[4]["rz_pbins_step1"]     = 20.;
  config[4]["rz_qbins_step1"]     = 20.;
  config[4]["rz_pmin_step1"]      = 40.;
  config[4]["rz_pmax_step1"]      = 140.;
  config[4]["rz_qmin_step1"]      = 0.;
  config[4]["rz_qmax_step1"]      = 120.;
  config[4]["rz_threshold_step1"] =  4.5;
  config[4]["rz_sigma1_step1"]    =  0.;
  config[4]["rz_sigma2_step1"]    =  0.;
  config[4]["rz_sigma3_step1"]    =  0.;
  config[4]["rz_sigma4_step1"]    =  0.;
  config[4]["rz_sigma2_step1"]    =  0.;
  config[4]["rz_sigma5_step1"]    =  0.;
  config[4]["rz_sigma6_step1"]    =  0.;
  config[4]["rz_sigma7_step1"]    =  0.;
  config[4]["rz_sigma8_step1"]    =  0.;
  config[4]["rz_pbins_step2"]     = 80.;
  config[4]["rz_qbins_step2"]     = 80.;
  config[4]["rz_zoom_step2"]      = 1.5;
  config[4]["rz_threshold_step2"] =  4.5;
  config[4]["rz_sigma1_step2"]    =  0.;
  config[4]["rz_sigma2_step2"]    =  0.;
  
  // --- eta range 6
  config[5]["xy_pbins_step1"]     = 40.;
  config[5]["xy_qbins_step1"]     = 40.;
  config[5]["xy_pmin_step1"]      = -10.05;
  config[5]["xy_pmax_step1"]      =  10.05;
  config[5]["xy_qmin_step1"]      = -10.05;
  config[5]["xy_qmax_step1"]      =  10.05;
  config[5]["xy_threshold_step1"] =  4.5;
  config[5]["xy_sigma1_step1"]    =  0.;
  config[5]["xy_sigma2_step1"]    =  0.;
  config[5]["xy_sigma3_step1"]    =  0.;
  config[5]["xy_sigma4_step1"]    =  0.;
  config[5]["xy_sigma2_step1"]    =  0.;
  config[5]["xy_sigma5_step1"]    =  0.;
  config[5]["xy_sigma6_step1"]    =  0.;
  config[5]["xy_sigma7_step1"]    =  0.;
  config[5]["xy_sigma8_step1"]    =  0.;
  config[5]["xy_pbins_step2"]     = 100.;
  config[5]["xy_qbins_step2"]     = 100.;
  config[5]["xy_zoom_step2"]      = 1.;
  config[5]["xy_threshold_step2"] =  4.5;
  config[5]["xy_sigma1_step2"]    =  0.;
  config[5]["xy_sigma2_step2"]    =  0.;

  config[5]["rz_pbins_step1"]     = 20.;
  config[5]["rz_qbins_step1"]     = 20.;
  config[5]["rz_pmin_step1"]      = 100.;
  config[5]["rz_pmax_step1"]      = 300.;
  config[5]["rz_qmin_step1"]      = 100.;
  config[5]["rz_qmax_step1"]      = 300.;
  config[5]["rz_threshold_step1"] =  4.5;
  config[5]["rz_sigma1_step1"]    =  0.;
  config[5]["rz_sigma2_step1"]    =  0.;
  config[5]["rz_sigma3_step1"]    =  0.;
  config[5]["rz_sigma4_step1"]    =  0.;
  config[5]["rz_sigma2_step1"]    =  0.;
  config[5]["rz_sigma5_step1"]    =  0.;
  config[5]["rz_sigma6_step1"]    =  0.;
  config[5]["rz_sigma7_step1"]    =  0.;
  config[5]["rz_sigma8_step1"]    =  0.;
  config[5]["rz_pbins_step2"]     = 80.;
  config[5]["rz_qbins_step2"]     = 80.;
  config[5]["rz_zoom_step2"]      = 1.5;
  config[5]["rz_threshold_step2"] = 4.5;
  config[5]["rz_sigma1_step2"]    = 0.;
  config[5]["rz_sigma2_step2"]    = 0.;
  

}
