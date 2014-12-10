/*! \class   TrackFitRetinaProducer
 *
 *  \author M Casarsa / L Martini (cut&paste from S.Viret and G.Baulieu's TrackFitHoughProducer) 
 *  \date   2014, May 20
 *
 */

#ifndef TRACK_FITTER_AM_H
#define TRACK_FITTER_AM_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

#include "L1Trigger/TrackFindingAM/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackFindingAM/interface/PatternFinder.h"
#include "L1Trigger/TrackFindingAM/interface/SectorTree.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>


//#ifndef __APPLE__
//BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer)
//#endif

class TrackFitRetinaProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TrackFitRetinaProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TrackFitRetinaProducer();

  private:
  
  /// Data members
  double                       mMagneticField;
  unsigned int                 nSectors;
  unsigned int                 nWedges;
  std::string                  nBKName;
  int                          nThresh;
  const StackedTrackerGeometry *theStackedTracker;
  edm::InputTag                TTStubsInputTag;
  edm::InputTag                TTPatternsInputTag;
  std::string                  TTTrackOutputTag;

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
TrackFitRetinaProducer::TrackFitRetinaProducer( const edm::ParameterSet& iConfig )
{
  TTStubsInputTag    = iConfig.getParameter< edm::InputTag >( "TTInputStubs" );
  TTPatternsInputTag = iConfig.getParameter< edm::InputTag >( "TTInputPatterns" );
  TTTrackOutputTag   = iConfig.getParameter< std::string >( "TTTrackName" );

  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTTrackOutputTag );
}

/// Destructor
TrackFitRetinaProducer::~TrackFitRetinaProducer() {}

/// Begin run
void TrackFitRetinaProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  edm::ESHandle< StackedTrackerGeometry > StackedTrackerGeomHandle;
  iSetup.get< StackedTrackerGeometryRecord >().get( StackedTrackerGeomHandle );
  theStackedTracker = StackedTrackerGeomHandle.product();

  /// Get magnetic field
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
  mMagneticField = (floor(mMagneticFieldStrength*10.0 + 0.5))/10.0;
}

/// End run
void TrackFitRetinaProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFitRetinaProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  // Get GEN particle collection
  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTPatternHandle;

  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );
  iEvent.getByLabel( TTPatternsInputTag, TTPatternHandle );

  /// STEP 0
  /// Prepare output
  TTTracksForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;

  /// STEP 1
  /// Loop over patterns

  //  std::cout << "Start the loop over pattern in order to recover the stubs" << std::endl;

  map<int,vector<Hit*>* > m_hitsPerSector;
  map<int,set<long>*> m_uniqueHitsPerSector;

  edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
  edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator stubIter;

  std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterTTTrack;

  /// Go on only if there are Patterns from PixelDigis
  if ( TTPatternHandle->size() > 0 )
  {

    cout << "-------------------------------------------------------------------" << endl;
    cout << " Generated particles:" << endl; 
    for (std::vector <reco::GenParticle>::const_iterator thepart  = genPart->begin(); 
	 thepart != genPart->end(); 
	 ++thepart) {
	  
      // curvature and helix radius:
      double c = thepart->charge()*0.003*mMagneticField/thepart->pt();
      double R = thepart->pt()/(0.003*mMagneticField);
	  
      // helix center:
      double x0 = thepart->vx() - thepart->charge()*R*thepart->py()/thepart->pt();
      double y0 = thepart->vy() + thepart->charge()*R*thepart->px()/thepart->pt();
	  
      // transverse and longitudinal impact parameters:
      double d0 = thepart->charge()*(sqrt(x0*x0+y0*y0)-R);
      double z0 = thepart->vz() - 2./c*thepart->pz()/thepart->pt()*
	asin(0.5*c*sqrt((thepart->vx()*thepart->vx()+thepart->vy()*thepart->vy()-d0*d0)/(1.+c*d0)));
    
      cout << "   " << iEvent.id().event() << "  -  "
	   << thepart->pdgId()
	   << ": c = " << c
	   << "  pt = " << thepart->pt()
	   << "  d0 = " << d0
	   << "  phi = " << thepart->phi()
	   << "  eta = " << thepart->eta()
	   << "  z0 = " << z0 
	   << endl;

    }

    //cout << "\n Number of roads = " << TTPatternHandle->size() << "\n" << endl;

    /// Loop over Patterns
    unsigned int tkCnt = 0;
    unsigned int j     = 0;
    unsigned int jreal = 0;

    std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMap;
    std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMapUsed;

    for ( iterTTTrack = TTPatternHandle->begin();
	  iterTTTrack != TTPatternHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > tempTrackPtr( TTPatternHandle, tkCnt++ );

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();
      //std::cout << "Pattern in sector " << seedSector << " with " 
      //		<< seedWedge << " active layers contains " 
      //		<< nStubs << " stubs" << std::endl;

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_  > >, TTStub< Ref_PixelDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      //get the hits list of this sector
      map<int,vector<Hit*>*>::iterator sec_it = m_hitsPerSector.find(seedSector);
      vector<Hit*> *m_hits;
      if(sec_it==m_hitsPerSector.end())
      {
	m_hits = new vector<Hit*>();
	m_hitsPerSector[seedSector]=m_hits;
      }
      else
      {
	m_hits = sec_it->second;
      }

      //get the hits set of this sector
      map<int,set<long>*>::iterator set_it = m_uniqueHitsPerSector.find(seedSector);
      set<long> *m_hitIDs;

      if(set_it==m_uniqueHitsPerSector.end())
      {
	m_hitIDs = new set<long>();
	m_uniqueHitsPerSector[seedSector]=m_hitIDs;
      }
      else
      {
	m_hitIDs = set_it->second;
      }

      // Loop over stubs contained in the pattern to recover the info

      vector<Hit*> road_hits;

      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	++j;

	edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = trackStubs.at(i);

	stubMap.insert( std::make_pair( j, tempStubRef ) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint mp0 = tempStubRef->getClusterRef(0)->findAverageLocalCoordinates();
	GlobalPoint posStub  = theStackedTracker->findGlobalPosition( &(*tempStubRef) );
	
	StackedTrackerDetId detIdStub( tempStubRef->getDetId() );
	//bool isPS = theStackedTracker->isPSModule( detIdStub );
	
	const GeomDetUnit* det0 = theStackedTracker->idToDetUnit( detIdStub, 0 );
	const GeomDetUnit* det1 = theStackedTracker->idToDetUnit( detIdStub, 1 );
	
	/// Find pixel pitch and topology related information
	const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelGeomDetUnit* pix1 = dynamic_cast< const PixelGeomDetUnit* >( det1 );
	const PixelTopology* top0    = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );
	const PixelTopology* top1    = dynamic_cast< const PixelTopology* >( &(pix1->specificTopology()) );
	
	/// Find the z-segment
	int cols0   = top0->ncolumns();
	int cols1   = top1->ncolumns();
	int ratio   = cols0/cols1; /// This assumes the ratio is integer!
	int segment = floor( mp0.y() / ratio );
	
	// Here we rearrange the number in order to be compatible with the AM emulator
	if ( detIdStub.isBarrel() )
	{
	  layer  = detIdStub.iLayer()+4;
	  ladder = detIdStub.iPhi()-1;
	  module = detIdStub.iZ()-1;
	}
	else if ( detIdStub.isEndcap() )
	{
	  layer  = 10+detIdStub.iZ()+abs((int)(detIdStub.iSide())-2)*7;
	  ladder = detIdStub.iRing()-1;
	  module = detIdStub.iPhi()-1;
	}

	//cout << layer << " / " << ladder << " / " << module << " / " << std::endl;

	int strip  =  mp0.x();
	int tp     = -1;
	float eta  = 0;
	float phi0 = 0;
	float spt  = 0;
	float x    = posStub.x();
	float y    = posStub.y();
	float z    = posStub.z();
	float x0   = 0.;
	float y0   = 0.;
	float z0   = 0.;
	float ip   = sqrt(x0*x0+y0*y0);
	
	//Check if the stub is already in the list
	long hit_id = (long)layer*10000000000+(long)ladder*100000000
	  +module*100000+segment*10000+strip;

	pair<set<long>::iterator,bool> result = m_hitIDs->insert(hit_id);

	if(result.second==true) //this is a new stub -> add it to the list
	{
	  ++jreal;
	  stubMapUsed.insert( std::make_pair( jreal, tempStubRef ) );

	  if (jreal>=16384)
	  {
	    cout << "Problem!!!" << endl;
	    continue;
	  }

	  Hit* h = new Hit(layer,ladder, module, segment, strip, 
			   jreal, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0);
	  m_hits->push_back(h);

	}

	Hit* h1 = new Hit(layer,ladder, module, segment, strip, 
			  j, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0);
	road_hits.push_back(h1);

      } /// End of loop over track stubs

      //cout << "\n road / number of stubs = " << j << " / " << road_hits.size() << "\n" << endl;

//      cout << " >>>>>>>>>>>>>>>>>> " <<  road_hits.size() << endl;
//      for (std::vector<Hit*>::iterator ihit=road_hits.begin(); ihit!=road_hits.end(); ++ihit)
//      	cout << "        " << (*ihit)->getX() << "  " << (*ihit)->getLayer() << endl;
     
      // Fit the road stubs:

      std::vector<Track*> tracks; 
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > tempVec;
      RetinaTrackFitter* fitter = new RetinaTrackFitter();

      fitter->setSectorID(seedSector);
      fitter->setEventID(iEvent.id());
     
      fitter->fit(road_hits);
      //fitter->mergeTracks();
      tracks.clear();
      tracks = fitter->getTracks();
      fitter->clean();

      // Store the tracks (no duplicate cleaning yet)
      if ( tracks.size() > 0 ) {
	cout << " Fitted tracks:" << endl;
	for(unsigned int tt=0;tt<tracks.size();tt++)
	  {	

	    ///// There is mismatch in the naming: pt <--> curvature to be fixed
	    //double pt_fit = 0.003*mMagneticField/ tracks[tt]->getCurve();
	    //tracks[tt]->setCurve(pt_fit);

	    cout << "   " << iEvent.id().event() << "  -  "
	    	 << tt << ": c = " << tracks[tt]->getCurve()
	    	 << "  pt = " << 0.003*mMagneticField/ tracks[tt]->getCurve()
	    	 << "  d0 = " << tracks[tt]->getD0()
	    	 << "  phi = " << tracks[tt]->getPhi0()
	    	 << "  eta = " << tracks[tt]->getEta0()
	    	 << "  z0 = " << tracks[tt]->getZ0() << endl;

	    tempVec.clear();

	    /////// Stubs used for the fit //////////
	    vector<int> stubs = tracks[tt]->getStubs();
	    for(unsigned int sti=0;sti<stubs.size();sti++)
	      {
		//cout<<stubs[sti]<<endl;
		tempVec.push_back( stubMapUsed[ stubs[sti] ] );
	      }
	    /////////////////////////////////////////


	    double pz = tracks[tt]->getCurve()/(tan(2.*atan(exp(-tracks[tt]->getEta0()))));

	    TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
	    GlobalPoint POCA(0.,0.,tracks[tt]->getZ0());
	    GlobalVector mom(tracks[tt]->getCurve()*cos(tracks[tt]->getPhi0()),
			     tracks[tt]->getCurve()*sin(tracks[tt]->getPhi0()),
			     pz);

	    //	std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << std::endl;

	    tempTrack.setSector( seedSector );
	    tempTrack.setWedge( -1 );
	    tempTrack.setMomentum( mom , 5);
	    tempTrack.setPOCA( POCA , 5);
	    //	std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << " / " << tempTrack.getPOCA().z() << std::endl;
	    TTTracksForOutput->push_back( tempTrack );
	    
	    delete tracks[tt];
	  }
    
      }

      // --- Clean-up memory:
      delete(fitter);

      // Clean-up the road stub vector:
      for(std::vector<Hit*>::iterator ihit=road_hits.begin(); ihit!=road_hits.end(); ++ihit)
	delete *ihit;

    } // End of loop over patterns

    
//    if ( TTTracksForOutput->size()>0 ){
//
//      cout << " Fitted tracks:" << endl;
//
//      for ( std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator itrk  = TTTracksForOutput->begin();
//	                                                             itrk != TTTracksForOutput->end(); 
//                                                                   ++itrk ){
//
//	//cout << itrk->getSector() << endl;
//	//GlobalVector mom = itrk->getMomentum(5);
//	
//	//cout << mom.
//
//	//double pt_fit = 0.003*mMagneticField/ tracks[tt]->getCurve();
//
////	    cout << "   " << iEvent.id().event() << "  -  "
////		 << tt << ": c = " << tracks[tt]->getCurve()
////		 << "  pt = " << 0.003*mMagneticField/ tracks[tt]->getCurve()
////		 << "  d0 = " << tracks[tt]->getD0()
////		 << "  phi = " << tracks[tt]->getPhi0()
////		 << "  eta = " << tracks[tt]->getEta0()
////		 << "  z0 = " << tracks[tt]->getZ0() << endl;
//
////
////    // --- Loop over the fitted tracks
////    unsigned int ntrk = 0;
////    edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksForOutputHandle;
////    for ( std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator itrk = TTTracksForOutputHandle->begin();
////	  itrk != TTTracksForOutputHandle->end(); ++itrk ){
////
////      edm::Ptr< TTTrack< Ref_PixelDigi_ > > trk( TTTracksForOutputHandle, ntrk++ );
////
////      cout << ntrk << " " << trk->getSector() << endl;
////
//
//      }
//    
//    }

    //free the map of sets
    for(map<int,set<long>*>::iterator set_it=m_uniqueHitsPerSector.begin();set_it!=m_uniqueHitsPerSector.end();set_it++)
      delete set_it->second;//delete the set*


    // cout<<"j value "<< j << " / " << jreal <<endl;
    
    /// STEP 2
    /// Passing the hits in the trackfit
    
//    RetinaTrackFitter* fitter = new RetinaTrackFitter();
//    vector<Track*> tracks; 
//    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > tempVec;
//
//    // Loop over the different sectors
//    for(map<int,vector<Hit*>*>::iterator sec_it=m_hitsPerSector.begin();sec_it!=m_hitsPerSector.end();sec_it++)
//    {
//      //cout<<"Sector "<<sec_it->first<<endl;
//      fitter->setSectorID(sec_it->first);
//      fitter->setEventID(iEvent.id());
//
//      // Do the fit
//      fitter->fit(*(sec_it->second));
////      fitter->mergeTracks();//remove some duplicates
//      tracks.clear();
//      tracks = fitter->getTracks();
//      fitter->clean();
//
//
//      if ( tracks.size() > 0 ) {
//	cout << endl;
//	cout << "Generated particles:" << endl; 
//	for (std::vector <reco::GenParticle>::const_iterator thepart  = genPart->begin(); 
//	     thepart != genPart->end(); 
//	     ++thepart) {
//
//	  // curvature and helix radius:
//	  double c = thepart->charge()*0.003*mMagneticField/thepart->pt();
//	  double R = thepart->pt()/(0.003*mMagneticField);
//	  
//	  // helix center:
//	  double x0 = thepart->vx() - thepart->charge()*R*thepart->py()/thepart->pt();
//	  double y0 = thepart->vy() + thepart->charge()*R*thepart->px()/thepart->pt();
//
//	  // transverse and longitudinal impact parameters:
//	  double d0 = thepart->charge()*(sqrt(x0*x0+y0*y0)-R);
//	  double z0 = thepart->vz() - 2./c*thepart->pz()/thepart->pt()*
//	    asin(0.5*c*sqrt((thepart->vx()*thepart->vx()+thepart->vy()*thepart->vy()-d0*d0)/(1+c*d0)));
//	  
//	  cout << iEvent.id().event() << "  -  "
//	       << thepart->pdgId()
//	       << ": c = " << c
//	       << "  pt = " << thepart->pt()
//	       << "  d0 = " << d0
//	       << "  phi = " << thepart->phi()
//	       << "  eta = " << thepart->eta()
//	       << "  z0 = " << z0 
//	       << endl;
//
//	}
//      }
// 
//
//
//      // Store the tracks (no duplicate cleaning yet)
//      cout << "Fitted tracks:" << endl;
//      for(unsigned int tt=0;tt<tracks.size();tt++)
//      {	
//
//	double pt_fit = 0.003*mMagneticField/ tracks[tt]->getCurve();
//
//	cout << iEvent.id().event() << "  -  "
//	     << tt << ": c = " << tracks[tt]->getCurve()
//	     << "  pt = " << 0.003*mMagneticField/ tracks[tt]->getCurve()
//	     << "  d0 = " << tracks[tt]->getD0()
//	     << "  phi = " << tracks[tt]->getPhi0()
//	     << "  eta = " << tracks[tt]->getEta0()
//	     << "  z0 = " << tracks[tt]->getZ0() << endl;
//
//
//	tracks[tt]->setCurve(pt_fit);
//
//	tempVec.clear();
//
//	/////// Stubs used for the fit //////////
//	vector<int> stubs = tracks[tt]->getStubs();
//	for(unsigned int sti=0;sti<stubs.size();sti++)
//	{
//	  //cout<<stubs[sti]<<endl;
//	  tempVec.push_back( stubMapUsed[ stubs[sti] ] );
//	}
//	/////////////////////////////////////////
//
//
//	double pz = tracks[tt]->getCurve()/(tan(2*atan(exp(-tracks[tt]->getEta0()))));
//
//	TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
//	GlobalPoint POCA(0.,0.,tracks[tt]->getZ0());
//	GlobalVector mom(tracks[tt]->getCurve()*cos(tracks[tt]->getPhi0()),
//			 tracks[tt]->getCurve()*sin(tracks[tt]->getPhi0()),
//			 pz);
//
//	//	std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << std::endl;
//
//	tempTrack.setSector( sec_it->first );
//	tempTrack.setWedge( -1 );
//	tempTrack.setMomentum( mom , 5);
//	tempTrack.setPOCA( POCA , 5);
//	//	std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << " / " << tempTrack.getPOCA().z() << std::endl;
//	TTTracksForOutput->push_back( tempTrack );
//
//	delete tracks[tt];
//      }
//
//      for(unsigned int i=0;i<sec_it->second->size();i++)
//      {
//	delete sec_it->second->at(i);//delete the Hit object
//      }
//
//      delete sec_it->second;//delete the vector*
//    }
//    
//    delete(fitter);
//    
  }



  /// Put in the event content
  iEvent.put( TTTracksForOutput, TTTrackOutputTag);
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFitRetinaProducer);

#endif

