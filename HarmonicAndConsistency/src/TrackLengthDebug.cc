#include "TrackLengthDebug.h"
#include "TrackLengthDebugUtils.h"

#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "MarlinTrk/Factory.h"
#include "EVENT/SimTrackerHit.h"
#include "marlinutil/DDMarlinCED.h"
#include "UTIL/TrackTools.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/MCParticle.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "marlinutil/MarlinUtil.h"
#include "IMPL/TrackImpl.h"
#include "MarlinTrk/MarlinTrkUtils.h"
using namespace TrackLengthDebugUtils;
using namespace EVENT;
using namespace UTIL;
using namespace IMPL;
using dd4hep::rec::Vector3D;
using std::pair;
using std::vector;
using std::string;
using MarlinTrk::IMarlinTrack;
using MarlinTrk::IMarlinTrkSystem;

TrackLengthDebug aTrackLengthDebug;


TrackLengthDebug::TrackLengthDebug() : marlin::Processor("TrackLengthDebug"){
    _description = "Processor that calculates track length with many different methods for comparison. Represents evolution of the track length calculation in iLCSoft";
}


void TrackLengthDebug::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);
    _bField = MarlinUtil::getBzAtOrigin();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();

    _file.reset( new TFile("results.root", "RECREATE") );
    _tree.reset( new TTree("treename", "treename") );

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("momentumIP", &_momentumIP);
    _tree->Branch("momentumCalo", &_momentumCalo);
    _tree->Branch("momentumHM", &_momentumHM);
    _tree->Branch("tof", &_tof);
    _tree->Branch("trackLength", &_trackLength);
    _tree->Branch("nBadPhi", &_nBadPhi);
    _tree->Branch("nBadZ", &_nBadZ);
    _tree->Branch("nSegments", &_nSegments);
}

void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    PIDHandler pidHandler( pfos );
    LCRelationNavigator nav ( evt->getCollection("RecoMCTruthLink") );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG7)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        if( nClusters != 1 || nTracks != 1) continue;

        const vector<LCObject*>& objects = nav.getRelatedToObjects(pfo);
        const std::vector<float>& weights = nav.getRelatedToWeights(pfo);
        int max_i = std::max_element(weights.begin(), weights.end(), [](float lhs, float rhs){return (int(lhs)%10000)/1000. < (int(rhs)%10000)/1000.;}) - weights.begin();
        if ( ( int(weights[max_i])%10000 )/1000. == 0 ){
            max_i = std::max_element(weights.begin(), weights.end(), [](float lhs, float rhs){return (int(lhs)/10000)/1000. < (int(rhs)/10000)/1000.;}) - weights.begin();
        }
        auto* mc = static_cast<MCParticle*> (objects[max_i]);
        Track* track = pfo->getTracks()[0];
        auto tsCalo = track->getTrackState(TrackState::AtCalorimeter);
        // store this to the TTree
        auto momAtCalo = getTrackMomentum(tsCalo, _bField);
        Vars results = getTrackLengthUsingZ(track, _bField, _trkSystem);
 
        _pdg = mc->getPDG();
        _momentumIP = std::hypot( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2] ); // in GeV
        _momentumCalo = std::hypot(momAtCalo[0], momAtCalo[1], momAtCalo[2]);
        _momentumHM = results.momentumHM;
        _tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "timeOfFlight"); // in ns
        _trackLength = results.trackLength;
        _nBadPhi = results.nBadPhi;
        _nBadZ = results.nBadZ;
        _nSegments = results.nSegments;
        _tree->Fill();
    }
}

void TrackLengthDebug::end(){
    _file->Write();
}

Vars TrackLengthDebug::getTrackLengthUsingZ(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    Vars vars;

    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
        covMatrix[0]  = 1e+06; //sigma_d0^2
        covMatrix[2]  = 100.; //sigma_phi0^2
        covMatrix[5]  = 0.00001; //sigma_omega^2
        covMatrix[9]  = 1e+06; //sigma_z0^2
        covMatrix[14] = 100.; //sigma_tanl^2
        double maxChi2PerHit = 100.;
        std::unique_ptr<MarlinTrk::IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        TrackStateImpl preFit = *track->getTrackState(TrackState::AtLastHit);
        preFit.setCovMatrix( covMatrix );
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::backward, &preFit, bField, maxChi2PerHit);
        //if fit fails, try also fit forward
        if (errorFit != 0){
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0) continue;

        //here hits are sorted by rho=(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        //Find which way to loop over the array of hits. We need to loop in the direction of the track.
        bool loopForward = true;
        double zFirst = std::abs( hitsInFit.front().first->getPosition()[2] );
        double zLast = std::abs( hitsInFit.back().first->getPosition()[2] );

        // OPTIMIZE: 10 mm is just a round number. With very small z difference it is more robust to use rho, to be sure z difference is not caused by tpc Z resolution or multiple scattering
        if ( std::abs(zLast - zFirst) > 10. ){
            if ( zLast < zFirst ) loopForward = false;
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }

    int nTrackStates = trackStates.size();
    vars.nSegments = nTrackStates - 1;
    if (nTrackStates <= 1) return vars;
    for( int j=1; j < nTrackStates; ++j ){
        auto ts1 = trackStates[j-1];
        auto ts2 = trackStates[j];
        double omega = ts1.getOmega();
        double tanL = ts1.getTanLambda();
        double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
        double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
        double phi1 = ts1.getPhi();
        double phi2 = ts2.getPhi();
        double dPhi = std::abs(phi2 - phi1);
        if (dPhi <= M_PI && (phi1 - phi2)/omega < 0) vars.nBadPhi++;
        else if (dPhi > M_PI && phi2/omega < 0) vars.nBadPhi++;

        if ( (z2-z1)/tanL < 0) vars.nBadZ++;

        double arcLength = std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
        auto mom = getTrackMomentum(&ts1, _bField);
        double AbsMomentum = std::hypot(mom[0], mom[1], mom[2]);
        vars.trackLength += arcLength;
        if(AbsMomentum == 0) continue;
        vars.momentumHM += arcLength/(AbsMomentum*AbsMomentum);
    }
    if (vars.momentumHM == 0) vars.momentumHM = 0.;
    else vars.momentumHM = std::sqrt(vars.trackLength/vars.momentumHM);

    return vars;
}
