#include "TrackLengthDebug.h"
#include "TrackLengthUtils.h"

#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "MarlinTrk/Factory.h"
#include "EVENT/SimTrackerHit.h"
#include "marlinutil/DDMarlinCED.h"

using namespace TrackLengthUtils;
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;
using EVENT::TrackerHit;
using EVENT::Track;
using EVENT::TrackState;
using EVENT::SimTrackerHit;
using EVENT::TrackState;
using EVENT::LCObject;
using dd4hep::rec::Vector3D;

TrackLengthDebug aTrackLengthDebug ;


TrackLengthDebug::TrackLengthDebug() : marlin::Processor("TrackLengthDebug") {
    _description = "TrackLengthDebug debugs track length";

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "ReconstructedParticleCollection",
                            "Name of the ReconstructedParticle collection",
                            _pfoCollectionName,
                            std::string("PandoraPFOs") );

    registerProcessorParameter("eventDisplay",
                             "eventDisplay",
                             _eventDisplay,
                             bool(false));
}


void TrackLengthDebug::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);

    _outputParNames = {"trackLengthToSET", "trackLengthToEcal", "momentumHMToSET", "momentumHMToEcal"};
    _bField = MarlinUtil::getBzAtOrigin();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();

    prepareRootTree();
    if (_eventDisplay) DDMarlinCED::init(this);

}


void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(DEBUG9)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection(_pfoCollectionName);

    PIDHandler pidHandler( pfos );
    int algoID = pidHandler.addAlgorithm( name(), _outputParNames );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG9)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        if( nClusters != 1 || nTracks != 1) continue;

        Track* track = pfo->getTracks()[0];
        const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
        //only barrel to check
        if ( std::abs(tsCalo->getReferencePoint()[2] + tsCalo->getZ0() ) > 2300. ) continue;
        vector<Track*> subTracks = getSubTracks(track);
        vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, _bField);

        _trackLengthDefault = getTrackLengthDefault(trackStates);
        _trackLengthTanL = getTrackLengthTanL(trackStates);
        _trackLengthZ = getTrackLengthZ(trackStates);
    
        double speedOfLight = 299.792458;
        _momentum = Vector3D( pfo->getMomentum() ).r();
        _tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "timeOfFlight"); // in ns
        _massDefault = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthDefault, 2) - 1 );
        _massTanL = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthTanL, 2) - 1 );
        _massZ = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthZ, 2) - 1 );

        std::cout<<"Final results for the "<<i+1<<" PFO"<<std::endl;
        std::cout<<"Track length using default: "<<_trackLengthDefault<<" mm"<<std::endl;
        std::cout<<"Track length using tanL: "<<_trackLengthTanL<<" mm"<<std::endl;
        std::cout<<"Track length using Z: "<<_trackLengthZ<<" mm"<<std::endl;
        std::cout<<std::endl<<std::endl;

        _tree->Fill();

        if(_eventDisplay){
            DDMarlinCED::newEvent(this, evt);
            DDMarlinCED::drawDD4hepDetector(_detector, false, vector<string>{""});
            DDCEDPickingHandler& pHandler= DDCEDPickingHandler::getInstance();
            pHandler.update(evt);
            drawPFO(pfo);
            DDMarlinCED::draw(this);            
        }


        // EVENT DISPLAY SHENANIGENS
        // if (trackLength == bullshit) {
        //     MyAmazingEventDisplay evDisplay;
        //     evDisplay.draw()
        // }

    

    }
}


double TrackLengthDebug::getTrackLengthDefault(std::vector<IMPL::TrackStateImpl> trackStates){
    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    //exclude last track state at the ECal
    for( int j=1; j < nTrackStates; ++j ){
        //we check which track length formula to use
        double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
        double arcLength;
        // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
        if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
        else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

        trackLength += arcLength;
    }
    return trackLength;
}

double TrackLengthDebug::getTrackLengthTanL(std::vector<IMPL::TrackStateImpl> trackStates){
    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    //exclude last track state at the ECal
    for( int j=1; j < nTrackStates; ++j ){
        //we check which track length formula to use
        double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
        double arcLength;
        // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
        if ( nTurns <= 0.5 ) arcLength = getHelixArcLengthTanL( trackStates[j-1], trackStates[j] );
        else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

        trackLength += arcLength;
    }
    return trackLength;
}



double TrackLengthDebug::getTrackLengthZ(std::vector<IMPL::TrackStateImpl> trackStates){
    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    //exclude last track state at the ECal
    for( int j=1; j < nTrackStates; ++j ){
        //we check which track length formula to use
        double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
        double arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

        trackLength += arcLength;
    }
    return trackLength;
}


float TrackLengthDebug::getParameterFromPID(ReconstructedParticle* pfo, PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
    int algorithmID = pidHandler.getAlgorithmID(algorithmName);
    const ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
    const std::vector<float>& parameters = pfoPID.getParameters();
    int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
    return parameters[parIdx]; 
}



void TrackLengthDebug::prepareRootTree(){
    _file.reset( new TFile("results.root", "RECREATE") );
    _tree.reset( new TTree("treename", "treename") );

    _tree->Branch("momentum", &_momentum);
    _tree->Branch("tof", &_tof);
    _tree->Branch("trackLengthDefault", &_trackLengthDefault);
    _tree->Branch("trackLengthTanL", &_trackLengthTanL);
    _tree->Branch("trackLengthZ", &_trackLengthZ);
    _tree->Branch("massDefault", &_massDefault);
    _tree->Branch("massTanL", &_massTanL);
    _tree->Branch("massZ", &_massZ);
    
}

void TrackLengthDebug::drawPFO(ReconstructedParticle* pfo){
    std::vector<Track*> tracks = pfo->getTracks();
    for(auto* track: tracks){
        auto hits = track->getTrackerHits();
        for (auto* hit: hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int layer = 1; // doesn't matter
            int size = 3; // larger point
            unsigned long color = 0x3232a8;
            ced_hit_ID(pos[0], pos[1], pos[2], type, layer, size, color, 0 ); // tracker hits
        }
    }
    std::vector<Cluster*> clusters = pfo->getClusters();
    for(auto* cluster: clusters){
        auto hits = cluster->getCalorimeterHits();
        for (auto* hit: hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int layer = 1; // doesn't matter
            int size = 6; // larger point
            unsigned long color = 0xbf2659;
            ced_hit_ID(pos[0], pos[1], pos[2], type, layer, size, color, 0 ); // tracker hits
        }
    }
}
