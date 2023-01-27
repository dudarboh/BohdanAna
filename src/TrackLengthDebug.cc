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
#include "UTIL/TrackTools.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/MCParticle.h"

using namespace TrackLengthUtils;
using namespace EVENT;
using namespace UTIL;
using dd4hep::rec::Vector3D;
using std::vector;
using std::string;

TrackLengthDebug aTrackLengthDebug ;


TrackLengthDebug::TrackLengthDebug() : marlin::Processor("TrackLengthDebug"), EventDisplayer(this) {
    _description = "TrackLengthDebug debugs track length";

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "ReconstructedParticleCollection",
                            "Name of the ReconstructedParticle collection",
                            _pfoCollectionName,
                            std::string("PandoraPFOs") );
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
}


void drawPFO(ReconstructedParticle* pfo){
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

void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(DEBUG8)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection(_pfoCollectionName);

    PIDHandler pidHandler( pfos );
    int algoID = pidHandler.addAlgorithm( name(), _outputParNames );

    LCRelationNavigator nav ( evt->getCollection("RecoMCTruthLink") );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG7)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        const MCParticle* mc = nullptr;
        double trackWeight = nav.getRelatedToMaxWeight(pfo, "track");
        if(trackWeight != 0) mc = static_cast<const MCParticle*> ( nav.getRelatedToMaxWeightObject(pfo, "track") );
        else mc = static_cast<const MCParticle*> ( nav.getRelatedToMaxWeightObject(pfo, "cluster") );

        _pdg = mc->getPDG();

        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        if( nClusters != 1 || nTracks != 1) continue;

        Track* track = pfo->getTracks()[0];

        _d0Ip = track->getD0();
        _z0Ip = track->getZ0();
        _omegaIp = track->getOmega();
        _tanLIp = track->getTanLambda();
        _phiIp = track->getPhi();

        const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
        _zCalo = tsCalo->getReferencePoint()[2] + tsCalo->getZ0();
        std::array<double, 3> momCalo = getTrackMomentum(tsCalo, _bField);
        _momCalo = Vector3D(momCalo[0], momCalo[1], momCalo[2]);
        _d0Calo = tsCalo->getD0();
        _z0Calo = tsCalo->getZ0();
        _omegaCalo = tsCalo->getOmega();
        _tanLCalo = tsCalo->getTanLambda();
        _phiCalo = tsCalo->getPhi();
        //only barrel to check
        // if ( std::abs(tsCalo->getReferencePoint()[2] + tsCalo->getZ0() ) > 2300. ) continue;
        vector<Track*> subTracks = getSubTracks(track);
        vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, _bField);

        _trackLengthDefault = getTrackLengthDefault(trackStates);
        _trackLengthTanL = getTrackLengthTanL(trackStates);
        _trackLengthZ = getTrackLengthZ(trackStates);
    
        double speedOfLight = 299.792458;
        _momIp = Vector3D( pfo->getMomentum() );
        _momentum = _momIp.r();
        _tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "timeOfFlight"); // in ns
        if ( std::pow(_tof*speedOfLight/_trackLengthDefault, 2) - 1 > 0 ) _massDefault = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthDefault, 2) - 1 );
        else _massDefault = 0;

        if ( std::pow(_tof*speedOfLight/_trackLengthTanL, 2) - 1 > 0 ) _massTanL = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthTanL, 2) - 1 );
        else _massTanL = 0;

        if ( std::pow(_tof*speedOfLight/_trackLengthZ, 2) - 1 > 0 ) _massZ = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthZ, 2) - 1 );
        else _massZ = 0;

        streamlog_out(DEBUG5)<<"Final results for the "<<i+1<<" PFO"<<std::endl;
        streamlog_out(DEBUG5)<<"Track length using default: "<<_trackLengthDefault<<" mm"<<std::endl;
        streamlog_out(DEBUG5)<<"Track length using tanL: "<<_trackLengthTanL<<" mm"<<std::endl;
        streamlog_out(DEBUG5)<<"Track length using Z: "<<_trackLengthZ<<" mm"<<std::endl;
        streamlog_out(DEBUG5)<<std::endl<<std::endl;

        drawDisplay(this, evt, drawPFO, pfo);
        _tree->Fill();


    }
}

void TrackLengthDebug::end(){
    _file->Write();
}

double TrackLengthDebug::getTrackLengthDefault(const std::vector<IMPL::TrackStateImpl>& trackStates){
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

double TrackLengthDebug::getTrackLengthTanL(const std::vector<IMPL::TrackStateImpl>& trackStates){
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



double TrackLengthDebug::getTrackLengthZ(const std::vector<IMPL::TrackStateImpl>& trackStates){
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

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("mom_ip", &_momIp);
    _tree->Branch("d0_ip", &_d0Ip);
    _tree->Branch("z0_ip", &_z0Ip);
    _tree->Branch("omega_ip", &_omegaIp);
    _tree->Branch("tanL_ip", &_tanLIp);
    _tree->Branch("phi_ip", &_phiIp);

    _tree->Branch("mom_calo", &_momCalo);
    _tree->Branch("d0_calo", &_d0Calo);
    _tree->Branch("z0_calo", &_z0Calo);
    _tree->Branch("omega_calo", &_omegaCalo);
    _tree->Branch("tanL_calo", &_tanLCalo);
    _tree->Branch("phi_calo", &_phiCalo);

    _tree->Branch("z_calo", &_zCalo);
    
}
