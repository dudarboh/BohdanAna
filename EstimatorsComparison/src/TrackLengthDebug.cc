#include "TrackLengthDebug.h"
#include "TrackLengthDebugUtils.h"
#include "TrackLengthEstimators.h"

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
#include "marlinutil/CalorimeterHitType.h"
using namespace TrackLengthDebugUtils;
using namespace EVENT;
using namespace UTIL;
using dd4hep::rec::Vector3D;
using std::vector;
using std::string;

TrackLengthDebug aTrackLengthDebug;


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVirtualMemoryUsage(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


int getPhysicalMemoryUsage(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


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
    _tree->Branch("momentum", &_momentum);
    _tree->Branch("tof", &_tof);
    _tree->Branch("tofv2", &_tofv2);
    _tree->Branch("hitLayer", &_hitLayer);
    // _tree->Branch("trackLengthIDR", &_trackLengthIDR);
    // _tree->Branch("trackLengthIDR2", &_trackLengthIDR2);
    // _tree->Branch("trackLengthIDR3", &_trackLengthIDR3);
    // _tree->Branch("trackLengthIDR4", &_trackLengthIDR4);

    _tree->Branch("trackLengthSHA1", &_trackLengthSHA1);
    _tree->Branch("trackLengthSHA2", &_trackLengthSHA2);
    _tree->Branch("trackLengthSHA3", &_trackLengthSHA3);
    _tree->Branch("trackLengthSHA4", &_trackLengthSHA4);
    _tree->Branch("trackLengthSHA5", &_trackLengthSHA5);
    _tree->Branch("trackLengthSHA6", &_trackLengthSHA6);

    // _tree->Branch("trackLengthWinni", &_trackLengthWinni);
    // _tree->Branch("trackLengthWinni2", &_trackLengthWinni2);
    // _tree->Branch("trackLengthUsingZ", &_trackLengthUsingZ);
    // _tree->Branch("trackLengthUsingZ2", &_trackLengthUsingZ2);
    // _tree->Branch("trackLengthUsingZ3", &_trackLengthUsingZ3);
    // _tree->Branch("trackLengthSimUsingZ", &_trackLengthSimUsingZ);

    _tree->Branch("trackLengthIKF1", &_trackLengthIKF1);
    _tree->Branch("trackLengthIKF2", &_trackLengthIKF2);
    _tree->Branch("trackLengthIKF3", &_trackLengthIKF3);
}

void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;
    
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
        Cluster* cluster = pfo->getClusters()[0];

        // store this to the TTree
        _pdg = mc->getPDG();
        _momentum = std::hypot( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2] ); // in GeV
        _tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "timeOfFlight"); // in ns
        _trackLengthIDR = getTrackLengthIDR(track);
        _trackLengthIDR2 = getTrackLengthIDR2(track);
        _trackLengthIDR3 = getTrackLengthIDR3(track);
        _trackLengthIDR4 = getTrackLengthIDR4(track);

        _trackLengthSHA1 = getTrackLengthSHA1(track);
        _trackLengthSHA2 = getTrackLengthSHA2(track);
        _trackLengthSHA3 = getTrackLengthSHA3(track);
        _trackLengthSHA4 = getTrackLengthSHA4(track);
        _trackLengthSHA5 = getTrackLengthSHA5(track);
        _trackLengthSHA6 = getTrackLengthSHA6(track);

        _trackLengthWinni = getTrackLengthWinni(track, _bField, _trkSystem);
        _trackLengthWinni2 = getTrackLengthWinni2(track, _bField, _trkSystem);
        _trackLengthUsingZ = getTrackLengthUsingZ(track, _bField, _trkSystem);
        _trackLengthUsingZ2 = getTrackLengthUsingZ2(pfo, _bField, _trkSystem);
        _trackLengthUsingZ3 = getTrackLengthUsingZ3(pfo, _bField, _trkSystem);
        _trackLengthSimUsingZ = getTrackLengthSimUsingZ(track, _bField, _trkSystem);

        _trackLengthIKF1 = getTrackLengthIKF1(pfo, _bField, _trkSystem);
        _trackLengthIKF2 = getTrackLengthIKF2(pfo, _bField, _trkSystem);
        _trackLengthIKF3 = getTrackLengthIKF3(pfo, _bField, _trkSystem);

        // FILL TIME OF FLIGHT with LATEST track state at Calorimeter
        const TrackState* tsEcal = nullptr;
        if (track->getTracks().size() <= 2) tsEcal = track->getTrackState( TrackState::AtCalorimeter );
        else{
            Track* lastSubTrack = track->getTracks().back();
            tsEcal = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
        }
        Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );

        double hitTime = std::numeric_limits<double>::max();
        double closestDistance = std::numeric_limits<double>::max();
        for( auto hit : cluster->getCalorimeterHits() ){
            CHT hitInfo( hit->getType() );
            bool isECALHit = ( hitInfo.caloID() == CHT::ecal );
            if (! isECALHit) continue;

            Vector3D hitPos( hit->getPosition() );
            double dToTrack = (hitPos - trackPosAtEcal).r();
            if( dToTrack < closestDistance ){
                closestDistance = dToTrack;
                hitTime = hit->getTime();
                _hitLayer = hitInfo.layer();
            }
        }
        if ( hitTime == std::numeric_limits<double>::max() ) _tofv2 = 0.;
        else  _tofv2 = hitTime - closestDistance/CLHEP::c_light;



        _tree->Fill();
    }
}

void TrackLengthDebug::end(){
    _file->Write();
}