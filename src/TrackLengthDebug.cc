#include "TrackLengthDebug.h"
#include "TrackLengthEstimators.h"

#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/MCParticle.h"
#include "UTIL/PIDHandler.h"
#include "UTIL/TrackTools.h"
#include "UTIL/LCRelationNavigator.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "MarlinTrk/Factory.h"
#include "marlinutil/DDMarlinCED.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "marlinutil/MarlinUtil.h"
#include "marlinutil/CalorimeterHitType.h"
#include "CLHEP/Random/Randomize.h"


#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"
#include "IMPL/SimTrackerHitImpl.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "UTIL/ILDConf.h"
#include "IMPL/TrackImpl.h"
#include <cmath>
#include <algorithm>
#include <limits>

// using namespace TrackLengthDebugUtils;
using namespace EVENT;
using namespace UTIL;
// using namespace IMPL;
using dd4hep::rec::Vector3D;
using CLHEP::RandGauss;
using std::vector;
using std::string;

TrackLengthDebug aTrackLengthDebug;

// might be useful at some point
float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
    int algorithmID = pidHandler.getAlgorithmID(algorithmName);
    const ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
    const std::vector<float>& parameters = pfoPID.getParameters();
    int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
    return parameters[parIdx]; 
}

//performance
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



double getHelixLength(double z_start, double z_end, Vector3D momentum, double bField){
    // theoretical formula to get the helix length
    double c_factor = 0.299792458;
    double pt = momentum.rho();
    double pz = momentum.z();
    return (z_end - z_start)/pz * std::sqrt(  std::pow( pt/(c_factor*bField), 2) + pz*pz  );
}

// drawing
void drawPFO(ReconstructedParticle* pfo, TrackStateImpl tsStdReco, TrackStateImpl tsEasy){
    vector<Track*> tracks;
    if ( not pfo->getTracks().empty() )
        tracks = getSubTracks(pfo->getTracks()[0]);
    for(auto* track: tracks){
        auto hits = track->getTrackerHits();
        for (auto* hit : hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int size = 4; // larger point
            unsigned long color = 0x000000;
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }
    //plot ecal state for fun
    const TrackState* tsCalo = tracks.back()->getTrackState( TrackState::AtCalorimeter );
    auto posCalo = tsCalo->getReferencePoint();
    ced_hit(posCalo[0], posCalo[1], posCalo[2] + tsCalo->getZ0(), 0, 4, 0x000000); // track state at calorimeterer


    std::vector<Cluster*> clusters = pfo->getClusters();
    for(auto* cluster: clusters){
        auto hits = cluster->getCalorimeterHits();
        for (auto* hit: hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int size = 4; // larger point
            unsigned long color = 0x000000;
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }

    auto printStateLong = [](TrackStateImpl ts){
        auto part = ts;
        std::stringstream tmp;
        //out << std::scientific << std::setprecision (2) << std::showpos;
        std::cout << std::noshowpos;
        std::cout << std::setw(41) << std::setfill('-') << std::right << "-- TrackState ---" << std::setfill('-') << std::setw(29) << "-" << std::endl;

        tmp.str("") ;
        switch( part.getLocation() ){
        case EVENT::TrackState::AtOther         :   tmp <<  "AtOther"        ;   break ;
        case EVENT::TrackState::AtIP            :   tmp <<  "AtIP"           ;   break ;
        case EVENT::TrackState::AtFirstHit      :   tmp <<  "AtFirstHit"     ;   break ;
        case EVENT::TrackState::AtLastHit       :   tmp <<  "AtLastHit"      ;   break ;
        case EVENT::TrackState::AtCalorimeter   :   tmp <<  "AtCalorimeter " ;   break ;
        case EVENT::TrackState::AtVertex        :   tmp <<  "AtVertex"       ;   break ;
        }
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Location" << std::right << std::setw(40) << tmp.str() << std::endl;
        tmp.str("") ;
        tmp << std::dec << std::setfill('0') << std::setw(8) << part.id();
        std::cout << std::scientific << std::setprecision(6) ;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Id"          << std::right << std::setw(40) << tmp.str() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "D0"          << std::right << std::setw(40) << part.getD0() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Phi"         << std::right << std::setw(40) << part.getPhi() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Omega"       << std::right << std::setw(40) << part.getOmega() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Z0"          << std::right << std::setw(40) << part.getZ0() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Tan Lambda"  << std::right << std::setw(40) << part.getTanLambda() << std::endl;
        tmp.str("");
        tmp  << std::dec << part.getReferencePoint()[0] << ", " << part.getReferencePoint()[1]  << ", " << part.getReferencePoint()[2]; 
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "ReferencePoint" << std::right << std::setw(40) << tmp.str() << std::endl;
        std::cout << "Cov matrix:" << std::showpos << std::scientific << std::setprecision(6) << std::setw(15) << std::setfill(' ')  ;
        // print cov matrix as lower triangle matrix 
        for( unsigned l=0 , N=part.getCovMatrix().size(), ncolumns = 1 , nele =1 ; l <N ; ++l , ++nele) {
        std::cout << part.getCovMatrix()[l];
        if(! ( (nele) % ncolumns ) ){ 
            nele = 0 ;
            ++ncolumns ;
            std::cout << std::endl << "             " ;
        } else {
            std::cout << ", ";
        } 
        }
        std::cout << std::noshowpos;
        std::cout<<std::endl;
    };

    auto plotHelixFromTrackState = [](TrackStateImpl ts, unsigned long color){
        double omega = ts.getOmega();
        double tanL = ts.getTanLambda();
        double phi = ts.getPhi();
        double d0 = ts.getD0();
        double z0 = ts.getZ0();
        auto ref = ts.getReferencePoint();
        std::array<double, 3> mom = UTIL::getTrackMomentum(&ts, 3.5);

        //
        int charge = omega > 0 ? 1 : -1;
        double x = -d0*std::sin(phi) + ref[0];
        double y = d0*std::cos(phi) + ref[1];
        double z = z0 + ref[2];


        int type = 0; // point
        int size = 8; // larger point
        ced_hit(x, y, z, type, size, color);
        size = 2;
        DDMarlinCED::drawHelix( 3.5, charge, x, y, z, mom[0], mom[1], mom[2], type, size, color, 0, 2000, 2600, 0);
        DDMarlinCED::drawHelix( -3.5, -charge, x, y, z, mom[0], mom[1], mom[2], type, size, color, 0, 2000, 2600, 0);
    };

    plotHelixFromTrackState(tsStdReco, 0xfc0505);
    printStateLong(tsStdReco);
    plotHelixFromTrackState(tsEasy, 0x0905fc);
    printStateLong(tsEasy);

}


TOFResults getTOFClosest(ReconstructedParticle* pfo){
    TOFResults results;
    if ( pfo->getTracks().empty() || pfo->getClusters().empty() ) return results;
    Track* track = pfo->getTracks()[0];
    Cluster* cluster = pfo->getClusters()[0];

    auto isTPCHit = [](TrackerHit* hit) -> bool {
        UTIL::BitField64 encoder( LCTrackerCellID::encoding_string() ) ;
        encoder.setValue( hit->getCellID0() ) ;
        int subdet = encoder[ LCTrackerCellID::subdet() ];
        return subdet == UTIL::ILDDetID::TPC;
    };

    int indexOfFirstTPCCurl = 0;
    int nSubTracks = track->getTracks().size();
    for(int i = 0; i < nSubTracks; ++i){
        Track* subTrack = track->getTracks()[i];
        auto hits = subTrack->getTrackerHits();
        if ( std::find_if(hits.begin(), hits.end(), isTPCHit) != hits.end() ){
            indexOfFirstTPCCurl = i;
            break;
        }
    }

    const TrackState* tsEcal = nullptr;
    if ( indexOfFirstTPCCurl == nSubTracks-1 ) tsEcal = track->getTrackState( TrackState::AtCalorimeter );
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
            results.layer = hitInfo.layer();
            hitTime = hit->getTime();
            closestDistance = dToTrack;
        }
    }
    if ( hitTime != std::numeric_limits<double>::max() ){
        results.dToTrack = closestDistance;
        results.tof0 = hitTime - closestDistance/CLHEP::c_light;
        results.tof10 = RandGauss::shoot(hitTime, 0.01) - closestDistance/CLHEP::c_light;
        results.tof20 = RandGauss::shoot(hitTime, 0.02) - closestDistance/CLHEP::c_light;
        results.tof30 = RandGauss::shoot(hitTime, 0.03) - closestDistance/CLHEP::c_light;
        results.tof40 = RandGauss::shoot(hitTime, 0.04) - closestDistance/CLHEP::c_light;
        results.tof50 = RandGauss::shoot(hitTime, 0.05) - closestDistance/CLHEP::c_light;
        results.tof60 = RandGauss::shoot(hitTime, 0.06) - closestDistance/CLHEP::c_light;
        results.tof70 = RandGauss::shoot(hitTime, 0.07) - closestDistance/CLHEP::c_light;
        results.tof80 = RandGauss::shoot(hitTime, 0.08) - closestDistance/CLHEP::c_light;
        results.tof90 = RandGauss::shoot(hitTime, 0.09) - closestDistance/CLHEP::c_light;
        results.tof100 = RandGauss::shoot(hitTime, 0.1) - closestDistance/CLHEP::c_light;
    }
    return results;
}

TrackLengthDebug::TrackLengthDebug() : marlin::Processor("TrackLengthDebug"), EventDisplayer(this){
    _description = "Main analysis of the track length and time-of-flight and momentum methods for time-of-flight pID";
}


void TrackLengthDebug::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);
    // DDMarlinCED::init(this);
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
    _tree->Branch("dToTrack", &_tofClosest.dToTrack);
    _tree->Branch("layer", &_tofClosest.layer);
    _tree->Branch("tofClosest0", &_tofClosest.tof0);
    _tree->Branch("tofClosest10", &_tofClosest.tof10);
    _tree->Branch("tofClosest20", &_tofClosest.tof20);
    _tree->Branch("tofClosest30", &_tofClosest.tof30);
    _tree->Branch("tofClosest40", &_tofClosest.tof40);
    _tree->Branch("tofClosest50", &_tofClosest.tof50);
    _tree->Branch("tofClosest60", &_tofClosest.tof60);
    _tree->Branch("tofClosest70", &_tofClosest.tof70);
    _tree->Branch("tofClosest80", &_tofClosest.tof80);
    _tree->Branch("tofClosest90", &_tofClosest.tof90);
    _tree->Branch("tofClosest100", &_tofClosest.tof100);


    _tree->Branch("trackLengthSHA1", &_trackLengthSHA1);
    _tree->Branch("trackLengthSHA2", &_trackLengthSHA2);
    _tree->Branch("trackLengthSHA3", &_trackLengthSHA3);
    _tree->Branch("trackLengthSHA4", &_trackLengthSHA4);
    _tree->Branch("trackLengthSHA5", &_trackLengthSHA5);
    _tree->Branch("trackLengthSHA6", &_trackLengthSHA6);
    _tree->Branch("trackLengthIKF1", &_trackLengthIKF1);
    _tree->Branch("trackLengthIKF2Bug", &_trackLengthIKF2Bug);
    _tree->Branch("trackLengthIKF2", &_trackLengthIKF2);
    _tree->Branch("trackLengthIKF3", &_trackLengthIKF3);

}

void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;

    //Do I need this !? I guess not!
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);


    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator nav ( evt->getCollection("RecoMCTruthLink") );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG7)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        // only simple cases for now
        if( nClusters != 1 || nTracks != 1) continue;
        Track* track = pfo->getTracks()[0];
        Cluster* cluster = pfo->getClusters()[0];

        auto getMC = [&pfo, &nav]() -> MCParticle* {
            const vector<LCObject*>& objects = nav.getRelatedToObjects(pfo);
            const std::vector<float>& weights = nav.getRelatedToWeights(pfo);

            auto getTrackWeight = [](float encodedWeight){ return float( int(encodedWeight) % 10000 ) / 1000.f;};
            auto getClusterWeight = [](float encodedWeight){ return float( int(encodedWeight) / 10000 ) / 1000.f;};

            int max_i = std::max_element(weights.begin(), weights.end(), [getTrackWeight](float lhs, float rhs){return getTrackWeight(lhs) < getTrackWeight(rhs);}) - weights.begin();
            if (getTrackWeight(max_i) == 0.f){
                max_i = std::max_element(weights.begin(), weights.end(), [getClusterWeight](float lhs, float rhs){return getClusterWeight(lhs) < getClusterWeight(rhs);}) - weights.begin();
            }
            return static_cast<MCParticle*> (objects[max_i]);
        };
        MCParticle* mc = getMC();

        // Fill all in the TTree
        _pdg = mc->getPDG();
        _tofClosest = getTOFClosest(pfo);
        _momentumIP = std::hypot( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2] ); // in GeV
        const TrackState* tsEcal = nullptr;
        if (track->getTracks().size() <= 2) tsEcal = track->getTrackState( TrackState::AtCalorimeter );
        else{
            Track* lastSubTrack = track->getTracks().back();
            tsEcal = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
        }
        std::array<double, 3> momAtCalo = getTrackMomentum(tsEcal, _bField);
        _momentumCalo = std::hypot(momAtCalo[0], momAtCalo[1], momAtCalo[2]);

        _trackLengthSHA1 = getTrackLengthSHA1(track);
        _trackLengthSHA2 = getTrackLengthSHA2(track);
        _trackLengthSHA3 = getTrackLengthSHA3(track);
        _trackLengthSHA4 = getTrackLengthSHA4(track);
        _trackLengthSHA5 = getTrackLengthSHA5(track);
        _trackLengthSHA6 = getTrackLengthSHA6(track);
        _trackLengthIKF2Bug = getTrackLengthIKF2Bug(pfo, _bField, _trkSystem);
        vector<TrackStateImpl> trackStates = getTrackStates(pfo, _bField, _trkSystem);
        _trackLengthIKF1 = getTrackLengthIKF1(trackStates);
        _trackLengthIKF2 = getTrackLengthIKF2(trackStates);
        _trackLengthIKF3 = getTrackLengthIKF3(trackStates);


        _momentumHM = 0.;
        int nTrackStates = trackStates.size();
        if (nTrackStates > 1){
            //IKF2 best and greatest
            auto getHelixLength = [](const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
                double tanL = ts1.getTanLambda();
                double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
                double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
                return std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
            };

            for( int j=1; j < nTrackStates; ++j ){
                double arcLength = getHelixLength( trackStates[j-1], trackStates[j] );
                std::array<double, 3> momArr = getTrackMomentum( &(trackStates[j-1]), _bField );
                Vector3D mom(momArr[0], momArr[1], momArr[2]);
                _momentumHM += arcLength/mom.r2();
            }
            _momentumHM = std::sqrt(_trackLengthIKF2/_momentumHM);
        }

        _tree->Fill();
        const TrackState* tsLastHit = nullptr;
        if (track->getTracks().size() <= 2) tsLastHit = track->getTrackState( TrackState::AtLastHit );
        else{
            Track* lastSubTrack = track->getTracks().back();
            tsLastHit = lastSubTrack->getTrackState( TrackState::AtLastHit );
        }


        TrackStateImpl tsStdReco = static_cast<TrackStateImpl> (*tsLastHit); // OUTERMOST!!!!!!! hit
        TrackStateImpl tsEasy = trackStates.rbegin()[1]; // last hit
        drawDisplay(this, evt, drawPFO, pfo, tsStdReco, tsEasy);
        // if(EventDisplay)
        //     printout();
    }
}

void TrackLengthDebug::end(){
    _file->Write();
}

void TrackLengthDebug::printout(){
    // std::cout<<"Event: "<<_nEvent<<std::endl;
    // std::cout<<"PFO: "<<i+1<<std::endl;
    // std::cout<<"PDG: "<<_pdg<<std::endl;
    // std::cout<<"Vertex: "<<mcPos.x()<<",  "<<mcPos.y()<<",  "<<mcPos.z()<<std::endl;
    // std::cout<<"Vertex radius: "<<mcPos.r()<<std::endl;
    // std::cout<<"MC momentum IP : ("<<mc->getMomentum()[0]<<",  "<<mc->getMomentum()[1]<<",  "<<mc->getMomentum()[2]<<")"<<std::endl;
    // std::cout<<"Reco momentum IP : ("<<pfo->getMomentum()[0]<<",  "<<pfo->getMomentum()[1]<<",  "<<pfo->getMomentum()[2]<<")"<<std::endl;
    // std::cout<<"MC abs p IP: "<<Vector3D(mc->getMomentum()).r()<<std::endl;
    // std::cout<<"Reco abs p IP: "<<Vector3D(pfo->getMomentum()).r()<<std::endl;
    // std::cout<<"HM momentumOld: "<<_momentumOld<<std::endl;
    // std::cout<<"HM momentumNew: "<<_momentumNew<<std::endl;
    // std::cout<<"TOF OLD: "<<_tofOld<<std::endl;
    // std::cout<<"TOF New: "<<_tofNew<<std::endl;
    // std::cout<<"trackLength Old: "<<_trackLengthOld<<std::endl;
    // std::cout<<"trackLength New: "<<_trackLengthNew<<std::endl;
    // std::cout<<"trackLength SHOULD BE for true mass: "<<trackLengthTrue<<std::endl;
    // std::cout<<"TOF SHOULD BE for true mass: "<<tofTrue<<std::endl;
    // std::cout<<"trackLength from helix formula reco: "<<getHelixLength(0., _tsFirstCaloPos.z(), Vector3D(momAtCalo[0], momAtCalo[1], momAtCalo[2]), 3.5)<<std::endl;
    // std::cout<<"trackLength from helix formula: "<<getHelixLength(mc->getVertex()[2], _tsFirstCaloPos.z(), mcMom, 3.5)<<std::endl;
    // std::cout<<"nSegments: "<<_nSegments<<std::endl;
    // std::cout<<"nBadZ: "<<_nBadZ<<std::endl;
    // std::cout<<"tsCaloOld: "<<_tsFirstCaloPos<<std::endl;
    // std::cout<<"tsCaloNew: "<<_tsLastCaloPos<<std::endl;
    // std::cout<<"Particle BORN time: "<<mc->getTime()<<std::endl;
    // if (_tofNew != 0 && _tofOld != 0){
    //     std::cout<<"NEW beta: "<<_trackLengthNew/(_tofNew*299.792458)<<std::endl;
    //     std::cout<<"OLD beta: "<<_trackLengthOld/(_tofOld*299.792458)<<std::endl;
    // }
    // if (_trackLengthNew != 0 && _trackLengthOld != 0){
    //     std::cout<<"NEW mass: "<<_momentumNew*std::sqrt( (_tofNew*299.792458*_tofNew*299.792458)/(_trackLengthNew*_trackLengthNew) - 1.)*1000<<std::endl;
    //     std::cout<<"OLD mass: "<<massOld<<std::endl;
    // }

    return;

}