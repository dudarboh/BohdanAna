#include "BohdanAna.h"
#include "BohdanDrawing.h"
#include "BohdanUtils.h"
#include "TOF.h"
#include "TrackLength.h"

#include "marlinutil/GeometryUtil.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/TrackTools.h"
#include "MarlinTrk/Factory.h"

using dd4hep::rec::Vector3D;

BohdanAna theBohdanAna;

BohdanAna::BohdanAna() : marlin::Processor("BohdanAna"), EventDisplayer(this){
    _description = "Main analysis of the track length and time-of-flight and momentum methods for time-of-flight pID";
}


void BohdanAna::init(){
    _bField = MarlinUtil::getBzAtOrigin();
    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();

    _file.reset( new TFile("results.root", "RECREATE") );
    _tree.reset( new TTree("treename", "treename") );

    //mc
    _tree->Branch("pdg", &_pdg);

    //momenta
    _tree->Branch("mcPx", &(_mcMom[0]) );
    _tree->Branch("mcPy", &(_mcMom[1]) );
    _tree->Branch("mcPz", &(_mcMom[2]) );
    _tree->Branch("recoIpPx", &(_recoIpMom[0]) );
    _tree->Branch("recoIpPy", &(_recoIpMom[1]) );
    _tree->Branch("recoIpPz", &(_recoIpMom[2]) );
    _tree->Branch("recoCaloPx", &(_recoCaloMom[0]) );
    _tree->Branch("recoCaloPy", &(_recoCaloMom[1]) );
    _tree->Branch("recoCaloPz", &(_recoCaloMom[2]) );
    _tree->Branch("harmonicMomToEcal", &_trackLength.harmonicMomToEcal);
    _tree->Branch("harmonicMomToSET", &_trackLength.harmonicMomToSET);

    //track lengths
    _tree->Branch("trackLengthToEcal", &_trackLength.trackLengthToEcal);
    _tree->Branch("trackLengthToSET", &_trackLength.trackLengthToSET);

    //tofs
    _tree->Branch("layerClosest", &_layerClosest);
    for (int i = 0; i < 11; i++){
        _tree->Branch(( "tofClosest"+std::to_string(i*10) ).c_str(), &( _tofClosest[i]) );
        _tree->Branch(( "tofAverage"+std::to_string(i*10) ).c_str(), &( _tofAverage[i]) );
        _tree->Branch(( "tofSET"+std::to_string(i*10) ).c_str(), &( _tofSET[i]) );
        _tree->Branch(( "tofFit"+std::to_string(i*10) ).c_str(), &( _tofFit[i]) );
    }

}

void BohdanAna::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<"==========Event========== "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG8)<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        resetVariables();
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nTracks = pfo->getTracks().size();
        int nClusters = pfo->getClusters().size();
        // only simple cases
        if( nTracks > 1 || nClusters != 1) continue;
        Cluster* cluster = pfo->getClusters().at(0);
        MCParticle* mc = getMC(pfo, pfo2mc);
        _pdg = mc->getPDG();
        for(int j=0; j<3; j++) _mcMom.at(j) = mc->getMomentum()[j];

        bool isHadron = std::abs(_pdg) == 211 || std::abs(_pdg) == 321 || std::abs(_pdg) == 2212;
        bool isPhoton = std::abs(_pdg) == 22;
        if (isHadron && nTracks == 1){
            Track* track = pfo->getTracks().at(0);
            auto tsCalo = getTrackStateAtCalorimeter( track );
            Vector3D trackPosAtCalo( tsCalo->getReferencePoint() );
            std::array<double, 3> mom = UTIL::getTrackMomentum(tsCalo, _bField);
            Vector3D trackMomAtCalo(mom[0], mom[1], mom[2]);

            for(int j=0; j<3; j++) _recoIpMom.at(j) = pfo->getMomentum()[j];
            for(int j=0; j<3; j++) _recoCaloMom.at(j) = trackMomAtCalo[j];

            std::vector<IMPL::TrackStateImpl> trackStates = getTrackStates(pfo, _bField, _trkSystem);
            _trackLength = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::zedLambda);

            _layerClosest = getTofClosest(cluster, trackPosAtCalo, 0.).first;

            auto selectedHits = selectFrankEcalHits(cluster, trackPosAtCalo, trackMomAtCalo, 10);
            for (int j = 0; j < 11; j++){
                _tofClosest.at(j) = getTofClosest(cluster, trackPosAtCalo, 0.01*j).second;
                _tofAverage.at(j) = getTofFrankAvg(selectedHits, trackPosAtCalo, 0.01*j);
                _tofFit.at(j) = getTofFrankFit(selectedHits, trackPosAtCalo, 0.01*j);
                _tofSET.at(j) = getTofSET(track, 0.01*j);
            }

            // plotCanvas(cluster, trackPosAtCalo, trackMomAtCalo, mc);

        }
        else if( isPhoton && nTracks == 0 && ( !mc->isDecayedInTracker() ) ) {
            Vector3D photonPosAtCalo = getPhotonAtCalorimeter(mc);
            Vector3D mom( mc->getMomentum() );

            _layerClosest = getTofClosest(cluster, photonPosAtCalo, 0.).first;
            auto selectedHits = selectFrankEcalHits(cluster, photonPosAtCalo, mom, 10);
            for (int j = 0; j < 11; j++){
                _tofClosest.at(j) = getTofClosest(cluster, photonPosAtCalo, 0.01*j).second;
                _tofAverage.at(j) = getTofFrankAvg(selectedHits, photonPosAtCalo, 0.01*j);
                _tofFit.at(j) = getTofFrankFit(selectedHits, photonPosAtCalo, 0.01*j);
                // no track - no SET!
            }
        }
        else{
            continue;
        }
        _tree->Fill();

        // Fill all in the TTree
        drawDisplay(this, evt, displayPFO, pfo);
    }
    // drawDisplay(this, evt, drawFTDSimHits, evt);
}

void BohdanAna::end(){
    _file->Write();
    _application.Run(true);
}

void BohdanAna::resetVariables(){
    _pdg = 0;
    _mcMom.fill(0.);
    _recoIpMom.fill(0.);
    _recoCaloMom.fill(0.);

    _trackLength = TrackLengthResult();

    _layerClosest = -1;
    _tofClosest.fill(0.);
    _tofAverage.fill(0.);
    _tofSET.fill(0.);
    _tofFit.fill(0.);
}