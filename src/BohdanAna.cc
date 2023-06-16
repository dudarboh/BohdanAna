#include "BohdanAna.h"
#include "BohdanDrawing.h"
#include "BohdanUtils.h"
#include "TOF.h"
#include "TrackLength.h"

#include "marlinutil/GeometryUtil.h"
#include "UTIL/LCRelationNavigator.h"
#include "MarlinTrk/Factory.h"

using dd4hep::rec::Vector3D;

BohdanAna theBohdanAna;

BohdanAna::BohdanAna() : marlin::Processor("BohdanAna"), EventDisplayer(this){
    _description = "Main analysis of the track length and time-of-flight and momentum methods for time-of-flight pID";
}


void BohdanAna::init(){
    // DDMarlinCED::init(this);
    _bField = MarlinUtil::getBzAtOrigin();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();

    _file.reset( new TFile("results.root", "RECREATE") );
    _tree.reset( new TTree("treename", "treename") );

    // _tree->Branch("pdg", &_pdg);
    // _tree->Branch("momentumIP", &_momentumIP);
    // _tree->Branch("momentumCalo", &_momentumCalo);
    // _tree->Branch("momentumHM", &_momentumHM);
    // _tree->Branch("dToTrack", &_tofClosest.dToTrack);
    // _tree->Branch("layer", &_tofClosest.layer);
    // _tree->Branch("tofClosest0", &_tofClosest.tof0);
    // _tree->Branch("tofClosest10", &_tofClosest.tof10);
    // _tree->Branch("tofClosest20", &_tofClosest.tof20);
    // _tree->Branch("tofClosest30", &_tofClosest.tof30);
    // _tree->Branch("tofClosest40", &_tofClosest.tof40);
    // _tree->Branch("tofClosest50", &_tofClosest.tof50);
    // _tree->Branch("tofClosest60", &_tofClosest.tof60);
    // _tree->Branch("tofClosest70", &_tofClosest.tof70);
    // _tree->Branch("tofClosest80", &_tofClosest.tof80);
    // _tree->Branch("tofClosest90", &_tofClosest.tof90);
    // _tree->Branch("tofClosest100", &_tofClosest.tof100);


    // _tree->Branch("trackLengthSHA1", &_trackLengthSHA1);
    // _tree->Branch("trackLengthSHA2", &_trackLengthSHA2);
    // _tree->Branch("trackLengthSHA3", &_trackLengthSHA3);
    // _tree->Branch("trackLengthSHA4", &_trackLengthSHA4);
    // _tree->Branch("trackLengthSHA5", &_trackLengthSHA5);
    // _tree->Branch("trackLengthSHA6", &_trackLengthSHA6);
    // _tree->Branch("trackLengthIKF1", &_trackLengthIKF1);
    // _tree->Branch("trackLengthIKF2Bug", &_trackLengthIKF2Bug);
    // _tree->Branch("trackLengthIKF2", &_trackLengthIKF2);
    // _tree->Branch("trackLengthIKF3", &_trackLengthIKF3);

}

void BohdanAna::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator nav ( evt->getCollection("RecoMCTruthLink") );

    auto getMC = [&nav](EVENT::ReconstructedParticle* pfo) -> MCParticle* {
        const std::vector<LCObject*>& objects = nav.getRelatedToObjects(pfo);
        const std::vector<float>& weights = nav.getRelatedToWeights(pfo);

        auto getTrackWeight = [](float encodedWeight){ return float( int(encodedWeight) % 10000 ) / 1000.f;};
        auto getClusterWeight = [](float encodedWeight){ return float( int(encodedWeight) / 10000 ) / 1000.f;};

        int max_i = std::max_element(weights.begin(), weights.end(), [getTrackWeight](float lhs, float rhs){return getTrackWeight(lhs) < getTrackWeight(rhs);}) - weights.begin();
        if (getTrackWeight(max_i) == 0.f){
            max_i = std::max_element(weights.begin(), weights.end(), [getClusterWeight](float lhs, float rhs){return getClusterWeight(lhs) < getClusterWeight(rhs);}) - weights.begin();
        }
        return static_cast<MCParticle*> (objects[max_i]);
    };


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG7)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        // only simple cases for now
        if( nClusters != 1 || nTracks != 1) continue;
        Track* track = pfo->getTracks()[0];
        Cluster* cluster = pfo->getClusters()[0];

        MCParticle* mc = getMC(pfo);

        // Fill all in the TTree
        auto tsCalo = getTrackStateAtCalorimeter(track);
        Vector3D trackPosAtcalo( tsCalo->getReferencePoint() );
        _pdg = mc->getPDG();
        _tofClosest = getTofClosest(cluster, trackPosAtcalo, 0.);
        // drawDisplay(this, evt, drawPFO, pfo, tsStdReco, tsEasy);
    }
}

void BohdanAna::end(){
    _file->Write();
}