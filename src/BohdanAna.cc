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
    streamlog_out(MESSAGE)<<"==========Event========== "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG7)<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        // only simple cases for now
        if( nTracks > 1 || nClusters > 1) continue;

        // Do everything mc particle related
        MCParticle* mc = getMC(pfo, pfo2mc);
        _pdg = mc->getPDG();

        // // Do everything track related
        // if (nTracks != 0){
        //     Track* track = pfo->getTracks()[0];        
        //     auto tsCalo = getTrackStateAtCalorimeter(track);
        //     Vector3D trackPosAtcalo( tsCalo->getReferencePoint() );

        // }

        // // Do everything cluster related
        // if ( nClusters != 0 ){
        //     Cluster* cluster = pfo->getClusters()[0];
        //     _tofClosest = getTofClosest(cluster, trackPosAtcalo, 0.);
        // }

        if (mc->getPDG() == 22){
            std::cout<<"Time: "<<getTofPhotonTrue(mc)<<std::endl;
            getchar();
        }

        // Fill all in the TTree
        // drawDisplay(this, evt, drawPFO, pfo, tsStdReco, tsEasy);
    }
}

void BohdanAna::end(){
    _file->Write();
}