#include "BohdanAna.h"
#include "BohdanDrawing.h"
#include "BohdanUtils.h"
#include "TOF.h"
#include "TrackLength.h"

#include "marlinutil/GeometryUtil.h"
#include "marlinutil/CalorimeterHitType.h"
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

    //mc
    auto _pdg = _model->MakeField<int>("pdg");

    //track parameters
    auto _dEdx = _model->MakeField<float>("dEdx");
    auto _omegaIP = _model->MakeField<float>("omegaIP");
    auto _omegaECAL = _model->MakeField<float>("omegaECAL");
    auto _tanLambdaIP = _model->MakeField<float>("tanLambdaIP");
    auto _tanLambdaECAL = _model->MakeField<float>("tanLambdaECAL");
    
    //momenta
    auto _mcPx = _model->MakeField<float>("mcPx");
    auto _mcPy = _model->MakeField<float>("mcPy");
    auto _mcPz = _model->MakeField<float>("mcPz");
    auto _recoIpPx = _model->MakeField<float>("recoIpPx");
    auto _recoIpPy = _model->MakeField<float>("recoIpPy");
    auto _recoIpPz = _model->MakeField<float>("recoIpPz");
    auto _recoCaloX = _model->MakeField<float>("recoCaloX");
    auto _recoCaloY = _model->MakeField<float>("recoCaloY");
    auto _recoCaloZ = _model->MakeField<float>("recoCaloZ");
    auto _recoCaloPx = _model->MakeField<float>("recoCaloPx");
    auto _recoCaloPy = _model->MakeField<float>("recoCaloPy");
    auto _recoCaloPz = _model->MakeField<float>("recoCaloPz");

    //track lengths
    auto _trackLength_IDR = _model->MakeField<float>("trackLength_IDR");
    auto _trackLength_SHA_phiLambda_IP = _model->MakeField<float>("trackLengthToEcal_SHA_phiLambda_IP");
    auto _trackLength_SHA_phiZed_IP = _model->MakeField<float>("trackLengthToEcal_SHA_phiZed_IP");
    auto _trackLength_SHA_zedLambda_IP = _model->MakeField<float>("trackLengthToEcal_SHA_phiZed_IP");
    auto _trackLength_SHA_phiLambda_ECAL = _model->MakeField<float>("trackLengthToEcal_SHA_phiLambda_ECAL");
    auto _trackLength_SHA_phiZed_ECAL = _model->MakeField<float>("trackLengthToEcal_SHA_phiZed_ECAL");
    auto _trackLength_SHA_zedLambda_ECAL = _model->MakeField<float>("trackLengthToEcal_SHA_zedLambda_ECAL");
    
    auto _trackLength_IKF_phiLambda_trackLengthToEcal = _model->MakeField<float>("trackLengthToEcal_IKF_phiLambda");
    auto _trackLength_IKF_phiLambda_trackLengthToSET = _model->MakeField<float>("trackLengthToSET_IKF_phiLambda");
    auto _trackLength_IKF_phiLambda_harmonicMomToEcal = _model->MakeField<float>("harmonicMomToEcal_IKF_phiLambda");
    auto _trackLength_IKF_phiLambda_harmonicMomToSET = _model->MakeField<float>("harmonicMomToSET_IKF_phiLambda");

    auto _trackLength_IKF_phiZed_trackLengthToEcal = _model->MakeField<float>("trackLengthToEcal_IKF_phiZed");
    auto _trackLength_IKF_phiZed_trackLengthToSET = _model->MakeField<float>("trackLengthToSET_IKF_phiZed");
    auto _trackLength_IKF_phiZed_harmonicMomToEcal = _model->MakeField<float>("harmonicMomToEcal_IKF_phiZed");
    auto _trackLength_IKF_phiZed_harmonicMomToSET = _model->MakeField<float>("harmonicMomToSET_IKF_phiZed");

    auto _trackLength_IKF_zedLambda_trackLengthToEcal = _model->MakeField<float>("trackLengthToEcal_IKF_zedLambda");
    auto _trackLength_IKF_zedLambda_trackLengthToSET = _model->MakeField<float>("trackLengthToSET_IKF_zedLambda");
    auto _trackLength_IKF_zedLambda_harmonicMomToEcal = _model->MakeField<float>("harmonicMomToEcal_IKF_zedLambda");
    auto _trackLength_IKF_zedLambda_harmonicMomToSET = _model->MakeField<float>("harmonicMomToSET_IKF_zedLambda");

    auto _cleanTrack = _model->MakeField<bool>("cleanTrack");

    //tofs
    auto _typeClosest = _model->MakeField<int>("typeClosest");
    auto _caloIDClosest = _model->MakeField<int>("caloIDClosest");
    auto _layoutClosest = _model->MakeField<int>("layoutClosest");
    auto _layerClosest = _model->MakeField<int>("layerClosest");
    auto _cleanClosestHit = _model->MakeField<bool>("cleanClosestHit");

    
    for (size_t i = 0; i < _resolutions.size(); i++){
        int res = int(_resolutions[i]);
        auto _tofClosest = _model->MakeField<float>( ( "tofClosest"+std::to_string(res) ).c_str() );


        _tree->Branch(( "tofClosest"+std::to_string(res) ).c_str(), &( _tofClosest[i]) );
        _tree->Branch(( "tofAverage"+std::to_string(res) ).c_str(), &( _tofAverage[i]) );
        _tree->Branch(( "tofSET"+std::to_string(res) ).c_str(), &( _tofSET[i]) );
        _tree->Branch(( "tofFit"+std::to_string(res) ).c_str(), &( _tofFit[i]) );
    }

    _tree->Branch("nHits", &_nHits);
    _tree->Branch("xHit", &_xHit);
    _tree->Branch("yHit", &_yHit);
    _tree->Branch("zHit", &_zHit);
    _tree->Branch("tHit", &_tHit);
    _tree->Branch("layerHit", &_layerHit);
    _tree->Branch("energyHit", &_energyHit);

}

void BohdanAna::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<"==========Event========== "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );
    LCRelationNavigator navToSimTrackerHits( evt->getCollection("TrackerHitsRelations") );
    LCRelationNavigator navToSimCalorimeterHits( evt->getCollection("CalorimeterHitsRelations") );

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
        _mcPx = mc->getMomentum()[0];
        _mcPy = mc->getMomentum()[1];
        _mcPz = mc->getMomentum()[2];

        bool isHadron = std::abs(_pdg) == 211 || std::abs(_pdg) == 321 || std::abs(_pdg) == 2212;
        bool isPhoton = std::abs(_pdg) == 22;

        for (const auto& hit:cluster->getCalorimeterHits()){
            //Count only ECAL hits
            CHT hitType( hit->getType() );
            bool isEcal = (hitType.caloID() == CHT::ecal);
            if (!isEcal) continue;
            _xHit.push_back(hit->getPosition()[0]);
            _yHit.push_back(hit->getPosition()[1]);
            _zHit.push_back(hit->getPosition()[2]);
            _tHit.push_back(hit->getTime());
            _layerHit.push_back( hitType.layer() );
            _energyHit.push_back( hit->getEnergy() );
        }
        _nHits = _tHit.size();

        if (isHadron && nTracks == 1){
            Track* track = pfo->getTracks().at(0);
            _dEdx = track->getdEdx();

            auto tsIP = track->getTrackState( TrackState::AtIP );
            _omegaIP = tsIP->getOmega();
            _tanLambdaIP = tsIP->getTanLambda();

            auto tsCalo = getTrackStateAtCalorimeter( track );
            _omegaECAL = tsCalo->getOmega();
            _tanLambdaECAL = tsCalo->getTanLambda();

            Vector3D trackPosAtCalo( tsCalo->getReferencePoint() );
            std::array<double, 3> mom = UTIL::getTrackMomentum(tsCalo, _bField);
            Vector3D trackMomAtCalo(mom[0], mom[1], mom[2]);

            _recoIpPx = pfo->getMomentum()[0];
            _recoIpPy = pfo->getMomentum()[1];
            _recoIpPz = pfo->getMomentum()[2];

            _recoCaloX = trackPosAtCalo[0];
            _recoCaloY = trackPosAtCalo[1];
            _recoCaloZ = trackPosAtCalo[2];

            _recoCaloPx = trackMomAtCalo[0];
            _recoCaloPy = trackMomAtCalo[1];
            _recoCaloPz = trackMomAtCalo[2];

            streamlog_out(DEBUG8)<<"getTrackStates()"<<std::endl;
            std::vector<HitState> trackHitStates = getTrackStates(pfo, _bField, _trkSystem, navToSimTrackerHits);
            std::vector<IMPL::TrackStateImpl> trackStates;
            for(auto& hitState: trackHitStates){
                trackStates.push_back(hitState.ts);

                if ( hitState.simHit != nullptr && hitState.simHit->getMCParticle() != mc ) _cleanTrack = false;
            }

            _trackLength_IDR = getTrackLengthIDR(track);
            streamlog_out(DEBUG8)<<"getTrackLengthSHA(AtIP)"<<std::endl;
            _trackLength_SHA_phiLambda_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::phiLambda);
            _trackLength_SHA_phiZed_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::phiZed);
            _trackLength_SHA_zedLambda_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::zedLambda);
            streamlog_out(DEBUG8)<<"getTrackLengthSHA(AtCalorimeter)"<<std::endl;
            _trackLength_SHA_phiLambda_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::phiLambda);
            _trackLength_SHA_phiZed_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::phiZed);
            _trackLength_SHA_zedLambda_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::zedLambda);
            streamlog_out(DEBUG8)<<"getTrackLengthIKF()"<<std::endl;

            _trackLength_IKF_phiLambda = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::phiLambda);
            _trackLength_IKF_phiLambda_trackLengthToEcal = _trackLength_IKF_phiLambda.trackLengthToEcal;
            _trackLength_IKF_phiLambda_trackLengthToSET = _trackLength_IKF_phiLambda.trackLengthToSET;
            _trackLength_IKF_phiLambda_harmonicMomToEcal = _trackLength_IKF_phiLambda.harmonicMomToEcal;
            _trackLength_IKF_phiLambda_harmonicMomToSET = _trackLength_IKF_phiLambda.harmonicMomToSET;

            _trackLength_IKF_phiZed = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::phiZed);
            _trackLength_IKF_phiZed_trackLengthToEcal = _trackLength_IKF_phiZed.trackLengthToEcal;
            _trackLength_IKF_phiZed_trackLengthToSET = _trackLength_IKF_phiZed.trackLengthToSET;
            _trackLength_IKF_phiZed_harmonicMomToEcal = _trackLength_IKF_phiZed.harmonicMomToEcal;
            _trackLength_IKF_phiZed_harmonicMomToSET = _trackLength_IKF_phiZed.harmonicMomToSET;

            _trackLength_IKF_zedLambda = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::zedLambda);
            _trackLength_IKF_zedLambda_trackLengthToEcal = _trackLength_IKF_zedLambda.trackLengthToEcal;
            _trackLength_IKF_zedLambda_trackLengthToSET = _trackLength_IKF_zedLambda.trackLengthToSET;
            _trackLength_IKF_zedLambda_harmonicMomToEcal = _trackLength_IKF_zedLambda.harmonicMomToEcal;
            _trackLength_IKF_zedLambda_harmonicMomToSET = _trackLength_IKF_zedLambda.harmonicMomToSET;


            streamlog_out(DEBUG8)<<"getTofClosest()"<<std::endl;
            CalorimeterHit* closestHit = getClosestHit(cluster, trackPosAtCalo);

            _typeClosest = getHitCaloType(closestHit);
            _caloIDClosest = getHitCaloID(closestHit);
            _layoutClosest = getHitCaloLayout(closestHit);
            _layerClosest = getHitCaloLayer(closestHit);
            _cleanClosestHit = getHitEarliestMC(closestHit, navToSimCalorimeterHits) == mc;

            streamlog_out(DEBUG8)<<"selectFrankEcalHits()"<<std::endl;
            auto selectedHits = selectFrankEcalHits(cluster, trackPosAtCalo, trackMomAtCalo, 10);
            streamlog_out(DEBUG8)<<"Calculating TOFs for all resolutions"<<std::endl;
            for (size_t j = 0; j < _resolutions.size(); j++){
                float res = _resolutions[j]/1000.; // in ns
                _tofClosest.at(j) = getHitTof(closestHit, trackPosAtCalo, res);
                _tofAverage.at(j) = getTofFrankAvg(selectedHits, trackPosAtCalo, res);
                _tofFit.at(j) = getTofFrankFit(selectedHits, trackPosAtCalo, res);
                _tofSET.at(j) = getTofSET(track, res);
            }

            // drawDisplay(this, evt, displayPFO, pfo, true);
            // plotCanvas(cluster, trackPosAtCalo, trackMomAtCalo, mc);
            // plotTrackParams(trackHitStates, pfo, mc, _bField);
        }
        else if( isPhoton && nTracks == 0 && ( !mc->isDecayedInTracker() ) ) {
            streamlog_out(DEBUG8)<<"Photon stuff"<<std::endl;

            Vector3D photonPosAtCalo = getPhotonAtCalorimeter(mc);
            Vector3D mom( mc->getMomentum() );

            CalorimeterHit* closestHit = getClosestHit(cluster, photonPosAtCalo);
            _typeClosest = getHitCaloType(closestHit);
            _caloIDClosest = getHitCaloID(closestHit);
            _layoutClosest = getHitCaloLayout(closestHit);
            _layerClosest = getHitCaloLayer(closestHit);
            _cleanClosestHit = getHitEarliestMC(closestHit, navToSimCalorimeterHits) == mc;


            auto selectedHits = selectFrankEcalHits(cluster, photonPosAtCalo, mom, 10);
            streamlog_out(DEBUG8)<<"Calculating TOFs for all resolutions"<<std::endl;
            for (size_t j = 0; j < _resolutions.size(); j++){
                float res = _resolutions[j]/1000.; // in ns
                _tofClosest.at(j) = getHitTof(closestHit, photonPosAtCalo, res);
                _tofAverage.at(j) = getTofFrankAvg(selectedHits, photonPosAtCalo, res);
                _tofFit.at(j) = getTofFrankFit(selectedHits, photonPosAtCalo, res);
                // no track - no SET!
            }
        }
        else{
            continue;
        }
        _tree->Fill();

        //DEBUGGING
        // float mom = _trackLength_IKF_zedLambda.harmonicMomToEcal;
        // float trackLength = _trackLength_IKF_zedLambda.trackLengthToEcal;
        // float beta = trackLength/(299.792458*_tofClosest[0]);
        // float m2 = mom*mom*(1./(beta*beta) - 1.);
        // if (m2 < -1. && _tofClosest[0] > 6. && mom < 2. && _layerClosest == 0 && _cleanClosestHit && _cleanTrack){
        //     streamlog_out(DEBUG8)<<" Momentum: "<<mom<<" GeV/c"<<std::endl;
        //     streamlog_out(DEBUG8)<<" Track length: "<<trackLength<<" mm"<<std::endl;
        //     streamlog_out(DEBUG8)<<" TOF: "<<_tofClosest[0]<<" ns"<<std::endl;
        //     streamlog_out(DEBUG8)<<" Beta: "<<beta<<std::endl;
        //     streamlog_out(DEBUG8)<<" m^2: "<<m2<<" GeV/c^2"<<std::endl;
        //     drawDisplay(this, evt, displayPFO, pfo, true);
        // }
    }
}

void BohdanAna::end(){
    _file->Write();
    // _application.Run(true);
}

void BohdanAna::resetVariables(){
    _pdg = 0;
    _dEdx = 0.f;
    _omegaIP = 0.f;
    _omegaECAL = 0.f;
    _tanLambdaIP = 0.f;
    _tanLambdaECAL = 0.f;
    _mcMom.fill(0.f);
    _recoIpMom.fill(0.f);
    _recoCaloPos.fill(0.f);
    _recoCaloMom.fill(0.f);

    _trackLength_IDR = 0.f;
    _trackLength_SHA_phiLambda_IP = 0.f;
    _trackLength_SHA_phiZed_IP = 0.f;
    _trackLength_SHA_zedLambda_IP = 0.f;
    _trackLength_SHA_phiLambda_ECAL = 0.f;
    _trackLength_SHA_phiZed_ECAL = 0.f;
    _trackLength_SHA_zedLambda_ECAL = 0.f;
    _trackLength_IKF_phiLambda = TrackLengthResult();
    _trackLength_IKF_phiZed = TrackLengthResult();
    _trackLength_IKF_zedLambda = TrackLengthResult();
    _cleanTrack = true;
    _typeClosest = -1;
    _caloIDClosest = -1;
    _layoutClosest = -1;
    _layerClosest = -1;
    _cleanClosestHit = false;
    _tofClosest.fill(0.f);
    _tofAverage.fill(0.f);
    _tofSET.fill(0.f);
    _tofFit.fill(0.f);

    _nHits = 0;
    _xHit.clear();
    _yHit.clear();
    _zHit.clear();
    _tHit.clear();
    _layerHit.clear();
    _energyHit.clear();
}