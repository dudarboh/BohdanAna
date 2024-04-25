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
#include "CLHEP/Random/Randomize.h"

using dd4hep::rec::Vector3D;

BohdanAna theBohdanAna;

BohdanAna::BohdanAna() : marlin::Processor("BohdanAna"), EventDisplayer(this){
    _description = "Main analysis of the track length and time-of-flight and momentum methods for time-of-flight pID";

    registerProcessorParameter( "produce_csv_output",
                                "Produce csv output file for Konrad or not",
                                _produce_csv_output,
                                false);
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

    //track parameters
    _tree->Branch("dEdx", &_dEdx);
    _tree->Branch("omegaIP", &_omegaIP);
    _tree->Branch("omegaECAL", &_omegaECAL);
    _tree->Branch("tanLambdaIP", &_tanLambdaIP);
    _tree->Branch("tanLambdaECAL", &_tanLambdaECAL);

    _tree->Branch("refittedOmegaIP", &_refittedOmegaIP);
    _tree->Branch("refittedOmegaECAL", &_refittedOmegaECAL);
    _tree->Branch("refittedTanLambdaIP", &_refittedTanLambdaIP);
    _tree->Branch("refittedTanLambdaECAL", &_refittedTanLambdaECAL);

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
    _tree->Branch("recoCaloX", &(_recoCaloPos[0]) );
    _tree->Branch("recoCaloY", &(_recoCaloPos[1]) );
    _tree->Branch("recoCaloZ", &(_recoCaloPos[2]) );

    _tree->Branch("refittedRecoIpPx", &(_refittedRecoIpMom[0]) );
    _tree->Branch("refittedRecoIpPy", &(_refittedRecoIpMom[1]) );
    _tree->Branch("refittedRecoIpPz", &(_refittedRecoIpMom[2]) );
    _tree->Branch("refittedRecoCaloPx", &(_refittedRecoCaloMom[0]) );
    _tree->Branch("refittedRecoCaloPy", &(_refittedRecoCaloMom[1]) );
    _tree->Branch("refittedRecoCaloPz", &(_refittedRecoCaloMom[2]) );
    _tree->Branch("refittedRecoCaloX", &(_refittedRecoCaloPos[0]) );
    _tree->Branch("refittedRecoCaloY", &(_refittedRecoCaloPos[1]) );
    _tree->Branch("refittedRecoCaloZ", &(_refittedRecoCaloPos[2]) );

    //track lengths
    _tree->Branch("trackLength_IDR", &_trackLength_IDR);
    _tree->Branch("trackLengthToEcal_SHA_phiLambda_IP", &_trackLength_SHA_phiLambda_IP);
    _tree->Branch("trackLengthToEcal_SHA_phiZed_IP", &_trackLength_SHA_phiZed_IP);
    _tree->Branch("trackLengthToEcal_SHA_zedLambda_IP", &_trackLength_SHA_zedLambda_IP);
    _tree->Branch("trackLengthToEcal_SHA_phiLambda_ECAL", &_trackLength_SHA_phiLambda_ECAL);
    _tree->Branch("trackLengthToEcal_SHA_phiZed_ECAL", &_trackLength_SHA_phiZed_ECAL);
    _tree->Branch("trackLengthToEcal_SHA_zedLambda_ECAL", &_trackLength_SHA_zedLambda_ECAL);

    _tree->Branch("trackLengthToEcal_IKF_phiLambda", &_trackLength_IKF_phiLambda);
    _tree->Branch("trackLengthToSET_IKF_phiLambda", &_trackLengthToSET_IKF_phiLambda);
    _tree->Branch("harmonicMomToEcal_IKF_phiLambda", &_harmonicMom_IKF_phiLambda);
    _tree->Branch("harmonicMomToSET_IKF_phiLambda", &_harmonicMomToSET_IKF_phiLambda);

    _tree->Branch("trackLengthToEcal_IKF_phiZed", &_trackLength_IKF_phiZed);
    _tree->Branch("trackLengthToSET_IKF_phiZed", &_trackLengthToSET_IKF_phiZed);
    _tree->Branch("harmonicMomToEcal_IKF_phiZed", &_harmonicMom_IKF_phiZed);
    _tree->Branch("harmonicMomToSET_IKF_phiZed", &_harmonicMomToSET_IKF_phiZed);

    _tree->Branch("trackLengthToEcal_IKF_zedLambda", &_trackLength_IKF_zedLambda);
    _tree->Branch("trackLengthToSET_IKF_zedLambda", &_trackLengthToSET_IKF_zedLambda);
    _tree->Branch("harmonicMomToEcal_IKF_zedLambda", &_harmonicMom_IKF_zedLambda);
    _tree->Branch("harmonicMomToSET_IKF_zedLambda", &_harmonicMomToSET_IKF_zedLambda);
    _tree->Branch("cleanTrack", &_cleanTrack);

    //tofs
    _tree->Branch("typeClosest", &_typeClosest);
    _tree->Branch("caloIDClosest", &_caloIDClosest);
    _tree->Branch("layoutClosest", &_layoutClosest);
    _tree->Branch("layerClosest", &_layerClosest);
    _tree->Branch("cleanClosestHit", &_cleanClosestHit);
    for (size_t i = 0; i < _resolutions.size(); i++){
        int res = int(_resolutions[i]);
        _tree->Branch(( "tofClosest"+std::to_string(res) ).c_str(), &( _tofClosest[i]) );
        _tree->Branch(( "tofAverage"+std::to_string(res) ).c_str(), &( _tofAverage[i]) );
        _tree->Branch(( "tofSETFront"+std::to_string(res) ).c_str(), &( _tofSETFront[i]) );
        _tree->Branch(( "tofSETBack"+std::to_string(res) ).c_str(), &( _tofSETBack[i]) );
        _tree->Branch(( "tofFit"+std::to_string(res) ).c_str(), &( _tofFit[i]) );
    }

    _tree->Branch("nHits", &_nHits);
    _tree->Branch("xHit", &_xHit);
    _tree->Branch("yHit", &_yHit);
    _tree->Branch("zHit", &_zHit);
    _tree->Branch("tHit", &_tHit);
    _tree->Branch("layerHit", &_layerHit);
    _tree->Branch("energyHit", &_energyHit);

    if(_produce_csv_output){
        _csv_output_file = std::ofstream("output.csv");
        _csv_output_file<<"PFO #,"
                    "PDG,"
                    "trk length (mm),"
                    "trk p (GeV),"
                    "trk pT (GeV),"
                    "trk px (GeV),"
                    "trk py (GeV),"
                    "trk pz (GeV),"
                    "trk Ecal x (mm),"
                    "trk Ecal y (mm),"
                    "trk Ecal z (mm),"
                    "true TOF (ns),"
                    "Hit type,"
                    "Hit Layout,"
                    "true hit time (ns),"
                    "hit time 50ps (ns),"
                    "hit time 100ps (ns),"
                    "hit energy (GeV),"
                    "hit layer,"
                    "x hit pos (mm),"
                    "y hit pos (mm),"
                    "z hit pos (mm)\n";
    }
}

void BohdanAna::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<"==================== Event: "<<_nEvent<<std::endl;
    // int vm = getVirtualMemoryUsage();
    // int rm = getPhysicalMemoryUsage();
    // streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCCollection* updatedPfos = evt->getCollection("updatedPandoraPFOs");
    LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );
    LCRelationNavigator navToSimTrackerHits( evt->getCollection("TrackerHitsRelations") );
    LCRelationNavigator navToSimCalorimeterHits( evt->getCollection("CalorimeterHitsRelations") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG8)<<"======== PFO: "<<i+1<<std::endl;
        resetVariables();
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nTracks = pfo->getTracks().size();
        int nClusters = pfo->getClusters().size();
        // only simple cases
        if( nTracks > 1 || nClusters != 1) continue;
        MCParticle* mc = getMC(pfo, pfo2mc);
        if (mc == nullptr) continue;
        _pdg = mc->getPDG();
        for(int j=0; j<3; j++) _mcMom.at(j) = mc->getMomentum()[j];

        bool isHadron = std::abs(_pdg) == 211 || std::abs(_pdg) == 321 || std::abs(_pdg) == 2212;
        bool isPhoton = std::abs(_pdg) == 22;
        streamlog_out(DEBUG8)<<"PDG: "<<_pdg<<"   isHadron: "<<isHadron<<"   isPhoton: "<<isPhoton<<std::endl;

        Cluster* cluster = pfo->getClusters().at(0);
        streamlog_out(DEBUG8)<<"Cluster with N hits: "<<cluster->getCalorimeterHits().size()<<std::endl;

        for (const auto& hit:cluster->getCalorimeterHits()){
            //Count only ECAL hits. No LumiCal, BeamCal, HCAL, Yoke hits are recorded!
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

            for(int j=0; j<3; j++) _recoIpMom.at(j) = pfo->getMomentum()[j];
            for(int j=0; j<3; j++) _recoCaloPos.at(j) = trackPosAtCalo[j];
            for(int j=0; j<3; j++) _recoCaloMom.at(j) = trackMomAtCalo[j];

            // Now fill the same for the refitted track.
            ReconstructedParticle* updatedPfo = static_cast <ReconstructedParticle*> ( updatedPfos->getElementAt(i) );
            Track* refittedTrack = updatedPfo->getTracks().at(0);
            if (refittedTrack->getOmega() == 0.f && refittedTrack->getTanLambda() == 0.f){
                //This is empty track! Fit has failed!?
                _refittedOmegaIP = 0.f;
                _refittedTanLambdaIP = 0.f;
                _refittedOmegaECAL = 0.f;
                _refittedTanLambdaECAL = 0.f;
                for(int j=0; j<3; j++) _refittedRecoIpMom.at(j) = 0.f;
                for(int j=0; j<3; j++) _refittedRecoCaloPos.at(j) = 0.f;
                for(int j=0; j<3; j++) _refittedRecoCaloMom.at(j) = 0.f;
            }
            else{
                auto refittedTsIP = refittedTrack->getTrackState( TrackState::AtIP );
                _refittedOmegaIP = refittedTsIP->getOmega();
                _refittedTanLambdaIP = refittedTsIP->getTanLambda();

                auto refittedTsCalo = getTrackStateAtCalorimeter( refittedTrack );
                _refittedOmegaECAL = refittedTsCalo->getOmega();
                _refittedTanLambdaECAL = refittedTsCalo->getTanLambda();

                Vector3D refittedTrackPosAtCalo( refittedTsCalo->getReferencePoint() );
                std::array<double, 3> refittedMom = UTIL::getTrackMomentum(refittedTsCalo, _bField);
                Vector3D refittedTrackMomAtCalo(refittedMom[0], refittedMom[1], refittedMom[2]);

                for(int j=0; j<3; j++) _refittedRecoIpMom.at(j) = updatedPfo->getMomentum()[j];
                for(int j=0; j<3; j++) _refittedRecoCaloPos.at(j) = refittedTrackPosAtCalo[j];
                for(int j=0; j<3; j++) _refittedRecoCaloMom.at(j) = refittedTrackMomAtCalo[j];
            }


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
            std::tie(_trackLength_IKF_phiLambda, _harmonicMom_IKF_phiLambda) = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::phiLambda);
            std::tie(_trackLength_IKF_phiZed, _harmonicMom_IKF_phiZed) = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::phiZed);
            std::tie(_trackLength_IKF_zedLambda, _harmonicMom_IKF_zedLambda) = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::zedLambda);

            auto it = std::find_if(trackHitStates.begin(), trackHitStates.end(), [](const HitState& hitState){return isSETHit(hitState.hit);});
            bool foundSETHit = it != trackHitStates.end();
            if (foundSETHit){
                int idx = it - trackHitStates.begin();
                std::vector<IMPL::TrackStateImpl> trackStatesToSET = std::vector<IMPL::TrackStateImpl>(trackStates.begin(), trackStates.begin() + idx + 1 );
                std::tie(_trackLengthToSET_IKF_phiLambda, _harmonicMomToSET_IKF_phiLambda) = getTrackLengthIKF(trackStatesToSET, _bField, TrackLengthOption::phiLambda);
                std::tie(_trackLengthToSET_IKF_phiZed, _harmonicMomToSET_IKF_phiZed) = getTrackLengthIKF(trackStatesToSET, _bField, TrackLengthOption::phiZed);
                std::tie(_trackLengthToSET_IKF_zedLambda, _harmonicMomToSET_IKF_zedLambda) = getTrackLengthIKF(trackStatesToSET, _bField, TrackLengthOption::zedLambda);
            }

            streamlog_out(DEBUG8)<<"getTofClosest()"<<std::endl;
            CalorimeterHit* closestHit = getClosestHit(cluster, trackPosAtCalo);
            //NOTE: We assume no time measurementin the LumiCal! This cut should be consistent with one in the loops over ECAL hits!
            // We ignore the particle if the closest hit is not in the ECAL, e.g. LumiCal.
            bool notEcalClosestHit = ! (CHT( closestHit->getType() ).caloID() == CHT::ecal);
            if ( notEcalClosestHit ) continue;

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
                std::tie(_tofSETFront.at(j), _tofSETBack.at(j)) = getTofSET(track, res);
            }

            if (_produce_csv_output){
                _global_pfo_number++;

                for (const auto& hit:cluster->getCalorimeterHits()){
                    //Count only ECAL hits. No LumiCal, BeamCal, HCAL, Yoke hits are recorded!
                    bool notEcalHit = ! (CHT( hit->getType() ).caloID() == CHT::ecal);
                    if ( notEcalHit ) continue;

                    auto hitCaloID = getHitCaloID(hit);
                    auto hitLayout = getHitCaloLayout(hit);
                    auto hitLayer = getHitCaloLayer(hit);

                    std::stringstream ss;
                    ss<<_global_pfo_number<<", ";
                    ss<<_pdg<<", ";
                    ss<<std::scientific<<std::setprecision(5)<<_trackLength_IKF_zedLambda<<", ";

                    ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.r()<<", ";
                    ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.trans()<<", ";
                    ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo[0]<<", ";
                    ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo[1]<<", ";
                    ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo[2]<<", ";

                    ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo[0]<<", ";
                    ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo[1]<<", ";
                    ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo[2]<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<_tofClosest[0]<<", ";
                    ss<<hitCaloID<<", "; // type of hit, i.e. ECal, LumiCal etc. hit
                    ss<<hitLayout<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<hit->getTime()<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<CLHEP::RandGauss::shoot(hit->getTime(), 0.05)<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<CLHEP::RandGauss::shoot(hit->getTime(), 0.1)<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<hit->getEnergy()<<", ";
                    ss<<hitLayer<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<hit->getPosition()[0]<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<hit->getPosition()[1]<<", ";
                    ss<<std::scientific<<std::setprecision(6)<<hit->getPosition()[2]<<"\n";

                    _csv_output_file << ss.str();
                }
            }

        }
        else if( isPhoton && nTracks == 0 && ( !mc->isDecayedInTracker() ) ) {
            streamlog_out(DEBUG8)<<"Photon stuff"<<std::endl;

            Vector3D photonPosAtCalo = getPhotonAtCalorimeter(mc);
            Vector3D mom( mc->getMomentum() );

            CalorimeterHit* closestHit = getClosestHit(cluster, photonPosAtCalo);
            CHT hitType( closestHit->getType() );
            bool isEcal = (hitType.caloID() == CHT::ecal);
            if (!isEcal) continue;

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

    }
}

void BohdanAna::end(){
    _file->Write();
    if (_produce_csv_output) _csv_output_file.close();
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
    _refittedOmegaIP = 0.f;
    _refittedOmegaECAL = 0.f;
    _refittedTanLambdaIP = 0.f;
    _refittedTanLambdaECAL = 0.f;
    _refittedRecoIpMom.fill(0.f);
    _refittedRecoCaloPos.fill(0.f);
    _refittedRecoCaloMom.fill(0.f);

    _trackLength_IDR = 0.f;
    _trackLength_SHA_phiLambda_IP = 0.f;
    _trackLength_SHA_phiZed_IP = 0.f;
    _trackLength_SHA_zedLambda_IP = 0.f;
    _trackLength_SHA_phiLambda_ECAL = 0.f;
    _trackLength_SHA_phiZed_ECAL = 0.f;
    _trackLength_SHA_zedLambda_ECAL = 0.f;
    _trackLength_IKF_phiLambda = 0.f;
    _trackLength_IKF_phiZed = 0.f;
    _trackLength_IKF_zedLambda = 0.f;
    _harmonicMom_IKF_phiLambda = 0.f;
    _harmonicMom_IKF_phiZed = 0.f;
    _harmonicMom_IKF_zedLambda = 0.f;
    _trackLengthToSET_IKF_phiLambda = 0.f;
    _trackLengthToSET_IKF_phiZed = 0.f;
    _trackLengthToSET_IKF_zedLambda = 0.f;
    _harmonicMomToSET_IKF_phiLambda = 0.f;
    _harmonicMomToSET_IKF_phiZed = 0.f;
    _harmonicMomToSET_IKF_zedLambda = 0.f;


    _cleanTrack = true;
    _typeClosest = -1;
    _caloIDClosest = -1;
    _layoutClosest = -1;
    _layerClosest = -1;
    _cleanClosestHit = false;
    _tofClosest.fill(0.f);
    _tofAverage.fill(0.f);
    _tofSETFront.fill(0.f);
    _tofSETBack.fill(0.f);
    _tofFit.fill(0.f);

    _nHits = 0;
    _xHit.clear();
    _yHit.clear();
    _zHit.clear();
    _tHit.clear();
    _layerHit.clear();
    _energyHit.clear();
}