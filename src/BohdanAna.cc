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

void BohdanAna::fillMCTrueInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo){
    _pdg = mc->getPDG();
    _mcVtx = Vector3D( mc->getVertex() );
    _mcMom = Vector3D( mc->getMomentum() );
    _omegaTrue = getOmegaTrue(mc);
    _tanLambdaTrue = getTanLTrue(mc);
    _d0True = getD0True(mc);
    _z0True = getZ0True(mc);
    _phiTrue = getPhiTrue(mc);
    _timeTrue = mc->getTime();
    _isOverlay = mc->isOverlay();
    _isSimulated = mc->isCreatedInSimulation();
    if ( not mc->getParents().empty() ){
        _imidiateParentPDG = mc->getParents()[0]->getPDG();
        MCParticle* firstStableParent = getFirstStableParent(mc);
        if ( firstStableParent != nullptr ) _firstStableParentPDG = firstStableParent->getPDG();
    }
    unsigned int quarkTypeDecay = getQuarkTypeDecay(mc);
    _isBottomQuarkDecay = quarkTypeDecay == 5;
    _isCharmQuarkDecay = quarkTypeDecay == 4;
    _isHadronisationDecay = not _isBottomQuarkDecay && not _isCharmQuarkDecay;

    _isV0DecayTrue = std::abs(_firstStableParentPDG) == 310  || std::abs(_firstStableParentPDG) == 3122;

    _isReconstructed = pfo != nullptr;

    if ( _isReconstructed ){
        // Have a WELL-DEFINED track/cluster. I do not store reco infromation of weird PFOs with 2 tracks/showers...
        _hasTrack = pfo->getTracks().size() == 1;
        _hasShower = pfo->getClusters().size() == 1;
    }

    std::cout<<"Checking MC particle information"<<std::endl;
    std::cout<<"_pdg: "<<_pdg<<std::endl;
    std::cout<<"_vtxPos: ("<<_mcVtx[0]<<", "<<_mcVtx[1]<<", "<<_mcVtx[2]<<")"<<std::endl;
    std::cout<<"_vtxMom: ("<<_mcMom[0]<<", "<<_mcMom[1]<<", "<<_mcMom[2]<<")"<<std::endl;
    std::cout<<"_d0True: "<<_d0True<<std::endl;
    std::cout<<"_z0True: "<<_z0True<<std::endl;
    std::cout<<"_omegaTrue: "<<_omegaTrue<<std::endl;
    std::cout<<"_tanLambdaTrue: "<<_tanLambdaTrue<<std::endl;
    std::cout<<"_timeTrue: "<<_timeTrue<<std::endl;
    std::cout<<"_isOverlay: "<<_isOverlay<<std::endl;
    std::cout<<"_isSimulated: "<<_isSimulated<<std::endl;
    std::cout<<"_imidiateParentPDG: "<<_imidiateParentPDG<<std::endl;
    std::cout<<"_firstStableParentPDG: "<<_firstStableParentPDG<<std::endl;
    std::cout<<"_isBottomQuarkDecay: "<<_isBottomQuarkDecay<<std::endl;
    std::cout<<"_isCharmQuarkDecay: "<<_isCharmQuarkDecay<<std::endl;
    std::cout<<"_isV0DecayTrue: "<<_isV0DecayTrue<<std::endl;
    std::cout<<"_isHadronisationDecay: "<<_isHadronisationDecay<<std::endl;
    std::cout<<"_isReconstructed: "<<_isReconstructed<<std::endl;
    std::cout<<"_hasTrack: "<<_hasTrack<<std::endl;
    std::cout<<"_hasShower: "<<_hasShower<<std::endl;

}

void BohdanAna::fillTOFInfo(EVENT::Cluster* cluster, EVENT::Track* track, EVENT::CalorimeterHit* closestHit, const dd4hep::rec::Vector3D& trackPosAtCalo, const dd4hep::rec::Vector3D& trackMomAtCalo){
    auto selectedFrankHits = selectFrankEcalHits(cluster, trackPosAtCalo, trackMomAtCalo, 10);

    for (size_t j = 0; j < _resolutions.size(); j++){
        float res = _resolutions[j]/1000.; // in ns
        _tofClosest.at(j) = getHitTof(closestHit, trackPosAtCalo, res);
        _tofAverage.at(j) = getTofFrankAvg(selectedFrankHits, trackPosAtCalo, res);
        _tofFit.at(j) = getTofFrankFit(selectedFrankHits, trackPosAtCalo, res);
        std::tie(_tofSETFront.at(j), _tofSETBack.at(j)) = getTofSET(track, res);
    }

    for (auto* hit: cluster->getCalorimeterHits()){
        //Count only ECAL hits. No LumiCal, BeamCal, HCAL, Yoke hits are recorded!
        if (not isEcalHit(hit) ) continue;
        _xHit.push_back(hit->getPosition()[0]);
        _yHit.push_back(hit->getPosition()[1]);
        _zHit.push_back(hit->getPosition()[2]);
        _tHit.push_back(hit->getTime());
        _layerHit.push_back( CHT( hit->getType() ).layer() );
        _energyHit.push_back( hit->getEnergy() );
    }
    _nHits = _tHit.size();
}

void BohdanAna::fillTrackLengthInfo(EVENT::ReconstructedParticle* pfo, EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& navToSimTrackerHits){
    if (pfo == nullptr || pfo->getTracks().size() != 1) return;
    auto track = pfo->getTracks()[0];
    _trackLength_IDR = getTrackLengthIDR(track);
    _trackLength_SHA_phiLambda_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::phiLambda);
    _trackLength_SHA_phiZed_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::phiZed);
    _trackLength_SHA_zedLambda_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::zedLambda);
    _trackLength_SHA_phiLambda_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::phiLambda);
    _trackLength_SHA_phiZed_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::phiZed);
    _trackLength_SHA_zedLambda_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::zedLambda);

    std::vector<HitState> trackHitStates = getTrackStates(pfo, _bField, _trkSystem, navToSimTrackerHits);
    std::vector<IMPL::TrackStateImpl> trackStates;
    for(auto& hitState: trackHitStates){
        trackStates.push_back(hitState.ts);

        if ( hitState.simHit != nullptr && hitState.simHit->getMCParticle() != mc ) _cleanTrack = false;
    }

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

}

void BohdanAna::fillCsvForKonrad(EVENT::Cluster* cluster, int pdg, double trueTOF, double trackLength, const dd4hep::rec::Vector3D& trackPosAtCalo, const dd4hep::rec::Vector3D& trackMomAtCalo){
    _global_pfo_number++;
    for (auto* hit:cluster->getCalorimeterHits()){
        //Count only ECAL hits. No LumiCal, BeamCal, HCAL, Yoke hits are recorded!
        if ( not isEcalHit(hit) ) continue;

        auto hitCaloID = getHitCaloID(hit);
        auto hitLayout = getHitCaloLayout(hit);
        auto hitLayer = getHitCaloLayer(hit);

        std::stringstream ss;
        ss<<_global_pfo_number<<", ";
        ss<<pdg<<", ";
        ss<<std::scientific<<std::setprecision(5)<<trackLength<<", ";

        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.r()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.trans()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.x()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.y()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.z()<<", ";

        ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo.x()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo.y()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo.z()<<", ";
        ss<<std::scientific<<std::setprecision(6)<<trueTOF<<", ";
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

void BohdanAna::fillTrackStates(EVENT::ReconstructedParticle* pfo, EVENT::ReconstructedParticle* refittedPfo){
    if (pfo == nullptr || pfo->getTracks().size() != 1) return;
    auto track = pfo->getTracks()[0];
    _dEdx = track->getdEdx();

    auto tsIP = track->getTrackState( TrackState::AtIP );
    _omegaIP = tsIP->getOmega();
    _tanLambdaIP = tsIP->getTanLambda();
    _d0IP = tsIP->getD0();
    _z0IP = tsIP->getZ0();
    _phiIP = tsIP->getPhi();
    _omegaErrIP = std::sqrt( tsIP->getCovMatrix()[5] );
    _tanLambdaErrIP = std::sqrt( tsIP->getCovMatrix()[14] );
    _d0ErrIP = std::sqrt( tsIP->getCovMatrix()[0] );
    _z0ErrIP = std::sqrt( tsIP->getCovMatrix()[9] );
    _phiErrIP = std::sqrt( tsIP->getCovMatrix()[2] );
    // _recoIpPos is always (0, 0, 0) so I do not store it...
    auto momIP = UTIL::getTrackMomentum(tsIP, _bField);
    _recoIpMom = Vector3D( momIP[0], momIP[1], momIP[2] ); // should be identical to pfo->getMomentum()

    auto tsECAL = getTrackStateAtCalorimeter( track );
    _omegaECAL = tsECAL->getOmega();
    _tanLambdaECAL = tsECAL->getTanLambda();
    _d0ECAL = tsECAL->getD0(); // NOTE: must be 0 by definition at ECAL
    _z0ECAL = tsECAL->getZ0(); // NOTE: must be 0 by definition at ECAL
    
    _phiECAL = tsECAL->getPhi();
    _omegaErrECAL = std::sqrt( tsECAL->getCovMatrix()[5] );
    _tanLambdaErrECAL = std::sqrt( tsECAL->getCovMatrix()[14] );
    _d0ErrECAL = std::sqrt( tsECAL->getCovMatrix()[0] );
    _z0ErrECAL = std::sqrt( tsECAL->getCovMatrix()[9] );
    _phiErrECAL = std::sqrt( tsECAL->getCovMatrix()[2] );
    _recoCaloPos = Vector3D( tsECAL->getReferencePoint() );
    auto momECAL = UTIL::getTrackMomentum(tsECAL, _bField);
    _recoCaloMom = Vector3D( momECAL[0], momECAL[1], momECAL[2] );

    if (refittedPfo == nullptr || refittedPfo->getTracks().size() != 1 ||  refittedPfo->getTracks()[0]->getOmega() == 0.f || refittedPfo->getTracks()[0]->getTanLambda() == 0.f) return;
    auto refittedTrack = refittedPfo->getTracks()[0];

    auto refittedTsIP = refittedTrack->getTrackState( TrackState::AtIP );
    _refittedOmegaIP = refittedTsIP->getOmega();
    _refittedTanLambdaIP = refittedTsIP->getTanLambda();
    _refittedD0IP = refittedTsIP->getD0();
    _refittedZ0IP = refittedTsIP->getZ0();
    _refittedPhiIP = refittedTsIP->getPhi();
    _refittedOmegaErrIP = std::sqrt( refittedTsIP->getCovMatrix()[5] );
    _refittedTanLambdaErrIP = std::sqrt( refittedTsIP->getCovMatrix()[14] );
    _refittedD0ErrIP = std::sqrt( refittedTsIP->getCovMatrix()[0] );
    _refittedZ0ErrIP = std::sqrt( refittedTsIP->getCovMatrix()[9] );
    _refittedPhiErrIP = std::sqrt( refittedTsIP->getCovMatrix()[2] );
    // _recoIpPos is always (0, 0, 0)
    auto refittedMomIP = UTIL::getTrackMomentum(refittedTsIP, _bField);
    _refittedRecoIpMom = Vector3D( refittedMomIP[0], refittedMomIP[1], refittedMomIP[2] ); // should be identical to refittedPfo->getMomentum();

    auto refittedTsECAL = getTrackStateAtCalorimeter( refittedTrack );
    _refittedOmegaECAL = refittedTsECAL->getOmega();
    _refittedTanLambdaECAL = refittedTsECAL->getTanLambda();
    _refittedD0ECAL = refittedTsECAL->getD0(); // NOTE: must be 0 by definition at ECAL
    _refittedZ0ECAL = refittedTsECAL->getZ0(); // NOTE: must be 0 by definition at ECAL
    _refittedPhiECAL = refittedTsECAL->getPhi();
    _refittedOmegaErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[5] );
    _refittedTanLambdaErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[14] );
    _refittedD0ErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[0] );
    _refittedZ0ErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[9] );
    _refittedPhiErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[2] );
    _refittedRecoCaloPos = Vector3D( refittedTsECAL->getReferencePoint() );
    auto refittedMomECAL = UTIL::getTrackMomentum(refittedTsECAL, _bField);
    _refittedRecoCaloMom = Vector3D( refittedMomECAL[0], refittedMomECAL[1], refittedMomECAL[2] );
}

void BohdanAna::fillRecoVertexInfo(EVENT::LCEvent* evt, EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& pfo2mc){
    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCCollection* updatedPfos = evt->getCollection("updatedPandoraPFOs");
    LCCollection* primVtxCol = evt->getCollection("PrimaryVertex_default");
    LCCollection* secondaryVtxCol = evt->getCollection("BuildUpVertex_default");
    LCCollection* secondaryV0VtxCol = evt->getCollection("BuildUpVertex_V0_default");
    LCCollection* primVtxRefitCol = evt->getCollection("PrimaryVertex_refit");
    LCCollection* secondaryVtxRefitCol = evt->getCollection("BuildUpVertex_refit");
    LCCollection* secondaryV0VtxRefitCol = evt->getCollection("BuildUpVertex_V0_refit");

    for(int i=0; i<primVtxCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (primVtxCol->getElementAt(i));
        auto prongs = vertex->getAssociatedParticle()->getParticles();
        for (auto* prong : prongs){
            auto prongMC = getMC(prong, pfo2mc);
            if (mc == prongMC) _isInRecoPrimaryVertex = true;
            break;
        }
    }
    for(int i=0; i<secondaryVtxCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryVtxCol->getElementAt(i));
        auto prongs = vertex->getAssociatedParticle()->getParticles();
        for (auto* prong : prongs){
            auto prongMC = getMC(prong, pfo2mc);
            if (mc == prongMC) _isInRecoSecondaryVertex = true;
            break;
        }
    }
    for(int i=0; i<secondaryV0VtxCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryV0VtxCol->getElementAt(i));
        auto prongs = vertex->getAssociatedParticle()->getParticles();
        for (auto* prong : prongs){
            auto prongMC = getMC(prong, pfo2mc);
            if (mc == prongMC){
                _isInRecoSecondaryVertex = true;
                _isV0DecayReco = true;
            }
            break;
        }
    }

    for(int i=0; i<primVtxRefitCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (primVtxRefitCol->getElementAt(i));
        auto prongs = vertex->getAssociatedParticle()->getParticles();
        for (auto* refittedProng : prongs){
            auto* prongObject = getMatchingElement(updatedPfos, refittedProng, pfos);
            if(prongObject == nullptr) continue;
            auto* prong = static_cast<ReconstructedParticle*> ( prongObject );
            auto prongMC = getMC(prong, pfo2mc);
            if (mc == prongMC) _isInRecoPrimaryRefitVertex = true;
            break;
        }
    }
    for(int i=0; i<secondaryVtxRefitCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryVtxRefitCol->getElementAt(i));
        auto prongs = vertex->getAssociatedParticle()->getParticles();
        for (auto* refittedProng : prongs){
            auto* prongObject = getMatchingElement(updatedPfos, refittedProng, pfos);
            if(prongObject == nullptr) continue;
            auto* prong = static_cast<ReconstructedParticle*> ( prongObject );
            auto prongMC = getMC(prong, pfo2mc);
            if (mc == prongMC) _isInRecoSecondaryRefitVertex = true;
            break;
        }
    }
    for(int i=0; i<secondaryV0VtxRefitCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryV0VtxRefitCol->getElementAt(i));
        auto prongs = vertex->getAssociatedParticle()->getParticles();
        for (auto* refittedProng : prongs){
            auto* prongObject = getMatchingElement(updatedPfos, refittedProng, pfos);
            if(prongObject == nullptr) continue;
            auto* prong = static_cast<ReconstructedParticle*> ( prongObject );
            auto prongMC = getMC(prong, pfo2mc);
            if (mc == prongMC){
                _isInRecoSecondaryRefitVertex = true;
                _isV0DecayRefitReco = true;
            }
            break;
        }
    }


}

void BohdanAna::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<"==================== Event: "<<_nEvent<<std::endl;

    LCCollection* mcs = evt->getCollection("MCParticle");
    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCCollection* updatedPfos = evt->getCollection("updatedPandoraPFOs");
    LCRelationNavigator mc2pfo ( evt->getCollection("MCTruthRecoLink") );
    LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );
    LCRelationNavigator navToSimTrackerHits( evt->getCollection("TrackerHitsRelations") );
    LCRelationNavigator navToSimCalorimeterHits( evt->getCollection("CalorimeterHitsRelations") );

    std::vector<VertexData> trueVertices = getReconstructableTrueVertices(evt);

    //Loop over MC particles
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        // Storing only charged hadrons - pi/K/p
        MCParticle* mc = static_cast <MCParticle*> ( mcs->getElementAt(i) );
        if (mc == nullptr) continue;
        bool isHadron = std::abs( mc->getPDG() ) == 211 || std::abs( mc->getPDG() ) == 321 || std::abs( mc->getPDG() ) == 2212;
        if ( not isHadron ) continue;

        resetVariables();

        _quarksToPythia = getQuarksToPythia(evt);

        // Is MC in the true vertex?
        for(auto& vtx : trueVertices) {
            auto it = std::find(vtx.mcs.begin(), vtx.mcs.end(), mc);
            if (it != vtx.mcs.end()){
                _isInTrueV0DecayTrue = vtx.isV0;
                if (vtx.isPrimary) _isInTruePrimaryVertex = true;
                else _isInTrueSecondaryVertex = true;
                break;
            }
        }

        auto pfo = getRelatedReconstructedParticle(mc, mc2pfo, pfo2mc);
        fillMCTrueInfo( mc, pfo );

        // Work only with simple reconstructed particles (1 track and 1 shower). Ignore rest ( ~ O( 0.1%) )
        if ( pfo == nullptr || pfo->getTracks().size() != 1 || pfo->getClusters().size() != 1){
            _tree->Fill();
            continue;
        }
        auto refittedPFO = static_cast<ReconstructedParticle* > ( getMatchingElement(pfos, pfo, updatedPfos) );
        fillRecoVertexInfo(evt, mc, pfo2mc);

        auto track = pfo->getTracks()[0];
        auto cluster = pfo->getClusters()[0];
        auto tsCalo = getTrackStateAtCalorimeter( track );
        Vector3D trackPosAtCalo( tsCalo->getReferencePoint() );
        auto closestHit = getClosestHit(cluster, trackPosAtCalo);

        auto mom = UTIL::getTrackMomentum(tsCalo, _bField);
        Vector3D trackMomAtCalo( mom[0], mom[1], mom[2] );

        fillTrackLengthInfo(pfo, mc, navToSimTrackerHits);
        fillTrackStates(pfo, refittedPFO);


        _typeClosest = getHitCaloType(closestHit);
        _caloIDClosest = getHitCaloID(closestHit);
        _layoutClosest = getHitCaloLayout(closestHit);
        _layerClosest = getHitCaloLayer(closestHit);
        _cleanClosestHit = getHitEarliestMC(closestHit, navToSimCalorimeterHits) == mc;

        //Do not fill TOF info when closest hit is in LumiCal
        if ( not isEcalHit(closestHit) ){
            _tree->Fill();
            continue;
        }
        fillTOFInfo(cluster, track, closestHit, trackPosAtCalo, trackMomAtCalo);
        if(_produce_csv_output) fillCsvForKonrad( cluster, _pdg, _tofClosest[0], _trackLength_IKF_zedLambda, trackPosAtCalo, trackMomAtCalo );

        // drawDisplay(this, evt, displayPFO, pfo);
        _tree->Fill();
    }


}

void BohdanAna::end(){
    _file->Write();
    if (_produce_csv_output) _csv_output_file.close();
    // _application.Run(true);
}

void BohdanAna::resetVariables(){
    _quarksToPythia = 0;
    _pdg = -1;
    _mcVtx = Vector3D();
    _mcMom = Vector3D();
    _omegaTrue = 0.f;
    _tanLambdaTrue = 0.f;
    _d0True = 0.f;
    _z0True = 0.f;
    _phiTrue = 0.f;
    _timeTrue = 0.f;
    _isOverlay = false;
    _isSimulated = false;
    _imidiateParentPDG = -1;
    _firstStableParentPDG = -1;
    _isBottomQuarkDecay = false;
    _isCharmQuarkDecay = false;
    _isHadronisationDecay = false;
    _isReconstructed = false;
    _hasTrack = false;
    _hasShower = false;
    _isInTruePrimaryVertex = false;
    _isInTrueSecondaryVertex = false;
    _isInRecoPrimaryVertex = false;
    _isInRecoSecondaryVertex = false;
    _isV0DecayTrue = false;
    _isV0DecayReco = false;

    _dEdx = 0.f;
    _omegaIP = 0.f;
    _tanLambdaIP = 0.f;
    _d0IP = 0.f;
    _z0IP = 0.f;
    _phiIP = 0.f;
    _omegaErrIP = 0.f;
    _tanLambdaErrIP = 0.f;
    _d0ErrIP = 0.f;
    _z0ErrIP = 0.f;
    _phiErrIP = 0.f;
    _recoIpMom = Vector3D();
    _omegaECAL = 0.f;
    _tanLambdaECAL = 0.f;
    _d0ECAL = 0.f;
    _z0ECAL = 0.f;
    _phiECAL = 0.f;
    _omegaErrECAL = 0.f;
    _tanLambdaErrECAL = 0.f;
    _d0ErrECAL = 0.f;
    _z0ErrECAL = 0.f;
    _phiErrECAL = 0.f;
    _recoCaloPos = Vector3D();
    _recoCaloMom = Vector3D();
    _refittedOmegaIP = 0.f;
    _refittedTanLambdaIP = 0.f;
    _refittedD0IP = 0.f;
    _refittedZ0IP = 0.f;
    _refittedPhiIP = 0.f;
    _refittedOmegaErrIP = 0.f;
    _refittedTanLambdaErrIP = 0.f;
    _refittedD0ErrIP = 0.f;
    _refittedZ0ErrIP = 0.f;
    _refittedPhiErrIP = 0.f;
    _refittedRecoIpMom = Vector3D();
    _refittedOmegaECAL = 0.f;
    _refittedTanLambdaECAL = 0.f;
    _refittedD0ECAL = 0.f;
    _refittedZ0ECAL = 0.f;
    _refittedPhiECAL = 0.f;
    _refittedOmegaErrECAL = 0.f;
    _refittedTanLambdaErrECAL = 0.f;
    _refittedD0ErrECAL = 0.f;
    _refittedZ0ErrECAL = 0.f;
    _refittedPhiErrECAL = 0.f;
    _refittedRecoCaloPos = Vector3D();
    _refittedRecoCaloMom = Vector3D();

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

    _xHit.clear();
    _yHit.clear();
    _zHit.clear();
    _tHit.clear();
    _layerHit.clear();
    _energyHit.clear();
    _nHits = 0;
}