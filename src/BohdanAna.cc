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

void BohdanAna::fillMCTrueInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo, bool isInReconstructablePrimaryVertex, bool isInReconstructableSecondaryVertex){
    _pdg = mc->getPDG();
    for(int j=0; j<3; j++) _mcVtx.at(j) = mc->getVertex()[j];
    for(int j=0; j<3; j++) _mcMom.at(j) = mc->getMomentum()[j];
    _d0True = getD0True(mc);
    _z0True = getZ0True(mc);
    _omegaTrue = getOmegaTrue(mc);
    _tanLambdaTrue = getTanLTrue(mc);
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

    _isInReconstructablePrimaryVertex = isInReconstructablePrimaryVertex;
    _isInReconstructableSecondaryVertex = isInReconstructableSecondaryVertex;
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
    std::cout<<"_isHadronisationDecay: "<<_isHadronisationDecay<<std::endl;
    std::cout<<"_isInReconstructablePrimaryVertex: "<<_isInReconstructablePrimaryVertex<<std::endl;
    std::cout<<"_isInReconstructableSecondaryVertex: "<<_isInReconstructableSecondaryVertex<<std::endl;
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

    _d0ErrBefore = std::sqrt( track->getCovMatrix()[0] );
    _d0ErrAfter = std::sqrt( trackAfter->getCovMatrix()[0] );
    _z0ErrBefore = std::sqrt( track->getCovMatrix()[9] );
    _z0ErrAfter = std::sqrt( trackAfter->getCovMatrix()[9] );

    auto tsECAL = getTrackStateAtCalorimeter( track );
    _omegaECAL = tsECAL->getOmega();
    _tanLambdaECAL = tsECAL->getTanLambda();
    _d0ECAL = tsECAL->getD0();
    _z0ECAL = tsECAL->getZ0();
    _phiECAL = tsECAL->getPhi();


    auto refittedTrack = refittedPfo->getTracks()[0];

    std::array<double, 3> mom = UTIL::getTrackMomentum(tsCalo, _bField);
    Vector3D trackMomAtCalo(mom[0], mom[1], mom[2]);

    for(int j=0; j<3; j++) _recoIpMom.at(j) = pfo->getMomentum()[j];
    for(int j=0; j<3; j++) _recoCaloPos.at(j) = trackPosAtCalo[j];
    for(int j=0; j<3; j++) _recoCaloMom.at(j) = trackMomAtCalo[j];

    // Now fill the same for the refitted track.
    ReconstructedParticle* updatedPfo = static_cast <ReconstructedParticle*> ( updatedPfos->getElementAt(i) );
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

    Vector3D truePrimaryVertex( static_cast<MCParticle*>(mcs->getElementAt(0))->getVertex() );
    std::map< Vector3D, std::vector<MCParticle*>, CompareVectors > reconstructableTrueVertices = getReconstructableTrueVertices(evt);

    // No primary vertex!? Throw away this event
    if ( reconstructableTrueVertices.size() == 0 ) return;
    // Split primary and secondary
    auto closestToTruePrimaryVertex = [&truePrimaryVertex](const std::pair<Vector3D, std::vector<MCParticle*> >& a, const std::pair<Vector3D, std::vector<MCParticle*> >& b){return (a.first - truePrimaryVertex).r() < (b.first - truePrimaryVertex).r();};
    auto primaryVertexIter = std::min_element(reconstructableTrueVertices.begin(), reconstructableTrueVertices.end(), closestToTruePrimaryVertex);
    Vector3D reconstructablePrimaryVertexPosition = primaryVertexIter->first;
    std::vector<MCParticle*> reconstructablePrimaryVertexMCs = primaryVertexIter->second;
    reconstructableTrueVertices.erase( primaryVertexIter );
    auto reconstructableSecondaryVertices = reconstructableTrueVertices;

    //LOOP OVER ALL MCParticles and store only pointers to pi/k/p/gamma
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        resetVariables();

        MCParticle* mc = static_cast <MCParticle*> ( mcs->getElementAt(i) );
        if (mc == nullptr) continue;
        bool isHadron = std::abs( mc->getPDG() ) == 211 || std::abs( mc->getPDG() ) == 321 || std::abs( mc->getPDG() ) == 2212;
        if ( not isHadron ) continue;
        ReconstructedParticle* pfo = getRelatedReconstructedParticle(mc, mc2pfo, pfo2mc);
        auto refittedPFO = getRefittedPFO(pfos, updatedPfos, pfo);

        bool isInRecoconstructablePrimaryVertex = std::find(reconstructablePrimaryVertexMCs.begin(), reconstructablePrimaryVertexMCs.end(), mc) != reconstructablePrimaryVertexMCs.end();
        bool isInRecoconstructableSecondaryVertex = false;
        for(auto const& [pos, reconstructableSecondaryVertexMCs] : reconstructableSecondaryVertices){
            if (std::find(reconstructableSecondaryVertexMCs.begin(), reconstructableSecondaryVertexMCs.end(), mc) != reconstructableSecondaryVertexMCs.end() ){
                isInRecoconstructableSecondaryVertex = true;
                break;
            }
        }

        _quarksToPythia = getQuarksToPythia(evt);
        fillMCTrueInfo( mc, pfo, isInRecoconstructablePrimaryVertex, isInRecoconstructableSecondaryVertex );
        std::cout<<"_quarksToPythia: "<<_quarksToPythia<<std::endl;

        if ( pfo == nullptr || pfo->getTracks().size() != 1 || pfo->getClusters().size() != 1){
            _tree->Fill();
            continue;
        }

        auto track = pfo->getTracks()[0];
        auto cluster = pfo->getClusters()[0];
        auto tsCalo = getTrackStateAtCalorimeter( track );
        Vector3D trackPosAtCalo( tsCalo->getReferencePoint() );
        auto closestHit = getClosestHit(cluster, trackPosAtCalo);

        std::array<double, 3> mom = UTIL::getTrackMomentum(tsCalo, _bField);
        Vector3D trackMomAtCalo(mom[0], mom[1], mom[2]);

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
    _mcVtx.fill(0.f);
    _mcMom.fill(0.f);
    _d0True = 0.f;
    _z0True = 0.f;
    _omegaTrue = 0.f;
    _tanLambdaTrue = 0.f;
    _timeTrue = 0.f;
    _isOverlay = false;
    _isSimulated = false;
    _imidiateParentPDG = -1;
    _firstStableParentPDG = -1;
    _isBottomQuarkDecay = false;
    _isCharmQuarkDecay = false;
    _isHadronisationDecay = false;
    _isInReconstructablePrimaryVertex = false;
    _isInReconstructableSecondaryVertex = false;
    _isReconstructed = false;
    _hasTrack = false;
    _hasShower = false;


    _dEdx = 0.f;
    _omegaIP = 0.f;
    _omegaECAL = 0.f;
    _tanLambdaIP = 0.f;
    _tanLambdaECAL = 0.f;
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