#include "TrackLengthDebug.h"
#include "TrackLengthDebugUtils.h"

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
#include "IMPL/TrackImpl.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "marlinutil/CalorimeterHitType.h"
#include <limits>
#include "marlinutil/DDMarlinCED.h"
using namespace TrackLengthDebugUtils;
using namespace EVENT;
using namespace UTIL;
using namespace IMPL;
using dd4hep::rec::Vector3D;
using std::pair;
using std::vector;
using std::string;
using MarlinTrk::IMarlinTrack;
using MarlinTrk::IMarlinTrkSystem;
using ROOT::Math::XYZVector;
TrackLengthDebug aTrackLengthDebug;


//TODOS:
// 1) Add event display for the whole event, only highlighting analysed particle
// 2) Check if the closest hit corresponds to the analysed MC particle
// 3) get distance between closest hit and tsCalo
// 4) Check if the track is pure or there are mixins
// 5) get pure TOF w/o distance correction
// 6) Figure out why the fuck track doesnt fit !?!!??!!?

double getHelixLength(double z_start, double z_end, Vector3D momentum, double bField){
    double c_factor = 0.299792458;
    double pt = momentum.rho();
    double pz = momentum.z();
    return (z_end - z_start)/pz * std::sqrt(  std::pow( pt/(c_factor*bField), 2) + pz*pz  );
}


void drawPFO(ReconstructedParticle* pfo, MCParticle* mc){
    Vector3D mcPos(mc->getVertex());
    Vector3D mcMom(mc->getMomentum());
    Vector3D recoPos(0, 0, 0);
    Vector3D recoMom(pfo->getMomentum());

    DDMarlinCED::drawHelix( 3.5, mc->getCharge(), recoPos.x(), recoPos.y(), recoPos.z(), recoMom.x(), recoMom.y(), recoMom.z(), 0, 3, 0x0062ff, 0, 2000, 2350, 0);
    DDMarlinCED::drawHelix( 3.5, mc->getCharge(), mcPos.x(), mcPos.y(), mcPos.z(), mcMom.x(), mcMom.y(), mcMom.z(), 0, 3, 0xff0026, 0, 2000, 2350, 0);

    vector<Track*> tracks;
    if ( not pfo->getTracks().empty() )
        tracks = getSubTracks(pfo->getTracks()[0]);
    int nHits = 0;
    for(auto* track: tracks){
        auto hits = track->getTrackerHits();
        for (auto* hit : hits){
            ++nHits;
            auto pos = hit->getPosition();
            int type = 0; // point
            int layer = 1; // doesn't matter
            int size = 4; // larger point

            unsigned long color = 0x3232a8;
            if ( (nHits) % 20 == 0){
                color = 0xfa0730;
                size = 8;
            }
            ced_hit_ID(pos[0], pos[1], pos[2], type, layer, size, color, 0 ); // tracker hits
        }
    }
    auto tsCalo = pfo->getTracks()[0]->getTrackState(TrackState::AtCalorimeter);
    auto posCalo = tsCalo->getReferencePoint();
    ced_hit_ID(posCalo[0], posCalo[1], posCalo[2] + tsCalo->getZ0(), 0, 1, 8, 0x00ffff, 0 ); // track state at calorimeter


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


TrackLengthDebug::TrackLengthDebug() : marlin::Processor("TrackLengthDebug"), EventDisplayer(this){
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

    _tree->Branch("isGoodTrack", &_isGoodTrack);
    _tree->Branch("nSegments", &_nSegments);
    _tree->Branch("nHitsFirstCurl", &_nHitsFirstCurl);
    _tree->Branch("nBadPhi", &_nBadPhi);
    _tree->Branch("nBadZ", &_nBadZ);
    _tree->Branch("momentumOld", &_momentumOld);
    _tree->Branch("momentumNew", &_momentumNew);
    _tree->Branch("trackLengthOld", &_trackLengthOld);
    _tree->Branch("trackLengthNew", &_trackLengthNew);

    _tree->Branch("tsFirstCaloPos", &_tsFirstCaloPos);
    _tree->Branch("tsLastCaloPos", &_tsLastCaloPos);

    _tree->Branch("isGoodClosestHitOld", &_isGoodClosestHitOld);
    _tree->Branch("isGoodClosestHitNew", &_isGoodClosestHitNew);
    _tree->Branch("hitClosestOldPos", &_hitClosestOldPos);
    _tree->Branch("hitClosestNewPos", &_hitClosestNewPos);
    _tree->Branch("hitTimeOld", &_hitTimeOld);
    _tree->Branch("hitTimeNew", &_hitTimeNew);
    _tree->Branch("tofOld", &_tofOld);
    _tree->Branch("tofNew", &_tofNew);
}

void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    PIDHandler pidHandler( pfos );
    LCRelationNavigator nav ( evt->getCollection("RecoMCTruthLink") );
    LCRelationNavigator navClusterHits ( evt->getCollection("ClusterHitsRelations") );


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
        _nHitsFirstCurl = track->getTrackerHits().size();
        Cluster* cluster = pfo->getClusters()[0];
        // store this to the TTree
        Vars resultsOld = getTrackLengthOld(pfo, _bField, _trkSystem);
        Vars resultsNew = getTrackLengthNew(pfo, _bField, _trkSystem);
 
        _pdg = mc->getPDG();

        _momentumOld = resultsOld.momentumHM;
        _momentumNew = resultsNew.momentumHM;

        _trackLengthOld = resultsOld.trackLength;
        _trackLengthNew = resultsNew.trackLength;

        auto tsCalo = track->getTrackState(TrackState::AtCalorimeter);
        _tsFirstCaloPos = Vector3D( tsCalo->getReferencePoint() );
        // get OLD TOF part
        CalorimeterHit* closestHitOld = nullptr;
        double hitTimeOld = std::numeric_limits<double>::max();
        double closestDistanceOld = std::numeric_limits<double>::max();
        for( auto hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (! isECALHit) continue;

            Vector3D hitPos( hit->getPosition() );
            double dToTrack = (hitPos - _tsFirstCaloPos).r();
            if( dToTrack < closestDistanceOld ){
                closestDistanceOld = dToTrack;
                hitTimeOld = hit->getTime();
                closestHitOld = hit;
            }
        }

        if ( hitTimeOld == std::numeric_limits<double>::max() ){
            _tofOld = 0.;
            _hitTimeOld = 0.;
        }
        else{
            _tofOld = hitTimeOld - closestDistanceOld/CLHEP::c_light;
            _hitTimeOld = hitTimeOld;
        }
        if (closestHitOld != nullptr){
            auto simObjects = navClusterHits.getRelatedToObjects(closestHitOld);
            if (simObjects.size() > 0){
                auto simHit = static_cast<SimCalorimeterHit*>(simObjects[0]);
                int fastestIdx = -1;
                double fastestTime = std::numeric_limits<double>::max();
                std::cout<<"Old closest hit has "<<simHit->getNMCContributions()<<" contributions:"<<std::endl;
                for(int j=0; j < simHit->getNMCContributions(); ++j){
                    std::cout<<"# "<<j+1<<"  time: "<<simHit->getTimeCont(j)<<"  energy: "<<simHit->getEnergyCont(j)<<"  mc: "<<simHit->getParticleCont(j)<<"  pdg: "<<simHit->getParticleCont(j)->getPDG()<<std::endl;
                    if ( simHit->getTimeCont(j) < fastestTime ){
                        fastestTime = simHit->getTimeCont(j);
                        fastestIdx = j;
                    }
                }
                MCParticle* mcClosestHit = simHit->getParticleCont(fastestIdx);
                std::cout<<"MC particle: "<<mc<<"  pdg: "<<_pdg<<" OLD closest hit MC:  "<<mcClosestHit<<"  pdg: "<<mcClosestHit->getPDG()<<std::endl;
                _isGoodClosestHitOld = (mc == mcClosestHit);
                _hitClosestOldPos = Vector3D( closestHitOld->getPosition() );
            }
            else{
                std::cout<<"WWATAFAAAAAAAAAACK It is emptY!? "<<simObjects.size()<<std::endl;
                _isGoodClosestHitOld = -1;
                _hitClosestOldPos = Vector3D();
            }
        }


        _tsLastCaloPos = Vector3D( (resultsNew.tsCalo).getReferencePoint() );
        // get NEW TOF part
        CalorimeterHit* closestHitNew = nullptr;
        double hitTimeNew = std::numeric_limits<double>::max();
        double closestDistanceNew = std::numeric_limits<double>::max();
        for( auto hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (! isECALHit) continue;

            Vector3D hitPos( hit->getPosition() );
            double dToTrack = (hitPos - _tsLastCaloPos).r();
            if( dToTrack < closestDistanceNew ){
                closestDistanceNew = dToTrack;
                hitTimeNew = hit->getTime();
                closestHitNew = hit;
            }
        }

        if ( hitTimeNew == std::numeric_limits<double>::max() ){
            _tofNew = 0.;
            _hitTimeNew = 0.;
        }
        else{
            _tofNew = hitTimeNew - closestDistanceNew/CLHEP::c_light;
            _hitTimeNew = hitTimeNew;
        }
        if (closestHitNew != nullptr){
            auto simObjects = navClusterHits.getRelatedToObjects(closestHitNew);
            if (simObjects.size() > 0){
                auto simHit = static_cast<SimCalorimeterHit*>(simObjects[0]);
                int fastestIdx = -1;
                double fastestTime = std::numeric_limits<double>::max();
                std::cout<<"New closest hit has "<<simHit->getNMCContributions()<<" contributions:"<<std::endl;
                for(int j=0; j < simHit->getNMCContributions(); ++j){
                    std::cout<<"# "<<j+1<<"  time: "<<simHit->getTimeCont(j)<<"  energy: "<<simHit->getEnergyCont(j)<<"  mc: "<<simHit->getParticleCont(j)<<"  pdg: "<<simHit->getParticleCont(j)->getPDG()<<std::endl;
                    if ( simHit->getTimeCont(j) < fastestTime ){
                        fastestTime = simHit->getTimeCont(j);
                        fastestIdx = j;
                    }
                }
                MCParticle* mcClosestHit = simHit->getParticleCont(fastestIdx);
                std::cout<<"MC particle: "<<mc<<"  pdg: "<<_pdg<<" NEW closest hit MC:  "<<mcClosestHit<<"  pdg: "<<mcClosestHit->getPDG()<<std::endl;
                _isGoodClosestHitNew = (mc == mcClosestHit);
                _hitClosestNewPos = Vector3D( closestHitNew->getPosition() );
            }
            else{
                std::cout<<"WWATAFAAAAAAAAAACK It is emptY!? "<<simObjects.size()<<std::endl;
                _isGoodClosestHitNew = -1;
                _hitClosestNewPos = Vector3D();
            }
        }



        _nBadPhi = resultsNew.nBadPhi;
        _nBadZ = resultsNew.nBadZ;
        _nSegments = resultsNew.nSegments;

        Vector3D mcPos = Vector3D(mc->getVertex());
        Vector3D mcMom = Vector3D(mc->getMomentum());
        Vector3D recoPos = Vector3D();
        Vector3D recoMom = Vector3D(pfo->getMomentum());
        auto momAtCalo = getTrackMomentum(tsCalo, _bField);
        // double momIP = Vector3D(pfo->getMomentum()).r();
        double massOld = _momentumOld*std::sqrt( (_tofOld*299.792458*_tofOld*299.792458)/(_trackLengthOld*_trackLengthOld) - 1.)*1000;
        double trackLengthTrue = _momentumOld*_tofOld*299.792458*std::sqrt(1./(0.13957039*0.13957039 + _momentumOld*_momentumOld));
        double tofTrue = _trackLengthOld / (_momentumOld*299.792458*std::sqrt(1./(0.13957039*0.13957039 + _momentumOld*_momentumOld)));
        if (std::abs(_pdg) == 211 && std::abs(_hitClosestOldPos.z()) < 2350. && std::abs(massOld - 139.57039) > 100.){
            std::cout<<"Event: "<<_nEvent<<std::endl;
            std::cout<<"PFO: "<<i+1<<std::endl;
            std::cout<<"PDG: "<<_pdg<<std::endl;
            std::cout<<"Vertex: "<<mcPos.x()<<",  "<<mcPos.y()<<",  "<<mcPos.z()<<std::endl;
            std::cout<<"Vertex radius: "<<mcPos.r()<<std::endl;
            std::cout<<"MC momentum IP : ("<<mc->getMomentum()[0]<<",  "<<mc->getMomentum()[1]<<",  "<<mc->getMomentum()[2]<<")"<<std::endl;
            std::cout<<"Reco momentum IP : ("<<pfo->getMomentum()[0]<<",  "<<pfo->getMomentum()[1]<<",  "<<pfo->getMomentum()[2]<<")"<<std::endl;
            std::cout<<"MC abs p IP: "<<Vector3D(mc->getMomentum()).r()<<std::endl;
            std::cout<<"Reco abs p IP: "<<Vector3D(pfo->getMomentum()).r()<<std::endl;
            std::cout<<"HM momentumOld: "<<_momentumOld<<std::endl;
            std::cout<<"HM momentumNew: "<<_momentumNew<<std::endl;
            std::cout<<"TOF OLD: "<<_tofOld<<std::endl;
            std::cout<<"TOF New: "<<_tofNew<<std::endl;
            std::cout<<"trackLength Old: "<<_trackLengthOld<<std::endl;
            std::cout<<"trackLength New: "<<_trackLengthNew<<std::endl;
            std::cout<<"trackLength SHOULD BE for true mass: "<<trackLengthTrue<<std::endl;
            std::cout<<"TOF SHOULD BE for true mass: "<<tofTrue<<std::endl;
            std::cout<<"trackLength from helix formula reco: "<<getHelixLength(0., _tsFirstCaloPos.z(), Vector3D(momAtCalo[0], momAtCalo[1], momAtCalo[2]), 3.5)<<std::endl;
            std::cout<<"trackLength from helix formula: "<<getHelixLength(mc->getVertex()[2], _tsFirstCaloPos.z(), mcMom, 3.5)<<std::endl;
            std::cout<<"nSegments: "<<_nSegments<<std::endl;
            std::cout<<"nBadZ: "<<_nBadZ<<std::endl;
            std::cout<<"tsCaloOld: "<<_tsFirstCaloPos<<std::endl;
            std::cout<<"tsCaloNew: "<<_tsLastCaloPos<<std::endl;
            std::cout<<"Particle BORN time: "<<mc->getTime()<<std::endl;
            if (_tofNew != 0 && _tofOld != 0){
                std::cout<<"NEW beta: "<<_trackLengthNew/(_tofNew*299.792458)<<std::endl;
                std::cout<<"OLD beta: "<<_trackLengthOld/(_tofOld*299.792458)<<std::endl;
            }
            if (_trackLengthNew != 0 && _trackLengthOld != 0){
                std::cout<<"NEW mass: "<<_momentumNew*std::sqrt( (_tofNew*299.792458*_tofNew*299.792458)/(_trackLengthNew*_trackLengthNew) - 1.)*1000<<std::endl;
                std::cout<<"OLD mass: "<<massOld<<std::endl;
            }
            auto hits = cluster->getCalorimeterHits();
            auto getEarliestHit = [this](CalorimeterHit* lhs, CalorimeterHit* rhs){
                Vector3D lhsPos( lhs->getPosition() );
                Vector3D rhsPos( rhs->getPosition() );
                double lhsTimeCorr = lhs->getTime() - (lhsPos - _tsFirstCaloPos).r()/299.792458;
                double rhsTimeCorr = rhs->getTime() - (rhsPos - _tsFirstCaloPos).r()/299.792458;
                return lhsTimeCorr < rhsTimeCorr;
            };

            std::sort(hits.begin(), hits.end(), getEarliestHit);
            Vector3D hitPos( hits[0]->getPosition() );
            double dToTrack = (hitPos - _tsFirstCaloPos).r();
            double earliestTof = hits[0]->getTime() - dToTrack/299.792458;
            std::cout<<"Mass using earliest corrected hit: "<<_momentumOld*std::sqrt( (earliestTof*299.792458*earliestTof*299.792458)/(_trackLengthOld*_trackLengthOld) - 1.)*1000<<std::endl;
            for(int h=0; h< hits.size(); ++h){
                auto hit = hits[h];
                Vector3D hitPos( hit->getPosition() );
                double dToTrack = (hitPos - _tsFirstCaloPos).r();
                std::cout<<"Hit #"<<h+1<<"  time: "<<hit->getTime()<<"  corrected time: "<<hit->getTime() - dToTrack/299.792458<<"   is closest:  "<<(hit == closestHitOld)<<std::endl;
            }

            drawDisplay(this, evt, drawPFO, pfo, mc);
        }
        _tree->Fill();
    }
}

void TrackLengthDebug::end(){
    _file->Write();
}

Vars TrackLengthDebug::getTrackLengthNew(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    Vars vars;
    if ( pfo->getTracks().empty() ) return vars;
    Track* track = pfo->getTracks()[0];

    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    TrackImpl lastGoodRefittedTrack;
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
        covMatrix[0]  = 1e+06; //sigma_d0^2
        covMatrix[2]  = 100.; //sigma_phi0^2
        covMatrix[5]  = 0.00001; //sigma_omega^2
        covMatrix[9]  = 1e+06; //sigma_z0^2
        covMatrix[14] = 100.; //sigma_tanl^2
        double maxChi2PerHit = 100.;
        std::unique_ptr<MarlinTrk::IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        TrackStateImpl preFit = *track->getTrackState(TrackState::AtLastHit);
        preFit.setCovMatrix( covMatrix );
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::backward, &preFit, bField, maxChi2PerHit);
        //if fit fails, try also fit forward
        if (errorFit != 0){
            std::cout<<" Fit failed, trying to fit forward"<<std::endl;
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0){
            std::cout<<" Fit failed completely, ignoring the subtrack"<<std::endl;
            continue;
        }
        lastGoodRefittedTrack = refittedTrack;
        //here hits are sorted by rho=sqrt(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        //Find which way to loop over the array of hits. We need to loop in the direction of the track.
        bool loopForward = true;
        double zFirst = std::abs( hitsInFit.front().first->getPosition()[2] );
        double zLast = std::abs( hitsInFit.back().first->getPosition()[2] );

        // OPTIMIZE: 10 mm is just a round number. With very small z difference it is more robust to use rho, to be sure z difference is not caused by tpc Z resolution or multiple scattering
        if ( std::abs(zLast - zFirst) > 10. ){
            if ( zLast < zFirst ) loopForward = false;
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            std::cout<<"Adding last and calorimeter hits"<<std::endl;
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
        }
    }

    const TrackStateImpl* tsCalo = static_cast<const TrackStateImpl*> (lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter) );
    if ( pfo->getClusters().size() > 0 && tsCalo != nullptr ){
        trackStates.push_back( *(tsCalo) );
        vars.tsCalo = *lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter);
    }



    int nTrackStates = trackStates.size();
    vars.nSegments = nTrackStates - 1;
    if (nTrackStates <= 1) return vars;
    std::cout<<"We found "<<nTrackStates<<" track states"<<std::endl;


    for( int j=1; j < nTrackStates; ++j ){
        auto ts1 = trackStates[j-1];
        auto ts2 = trackStates[j];
        double omega = ts1.getOmega();
        double tanL = ts1.getTanLambda();
        double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
        double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
        double phi1 = ts1.getPhi();
        double phi2 = ts2.getPhi();
        double dPhi = std::abs(phi2 - phi1);
        if (dPhi <= M_PI && (phi1 - phi2)/omega < 0) vars.nBadPhi++;
        else if (dPhi > M_PI && phi2/omega < 0) vars.nBadPhi++;

        if ( (z2-z1)/tanL < 0) vars.nBadZ++;

        double arcLength = std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
        std::cout<<"Segment "<<j<<":  track length between z1 = "<<z1<<"  and z2 = "<<z2<<"   is   "<<arcLength<<std::endl;
        auto mom = getTrackMomentum(&ts1, _bField);
        double AbsMomentum = std::hypot(mom[0], mom[1], mom[2]);

        vars.trackLength += arcLength;
        if(AbsMomentum == 0) continue;
        vars.momentumHM += arcLength/(AbsMomentum*AbsMomentum);
    }
    if (vars.momentumHM != 0) vars.momentumHM = std::sqrt(vars.trackLength/vars.momentumHM);

    return vars;
}


Vars TrackLengthDebug::getTrackLengthOld(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    Vars vars;
    if ( pfo->getTracks().empty() ) return vars;
    Track* track = pfo->getTracks()[0];

    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
        covMatrix[0]  = 1e+06; //sigma_d0^2
        covMatrix[2]  = 100.; //sigma_phi0^2
        covMatrix[5]  = 0.00001; //sigma_omega^2
        covMatrix[9]  = 1e+06; //sigma_z0^2
        covMatrix[14] = 100.; //sigma_tanl^2
        double maxChi2PerHit = 100.;
        std::unique_ptr<MarlinTrk::IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        TrackStateImpl preFit = *track->getTrackState(TrackState::AtLastHit);
        preFit.setCovMatrix( covMatrix );
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::backward, &preFit, bField, maxChi2PerHit);
        //if fit fails, try also fit forward
        if (errorFit != 0){
            std::cout<<" Fit failed, trying to fit forward"<<std::endl;
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0){
            std::cout<<" Fit failed completely, ignoring the subtrack"<<std::endl;
            continue;
        }
        //here hits are sorted by rho=sqrt(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        //Find which way to loop over the array of hits. We need to loop in the direction of the track.
        bool loopForward = true;
        double zFirst = std::abs( hitsInFit.front().first->getPosition()[2] );
        double zLast = std::abs( hitsInFit.back().first->getPosition()[2] );

        // OPTIMIZE: 10 mm is just a round number. With very small z difference it is more robust to use rho, to be sure z difference is not caused by tpc Z resolution or multiple scattering
        if ( std::abs(zLast - zFirst) > 10. ){
            if ( zLast < zFirst ) loopForward = false;
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            std::cout<<"Adding last and calorimeter hits"<<std::endl;
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }

    int nTrackStates = trackStates.size();
    vars.nSegments = nTrackStates - 1;
    if (nTrackStates <= 1) return vars;
    std::cout<<"We found "<<nTrackStates<<" track states"<<std::endl;
    for( int j=1; j < nTrackStates; ++j ){
        auto ts1 = trackStates[j-1];
        auto ts2 = trackStates[j];
        double omega = ts1.getOmega();
        double tanL = ts1.getTanLambda();
        double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
        double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
        double phi1 = ts1.getPhi();
        double phi2 = ts2.getPhi();
        double dPhi = std::abs(phi2 - phi1);
        if (dPhi <= M_PI && (phi1 - phi2)/omega < 0) vars.nBadPhi++;
        else if (dPhi > M_PI && phi2/omega < 0) vars.nBadPhi++;

        if ( (z2-z1)/tanL < 0) vars.nBadZ++;

        double arcLength = std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
        std::cout<<"Segment "<<j<<":  track length between z1 = "<<z1<<"  and z2 = "<<z2<<"   is   "<<arcLength<<std::endl;
        auto mom = getTrackMomentum(&ts1, _bField);
        double AbsMomentum = std::hypot(mom[0], mom[1], mom[2]);

        vars.trackLength += arcLength;
        if(AbsMomentum == 0) continue;
        vars.momentumHM += arcLength/(AbsMomentum*AbsMomentum);
    }
    if (vars.momentumHM == 0) vars.momentumHM = 0.;
    else vars.momentumHM = std::sqrt(vars.trackLength/vars.momentumHM);

    return vars;
}