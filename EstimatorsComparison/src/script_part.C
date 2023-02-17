// #include "TrackLengthDebug.h"
// #include "TrackLengthUtils.h"

// #include "EVENT/LCCollection.h"
// #include "UTIL/PIDHandler.h"

// #include "marlin/Global.h"
// #include "marlin/ProcessorEventSeeder.h"
// #include "marlin/VerbosityLevels.h"
// #include "marlinutil/GeometryUtil.h"
// #include "MarlinTrk/Factory.h"
// #include "EVENT/SimTrackerHit.h"
// #include "marlinutil/DDMarlinCED.h"
// #include "UTIL/TrackTools.h"
// #include "UTIL/LCRelationNavigator.h"
// #include "EVENT/MCParticle.h"
// #include "TCanvas.h"
// #include "TSystem.h"
// #include "marlinutil/MarlinUtil.h"

// using namespace TrackLengthUtils;
// using namespace EVENT;
// using namespace UTIL;
// using dd4hep::rec::Vector3D;
// using std::vector;
// using std::string;

// TrackLengthDebug aTrackLengthDebug;

// // header
// // class TrackLengthDebug : public marlin::Processor, EventDisplayer {
// //     friend class EventDisplayer;
// TrackLengthDebug::TrackLengthDebug() : marlin::Processor("TrackLengthDebug"), EventDisplayer(this) {
//     _description = "TrackLengthDebug debugs track length";

//     registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
//                             "ReconstructedParticleCollection",
//                             "Name of the ReconstructedParticle collection",
//                             _pfoCollectionName,
//                             std::string("PandoraPFOs") );
// }


// void TrackLengthDebug::init(){
//     marlin::Global::EVENTSEEDER->registerProcessor(this);

//     _outputParNames = {"trackLengthToSET", "trackLengthToEcal", "momentumHMToSET", "momentumHMToEcal"};
//     _bField = MarlinUtil::getBzAtOrigin();

//     _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
//     _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
//     _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
//     _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
//     _trkSystem->init();

//     prepareRootTree();

// }


// void drawPFO(ReconstructedParticle* pfo){
//     Track* track = pfo->getTracks()[0];
//     vector<Track*> tracks = getSubTracks(track);
//     int nHits = 0;
//     for(auto* track: tracks){
//         auto hits = track->getTrackerHits();
//         for (auto* hit : hits){
//             ++nHits;
//             auto pos = hit->getPosition();
//             int type = 0; // point
//             int layer = 1; // doesn't matter
//             int size = 4; // larger point

//             unsigned long color = 0x3232a8;
//             if ( (nHits) % 20 == 0){
//                 color = 0xfa0730;
//                 size = 8;
//             }
//             ced_hit_ID(pos[0], pos[1], pos[2], type, layer, size, color, 0 ); // tracker hits
//         }
//     }
//     std::vector<Cluster*> clusters = pfo->getClusters();
//     for(auto* cluster: clusters){
//         auto hits = cluster->getCalorimeterHits();
//         for (auto* hit: hits){
//             auto pos = hit->getPosition();
//             int type = 0; // point
//             int layer = 1; // doesn't matter
//             int size = 6; // larger point
//             unsigned long color = 0xbf2659;
//             ced_hit_ID(pos[0], pos[1], pos[2], type, layer, size, color, 0 ); // tracker hits
//         }
//     }
// }

// void TrackLengthDebug::processEvent(EVENT::LCEvent * evt){
//     ++_nEvent;
//     streamlog_out(DEBUG8)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

//     LCCollection* pfos = evt->getCollection(_pfoCollectionName);

//     PIDHandler pidHandler( pfos );
//     int algoID = pidHandler.addAlgorithm( name(), _outputParNames );

//     LCRelationNavigator nav ( evt->getCollection("RecoMCTruthLink") );

//     auto getArcLenDefault = [](double dPhi, double omega, double z1, double z2){ return std::sqrt( std::pow(dPhi/omega, 2) + std::pow(z2-z1, 2) );};
//     auto getArcLenTanL = [](double dPhi, double omega, double tanL){ return std::abs(dPhi/omega) * std::sqrt( 1. + std::pow(tanL, 2) );};
//     auto getArcLenZ = [](double z1, double z2, double tanL){ return std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );};

//     for (int i=0; i<pfos->getNumberOfElements(); ++i){
//         streamlog_out(DEBUG7)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
//         ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

//         const vector<LCObject*>& objects = nav.getRelatedToObjects(pfo);
//         const std::vector<float>& weights = nav.getRelatedToWeights(pfo);
//         // int max_i = std::max_element(weights.begin(), weights.end(), [](float lhs, float rhs){return (int(lhs)%10000)/1000. < (int(rhs)%10000)/1000.;}) - weights.begin();
//         int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
//         if ( ( int(weights[max_i])%10000 )/1000. == 0 ){
//             max_i = std::max_element(weights.begin(), weights.end(), [](float lhs, float rhs){return (int(lhs)/10000)/1000. < (int(rhs)/10000)/1000.;}) - weights.begin();
//         }
//         auto* mc = static_cast<MCParticle*> (objects[max_i]);
//         std::cout<<"MC1:"<<mc<<std::endl;

//         // const auto* mc2 = static_cast<const MCParticle*> (nav.getRelatedToMaxWeightObject(pfo,  MarlinUtil::getTrackWeight));
//         // float test = nav.getRelatedToMaxWeight(pfo);
//         // float test2 = nav.getRelatedToMaxWeight(pfo, MarlinUtil::getClusterWeight);
//         // std::cout<<"MC2:"<<mc2<<"    weight track: "<<test<<"          cluster: "<<test2<<std::endl;

//         _pdg = mc->getPDG();

//         int nClusters = pfo->getClusters().size();
//         int nTracks = pfo->getTracks().size();

//         if( nClusters != 1 || nTracks != 1) continue;
//         Track* track = pfo->getTracks()[0];

//         _d0Ip = track->getD0();
//         _z0Ip = track->getZ0();
//         _omegaIp = track->getOmega();
//         _tanLIp = track->getTanLambda();
//         _phiIp = track->getPhi();

//         const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
//         _zCalo = tsCalo->getReferencePoint()[2] + tsCalo->getZ0();
//         std::array<double, 3> momCalo = getTrackMomentum(tsCalo, _bField);
//         _momCalo = Vector3D(momCalo[0], momCalo[1], momCalo[2]);
//         _d0Calo = tsCalo->getD0();
//         _z0Calo = tsCalo->getZ0();
//         _omegaCalo = tsCalo->getOmega();
//         _tanLCalo = tsCalo->getTanLambda();
//         _phiCalo = tsCalo->getPhi();
//         //only barrel to check
//         // if ( std::abs(tsCalo->getReferencePoint()[2] + tsCalo->getZ0() ) > 2300. ) continue;
//         vector<Track*> subTracks = getSubTracks(track);
//         vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, _bField);
//         std::vector<EVENT::SimTrackerHit*> simHits = convertHitsToSimHits(evt, subTracks, _trkSystem, _bField);


//         _trackLengthDefault = getTrackLengthDefault(trackStates);
//         _trackLengthTanL = getTrackLengthTanL(trackStates);
//         _trackLengthZ = getTrackLengthZ(trackStates);
    
//         double speedOfLight = 299.792458;
//         _momIp = Vector3D( pfo->getMomentum() );
//         _momentum = _momIp.r();
//         _tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "timeOfFlight"); // in ns
//         if ( std::pow(_tof*speedOfLight/_trackLengthDefault, 2) - 1 > 0 ) _massDefault = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthDefault, 2) - 1 );
//         else _massDefault = 0;

//         if ( std::pow(_tof*speedOfLight/_trackLengthTanL, 2) - 1 > 0 ) _massTanL = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthTanL, 2) - 1 );
//         else _massTanL = 0;

//         if ( std::pow(_tof*speedOfLight/_trackLengthZ, 2) - 1 > 0 ) _massZ = _momentum * std::sqrt( std::pow(_tof*speedOfLight/_trackLengthZ, 2) - 1 );
//         else _massZ = 0;

//         streamlog_out(DEBUG5)<<"Final results for the "<<i+1<<" PFO"<<std::endl;
//         streamlog_out(DEBUG5)<<"Track length using default: "<<_trackLengthDefault<<" mm"<<std::endl;
//         streamlog_out(DEBUG5)<<"Track length using tanL: "<<_trackLengthTanL<<" mm"<<std::endl;
//         streamlog_out(DEBUG5)<<"Track length using Z: "<<_trackLengthZ<<" mm"<<std::endl;
//         streamlog_out(DEBUG5)<<std::endl<<std::endl;


//         if(std::abs(_trackLengthDefault - _trackLengthZ) > 5. ){
//             double cumTrackLengthDefault = 0;
//             double cumTrackLengthTanL = 0;
//             double cumTrackLengthZ = 0;
//             double cumMcTrackLengthDefault = 0;
//             double cumMcTrackLengthTanL = 0;
//             double cumMcTrackLengthZ = 0;


//             _hD0 = new TGraph();
//             _hD0->SetTitle("d0; Hit number; d0 (mm)");
//             _hZ0 = new TGraph();
//             _hZ0->SetTitle("z0; Hit number; z0 (mm)");
//             _hOmega = new TGraph();
//             _hOmega->SetTitle("#Omega; Hit number; #Omega (1/mm)");
//             _hTanL = new TGraph();
//             _hTanL->SetTitle("tan #lambda; Hit number; tan #lambda");
//             _hPhi = new TGraph();
//             _hPhi->SetTitle("#varphi; Hit number; #varphi");
//             _hZ = new TGraph();
//             _hZ->SetTitle("z; Hit number; z (mm)");
//             _hPt = new TGraph();
//             _hPt->SetTitle("pt; Hit number; pt (GeV)");
//             _hPz = new TGraph();
//             _hPz->SetTitle("pz; Hit number; pz (GeV)");
//             _hP = new TGraph();
//             _hP->SetTitle("p; Hit number; p (GeV)");
//             _hTrkLenDefault = new TGraph();
//             _hTrkLenDefault->SetTitle("Cumulative track length DEFAULT; Hit number; #ell_{track} (mm)");
//             _hTrkLenTanL = new TGraph();
//             _hTrkLenTanL->SetTitle("Cumulative track length TANL; Hit number; #ell_{track} (mm)");
//             _hTrkLenZ = new TGraph();
//             _hTrkLenZ->SetTitle("Cumulative track length Z; Hit number; #ell_{track} (mm)");

//             _hMcOmega = new TGraph();
//             _hMcOmega->SetTitle("MC #Omega; Hit number; #Omega (1/mm)");
//             _hMcOmega->SetLineColor(2);
//             _hMcTanL = new TGraph();
//             _hMcTanL->SetTitle("MC tan #lambda; Hit number; tan #lambda");
//             _hMcTanL->SetLineColor(2);
//             _hMcPhi = new TGraph();
//             _hMcPhi->SetTitle("MC #varphi; Hit number; #varphi");
//             _hMcPhi->SetLineColor(2);
//             _hMcZ = new TGraph();
//             _hMcZ->SetTitle("MC z; Hit number; z (mm)");
//             _hMcZ->SetLineColor(2);
//             _hMcPt = new TGraph();
//             _hMcPt->SetTitle("MC pt; Hit number; pt (GeV)");
//             _hMcPt->SetLineColor(2);
//             _hMcPz = new TGraph();
//             _hMcPz->SetTitle("MC pz; Hit number; pz (GeV)");
//             _hMcPz->SetLineColor(2);
//             _hMcP = new TGraph();
//             _hMcP->SetTitle("MC p; Hit number; p (GeV)");
//             _hMcP->SetLineColor(2);
//             _hMcTrkLenDefault = new TGraph();
//             _hMcTrkLenDefault->SetTitle("MC Cumulative track length DEFAULT; Hit number; #ell_{track} (mm)");
//             _hMcTrkLenDefault->SetLineColor(2);
//             _hMcTrkLenTanL = new TGraph();
//             _hMcTrkLenTanL->SetTitle("MC Cumulative track length TANL; Hit number; #ell_{track} (mm)");
//             _hMcTrkLenTanL->SetLineColor(2);
//             _hMcTrkLenZ = new TGraph();
//             _hMcTrkLenZ->SetTitle("MC Cumulative track length Z; Hit number; #ell_{track} (mm)");
//             _hMcTrkLenZ->SetLineColor(2);

//             TFile plot(Form("evt_%d_pfo_%d.root", _nEvent, i+1), "RECREATE");
//             plot.cd();
//             for(int j=0; j < trackStates.size(); ++j){
//                 //reco
//                 auto ts = trackStates[j];
//                 _hD0->SetPoint(j, j+1, ts.getD0());
//                 _hZ0->SetPoint(j, j+1, ts.getZ0());
//                 _hOmega->SetPoint(j, j+1, ts.getOmega());
//                 _hTanL->SetPoint(j, j+1, ts.getTanLambda());
//                 _hPhi->SetPoint(j, j+1, ts.getPhi());
//                 _hZ->SetPoint(j, j+1, ts.getReferencePoint()[2] + ts.getZ0());

//                 std::array<double, 3> mom = getTrackMomentum(&ts, _bField);
//                 _hPt->SetPoint(j, j+1, std::hypot(mom[0], mom[1]) );
//                 _hPz->SetPoint(j, j+1, mom[2]);
//                 _hP->SetPoint(j, j+1, std::hypot(mom[0], mom[1], mom[2]));

//                 if(j != 0){
//                     double omega = trackStates[j-1].getOmega();
//                     double tanL = trackStates[j-1].getTanLambda();
//                     double z1 = trackStates[j-1].getReferencePoint()[2] + trackStates[j-1].getZ0();
//                     double z2 = trackStates[j].getReferencePoint()[2] + trackStates[j].getZ0();
//                     double dPhi = std::abs( trackStates[j].getPhi() - trackStates[j-1].getPhi() );
//                     if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

//                     cumTrackLengthDefault += getArcLenDefault(dPhi, omega, z1, z2);
//                     cumTrackLengthTanL += getArcLenTanL(dPhi, omega, tanL);
//                     cumTrackLengthZ += getArcLenZ(z1, z2, tanL);
//                 }
//                 _hTrkLenDefault->SetPoint(j, j+1, cumTrackLengthDefault);
//                 _hTrkLenTanL->SetPoint(j, j+1, cumTrackLengthTanL);
//                 _hTrkLenZ->SetPoint(j, j+1, cumTrackLengthZ);


//                 //true
//                 auto hit = simHits[j];
//                 double c_light = 299.792458; // mm/ns
//                 double pt = std::hypot(hit->getMomentum()[0], hit->getMomentum()[1]);
//                 double omega = (1e-6 * c_light * _bField) / pt;
//                 double tanL = hit->getMomentum()[2] / pt;
//                 _hMcOmega->SetPoint(j, j+1, omega);
//                 _hMcTanL->SetPoint(j, j+1, tanL);
//                 _hMcPhi->SetPoint(j, j+1, std::atan2(hit->getMomentum()[1], hit->getMomentum()[0]));
//                 _hMcZ->SetPoint(j, j+1, hit->getPosition()[2]);
//                 _hMcPt->SetPoint(j, j+1, pt);
//                 _hMcPz->SetPoint(j, j+1, hit->getMomentum()[2]);
//                 _hMcP->SetPoint(j, j+1, std::hypot(pt, hit->getMomentum()[2]));

//                 if(j != 0){
//                     auto prevHit = simHits[j-1];

//                     double pt = std::hypot(prevHit->getMomentum()[0], prevHit->getMomentum()[1]);
//                     double omega = (1e-6 * c_light * _bField) / pt;
//                     double tanL = prevHit->getMomentum()[2] / pt;
//                     double z1 = prevHit->getPosition()[2];
//                     double z2 = hit->getPosition()[2];
//                     double phi1 = std::atan2(prevHit->getMomentum()[1], prevHit->getMomentum()[0]);
//                     double phi2 = std::atan2(hit->getMomentum()[1], hit->getMomentum()[0]);
//                     double dPhi = std::abs( phi2 - phi1 );
//                     if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

//                     cumMcTrackLengthDefault += getArcLenDefault(dPhi, omega, z1, z2);
//                     cumMcTrackLengthTanL += getArcLenTanL(dPhi, omega, tanL);
//                     cumMcTrackLengthZ +=getArcLenZ(z1, z2, tanL);
//                 }
//                 _hMcTrkLenDefault->SetPoint(j, j+1, cumMcTrackLengthDefault);
//                 _hMcTrkLenTanL->SetPoint(j, j+1, cumMcTrackLengthTanL);
//                 _hMcTrkLenZ->SetPoint(j, j+1, cumMcTrackLengthZ);

//             }
//             _hD0->Write("d0");
//             _hZ0->Write("z0");
//             _hOmega->Write("omega");
//             _hTanL->Write("tanl");
//             _hPhi->Write("phi");
//             _hZ->Write("z");
//             _hPt->Write("pt");
//             _hPz->Write("pz");
//             _hP->Write("p");
//             _hTrkLenDefault->Write("trk_len_default");
//             _hTrkLenTanL->Write("trk_len_tanl");
//             _hTrkLenZ->Write("trk_len_z");

//             _hMcOmega->Write("mc_omega");
//             _hMcTanL->Write("mc_tanl");
//             _hMcPhi->Write("mc_phi");
//             _hMcZ->Write("mc_z");
//             _hMcPt->Write("mc_pt");
//             _hMcPz->Write("mc_pz");
//             _hMcP->Write("mc_p");
//             _hMcTrkLenDefault->Write("mc_trk_len_default");
//             _hMcTrkLenTanL->Write("mc_trk_len_tanl");
//             _hMcTrkLenZ->Write("mc_trk_len_z");
//             // plot.Write();
//             plot.Close();
//             std::cout<<"EVENT: "<<_nEvent<<"    PFO: "<<i+1<<"    len diff: "<<std::abs(_trackLengthDefault - _trackLengthZ)<<std::endl;
//             std::cout<<"PDG: "<<_pdg<<"    mom: "<<_momIp.r()<<std::endl;
//             std::cout<<"pt: "<<_momIp.rho()<<"    pz: "<<_momIp.z()<<std::endl;
//             std::cout<<"mass_default: "<<_massDefault<<"    mass_tanl: "<<_massTanL<<"    mass_z: "<<_massZ<<std::endl;
//             drawDisplay(this, evt, drawPFO, pfo);
//         }
//         _tree->Fill();


//     }
// }

// void TrackLengthDebug::end(){
//     _file->Write();
// }


// double TrackLengthDebug::getTrackLengthTanL(const std::vector<IMPL::TrackStateImpl>& trackStates){
//     double trackLength = 0.;
//     int nTrackStates = trackStates.size();
//     //exclude last track state at the ECal
//     for( int j=1; j < nTrackStates; ++j ){
//         //we check which track length formula to use
//         double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
//         double arcLength;
//         // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
//         if ( nTurns <= 0.5 ) arcLength = getHelixArcLengthTanL( trackStates[j-1], trackStates[j] );
//         else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

//         trackLength += arcLength;
//     }
//     return trackLength;
// }



// double TrackLengthDebug::getTrackLengthZ(const std::vector<IMPL::TrackStateImpl>& trackStates){
//     double trackLength = 0.;
//     int nTrackStates = trackStates.size();
//     //exclude last track state at the ECal
//     for( int j=1; j < nTrackStates; ++j ){
//         //we check which track length formula to use
//         double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
//         double arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

//         trackLength += arcLength;
//     }
//     return trackLength;
// }


// float TrackLengthDebug::getParameterFromPID(ReconstructedParticle* pfo, PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
//     int algorithmID = pidHandler.getAlgorithmID(algorithmName);
//     const ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
//     const std::vector<float>& parameters = pfoPID.getParameters();
//     int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
//     return parameters[parIdx]; 
// }


// void TrackLengthDebug::prepareRootTree(){
//     _file.reset( new TFile("results.root", "RECREATE") );
//     _tree.reset( new TTree("treename", "treename") );

//     _tree->Branch("momentum", &_momentum);
//     _tree->Branch("tof", &_tof);
//     _tree->Branch("trackLengthDefault", &_trackLengthDefault);
//     _tree->Branch("trackLengthTanL", &_trackLengthTanL);
//     _tree->Branch("trackLengthZ", &_trackLengthZ);
//     _tree->Branch("massDefault", &_massDefault);
//     _tree->Branch("massTanL", &_massTanL);
//     _tree->Branch("massZ", &_massZ);

//     _tree->Branch("pdg", &_pdg);
//     _tree->Branch("mom_ip", &_momIp);
//     _tree->Branch("d0_ip", &_d0Ip);
//     _tree->Branch("z0_ip", &_z0Ip);
//     _tree->Branch("omega_ip", &_omegaIp);
//     _tree->Branch("tanL_ip", &_tanLIp);
//     _tree->Branch("phi_ip", &_phiIp);

//     _tree->Branch("mom_calo", &_momCalo);
//     _tree->Branch("d0_calo", &_d0Calo);
//     _tree->Branch("z0_calo", &_z0Calo);
//     _tree->Branch("omega_calo", &_omegaCalo);
//     _tree->Branch("tanL_calo", &_tanLCalo);
//     _tree->Branch("phi_calo", &_phiCalo);

//     _tree->Branch("z_calo", &_zCalo);
    
// }
