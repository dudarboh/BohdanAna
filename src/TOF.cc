#include "TOF.h"
#include "BohdanUtils.h"
#include "marlinutil/CalorimeterHitType.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TrackState.h"
#include "IMPL/TrackStateImpl.h"
#include "UTIL/ILDConf.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TGraph.h"
#include "TF1.h"
#include "marlin/VerbosityLevels.h"
#include <limits>

using namespace EVENT;
using dd4hep::rec::Vector3D;
using CLHEP::RandGauss;



std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, Vector3D posAtEcal, Vector3D momAtEcal, int maxEcalLayer ){
    std::vector<CalorimeterHit*> selectedHits(maxEcalLayer, nullptr);
    std::vector<double> minDistances(maxEcalLayer, std::numeric_limits<double>::max());

    for ( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        int layer = hitType.layer();
        if ( (!isECALHit) || ( layer >= maxEcalLayer) ) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToLine = (hitPos - posAtEcal).cross(momAtEcal.unit()).r();
        if ( dToLine < minDistances[layer] ){
            minDistances[layer] = dToLine;
            selectedHits[layer] = hit;
        }
    }
    selectedHits.erase( std::remove_if( selectedHits.begin(), selectedHits.end(), [](CalorimeterHit* h) { return h == nullptr; } ), selectedHits.end() );

    return selectedHits;
}


std::pair<int, double> getTofClosest( EVENT::Cluster* cluster, dd4hep::rec::Vector3D posAtEcal, double timeResolution){
    double hitTime = std::numeric_limits<double>::max();
    double closestDistance = std::numeric_limits<double>::max();
    int layer = -1;
    for( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - posAtEcal).r();
        if( dToTrack < closestDistance ){
            closestDistance = dToTrack;
            hitTime = hit->getTime();
            layer = hitType.layer();
        }
    }
    if ( hitTime == std::numeric_limits<double>::max() ) return {-1, 0.};
    return {layer, RandGauss::shoot(hitTime, timeResolution) - closestDistance/CLHEP::c_light};
}



double getTofFrankAvg( std::vector<EVENT::CalorimeterHit*> selectedHits, Vector3D posAtEcal, double timeResolution){
    double tof = 0.;
    if ( selectedHits.empty() ) return tof;

    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - posAtEcal).r();
        tof += RandGauss::shoot(hit->getTime(), timeResolution) - dToTrack/CLHEP::c_light;
    }
    return tof/selectedHits.size();
}


double getTofFrankFit( std::vector<EVENT::CalorimeterHit*> selectedHits, Vector3D posAtEcal, double timeResolution){
    double tof = 0.;
    if ( selectedHits.empty() ) return tof;
    else if ( selectedHits.size() == 1 ){
        //we can't fit 1 point, but lets return something reasonable
        Vector3D hitPos( selectedHits[0]->getPosition() );
        double dToTrack = (hitPos - posAtEcal).r();
        return RandGauss::shoot(selectedHits[0]->getTime(), timeResolution) - dToTrack/CLHEP::c_light;
    }

    std::vector <double> x, y;
    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - posAtEcal).r();
        x.push_back(dToTrack);
        double time = RandGauss::shoot(hit->getTime(), timeResolution);
        y.push_back(time);
    }

    TGraph gr(x.size(), x.data(), y.data());
    gr.Fit("pol1", "Q");
    return gr.GetFunction("pol1")->GetParameter(0);
}


double getTofSET(EVENT::Track* track, double timeResolution){
    EVENT::TrackerHit* setHit = getSETHit(track);
    if (setHit == nullptr) return 0.;
    auto stripObjects = setHit->getRawHits();

    if ( stripObjects.empty() ) return 0.;
    else if (stripObjects.size() == 1){
        streamlog_out(WARNING)<<"Found only one SET strip hit, how is this possible!? Writing TOF from a single strip."<<std::endl;
        auto strip = static_cast<EVENT::TrackerHitPlane*> (stripObjects[0]);
        return RandGauss::shoot(strip->getTime(), timeResolution);
    }
    if (stripObjects.size() > 2) streamlog_out(WARNING)<<"Found more than two SET strip hits, how is this possible!? TOF is an average of the first two elements."<<std::endl;
    auto stripFront = static_cast<EVENT::TrackerHitPlane*> (stripObjects[0]);
    auto stripBack = static_cast<EVENT::TrackerHitPlane*> (stripObjects[1]);
    double timeFront = RandGauss::shoot(stripFront->getTime(), timeResolution);
    double timeBack = RandGauss::shoot(stripBack->getTime(), timeResolution);
    return (timeFront + timeBack)/2.;
}


double getTofPhotonTrue(EVENT::MCParticle* mc){
    // ignore non photons and those who doesn't have a straight path to the calorimeter.
    if (mc->getPDG() != 22 || mc->isDecayedInTracker()) return 0.;

    Vector3D startPos( mc->getVertex() );
    //find intersection points between photon momentum line and ECAL surface planes
    Vector3D finishPos = getPhotonAtCalorimeter(mc);

    //NOTE: photon might not be created at event_time = 0
    return mc->getTime() + (finishPos - startPos).r()/CLHEP::c_light; // in ns
}
