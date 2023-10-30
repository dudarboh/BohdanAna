#include "BohdanUtils.h"
#include "UTIL/ILDConf.h"
#include "marlinutil/GeometryUtil.h"

#include <cstring>
#include <iomanip>

using dd4hep::rec::LayeredCalorimeterData;
using dd4hep::DetType;
using dd4hep::rec::Vector3D;

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVirtualMemoryUsage(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPhysicalMemoryUsage(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track){
    std::vector<EVENT::Track*> subTracks;
    // add track itself, which contains VXD+FTD+SIT+TPC hits of the first curl.
    subTracks.push_back(track);

    int nSubTracks = track->getTracks().size();
    if (nSubTracks <= 1) return subTracks;

    UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
    auto isTPCHit = [&encoder](EVENT::TrackerHit* hit) -> bool {
        encoder.setValue( hit->getCellID0() ) ;
        int subdet = encoder[ UTIL::LCTrackerCellID::subdet() ];
        return subdet == UTIL::ILDDetID::TPC;
    };

    int indexOfFirstTPCCurl = 0;
    for(int i = 0; i < nSubTracks; ++i){
        EVENT::Track* subTrack = track->getTracks()[i];
        auto hits = subTrack->getTrackerHits();
        if ( std::find_if(hits.begin(), hits.end(), isTPCHit) != hits.end() ){
            indexOfFirstTPCCurl = i;
            break;
        }
    }

    for(int j=indexOfFirstTPCCurl+1; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
    return subTracks;
}

float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
    int algorithmID = pidHandler.getAlgorithmID(algorithmName);
    const EVENT::ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
    const std::vector<float>& parameters = pfoPID.getParameters();
    int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
    return parameters[parIdx]; 
}


EVENT::TrackerHit* getSETHit(EVENT::Track* track){
    std::vector<EVENT::TrackerHit*> hits = track->getTrackerHits();
    UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
    auto isSETHit = [&encoder](EVENT::TrackerHit* hit) -> bool {
        encoder.setValue( hit->getCellID0() ) ;
        int subdet = encoder[ UTIL::LCTrackerCellID::subdet() ];
        return subdet == UTIL::ILDDetID::SET;
    };
    auto it = std::find_if(hits.begin(), hits.end(), isSETHit);
    if ( it == hits.end() ) return nullptr;
    return *it;
}

const EVENT::TrackState* getTrackStateAtCalorimeter(EVENT::Track* track){
    return getSubTracks(track).back()->getTrackState( EVENT::TrackState::AtCalorimeter );
}

IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrack, EVENT::TrackerHit* hit){
    IMPL::TrackStateImpl ts;
    double chi2Dummy;
    int ndfDummy;
    marlinTrack->getTrackState(hit, ts, chi2Dummy, ndfDummy);
    return ts;
}


EVENT::MCParticle* getMC(EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& pfo2mc){
    const std::vector<EVENT::LCObject*>& objects = pfo2mc.getRelatedToObjects(pfo);
    const std::vector<float>& weights = pfo2mc.getRelatedToWeights(pfo);
    if ( objects.empty() ) return nullptr;

    auto getTrackWeight = [](float encodedWeight){ return float( int(encodedWeight) % 10000 ) / 1000.f;};
    auto getClusterWeight = [](float encodedWeight){ return float( int(encodedWeight) / 10000 ) / 1000.f;};

    int max_i = std::max_element(weights.begin(), weights.end(), [getTrackWeight](float lhs, float rhs){return getTrackWeight(lhs) < getTrackWeight(rhs);}) - weights.begin();
    if (getTrackWeight(max_i) == 0.f){
        max_i = std::max_element(weights.begin(), weights.end(), [getClusterWeight](float lhs, float rhs){return getClusterWeight(lhs) < getClusterWeight(rhs);}) - weights.begin();
    }
    return static_cast<EVENT::MCParticle*> (objects[max_i]);
}

double getECALBarelRMin(){
    double cm2mm = 10.;
    auto ecalBarrelData = MarlinUtil::getLayeredCalorimeterData(( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::BARREL),
                                                                ( DetType::AUXILIARY | DetType::FORWARD ) );

    return ecalBarrelData->extent[0]*cm2mm; // rmin in mm
}

double getECALEndcapZMin(){
    double cm2mm = 10.;
    auto ecalEndcapData = MarlinUtil::getLayeredCalorimeterData( ( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::ENDCAP),
                                                                 ( DetType::AUXILIARY | DetType::FORWARD ) );
    return ecalEndcapData->extent[2]*cm2mm; // zmin in mm
}

dd4hep::rec::Vector3D getPhotonAtCalorimeter(EVENT::MCParticle* mc){
    // find intersection point between photon momentum line and ECAL surface planes
    // https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection

    // startPos - starting position of the photon (l0)
    // mom - momentum of the photon (l)
    Vector3D startPos( mc->getVertex() );
    Vector3D mom( mc->getMomentum() );

    // rMin - minimal radial distance to the barrel ECAL surface to deduce normal vector n and point at the surface p0.
    // zMin - minimal longitudinal distance to the endcap ECAL surface to deduce normal vector n and point at the surface p0.
    double rMin = getECALBarelRMin();
    double zMin = getECALEndcapZMin();


    // ENDCAP plane parameters:
    int direction = (mom.z() > 0) ? 1 : -1;
    Vector3D p0Endcap(0, 0, direction*zMin);
    Vector3D nEndcap = p0Endcap.unit();

    // BARREL plane parameters:
    Vector3D p0Barrel, nBarrel;
    int nSides = 8;
    double step = M_PI/nSides;
    double phi = mom.phi();
    // phi is in the range of singularity point [pi, -pi]. Check this first.

    if( phi < (- M_PI + step) || phi > (M_PI - step) ){
        p0Barrel = Vector3D(rMin, M_PI, M_PI/2., Vector3D::spherical);
        nBarrel = p0Barrel.unit();
    }
    else{
        double ecalPhi = -M_PI + 2*step;
        for ( int i=0; i < nSides-1; ++i ){
            if ( ecalPhi-step <= phi && phi < ecalPhi+step ){
                p0Barrel = Vector3D(rMin, ecalPhi, M_PI/2., Vector3D::spherical);
                nBarrel = p0Barrel.unit();
                break;
            }
            else ecalPhi += 2*step;
        }
    }

    //find intersection point, but don't divide by zero
    if ( mom.z() == 0 ){
        double d = (p0Barrel - startPos).dot(nBarrel)/(mom.dot(nBarrel));
        return startPos + d*mom;
    }
    else if( mom.rho() == 0 ){
        double d = (p0Endcap - startPos).dot(nEndcap)/(mom.dot(nEndcap));
        return startPos + d*mom;
    }
    //choose closest intersection point to the 0,0,0
    double dBarrel = (p0Barrel - startPos).dot(nBarrel)/(mom.dot(nBarrel));
    Vector3D intersectionBarrel = startPos + dBarrel*mom;
    double dEndcap = (p0Endcap - startPos).dot(nEndcap)/(mom.dot(nEndcap));
    Vector3D intersectionEndcap = startPos + dEndcap*mom;
    if ( intersectionBarrel.r() <= intersectionEndcap.r() ) return intersectionBarrel;
    return intersectionEndcap;
};

EVENT::SimTrackerHit* getSimTrackerHit(EVENT::TrackerHit* hit, const UTIL::LCRelationNavigator& navToSimTrackerHits){
    // I merge all tracker hit relation collections in the steering file. ENSURE this happens!
    // Otherwise I need to check every possible tracker hit relation collection, which makes this code x10 longer.
    // In case collection doesn't exist, merging is still happens (I think..) with a warning, which is good.
    if (navToSimTrackerHits.getRelatedToObjects(hit).empty()) return nullptr;
    
    const std::vector<float>& weights = navToSimTrackerHits.getRelatedToWeights(hit);
    int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
    EVENT::SimTrackerHit* simHit = static_cast<EVENT::SimTrackerHit*> (navToSimTrackerHits.getRelatedToObjects(hit)[max_i]);
    return simHit;
}

unsigned long interpolateHexColor(unsigned long startColor, unsigned long endColor, float ratio) {
    //WARNING: chatGPT's work
    unsigned char startR = (startColor >> 16) & 0xFF;
    unsigned char startG = (startColor >> 8) & 0xFF;
    unsigned char startB = startColor & 0xFF;

    unsigned char endR = (endColor >> 16) & 0xFF;
    unsigned char endG = (endColor >> 8) & 0xFF;
    unsigned char endB = endColor & 0xFF;

    unsigned char interpolatedR = static_cast<unsigned char>(startR + ratio * (endR - startR));
    unsigned char interpolatedG = static_cast<unsigned char>(startG + ratio * (endG - startG));
    unsigned char interpolatedB = static_cast<unsigned char>(startB + ratio * (endB - startB));

    return (interpolatedR << 16) | (interpolatedG << 8) | interpolatedB;
}