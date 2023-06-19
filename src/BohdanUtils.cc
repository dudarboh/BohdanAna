#include "BohdanUtils.h"
#include "UTIL/ILDConf.h"
#include "marlinutil/GeometryUtil.h"

#include <cstring>
#include <iomanip>

using dd4hep::rec::LayeredCalorimeterData;
using dd4hep::DetType;

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


EVENT::MCParticle* getMC(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator pfo2mc){
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