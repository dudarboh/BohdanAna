#include "BohdanUtils.h"
#include "UTIL/ILDConf.h"

#include <cstring>
#include <iomanip>

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
    if ( it != hits.end() ) return *it;
    return nullptr;
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


void printTrackStateLong(IMPL::TrackStateImpl ts){
    std::stringstream tmp;
    //out << std::scientific << std::setprecision (2) << std::showpos;
    std::cout << std::noshowpos;
    std::cout << std::setw(41) << std::setfill('-') << std::right << "-- TrackState ---" << std::setfill('-') << std::setw(29) << "-" << std::endl;

    tmp.str("") ;
    switch( ts.getLocation() ){
    case EVENT::TrackState::AtOther         :   tmp <<  "AtOther"        ;   break ;
    case EVENT::TrackState::AtIP            :   tmp <<  "AtIP"           ;   break ;
    case EVENT::TrackState::AtFirstHit      :   tmp <<  "AtFirstHit"     ;   break ;
    case EVENT::TrackState::AtLastHit       :   tmp <<  "AtLastHit"      ;   break ;
    case EVENT::TrackState::AtCalorimeter   :   tmp <<  "AtCalorimeter " ;   break ;
    case EVENT::TrackState::AtVertex        :   tmp <<  "AtVertex"       ;   break ;
    }
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Location" << std::right << std::setw(40) << tmp.str() << std::endl;
    tmp.str("") ;
    tmp << std::dec << std::setfill('0') << std::setw(8) << ts.id();
    std::cout << std::scientific << std::setprecision(6) ;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Id"          << std::right << std::setw(40) << tmp.str() << std::endl;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "D0"          << std::right << std::setw(40) << ts.getD0() << std::endl;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Phi"         << std::right << std::setw(40) << ts.getPhi() << std::endl;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Omega"       << std::right << std::setw(40) << ts.getOmega() << std::endl;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Z0"          << std::right << std::setw(40) << ts.getZ0() << std::endl;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Tan Lambda"  << std::right << std::setw(40) << ts.getTanLambda() << std::endl;
    tmp.str("");
    tmp  << std::dec << ts.getReferencePoint()[0] << ", " << ts.getReferencePoint()[1]  << ", " << ts.getReferencePoint()[2]; 
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "ReferencePoint" << std::right << std::setw(40) << tmp.str() << std::endl;
    std::cout << "Cov matrix:" << std::showpos << std::scientific << std::setprecision(6) << std::setw(15) << std::setfill(' ')  ;
    // print cov matrix as lower triangle matrix 
    for( unsigned l=0 , N=ts.getCovMatrix().size(), ncolumns = 1 , nele =1 ; l <N ; ++l , ++nele) {
    std::cout << ts.getCovMatrix()[l];
    if(! ( (nele) % ncolumns ) ){ 
        nele = 0 ;
        ++ncolumns ;
        std::cout << std::endl << "             " ;
    } else {
        std::cout << ", ";
    } 
    }
    std::cout << std::noshowpos;
    std::cout<<std::endl;
}


void printTrackStateShort(IMPL::TrackStateImpl ts){
    std::cout<<std::endl;
    std::cout<<"Location:    "<<ts.getLocation()<<std::endl;
    std::cout<<"D0 Z0:    "<<ts.getD0()<<"  "<<ts.getZ0()<<std::endl;
    std::cout<<"Phi:    "<<ts.getPhi()<<std::endl;
    std::cout<<"Omega:    "<<ts.getOmega()<<std::endl;
    std::cout<<"TanLambda:    "<<ts.getTanLambda()<<std::endl;
    std::cout<<"Ref:    "<<ts.getReferencePoint()[0]<<"  "<<ts.getReferencePoint()[1]<<"  "<<ts.getReferencePoint()[2]<<std::endl;
}
