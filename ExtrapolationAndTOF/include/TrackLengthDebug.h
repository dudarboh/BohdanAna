#ifndef TrackLengthDebug_h
#define TrackLengthDebug_h 1

#include "EventDisplayer.h"
#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "EVENT/Track.h"
#include "IMPL/TrackStateImpl.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "DD4hep/Detector.h"
#include "EVENT/SimTrackerHit.h"
#include <string>
#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "DDRec/Vector3D.h"
#include "EVENT/Track.h"
#include <memory>
#include "Math/Vector3D.h"

struct Vars{
    double trackLength = 0.;
    double momentumHM = 0.;
    int nBadPhi = 0;
    int nBadZ = 0;
    int nSegments = 0;
    IMPL::TrackStateImpl tsCalo;
};

class TrackLengthDebug : public marlin::Processor, EventDisplayer{
    friend class EventDisplayer;
    public:
        TrackLengthDebug(const TrackLengthDebug&) = delete;
        TrackLengthDebug& operator=(const TrackLengthDebug&) = delete;

        marlin::Processor* newProcessor() { return new TrackLengthDebug; }

        TrackLengthDebug();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();
        Vars getTrackLengthOld(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);
        Vars getTrackLengthNew(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);

    private:
        int _nEvent{};
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        int _pdg;
        bool _isGoodTrack;
        int _isGoodClosestHitOld, _isGoodClosestHitNew;
        double _tofOld, _tofNew;
        double _momentumOld;
        double _momentumNew;
        double _trackLengthOld;
        double _trackLengthNew;
        int _nBadPhi;
        int _nBadZ;
        int _nSegments;
        int _nHitsFirstCurl;
        dd4hep::rec::Vector3D _tsLastCaloPos, _tsFirstCaloPos, _hitClosestOldPos, _hitClosestNewPos;
        double _hitTimeOld, _hitTimeNew;
        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
};

#endif
