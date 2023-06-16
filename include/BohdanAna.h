#ifndef BohdanAna_h
#define BohdanAna_h 1

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "EventDisplayer.h"

class BohdanAna : public marlin::Processor, EventDisplayer {
    friend class EventDisplayer;
    public:
        BohdanAna(const BohdanAna&) = delete;
        BohdanAna& operator=(const BohdanAna&) = delete;

        marlin::Processor* newProcessor() { return new BohdanAna; }

        BohdanAna();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();
    private:
        int _nEvent{};
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        int _pdg;
        double _momentumIP;
        double _momentumCalo;
        double _momentumHM;

        std::pair<int, double> _tofClosest;

        double _trackLengthSHA1;
        double _trackLengthSHA2;
        double _trackLengthSHA3;
        double _trackLengthSHA4;
        double _trackLengthSHA5;
        double _trackLengthSHA6;

        double _trackLengthIKF1;
        double _trackLengthIKF2Bug;
        double _trackLengthIKF2;
        double _trackLengthIKF3;

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
};

#endif
