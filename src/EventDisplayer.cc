#include "EventDisplayer.h"
#include "marlin/Processor.h"
#include "marlinutil/DDMarlinCED.h"

void EventDisplayer::initDisplay(marlin::Processor* proc){
    if (_eventDisplay && !_isInit) {
        DDMarlinCED::init(proc);
        _isInit = true;
    }
}
