#ifndef BohdanDrawing_h
#define BohdanDrawing_h 1

#include "BohdanUtils.h"

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"

void drawPFO(EVENT::ReconstructedParticle* pfo, IMPL::TrackStateImpl tsStdReco, IMPL::TrackStateImpl tsEasy);
void drawCanvas();


#endif