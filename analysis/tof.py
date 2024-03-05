import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
ROOT.EnableImplicitMT()

def fit90(x): 
    x = np.sort(x)
    n10percent = int(round(len(x)*0.1))
    n90percent = len(x) - n10percent
    
    def interval_quantile_(x, quant=0.9):
        '''Calculate the shortest interval that contains the desired quantile'''
        # the number of possible starting points
        n_low = int(len(x) * (1. - quant))
        # the number of events contained in the quantil
        n_quant = len(x) - n_low

        # Calculate all distances in one go
        distances = x[-n_low:] - x[:n_low]
        i_start = np.argmin(distances)

        return i_start, i_start + n_quant

    start, end = interval_quantile_(x, quant=0.9)
    
    rms90 = np.std(x[start:end])
    mean90 = np.mean(x[start:end])
    mean90_err = rms90/np.sqrt(n90percent)
    rms90_err = rms90/np.sqrt(2*n90percent)   # estimator in root
    return mean90, rms90, mean90_err, rms90_err


code='''
#include "Math/Vector3D.h"
#include <random>
using namespace ROOT::Math;
using namespace ROOT::VecOps;
using std::cout, std::endl, std::stringstream, std::vector, std::string, std::runtime_error;

// in mm/ns
#define SPEED_OF_LIGHT 299.792458
// const float rInner = 1804.8;

// smearing generator and distributions
std::default_random_engine gen;
std::normal_distribution<float> gaus50(0., 0.05);
std::normal_distribution<float> gaus100(0., 0.1);

RVec <XYZVector> hitPos(const RVec<float>& x, const RVec<float>& y, const RVec<float>& z){
    auto constructHit = [](float x, float y, float z) { return XYZVector(x, y, z); };
    return Map(x, y, z, constructHit);
}


RVec <float> dToImpact(const RVec<XYZVector>& hits, const XYZVector& rImpact){
    auto getDistance = [&](const XYZVector& hit) { return (hit - rImpact).R(); };
    return Map(hits, getDistance);
}

RVec <float> getTimeAtSurface(const RVec<float>& tHit, const RVec<float>& dToImpact){
    auto correctForDistance = [](float tHit, float dToImpact) { return tHit - dToImpact/SPEED_OF_LIGHT; };
    return Map(tHit, dToImpact, correctForDistance);
}



RVec <float> dToLine(const RVec <XYZVector>& hits, const XYZVector& p0, const XYZVector& p){
    RVec <float> distance;
    for (const auto& hit:hits){
        float d = (hit - p0).Cross(p.Unit()).R();
        distance.push_back(d);
    }
    return distance;
}

RVec <float> smear(const RVec <float>& times, std::normal_distribution<float>& gaus){
    auto smearWithGaussian = [&](float time){ return time + gaus(gen); };
    return Map(times, smearWithGaussian);
}

RVec <bool> selectHits(const RVec<float>& dToLine, const RVec<int>& layer_hit, bool only_closest=true, int n_layers=10, float cyl_cut=9999.){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    if(only_closest){
        std::map<int, std::vector< std::pair<int, float> > > layer2hit;
        for (int i=0; i < nHits; ++i){
            if( dToLine[i] < cyl_cut ){
                layer2hit[ layer_hit[i] ].push_back( { i, dToLine[i] } );
            }
        }

        for (int layer=0; layer<n_layers; ++layer){
            if ( layer2hit[layer].empty() ) continue;
            int idx = (*std::min_element( layer2hit[layer].begin(), layer2hit[layer].end(), [](const auto& l, const auto& r) { return l.second < r.second; })).first;
            selected[idx] = true;
        }
    }
    else{
        for (int i=0; i < nHits; ++i){
            selected[i] = dToLine[i] < cyl_cut && layer_hit[i] < n_layers;
        }
    }
    return selected;
}

RVec <bool> selectCylinderHits(const RVec<float>& dToLine, const RVec<int>& layer_hit, float start_radii=6., int n_layers=10){
    // WARNING: ENSURE AT LEAST ONE HIT IS PRESENT WITHIN GIVEN n_layers! Otherwise infinite loop
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);
    if (nHits == 0) return selected;

    float min_hit_radii = *std::min_element(dToLine.begin(), dToLine.end());

    float r = max(start_radii, min_hit_radii);
    bool foundHit = false;
    while (not foundHit){
        for (int i=0; i < nHits; ++i){
            if ( dToLine[i] < r && layer_hit[i] < n_layers){
                selected[i] = true;
                foundHit = true;
            }
        }
        if (foundHit) return selected;
        r += 0.1;
    }
    return selected;
}


//Select hits, if their TimeAtSurface is within +/- sigma of the time resolution.
RVec <bool> selectReasonableHits(const RVec<float> timeAtSurface, const RVec<float>& dToImpact, const RVec<float>& dToLine, const RVec<int>& layerHit, float timeResolution, float nSigma=3, int nLayers=10){
    int nHits = timeAtSurface.size();

    //use the closest hit as the initial reference time
    int closestHitIdx = std::min_element(dToImpact.begin(), dToImpact.end()) - dToImpact.begin();
    float referenceTime = timeAtSurface[closestHitIdx];

    //Find the list of all "good" indicies/hits
    while(true){
        RVec <bool> selected(nHits);
        RVec <float> selectedTimes;
        for (int i=0; i < nHits; ++i){
            if ( layerHit[i] >= nLayers || dToLine[i] > 7.) continue;
            float timeDiffInPicoSeconds = (timeAtSurface[i] - referenceTime)*1000.;
            if ( std::abs(timeDiffInPicoSeconds) < nSigma*timeResolution ){
                selectedTimes.push_back( timeAtSurface[i] );
                selected[i] = true;
            }
        }
        if( selectedTimes.empty() ) return selected;
        auto const count = static_cast<float>( selectedTimes.size() );
        float tmpReferenceTime = std::reduce(selectedTimes.begin(), selectedTimes.end()) / count;
        if (tmpReferenceTime == referenceTime) return selected;
        else referenceTime = tmpReferenceTime;
    }
}

/////////////////////////////////TOFS////////////////////////////////

float tofClosest(const RVec<float>& tHit, const RVec<float>& dToImpact){
    int min_idx = std::min_element(dToImpact.begin(), dToImpact.end()) - dToImpact.begin();
    return tHit[min_idx] - dToImpact[min_idx]/SPEED_OF_LIGHT;
}

float getAverage(const RVec<float>& v){
    if( v.empty() ) return 0;
    auto const count = static_cast<float>( v.size() );
    return std::reduce(v.begin(), v.end()) / count;
}

int getNHitsInLayers(const RVec<int>& layer_hit, int n_layers=10){
    return std::count_if(layer_hit.begin(), layer_hit.end(), [&](int layer){return layer < n_layers;});
}


'''
ROOT.gInterpreter.Declare(code)


colors = [ROOT.TColor.GetColor('#1b9e77'), ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#e7298a')]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root", ["pdg", "tofClosest0", "layerHit", "recoCaloX", "recoCaloPx", "recoCaloY", "recoCaloPy", "recoCaloZ", "recoCaloPz", "xHit", "yHit", "zHit", "tHit"] )
df = df.Filter("if (rdfentry_ % 1000000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")\
        .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
        .Filter("tofClosest0 > 6. && getNHitsInLayers(layerHit, 10) > 0")\
        .Define("rImpact", "ROOT::Math::XYZVector(recoCaloX, recoCaloY, recoCaloZ)")\
        .Define("momImpact", "ROOT::Math::XYZVector(recoCaloPx, recoCaloPy, recoCaloPz)")\
        .Define("hitPos", "hitPos(xHit, yHit, zHit)")\
        .Define("dToImpact", "dToImpact(hitPos, rImpact)")\
        .Define("dToLine", "dToLine(hitPos, rImpact, momImpact)")\
        .Define("trueTOF", "tofClosest(tHit, dToImpact)")\
        .Define("tHit50ps", "smear(tHit, gaus50)")\
        .Define("tSurface50ps", "getTimeAtSurface(tHit50ps, dToImpact)")\
        .Define("tHit100ps", "smear(tHit, gaus100)")\
        .Define("tSurface100ps", "getTimeAtSurface(tHit100ps, dToImpact)")\
        .Define("frankHits", "selectHits(dToLine, layerHit, true, 10, 999999.)")\
        .Define("frankTof50ps", "getAverage(tSurface50ps[frankHits])")\
        .Define("diff_frankTof50ps", "(frankTof50ps - trueTOF) * 1000")

names = ["diff_frankTof50ps"]
for r in [6.]:
    df = df.Define(f"cylinder{int(r*10)}Hits", f"selectCylinderHits(dToLine, layerHit, {r}, 10)")\
            .Define(f"cyl{int(r*10)}Tof50ps", f"getAverage(tSurface50ps[cylinder{int(r*10)}Hits])")\
            .Define(f"diff_cyl{int(r*10)}Tof50ps", f"(cyl{int(r*10)}Tof50ps - trueTOF) * 1000")
    names.append(f"diff_cyl{int(r*10)}Tof50ps")

histos= []
for name in names:
    h = df.Histo1D((f"{name}", ";#Delta T (ps);N entries", 2000, -200, 200), f"{name}")
    histos.append(h)

# for i, nSigma in enumerate([2.5]):
#     df = df.Define(f"HitSelection50ps{i}", f"selectReasonableHits(tSurface50ps, dToImpact, dToLine, layerHit, 50, {nSigma}, 10)")\
#            .Define(f"Tof50ps{i}", f"getAverage(tSurface50ps[HitSelection50ps{i}])")
#     names.append(f"Tof50ps{i}")


arrays = df.AsNumpy(columns=names)

canvas = create_canvas()

for i, h in enumerate(histos):
    h.Draw("" if i == 0 else "same")
    h.SetLineColor(i+1)

for name,arr in arrays.items():
    mean90, rms90, mean90_err, rms90_err = fit90(arr)
    print(f"Method: {name}    RMS90: {rms90:.3f}    RMS: {np.std(arr):.3f}")

canvas.Update()
input("wait")