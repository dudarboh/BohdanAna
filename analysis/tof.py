import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
ROOT.EnableImplicitMT()


code='''
#include "Math/Vector3D.h"
#include <random>
using namespace ROOT::Math;
using namespace ROOT::VecOps;
using std::cout, std::endl, std::stringstream, std::vector, std::string, std::runtime_error;

// in mm/ns
#define SPEED_OF_LIGHT 299.792458
// const double rInner = 1804.8;

// smearing generator and distributions
std::default_random_engine gen;
std::normal_distribution gaus100(0., 0.1);

RVec <XYZVector> hitPos(const RVec<double>& x, const RVec<double>& y, const RVec<double>& z){
    auto constructHit = [](double x, double y, double z) { return XYZVector(x, y, z); };
    return Map(x, y, z, constructHit);
}


RVec <double> dToImpact(const RVec<XYZVector>& hits, const XYZVector& rImpact){
    auto getDistance = [&](const XYZVector& hit) { return (hit - rImpact).R(); };
    return Map(hits, getDistance);
}


RVec <double> dToLine(const RVec <XYZVector>& hits, const XYZVector& p0, const XYZVector& p){
    RVec <double> distance;
    for (const auto& hit:hits){
        double d = (hit - p0).Cross(p.Unit()).R();
        distance.push_back(d);
    }
    return distance;
}

RVec <bool> selectHits(const RVec<double>& dToLine, const RVec<int>& layer_hit, bool only_closest=true, int n_layers=10, double cyl_cut=9999.){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    if(only_closest){
        std::map<int, std::vector< std::pair<int, double> > > layer2hit;
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
        for (int i=0; i < nHits; ++i) selected[i] = dToLine[i] < cyl_cut && layer_hit[i] < n_layers;
    }
    return selected;
}

/////////////////////////////////TOFS////////////////////////////////
RVec <double> smearTime(const RVec <double>& times){
    auto smear = [&](double time){ return time + gaus100(gen); };
    return Map(times, smear);
}

double tofClosest(const RVec<double>& tHit, const RVec<double>& dToImpact){
    int min_idx = std::min_element(dToImpact.begin(), dToImpact.end()) - dToImpact.begin();
    return tHit[min_idx] - dToImpact[min_idx]/SPEED_OF_LIGHT;
}

double tofAvg(const RVec <double>& tofHit, const RVec<double>& dToImpact){
    int nHits = tofHit.size();
    double tofSum = 0.;
    for(int i=0; i < nHits; ++i) tofSum += tofHit[i] - dToImpact[i]/SPEED_OF_LIGHT;
    return tofSum/nHits;
}

'''
ROOT.gInterpreter.Declare(code)


colors = [ROOT.TColor.GetColor('#1b9e77'), ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#e7298a')]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("if (rdfentry_ % 1000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")
df = df.Filter("tofClosest0 > 6. && nHits > 0").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Define("rImpact", "ROOT::Math::XYZVector(recoCaloX, recoCaloY, recoCaloZ)")\
       .Define("momImpact", "ROOT::Math::XYZVector(recoCaloPx, recoCaloPy, recoCaloPz)")\
       .Define("hitPos", "hitPos(xHit, yHit, zHit)")\
       .Define("dToImpact", "dToImpact(hitPos, rImpact)")\
       .Define("dToLine", "dToLine(hitPos, rImpact, momImpact)")\
       .Define("tHitSmeared", "smearTime(tHit)")\
       .Define("trueTOF", "tofClosest(tHit, dToImpact)")

df = df.Define("defaultSelect", "selectHits(dToLine, layerHit, true, 10, 9999.)")\
       .Define("tof_default", "tofAvg(tHitSmeared[defaultSelect], dToImpact[defaultSelect])")\
       .Define("select_cut5", "selectHits(dToLine, layerHit, true, 10, 5.)")\
       .Define("tof_cut5", "tofAvg(tHitSmeared[select_cut5], dToImpact[select_cut5])")\
       .Define("select_cut10", "selectHits(dToLine, layerHit, true, 10, 10.)")\
       .Define("tof_cut10", "tofAvg(tHitSmeared[select_cut10], dToImpact[select_cut10])")\
       .Define("select_cut15", "selectHits(dToLine, layerHit, true, 10, 15.)")\
       .Define("tof_cut15", "tofAvg(tHitSmeared[select_cut15], dToImpact[select_cut15])")\
       .Define("select_cyl_all", "selectHits(dToLine, layerHit, false, 10, 9999.)")\
       .Define("tof_cyl", "tofAvg(tHitSmeared[select_cyl_all], dToImpact[select_cyl_all])")\
       .Define("select_cyl_cut5", "selectHits(dToLine, layerHit, false, 10, 5.)")\
       .Define("tof_cyl_cut5", "tofAvg(tHitSmeared[select_cyl_cut5], dToImpact[select_cyl_cut5])")\
       .Define("select_cyl_cut10", "selectHits(dToLine, layerHit, false, 10, 10.)")\
       .Define("tof_cyl_cut10", "tofAvg(tHitSmeared[select_cyl_cut10], dToImpact[select_cyl_cut10])")\
       .Define("select_cyl_cut15", "selectHits(dToLine, layerHit, false, 10, 15.)")\
       .Define("tof_cyl_cut15", "tofAvg(tHitSmeared[select_cyl_cut15], dToImpact[select_cyl_cut15])")

names = ["tof_default", "tof_cut5", "tof_cut10", "tof_cut15", "tof_cyl", "tof_cyl_cut5", "tof_cyl_cut10", "tof_cyl_cut15"]
histos= []

for n in names:
    h = df.Define("diff", f"({n} - trueTOF) * 1000").Histo1D((f"h_{n}", ";#Delta T (ps);N entries", 300, -100, 100), "diff")
    histos.append(h)

hs = ROOT.THStack()
for i, h in enumerate(histos):
    hs.Add(h.GetPtr())
    h.SetLineColor(i+1)

hs.Draw("nostack")
input("wait")