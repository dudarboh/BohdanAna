import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
ROOT.EnableImplicitMT()

ROOT.gInterpreter.Declare('#include "tof.hpp"')

colors = [ROOT.TColor.GetColor('#1b9e77'), ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#e7298a')]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root", ["pdg", "tofClosest0", "layerHit", "recoCaloX", "recoCaloPx", "recoCaloY", "recoCaloPy", "recoCaloZ", "recoCaloPz", "xHit", "yHit", "zHit", "tHit", "cleanTrack", "cleanShower"] )
df = df.Filter("if (rdfentry_ % 1000000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")\
        .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
        .Define("nHitsIn10Layers", "getNHitsInLayers(layerHit, 10)")\
        .Filter("tofClosest0 > 6. && nHitsIn10Layers > 0")\
        .Define("rImpact", "ROOT::Math::XYZVector(recoCaloX, recoCaloY, recoCaloZ)")\
        .Define("momImpact", "ROOT::Math::XYZVector(recoCaloPx, recoCaloPy, recoCaloPz)")\
        .Define("hitPos", "hitPos(xHit, yHit, zHit)")\
        .Define("dToImpact", "dToImpact(hitPos, rImpact)")\
        .Define("dToLine", "dToLine(hitPos, rImpact, momImpact)")\
        .Define("trueTOF", "tofClosest(tHit, dToImpact)")\
        .Define("tHit50ps", "smear(rdfslot_, tHit, gaus50)")\
        .Define("tSurface50ps", "getTimeAtSurface(tHit50ps, dToImpact)")\
        .Define("tHit100ps", "smear(rdfslot_, tHit, gaus100)")\
        .Define("tSurface100ps", "getTimeAtSurface(tHit100ps, dToImpact)")\
        .Define("frankHits", "selectHits_OLD(dToLine, layerHit, true, 10, 999999.)")\
        .Define("frankTof50ps", "getAverage(tSurface50ps[frankHits])")\
        .Define("diff_frankTof50ps", "(frankTof50ps - trueTOF) * 1000")

names = ["diff_frankTof50ps"]

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
    h.SetMarkerStyle(0)
    h.GetXaxis().SetNdivisions(506)
    h.GetYaxis().SetMaxDigits(3)

for name,arr in arrays.items():
    mean90, rms90, mean90_err, rms90_err = fit90(arr)
    print(f"Method: {name}    RMS90: {rms90:.3f}    RMS: {np.std(arr):.3f}")

canvas.Update()
input("wait")