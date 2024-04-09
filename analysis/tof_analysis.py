import ROOT
import numpy as np
from sklearn import linear_model
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
ROOT.EnableImplicitMT()

ROOT.gInterpreter.Declare('#include "tof.hpp"')
df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/tof_studies.root")\
    .Filter("if (rdfentry_ % 1000000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")\
    .Filter("sqrt(recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy + recoCaloPz*recoCaloPz) < 10.")\
    .Filter("rdfentry_ < 100000")

def filter_n_layers(df, n_layers):
    '''Return df with hits only within first n_layers'''
    df = df.Redefine("dl", f"dl[layerHit<{n_layers}]")\
            .Redefine("dToImpact", f"dToImpact[layerHit<{n_layers}]")\
            .Redefine("dToLine", f"dToLine[layerHit<{n_layers}]")\
            .Redefine("tHit", f"tHit[layerHit<{n_layers}]")\
            .Redefine("tSurface", f"tSurface[layerHit<{n_layers}]")\
            .Redefine("layerHit", f"layerHit[layerHit<{n_layers}]")
    return df

def smear_time(df):
    '''Return df with time smeared with time_resolution (in ps)'''
    # Only with 50 ps for now
    df = df.Redefine("tHit", "smear(rdfslot_, tHit, gaus50)")\
            .Redefine("tSurface", "getTimeAtSurface(tHit, dToImpact)")
    return df

def scan_cuts(df):
    '''Return numpy array of the dt results using scan cuts'''
    d_perp_cuts = np.linspace(7, 11, 20)
    dt_cuts = np.linspace(120, 300, 20)
    results = []
    for i, d_perp_cut in enumerate(d_perp_cuts):
        for j, dt_cut in enumerate(dt_cuts):
            results.append(f"dt_{i}_{j}")
            df = df.Define(f"tSurface_cyl_mask_{i}_{j}", f"tSurface[selectCylinderHits(dToLine, {d_perp_cut})]")\
            .Define(f"tSurface_both_masks_{i}_{j}", f"tSurface_cyl_mask_{i}_{j}[selectMedianHits(tSurface_cyl_mask_{i}_{j}, {dt_cut})]")\
            .Define(f"dt_{i}_{j}", f"1000*(Mean(tSurface_both_masks_{i}_{j}) - tofClosest0)")
    return d_perp_cuts, dt_cuts, df.AsNumpy(results)

df = filter_n_layers(df, 10)
df = smear_time(df)
print("Scanning...")
d_perp_cuts, dt_cuts, results = scan_cuts(df)

print("Plotting...")
gr = ROOT.TGraph2D()
for i, d_perp_cut in enumerate(d_perp_cuts):
    for j, dt_cut in enumerate(dt_cuts):
        _, rms90, _, _ = fit90( results[f"dt_{i}_{j}"] )
        gr.SetPoint( i*len(d_perp_cuts) + j, d_perp_cut, dt_cut, rms90 )
gr.Draw("COLZ")

input("wait")