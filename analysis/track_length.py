import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

def draw_lines():
    lines = {}
    for p in particles:
        lines[p] = ROOT.TLine(0., p.mass2, 15., p.mass2)
        lines[p].SetLineColor(2)
        lines[p].SetLineWidth(1)
        lines[p].SetLineStyle(9)
        lines[p].Draw()
    return lines

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Define("mom2", "recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz")
df = df.Define("mom", "sqrt(mom2)")
trk_len_methods = ["trackLengthToEcal_SHA_phiLambda_IP", "trackLengthToEcal_SHA_phiZed_IP", "trackLengthToEcal_SHA_zedLambda_IP"]

canvases = []
histos = []
x_title, y_title = "Momentum (GeV/c)", "Mass^{2} (GeV^{2}/c^{4})"
for trk_len in trk_len_methods:
    df_new = df.Define("beta", f"{trk_len}/(tofClosest0*299.792458)")\
               .Define("mass2", "mom2*(1./(beta*beta) - 1.)")
    h = df_new.Histo2D( (f"h_{trk_len}", "", 1000, 0, 10, 1000, -0.3, 1.2), "mom","mass2" )
    print( h.GetEntries() )
    h.GetXaxis().SetTitle(x_title)
    h.GetYaxis().SetTitle(y_title)
    # print( h.Inspect())
    histos.append(h)
    canvases.append( draw_2d_plot(h) )
    input("wait")