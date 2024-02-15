import ROOT
import numpy as np
from utils import *
ROOT.EnableImplicitMT()

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
                .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
                .Filter("tofClosest0 > 6.")
df = df.Define("pzvspt", "abs(mcPz)/sqrt(mcPx*mcPx + mcPy*mcPy)")

canvas = create_canvas()

h = df.Histo1D(("h", "",1000, 0, 10), "pzvspt")
h.Draw()

print("Total particles: ", df.Count().GetValue())
print("Super perp particles: ", df.Filter("pzvspt < 0.0055431805").Count().GetValue())
print(" % ", df.Filter("pzvspt < 0.0055431805").Count().GetValue()/df.Count().GetValue() * 100)

canvas.Update()
input("wait")

df_pion = df.Filter("abs(pdg) == 211")
df_kaon = df.Filter("abs(pdg) == 321")
df_proton = df.Filter("abs(pdg) == 2212")

# # MOMENTUM PLOT
h_mom_pions = df_pion.Histo1D(("h_mom_pions", "Pions; momentum (GeV); N particles",1000, 0, 1), "pzvspt")
h_mom_kaons = df_kaon.Histo1D(("h_mom_kaons", "Kaons; momentum (GeV); N particles",1000, 0, 1), "pzvspt")
h_mom_protons = df_proton.Histo1D(("h_mom_protons", "Protons; momentum (GeV); N particles",1000, 0, 1), "pzvspt")
h_mom_pions.Scale( 1. / h_mom_pions.GetEntries() )
h_mom_kaons.Scale( 1. / h_mom_kaons.GetEntries() )
h_mom_protons.Scale( 1. / h_mom_protons.GetEntries() )

h_mom_pions.Draw("histo ")
h_mom_pions.SetLineColor(colors[0])
h_mom_pions.SetLineWidth(5)
h_mom_protons.Draw("histo same")
h_mom_protons.SetLineColor(colors[2])
h_mom_protons.SetLineWidth(5)
h_mom_kaons.Draw("histo same")
h_mom_kaons.SetLineColor(colors[1])
h_mom_kaons.SetLineWidth(5)

print("Bin number for 3 GeV: ", h_mom_pions.FindBin(3.))
print("Pion integral below 3 GeV: ", h_mom_pions.Integral(0, h_mom_pions.FindBin(3.)))
print("Kaon integral below 3 GeV: ", h_mom_kaons.Integral(0, h_mom_kaons.FindBin(3.)))
print("Proton integral below 3 GeV: ", h_mom_protons.Integral(0, h_mom_protons.FindBin(3.)))
canvas.Update()
input("wait")
# # END OF MOMENTUM PLOT



