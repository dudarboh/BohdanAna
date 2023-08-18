import ROOT
import numpy as np
ROOT.gStyle.SetOptTitle(1)
ROOT.EnableImplicitMT()

colors = ['#1b9e77', '#d95f02', '#7570b3']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

df = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/refit_default_allparticles.root")
df = df.Filter("mom_true.r() < 1")\
       .Define("diff_default", "omega_before - omega_true")\
       .Define("diff_refit", "omega_after - omega_true")\

h_pion_default = df.Filter("abs(pdg) == 211").Histo1D(("h_pion_default", "#pi; #Omega_{reco} - #Omega_{true} [1/mm]; N entries",300, -6e-5, 6e-5), "diff_default")
h_kaon_default = df.Filter("abs(pdg) == 321").Histo1D(("h_kaon_default", "K; #Omega_{reco} - #Omega_{true} [1/mm]; N entries",300, -6e-5, 6e-5), "diff_default")
h_proton_default = df.Filter("abs(pdg) == 2212").Histo1D(("h_proton_default", "p; #Omega_{reco} - #Omega_{true} [1/mm]; N entries",300, -6e-5, 6e-5), "diff_default")

h_pion_refit = df.Filter("abs(pdg) == 211").Histo1D(("h_pion_refit", "#pi; #Omega_{reco} - #Omega_{true} [1/mm]; N entries",300, -6e-5, 6e-5), "diff_refit")
h_kaon_refit = df.Filter("abs(pdg) == 321").Histo1D(("h_pion_refit", "K; #Omega_{reco} - #Omega_{true} [1/mm]; N entries",300, -6e-5, 6e-5), "diff_refit")
h_proton_refit = df.Filter("abs(pdg) == 2212").Histo1D(("h_pion_refit", "p; #Omega_{reco} - #Omega_{true} [1/mm]; N entries",300, -6e-5, 6e-5), "diff_refit")


ROOT.gStyle.SetPadLeftMargin(0.18)
ROOT.gStyle.SetPadRightMargin(0.1)

canvas = ROOT.TCanvas("omega_change_distribution",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )

canvas.SetGridx()
canvas.SetGridy()

h_pion_default.Scale( 1./h_pion_default.GetEntries() )
h_pion_default.Draw("histo ")
h_pion_default.SetLineColor(colors[0])
h_pion_default.SetLineWidth(3)

h_kaon_default.Scale( 1./h_kaon_default.GetEntries() )
h_kaon_default.Draw("histo same")
h_kaon_default.SetLineColor(colors[1])
h_kaon_default.SetLineWidth(3)

h_proton_default.Scale( 1./h_proton_default.GetEntries() )
h_proton_default.Draw("histo same")
h_proton_default.SetLineColor(colors[2])
h_proton_default.SetLineWidth(3)
canvas.BuildLegend()
canvas.Update()
input("wait")


h_pion_refit.Scale( 1./h_pion_refit.GetEntries() )
h_pion_refit.Draw("histo ")
h_pion_refit.SetLineColor(colors[0])
h_pion_refit.SetLineWidth(3)

h_kaon_refit.Scale( 1./h_kaon_refit.GetEntries() )
h_kaon_refit.Draw("histo same")
h_kaon_refit.SetLineColor(colors[1])
h_kaon_refit.SetLineWidth(3)

h_proton_refit.Scale( 1./h_proton_refit.GetEntries() )
h_proton_refit.Draw("histo same")
h_proton_refit.SetLineColor(colors[2])
h_proton_refit.SetLineWidth(3)
canvas.Update()
input("wait")

