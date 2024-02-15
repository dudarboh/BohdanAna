import ROOT
import numpy as np
from utils import *
ROOT.EnableImplicitMT()

particles = create_list_of_particles(1., 1) # just for colours
(pion, kaon, proton) = particles

VAR_NAME = "tanL"
if VAR_NAME == "omega":
    x_title, y_title = "#Omega_{reco} - #Omega_{true} (1/mm)", "Normalised N entries"
    n_bins, x_min, x_max = 100, -20e-6, 20e-6
    x_text_pos, y_text_pos = 2.58e-6, 13.34e-3
elif VAR_NAME == "tanL":
    x_title, y_title = "tan(#lambda)_{reco} - tan(#lambda)_{true}", "Normalised N entries"
    n_bins, x_min, x_max = 300, -10e-3, 10e-3
    x_text_pos, y_text_pos = 3e-3, 10e-3



df = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/refit_default_allparticles.root")
df = df.Filter("mom_true.r() < 1")\
       .Define("diff_default", f"{VAR_NAME}_before - {VAR_NAME}_true")\
       .Define("diff_refit", f"{VAR_NAME}_after - {VAR_NAME}_true")\

h_pion_default = df.Filter("abs(pdg) == 211").Histo1D(("h_pion_default", pion.legend, n_bins, x_min, x_max), "diff_default")
h_kaon_default = df.Filter("abs(pdg) == 321").Histo1D(("h_kaon_default", kaon.legend, n_bins, x_min, x_max), "diff_default")
h_proton_default = df.Filter("abs(pdg) == 2212").Histo1D(("h_proton_default", proton.legend, n_bins, x_min, x_max), "diff_default")

h_pion_refit = df.Filter("abs(pdg) == 211").Histo1D(("h_pion_refit", pion.legend, n_bins, x_min, x_max), "diff_refit")
h_kaon_refit = df.Filter("abs(pdg) == 321").Histo1D(("h_pion_refit", kaon.legend, n_bins, x_min, x_max), "diff_refit")
h_proton_refit = df.Filter("abs(pdg) == 2212").Histo1D(("h_pion_refit", proton.legend, n_bins, x_min, x_max), "diff_refit")

margin = 0.3
ROOT.gStyle.SetPadLeftMargin(0.55*margin)
ROOT.gStyle.SetPadRightMargin(0.45*margin)
ROOT.gStyle.SetPadTopMargin(0.3*margin)
ROOT.gStyle.SetPadBottomMargin(0.7*margin)
canvas = ROOT.TCanvas(get_rand_string(),"", 600, 600)
legend =create_legend(0.2, 0.65, 0.62, 0.87)

histos = [h_pion_default, h_kaon_default, h_proton_default]

for p, h in zip(particles, histos):
    h.Scale(1./h.GetEntries())
    h.SetLineColor(p.color)
    h.SetLineWidth(3)
    legend.AddEntry(p.legend_graph, p.legend, "f")
    if p.name == "pion":
        h.Draw("histo")
        h.GetXaxis().SetTitle(x_title)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetXaxis().SetMaxDigits(3)
        h.GetXaxis().SetNdivisions(506)
        h.GetYaxis().SetTitle(y_title)
        h.GetYaxis().SetMaxDigits(3)
    else:
        h.Draw("histo same")

latex.DrawLatex(x_text_pos, y_text_pos, "p < 1 GeV/c")

legend.Draw()
canvas.Update()
input("wait")