import ROOT
import numpy as np
# ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
# ROOT.EnableImplicitMT()


colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]

def draw_lines():
    lines = {}
    pdgs = [211, 321, 2212]
    m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
    for pdg in pdgs:
        lines[pdg] = ROOT.TLine(0., m_pdg[pdg], 15., m_pdg[pdg])
        lines[pdg].SetLineColor(2)
        lines[pdg].SetLineWidth(1)
        lines[pdg].SetLineStyle(9)
        lines[pdg].Draw()
    return lines

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_lengths.root")
df = df.Filter("tofv2 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")

# algorithms = ["trackLengthIDR", "trackLengthIDR2", "trackLengthIDR3",
            #   "trackLengthIDR4", "trackLengthWinni", "trackLengthWinni2",
            #   "trackLengthUsingZ", "trackLengthUsingZ2", "trackLengthUsingZ3", "trackLengthSimUsingZ"]
# algorithms = ["trackLengthUsingZ2", "trackLengthUsingZ3"]
# algorithms = ["trackLengthSHA1" ,"trackLengthSHA2" ,"trackLengthSHA3" ,"trackLengthSHA4" ,"trackLengthSHA5" ,"trackLengthSHA6", "trackLengthIKF1" ,"trackLengthIKF2" ,"trackLengthIKF3"]
algorithms = ["trackLengthIKF2"]


#PLOT 2D
histos = []
for name in algorithms:
    df_beta = df.Define("beta", f"{name}/(tofv2*299.792458)").Filter("beta >= 0 && beta <= 1")
    df_mass = df_beta.Define("mass", "momentum*sqrt( 1./(beta*beta) - 1.)*1000")
    h = df_mass.Histo2D((f"h_{name}", f"{name}; momentum [GeV]; Mass [MeV]", 500, 0, 15, 500, -100, 1300.), "momentum","mass")
    histos.append(h)

for name, h in zip(algorithms, histos):
    canvas = ROOT.TCanvas(f"{name}")
    h.Draw("colz")
    lines = draw_lines()
    h.SetMinimum(0.1)
    h.SetMaximum(3000)
    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.Update()
    # h.SetStats(0)
    stats = h.FindObject("stats")
    stats.SetOptStat(110010)
    stats.SetX1NDC(0.74)
    stats.SetX2NDC(0.92)
    stats.SetY1NDC(0.89)
    stats.SetY2NDC(0.99)
    palette = h.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.93)
    palette.SetX2NDC(0.95)
    canvas.Modified()
    canvas.Update()
    input("wait")
