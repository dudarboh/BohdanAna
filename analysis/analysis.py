import ROOT
import numpy as np
ROOT.EnableImplicitMT()

class Algorithm:
    def __init__(self, name, momentum, track_length, tof):
        self.name = name
        self.mom = momentum
        self.len = track_length
        self.tof = tof

colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]

algoClosest0 = Algorithm("Closest 0 ps", "harmonicMomToEcal", "trackLengthToEcal", "tofClosest0")
algoClosest30 = Algorithm("Closest 30 ps", "harmonicMomToEcal", "trackLengthToEcal", "tofClosest30")
algoSET50 = Algorithm("SET 50 ps", "harmonicMomToSET", "trackLengthToSET", "tofSET50")
algoAverage100 = Algorithm("Average 100 ps", "harmonicMomToEcal", "trackLengthToEcal", "tofAverage100")
# algorithms = [algoClosest0, algoClosest30, algoSET50, algoAverage100]
algorithms = [algoClosest30, algoSET50, algoAverage100]


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


df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")

histos1d = []
histos2d = []
for alg in algorithms:
    df_beta = df.Define("beta", f"{alg.len}/({alg.tof}*299.792458)").Filter("beta >= 0 && beta <= 1")
    df_mass = df_beta.Define("mass", f"{alg.mom}*sqrt( 1./(beta*beta) - 1.)*1000")

    h1 = df_mass.Histo1D((f"{alg.name}_1d", f"{alg.name}; mass [MeV]; N entries", 2000, 0, 1100), "mass")
    histos1d.append(h1)

    h1 = df_mass.Histo1D((f"{alg.name}_1d", f"{alg.name}; mass [MeV]; N entries", 2000, 0, 1100), "mass")
    histos2d.append(h2)

    h3 = df_mass.Histo2D((f"{alg.name}_2d", f"{alg.name}; momentum [GeV]; Mass [MeV]", 500, 0, 15, 500, -100, 1300.), f"{alg.mom}","mass")


# h.Draw("colz")
# ROOT.gStyle.SetPalette(ROOT.kBird)
# lines = draw_lines()
# h.SetMinimum(0.1)
# h.SetMaximum(10000)
# canvas.SetLogz()
# canvas.SetRightMargin(0.16)
# canvas.Update()
# # # h.SetStats(0)
# # stats = h.FindObject("stats")
# # stats.SetOptStat(110010)
# # stats.SetX1NDC(0.74)
# # stats.SetX2NDC(0.92)
# # stats.SetY1NDC(0.89)
# # stats.SetY2NDC(0.99)
# palette = h.GetListOfFunctions().FindObject("palette")
# palette.SetX1NDC(0.86)
# palette.SetX2NDC(0.9)
# canvas.Modified()
# canvas.Update()
# input("wait")



def plot_ptpz():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/trk_len.root")
    df = df.Filter("massDefault != 0 && massTanL != 0 && massZ != 0 && abs(pdg) == 321")
    df = df.Define("pt", "hypot(mom_ip._x,mom_ip._y)")
    # Make binning consistent with old plots
    # df.AsNumpy(["pt", "mom_ip._z", "m"])

    h_default = df.Define("m", "massDefault*1000 - 0.493677*1000").Profile2D(("h_default", "Default; pt [GeV]; pz [GeV]", 600, 0, 15, 600, -15, 15, -300, 300), "pt","mom_ip._z", "m")
    h_tanl = df.Define("m", "massTanL*1000 - 0.493677*1000").Profile2D(("h_tanL", "TanL; pt [GeV]; pz [GeV]", 600, 0, 15, 600, -15, 15, -300, 300), "pt","mom_ip._z", "m")
    h_z = df.Define("m", "massZ*1000 - 0.493677*1000").Profile2D(("h_z", "Z; pt [GeV]; pz [GeV]", 600, 0, 15, 600, -15, 15, -300, 300), "pt","mom_ip._z", "m")


    canvas = ROOT.TCanvas()

    h_default.Draw("colz")
    h_default.SetStats(0)
    h_default.SetMinimum(-200)
    h_default.SetMaximum(200)
    canvas.Update()
    input("wait")

    canvas = ROOT.TCanvas()
    h_tanl.Draw("colz")
    h_tanl.SetStats(0)
    h_tanl.SetMinimum(-200)
    h_tanl.SetMaximum(200)
    canvas.Update()
    input("wait")

    canvas = ROOT.TCanvas()
    h_z.Draw("colz")
    h_z.SetStats(0)
    h_z.SetMinimum(-200)
    h_z.SetMaximum(200)
    canvas.Update()
    input("wait")

