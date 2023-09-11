import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
# ROOT.EnableImplicitMT()


# colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]
colors = [ROOT.TColor.GetColor('#1b9e77'), ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#e7298a')]

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

#PLOT 2D
def plot_2d(track_length_column="trackLengthToEcal_IKF_zedLambda", tof_column="tofClosest0"):
    print("1")
    df_beta = df.Define("mom", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)").\
                Define("beta", f"{track_length_column}/({tof_column}*299.792458)").Filter("beta >= 0 && beta <= 1")
    df_mass = df_beta.Define("mass", "mom*sqrt( 1./(beta*beta) - 1.)*1000")
    print("2")
    h = df_mass.Histo2D((f"h_{track_length_column}", f"{track_length_column}; momentum [GeV]; Mass [MeV]", 500, 0, 15, 500, -100, 1300.), "mom","mass")
    print("3")

    canvas = ROOT.TCanvas(f"{track_length_column}",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    print("4")

    h.Draw("colz")
    print("5")
    h.SetMinimum(1)
    h.SetMaximum(10000)
    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.SetRightMargin(0.12)
    canvas.Update()
    print("6")
    palette = h.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.89)
    palette.SetX2NDC(0.91)
    canvas.Modified()
    canvas.Update()
    print("7")
    input("wait")
    canvas.Close()
    print("8")


#PLOT 1D
def plot_1d():
    algorithms = ["trackLengthToEcal_IKF_zedLambda",
                "trackLengthToEcal_SHA_zedLambda_ECAL",
                "trackLengthToEcal_SHA_phiLambda_IP"]

    histos = []
    for name in algorithms:
        df_beta = df.Define("mom", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)").\
                    Define("beta", f"{name}/(tofClosest0*299.792458)").Filter("beta >= 0 && beta <= 1")
        df_mass = df_beta.Define("mass", "mom*sqrt( 1./(beta*beta) - 1.)*1000")
        h = df_mass.Histo1D((f"h_{name}", f"{name}; mass [MeV]; N entries", 2000, 0, 1100), "mass")
        histos.append(h)

    ROOT.gStyle.SetPadLeftMargin(0.18)
    hs = ROOT.THStack()
    canvas = ROOT.TCanvas("mass_1d_track_lenghts",
                            "",
                            int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                            int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    legend = ROOT.TLegend()
    for i, (name, h) in enumerate(zip(algorithms, histos)):
        h.Draw("" if i == 0 else "Lsame")
        h.SetLineColor(colors[i])
        h.SetLineWidth(3)
        legend.AddEntry(h.GetPtr(),name,"l")

    legend.Draw()
    canvas.Modified()
    canvas.Update()
    input("wait")


plot_1d()
# plot_2d("trackLengthToEcal_SHA_phiLambda_IP", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_phiZed_IP", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_zedLambda_IP", "tofClosest0")

# plot_2d("trackLengthToEcal_SHA_phiLambda_ECAL", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_phiZed_ECAL", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_zedLambda_ECAL", "tofClosest0")

# plot_2d("trackLengthToEcal_IKF_phiLambda", "tofClosest0")
# plot_2d("trackLengthToEcal_IKF_phiZed", "tofClosest0")
# plot_2d("trackLengthToEcal_IKF_zedLambda", "tofClosest0")