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

#PLOT 1D
histos = []
df = df.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)").Filter("beta >= 0 && beta <= 1").\
        Define("mom_ip", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)").\
        Define("mom_ecal", "sqrt(recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy + recoCaloPz*recoCaloPz)").\
        Define("mom_hm", "harmonicMomToEcal_IKF_zedLambda")
        
h_ip = df.Define("mass", "mom_ip*sqrt( 1./(beta*beta) - 1.)*1000").Histo1D(("h_ip", "; mass [MeV]; N entries", 500, 130, 150), "mass")
h_ecal = df.Define("mass", "mom_ecal*sqrt( 1./(beta*beta) - 1.)*1000").Histo1D(("h_ecal", "; mass [MeV]; N entries", 500, 130, 150), "mass")
h_hm = df.Define("mass", "mom_hm*sqrt( 1./(beta*beta) - 1.)*1000").Histo1D(("h_hm", "; mass [MeV]; N entries", 500, 130, 150), "mass")
histos.append(h_ip)
histos.append(h_ecal)
histos.append(h_hm)

ROOT.gStyle.SetPadLeftMargin(0.18)
ROOT.gStyle.SetPadRightMargin(0.08)
canvas = ROOT.TCanvas("mass_1d_mom_options",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
canvas.SetGridx()
canvas.SetGridy()

legend = ROOT.TLegend()

names = ["p_{IP}", "p_{ECAL}", "p_{HM}"]
for i, (name, h) in enumerate(zip(names, histos)):
    h.Draw("L" if i == 0 else "Lsame")
    h.SetLineColor(colors[i])
    h.SetLineWidth(3)
    legend.AddEntry(h.GetPtr(),name,"l")

legend.Draw()
histos[0].GetYaxis().SetTitleOffset(1.45)
histos[0].SetMinimum(0)
canvas.Modified()
canvas.Update()
input("wait")