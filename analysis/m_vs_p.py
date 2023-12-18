import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

colors = ['#1b9e77', '#d95f02', '#7570b3']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Filter("tofClosest0 > 6.")
# ROOT.Experimental.AddProgressBar(df)
n_mom_bins, mom_min, mom_max = 2000, 0, 20
n_mass_bins, mass_min, mass_max = 2000, -3, 3

def get_2d_histo(df, tof_column="tofClosest0"):
    df = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")
    df_total = df.Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")
    h_2d_total = df_total.Histo2D(("h_total", "All; Momentum (GeV/c); Mass^{2} (GeV^{2}/c^{4})", n_mom_bins, mom_min, mom_max, n_mass_bins, mass_min, mass_max), "harmonicMomToEcal_IKF_zedLambda","mass2")

    n_pions = df.Filter("abs(pdg) == 211").Count().GetValue()
    n_kaons = df.Filter("abs(pdg) == 321").Count().GetValue()
    n_protons = df.Filter("abs(pdg) == 2212").Count().GetValue()

    n_pions_passed = df_total.Filter("abs(pdg) == 211").Count().GetValue()
    n_kaons_passed = df_total.Filter("abs(pdg) == 321").Count().GetValue()
    n_protons_passed = df_total.Filter("abs(pdg) == 2212").Count().GetValue()
    print(f"N pions: total: {n_pions} / passed {n_pions_passed} ({n_pions_passed/n_pions*100:.0f} %) / ignored {n_pions-n_pions_passed} ({((n_pions - n_pions_passed)/n_pions*100):.0f} %)")

    print(f"N kaons: total: {n_kaons} / passed {n_kaons_passed} ({n_kaons_passed/n_kaons*100:.0f} %) / ignored {n_kaons-n_kaons_passed} ({((n_kaons - n_kaons_passed)/n_kaons*100):.0f} %)")

    print(f"N protons: total: {n_protons} / passed {n_protons_passed} ({n_protons_passed/n_protons*100:.0f} %) / ignored {n_protons-n_protons_passed} ({((n_protons - n_protons_passed)/n_protons*100):.0f} %)")

    # h_2d_pion = df_total.Filter("abs(pdg) == 211")\
    #                     .Histo2D(("h_pion", "#pi; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    # h_2d_kaon = df_total.Filter("abs(pdg) == 321")\
    #                     .Histo2D(("h_kaon", "K; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    # h_2d_proton = df_total.Filter("abs(pdg) == 2212")\
    #                       .Histo2D(("h_proton", "proton; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    # DRAW FANCY 2D plot
    ROOT.gStyle.SetPadLeftMargin(0.17)
    ROOT.gStyle.SetPadRightMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.04)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    canvas = ROOT.TCanvas("c_2d_m_vs_p_total","",600,600)
    h_2d_total.Draw("colz")
    h_2d_total.GetYaxis().SetTitleOffset(1.1)
    print(h_2d_total.GetMaximum())
    h_2d_total.SetMinimum(1)
    h_2d_total.SetMaximum(10000)
    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.Update()
    palette = h_2d_total.GetListOfFunctions().FindObject("palette")
    palette.SetX1(mom_min + 1.01*(mom_max-mom_min))
    palette.SetX2(mom_min + 1.05*(mom_max-mom_min))
    palette.SetLabelOffset(0.001)

    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetTextSize(0.025)

    
    
    


    # latex.DrawLatex(0, -0.08, f"{'TOTAL':^28}" + "/" + f"{'PASSED (0 < #beta < 1)':^28}" + "/" + f"{'IGNORED ( #beta > 1 )':^28}")
                                        

    # latex.DrawLatex(0, -0.14, f"{f'pions: {n_pions}':^28}" + "/" + f"{f'{n_pions_passed} ({n_pions_passed/n_pions*100:.0f}%)':^28}" + "/" + f"{f'{n_pions-n_pions_passed} ({(n_pions - n_pions_passed)/n_pions*100:.0f}%)':^28}")
    # latex.DrawLatex(0, -0.2, f"{f'kaons: {n_kaons}':^28}" + "/" + f"{f'{n_kaons_passed} ({n_kaons_passed/n_kaons*100:.0f}%)':^28}" + "/" + f"{f'{n_kaons-n_kaons_passed} ({(n_kaons - n_kaons_passed)/n_kaons*100:.0f}%)':^28}")
    # latex.DrawLatex(0, -0.26, f"{f'protons: {n_protons}':^28}" + "/" + f"{f'{n_protons_passed} ({n_protons_passed/n_protons*100:.0f}%)':^28}" + "/" + f"{f'{n_protons-n_protons_passed} ({(n_protons - n_protons_passed)/n_protons*100:.0f}%)':^28}")

    canvas.Modified()
    canvas.Update()
    input("wait")

get_2d_histo(df, "tofClosest0")
