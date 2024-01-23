import ROOT
import numpy as np
import my_utils

ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

# NOTE: I use equal amount of bins in x and y to preserve each bin to be a square, like a pixel.
N_BINS, MOM_MIN, MOM_MAX = 1000, 0, 10

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
         .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Filter("tofClosest0 > 6.")
# ROOT.Experimental.AddProgressBar(df)

def get_2d_histo(df, tof_column="tofClosest0"):
    df = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")
    df = df.Define("inverseBeta", "1/beta")
    df = df.Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")
    df_mass = df.Filter("mass2 > 0.").Define("mass", "sqrt(mass2)")

    h_beta = df.Histo2D(("h_beta", "; Momentum (GeV/c); #beta ", N_BINS, MOM_MIN, MOM_MAX, N_BINS, 0.5, 1.1), "harmonicMomToEcal_IKF_zedLambda","beta")
    h_inverseBeta = df.Histo2D(("h_inverseBeta", "; Momentum (GeV/c); 1/#beta ", N_BINS, MOM_MIN, MOM_MAX, N_BINS, 0.9, 1.5), "harmonicMomToEcal_IKF_zedLambda","inverseBeta")
    h_mass2 = df.Histo2D(("h_mass2", "; Momentum (GeV/c); Mass^{2} (GeV^{2}/c^{4})", N_BINS, MOM_MIN, MOM_MAX, N_BINS, -0.3, 1.2), "harmonicMomToEcal_IKF_zedLambda","mass2")
    h_mass = df_mass.Histo2D(("h_mass", "; Momentum (GeV/c); Mass (GeV/c^{2})", N_BINS, MOM_MIN, MOM_MAX, N_BINS, -0.3, 1.2), "harmonicMomToEcal_IKF_zedLambda","mass")

    n_pions = df.Filter("abs(pdg) == 211").Count().GetValue()
    n_kaons = df.Filter("abs(pdg) == 321").Count().GetValue()
    n_protons = df.Filter("abs(pdg) == 2212").Count().GetValue()

    n_pions_passed = df_mass.Filter("abs(pdg) == 211").Count().GetValue()
    n_kaons_passed = df_mass.Filter("abs(pdg) == 321").Count().GetValue()
    n_protons_passed = df_mass.Filter("abs(pdg) == 2212").Count().GetValue()

    print(f"N pions: total: {n_pions} / passed {n_pions_passed} ({n_pions_passed/n_pions*100:.0f} %) / ignored {n_pions-n_pions_passed} ({((n_pions - n_pions_passed)/n_pions*100):.0f} %)")
    print(f"N kaons: total: {n_kaons} / passed {n_kaons_passed} ({n_kaons_passed/n_kaons*100:.0f} %) / ignored {n_kaons-n_kaons_passed} ({((n_kaons - n_kaons_passed)/n_kaons*100):.0f} %)")
    print(f"N protons: total: {n_protons} / passed {n_protons_passed} ({n_protons_passed/n_protons*100:.0f} %) / ignored {n_protons-n_protons_passed} ({((n_protons - n_protons_passed)/n_protons*100):.0f} %)")

    # DRAW FANCY 2D plot
    histos = [h_beta, h_inverseBeta, h_mass2, h_mass]
    # margin = 0.33
    # ROOT.gStyle.SetPadLeftMargin(0.6*margin)
    # ROOT.gStyle.SetPadRightMargin(0.4*margin)
    # ROOT.gStyle.SetPadTopMargin(0.35*margin)
    # ROOT.gStyle.SetPadBottomMargin(0.65*margin)
    # canvas = ROOT.TCanvas("c_2d_m_vs_p_total","",600,600)
    for h in histos:
        draw_2d_plot(h)
        # h.Draw("colz")
        # h.GetYaxis().SetTitleOffset(1.4)
        # h.GetXaxis().SetTitleOffset(1.1)
        # print(h.GetMaximum())
        # h.SetMinimum(1)
        # h.SetMaximum(10000)
        # canvas.SetLogz()
        # canvas.SetGridx(0)
        # canvas.SetGridy(0)
        # canvas.Update()
        # palette = h.GetListOfFunctions().FindObject("palette")
        # palette.SetX1(MOM_MIN + 1.01*(MOM_MAX-MOM_MIN))
        # palette.SetX2(MOM_MIN + 1.05*(MOM_MAX-MOM_MIN))
        # palette.SetLabelOffset(0.001)
        # # if h.GetName() == "h_mass":
        # #     latex = ROOT.TLatex()
        # #     latex.SetTextFont(62)
        # #     latex.SetTextSize(0.02)
        # #     latex.DrawLatex(0, -0.08, f"{'TOTAL':^28}" + "/" + f"{'PASSED (0 < #beta < 1)':^28}" + "/" + f"{'IGNORED':^28}")
        # #     latex.DrawLatex(0, -0.14, f"{f'pions: {n_pions}':^28}" + "/" + f"{f'{n_pions_passed} ({n_pions_passed/n_pions*100:.0f}%)':^28}" + "/" + f"{f'{n_pions-n_pions_passed} ({(n_pions - n_pions_passed)/n_pions*100:.0f}%)':^28}")
        # #     latex.DrawLatex(0, -0.2, f"{f'kaons: {n_kaons}':^28}" + "/" + f"{f'{n_kaons_passed} ({n_kaons_passed/n_kaons*100:.0f}%)':^28}" + "/" + f"{f'{n_kaons-n_kaons_passed} ({(n_kaons - n_kaons_passed)/n_kaons*100:.0f}%)':^28}")
        # #     latex.DrawLatex(0, -0.26, f"{f'protons: {n_protons}':^28}" + "/" + f"{f'{n_protons_passed} ({n_protons_passed/n_protons*100:.0f}%)':^28}" + "/" + f"{f'{n_protons-n_protons_passed} ({(n_protons - n_protons_passed)/n_protons*100:.0f}%)':^28}")

        # canvas.Modified()
        # canvas.Update()
        # input("wait")

get_2d_histo(df, "tofClosest30")
