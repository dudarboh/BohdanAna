import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

def plots_2d():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
            .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212").Filter("tofClosest0 > 6.").Filter("tofSETFront0 > 0.1")

    df_latest = df.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)")\
            .Define("mom2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda")\
            .Define("mom", "sqrt(mom2)")\
            .Define("mass2", "mom2*(1./(beta*beta) - 1.)")\
            .Define("mass", "sqrt(mass2)*1000.")
    h2d1 = df_latest.Histo2D( (get_rand_string(), ";Momentum (GeV/c);Mass^{2} (GeV^{2}/c^{4})", 1000, 0, 10, 1000, -0.3, 1.2), "mom","mass2" )


    df_latest_to_set = df.Define("beta", "trackLengthToSET_IKF_zedLambda/(tofSETFront0*299.792458)")\
            .Define("mom2", "harmonicMomToSET_IKF_zedLambda*harmonicMomToSET_IKF_zedLambda")\
            .Define("mom", "sqrt(mom2)")\
            .Define("mass2", "mom2*(1./(beta*beta) - 1.)")\
            .Define("mass", "sqrt(mass2)*1000.")
    h2d2 = df_latest_to_set.Histo2D( (get_rand_string(), ";Momentum (GeV/c);Mass^{2} (GeV^{2}/c^{4})", 1000, 0, 10, 1000, -0.3, 1.2), "mom","mass2" )

    c1 = draw_2d_plot(h2d1)
    c2 = draw_2d_plot(h2d2)
    input("wait111111")

def check_set_bias():

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
            .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212").Filter("tofClosest0 > 6.")\
            .Filter("tofSETFront0 > 0.1")\
            .Define("beta", "trackLengthToSET_IKF_zedLambda/(tofSETFront0*299.792458)")\
            .Define("mom2", "harmonicMomToSET_IKF_zedLambda*harmonicMomToSET_IKF_zedLambda")\
            .Define("mom", "sqrt(mom2)")\
            .Define("mass2", "mom2*(1./(beta*beta) - 1.)")\
            .Define("mass", "sqrt(mass2)*1000.")

    histos_2d = {}
    histos = {}
    for p in particles:
        histos_2d[p.name] = df.Histo2D( (get_rand_string(), ";Momentum (GeV/c);Mass (MeV/c^{2})", 100, 0, 5, 100, 1000*p.mass-15, 1000*p.mass+15), "mom","mass" )

        histos[p.name] = df.Histo1D( (get_rand_string(), ";Mass (MeV/c^{2});N entries", 1000, 1000*p.mass-15, 1000*p.mass+15), "mass" )

    canvases = []
    lines = []

    for p, h in zip(particles, histos.values()):
        c2d = draw_2d_plot(histos_2d[p.name], histos_2d[p.name].GetMaximum())
        c2d.SetLogz(False)
        histos_2d[p.name].SetMinimum(0)
        line = ROOT.TLine(0., p.mass*1000., 5., p.mass*1000.)
        line.SetLineStyle(9)
        line.SetLineWidth(3)
        line.SetLineColor(ROOT.kRed+1)
        line.Draw()
        lines.append(line)
        c2d.Update()
        canvases.append(c2d)
        canvas = create_canvas(margin=0.33, left_margin_fraction=0.55, bottom_margin_fraction=0.65)
        h.GetXaxis().SetMaxDigits(3)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetYaxis().SetTitleOffset(1.4)
        h.GetYaxis().SetMaxDigits(3)
        h.SetLineWidth(3)
        h.SetMinimum( 0. )
        h.Draw("L")
        peak_mass = h.GetXaxis().GetBinCenter( h.GetMaximumBin() )
        fit = ROOT.TF1(f"fit{p.name}", "gaus", peak_mass-1., peak_mass+1.)
        fit.SetLineWidth(3)
        h.Fit(fit, "RQ")
        print(f"Mean for {p.name}", h.GetFunction(f"fit{p.name}").GetParameter(1), "+-",  h.GetFunction(f"fit{p.name}").GetParameter(2) )
        line = ROOT.TLine(p.mass*1000, 0., p.mass*1000, h.GetMaximum())
        line.SetLineStyle(9)
        line.SetLineWidth(3)
        line.SetLineColor(ROOT.kRed+1)
        line.Draw()
        lines.append(line)
        canvas.Update()
        canvases.append(canvas)

    input("wait")

def refit_impact():
    ROOT.gStyle.SetOptTitle(True)
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAnaNewUncomplete.root")\
            .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
            .Filter("tofClosest0 > 6.")\
            .Filter("caloIDClosest == 1")

    df = df.Define("pTrue", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df = df.Define("ptTrue", "sqrt(mcPx*mcPx + mcPy*mcPy)")
    df = df.Define("ptIp", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy)")
    df = df.Define("refittedPtIp", "sqrt(refittedRecoIpPx*refittedRecoIpPx + refittedRecoIpPy*refittedRecoIpPy)")
    df = df.Define("omegaTrue", "1e-6 * 299.792458 * 3.5 / ptTrue")
    df = df.Define("tanLambdaTrue", "1e-6 * 299.792458 * 3.5 / ptTrue")

    # Do omega vs omega plots
    # h_omega_k = df.Filter("abs(pdg) == 321").Define("diff_y", "refittedOmegaIP - omegaIP")\
    #             .Define("diff_x", "omegaIP - omegaTrue")\
    #             .Histo2D(("h_omega_k", "Kaons;#Omega_{default} - #Omega_{true} (#frac{1}{mm});#Omega_{refitted} - #Omega_{default} (#frac{1}{mm})", 500, -0.0001, 0.0001, 500, -0.0001, 0.0001), "diff_x", "diff_y")
    # h_omega_p = df.Filter("abs(pdg) == 2212").Define("diff_y", "refittedOmegaIP - omegaIP")\
    #             .Define("diff_x", "omegaIP - omegaTrue")\
    #             .Histo2D(("h_omega_p", "Protons;#Omega_{default} - #Omega_{true} (#frac{1}{mm});#Omega_{refitted} - #Omega_{default} (#frac{1}{mm})", 500, -0.0001, 0.0001, 500, -0.0001, 0.0001), "diff_x", "diff_y")

    # c2 = draw_2d_plot(h_omega_k)
    # c3 = draw_2d_plot(h_omega_p)
    # input("wait")

    VARIABLE = "omega"
    if VARIABLE == "omega":
        title_default = ";Momentum (GeV/c);#Omega_{default} - #Omega_{true} (1/mm)"
        title_refitted = ";Momentum (GeV/c);#Omega_{refitted} - #Omega_{true} (1/mm)"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.0001, 0.0001
        y_column_default = "omegaIP - omegaTrue"
        y_column_refitted = "refittedOmegaIP - omegaTrue"
    elif VARIABLE == "tanLambda":
        title_default = ";Momentum (GeV/c);tan#lambda_{default} - tan#lambda_{true} "
        title_refitted = ";Momentum (GeV/c);tan#lambda_{refitted} - tan#lambda_{true}"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.01, 0.01
        y_column_default = "tanLambdaIP - omegaTrue"
        y_column_refitted = "refittedTanLambdaIP - omegaTrue"

    # Omega vs momentum
    h_k_default = df.Filter("abs(pdg) == 321").Define("diff", y_column_default)\
                .Histo2D((get_rand_string(), title_default, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_p_default = df.Filter("abs(pdg) == 2212").Define("diff", y_column_default)\
                .Histo2D((get_rand_string(), title_default, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_k_refitted = df.Filter("abs(pdg) == 321").Define("diff", y_column_refitted)\
                .Histo2D((get_rand_string(), title_refitted, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_p_refitted = df.Filter("abs(pdg) == 2212").Define("diff", y_column_refitted)\
                .Histo2D((get_rand_string(), title_refitted, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")

    c1 = draw_2d_plot(h_k_default, 0.4, 0.65, 0.65)
    h_k_default.GetYaxis().SetTitleOffset(1.8)
    h_p_default.GetYaxis().SetTitleOffset(1.8)
    h_k_refitted.GetYaxis().SetTitleOffset(1.8)
    h_p_refitted.GetYaxis().SetTitleOffset(1.8)

    c1.Print("kaons_default.png")
    c2 = draw_2d_plot(h_p_default, 0.4, 0.65, 0.65)
    c2.Print("protons_default.png")
    c3 = draw_2d_plot(h_k_refitted, 0.4, 0.65, 0.65)
    c3.Print("kaons_refit.png")
    c4 = draw_2d_plot(h_p_refitted, 0.4, 0.65, 0.65)
    c4.Print("protons_refit.png")

    input("wait")

    # pt vs momentum
    # h_k1 = df.Filter("abs(pdg) == 321").Define("diff", "ptIp - ptTrue")\
    #             .Histo2D((get_rand_string(), ";Momentum (GeV/c);p_{default} - p_{true} (GeV/c)", 200, 0, 2, 200, -0.1, 0.1), "pTrue", "diff")
    # h_p1 = df.Filter("abs(pdg) == 2212").Define("diff", "ptIp - ptTrue")\
    #             .Histo2D((get_rand_string(), ";Momentum (GeV/c);p_{default} - p_{true} (GeV/c)", 200, 0, 2, 200, -0.1, 0.1), "pTrue", "diff")

    # h_k2 = df.Filter("abs(pdg) == 321").Define("diff", "refittedPtIp - ptTrue")\
    #             .Histo2D((get_rand_string(), ";Momentum (GeV/c);p_{refitted} - p_{true} (GeV/c)", 200, 0, 2, 200, -0.1, 0.1), "pTrue", "diff")
    # h_p2 = df.Filter("abs(pdg) == 2212").Define("diff", "refittedPtIp - ptTrue")\
    #             .Histo2D((get_rand_string(), ";Momentum (GeV/c);p_{refitted} - p_{true} (GeV/c)", 200, 0, 2, 200, -0.1, 0.1), "pTrue", "diff")

    # c5 = draw_2d_plot(h_k1, 0.4, 0.65, 0.65)
    # c5.Print("kaons_default.png")
    # c6 = draw_2d_plot(h_k2, 0.4, 0.65, 0.65)
    # c6.Print("protons_default.png")
    # c7 = draw_2d_plot(h_p1, 0.4, 0.65, 0.65)
    # c7.Print("kaons_refit.png")
    # c8 = draw_2d_plot(h_p2, 0.4, 0.65, 0.65)
    # c8.Print("protons_refit.png")

    # input("wait")

    # df = df.Define("ptRefitDiff", "refittedPtIp - ptTrue").Define("diffOfdiffs", "ptRefitDiff - ptDefaultDiff")

    # h = df.Histo2D(("h", "title;reco default;reco perfect pid", 500, -0.01, 0.01, 500, -0.02, 0.02), "ptDefaultDiff", "diffOfdiffs")
    # c1 = draw_2d_plot(h)

    # c2 = create_canvas(0.33, 0.58, 0.65)
    # h2 = df.Histo2D(("h2", ";p_{T}^{true} (GeV/c);p_{T}^{reco, default} (GeV/c)", 500, 0, 1, 500, 0, 1), "ptTrue", "ptIp")
    # h2.Draw("colz")
    # c2.SetLogz()
    # c2.Update()

    # c3 = create_canvas(0.33, 0.58, 0.65)
    # h3 = df.Histo2D(("h3", ";p_{T}^{true} (GeV/c);p_{T}^{reco, refitted} (GeV/c)", 500, 0, 1, 500, 0, 1), "ptTrue", "refittedPtIp")
    # h3.Draw("colz")
    # c3.SetLogz()
    # c3.Update()

    input("wait")


def tof_vs_dedx():
    resolutions = np.array([0, 1, 5, 10, 17, 30, 50, 100]) # ps
    pik_mom_low = np.array([0.75, 0.75])
    pik_mom_high = np.array([6.2, 6.14, 5.3])
    kp_mom_low = np.array([1.26, 1.26])
    kp_mom_high = np.array([23, 23])

    canvas = create_canvas()
    pass


# refit_impact()
tof_vs_dedx()