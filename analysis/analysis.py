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

check_set_bias()