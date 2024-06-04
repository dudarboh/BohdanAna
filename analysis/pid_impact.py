#!/usr/bin/env python3

"""Produces momentum distribution plots for chapter 8"""

import numpy as np
import ROOT
from utils import *

ROOT.EnableImplicitMT()

dark24 = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']
dark24 = [ROOT.TColor.GetColor(c) for c in dark24]

#Linear momentum bins
N_MOMENTUM_BINS, MIN_MOMENTUM, MAX_MOMENTUM = 70, 0, 20 # GeV/c
#Log momentum bins for dE/dx
N_LOG_MOMENTUM_BINS, MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM = 30, -0.3, 1.3 # (0.5 to ~20) GeV/c
# LOG_MOMENTUM_BINS = np.array([ 10**(MIN_LOG_MOMENTUM + (MAX_LOG_MOMENTUM-MIN_LOG_MOMENTUM)*i/N_LOG_MOMENTUM_BINS ) for i in range(N_LOG_MOMENTUM_BINS+1) ])
LOG_MOMENTUM_BINS = np.logspace(MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM, N_LOG_MOMENTUM_BINS+1)
N_DEDX_BINS, MIN_DEDX, MAX_DEDX = 3000, 0, 1e-6
DEDX_BINS = np.linspace(MIN_DEDX, MAX_DEDX, N_DEDX_BINS+1)
N_MASS2_BINS, MIN_MASS2, MAX_MASS2 = 3000, -3, 3  # GeV^2/c^4
# MASS2_BINS = np.array([ MIN_MASS2 + (MAX_MASS2 - MIN_MASS2)*i/N_MASS2_BINS for i in range(N_MASS2_BINS+1) ])
MASS2_BINS = np.linspace(MIN_MASS2, MAX_MASS2, N_MASS2_BINS+1)
MOMENTUM_COLUMN = "harmonicMomToEcal_IKF_zedLambda" 
TRACK_LENGTH_COLUMN = "trackLengthToEcal_IKF_zedLambda"
RESOLUTIONS = [0, 1, 5, 10, 30, 50, 100, 300] # ps
COLORS_RESOLUTION = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#0091ea", "#0079d3", "#0061bd", "#004aa5", "#00348d", "#001d75", "#00045c"] ]
COLORS_DEDX = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#cd54b9", "#c52700"] ]
MIN_SEP_POWER, MAX_SEP_POWER = 0, 6

def apply_efficiency(h_mom, gr_eff):
    '''Return histogram where each bin is multiplied by the interpolated graph value'''
    h = h_mom.Clone()
    for i in range( 1, h_mom.GetXaxis().GetNbins() + 1 ):
        x = h_mom.GetXaxis().GetBinCenter(i)
        eff = gr_eff.Eval(x)
        h.SetBinContent(i, h_mom.GetBinContent(i)*eff )
    return h

def apply_misid(h_mom, gr_eff):
    '''Return histogram where each bin is multiplied by the interpolated 1-graph value'''
    #NOTE: this assumes misid = 1-eff
    h = h_mom.Clone()
    for i in range( 1, h_mom.GetXaxis().GetNbins() + 1 ):
        x = h_mom.GetXaxis().GetBinCenter(i)
        eff = gr_eff.Eval(x)
        h.SetBinContent(i, h_mom.GetBinContent(i)*(1-eff) )
    return h


def get_efficiency_graphs(tof_only=False):
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")

    h_pi_dedx = df.Filter("dEdx > 0.")\
            .Filter("abs(pdg) == 211")\
            .Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), MOMENTUM_COLUMN, "dEdx")
    h_k_dedx = df.Filter("dEdx > 0.")\
            .Filter("abs(pdg) == 321")\
            .Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), MOMENTUM_COLUMN, "dEdx")

    h_pi_tof30 = df.Filter("tofClosest0 > 6.")\
                .Filter("abs(pdg) == 211")\
                .Define("beta", f"({TRACK_LENGTH_COLUMN})/(tofClosest30*299.792458)" )\
                .Define("mass2", f"{MOMENTUM_COLUMN}*{MOMENTUM_COLUMN}*( 1./(beta*beta) - 1.)")\
                .Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")

    h_k_tof30 = df.Filter("tofClosest0 > 6.")\
                .Filter("abs(pdg) == 321")\
                .Define("beta", f"{TRACK_LENGTH_COLUMN}/(tofClosest30*299.792458)")\
                .Define("mass2", f"{MOMENTUM_COLUMN}*{MOMENTUM_COLUMN}*( 1./(beta*beta) - 1.)")\
                .Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")


    h_pi_tof10 = df.Filter("tofClosest0 > 6.")\
                .Filter("abs(pdg) == 211")\
                .Define("beta", f"{TRACK_LENGTH_COLUMN}/(tofClosest10*299.792458)")\
                .Define("mass2", f"{MOMENTUM_COLUMN}*{MOMENTUM_COLUMN}*( 1./(beta*beta) - 1.)")\
                .Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")

    h_k_tof10 = df.Filter("tofClosest0 > 6.")\
                .Filter("abs(pdg) == 321")\
                .Define("beta", f"{TRACK_LENGTH_COLUMN}/(tofClosest10*299.792458)")\
                .Define("mass2", f"{MOMENTUM_COLUMN}*{MOMENTUM_COLUMN}*( 1./(beta*beta) - 1.)")\
                .Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")

    _, gr_eff_dedx, _, gr_dedx_sp = calculate_pid_graphs(h_pi_dedx, h_k_dedx)
    _, gr_eff_tof30, _, gr_tof30_sp = calculate_pid_graphs(h_pi_tof30, h_k_tof30)
    _, gr_eff_tof10, _, gr_tof10_sp = calculate_pid_graphs(h_pi_tof10, h_k_tof10)

    gr_sp_dedx_tof30 = combine_two_graphs(gr_dedx_sp, gr_tof30_sp)
    gr_sp_dedx_tof10 = combine_two_graphs(gr_dedx_sp, gr_tof10_sp)

    gr_eff_dedx_tof30 = convert_graph_sp_to_eff ( gr_sp_dedx_tof30 )
    gr_eff_dedx_tof10 = convert_graph_sp_to_eff ( gr_sp_dedx_tof10 )
    gr_eff_random = ROOT.TGraph()

    for i in range( gr_eff_dedx.GetN() ):
        x = gr_eff_dedx.GetPointX(i)
        gr_eff_random.SetPoint(i, x, 0.5)

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")\
            .Filter("!isSimulated && !isOverlay")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Filter("hasTrack")

    n_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root").Count()


    df_k = df.Filter("abs(pdg) == 321")
    h_k_total = df_k.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);N entries", 400, 0, 8), "mom" )
    h_k_total_shower = df_k.Filter("hasShower").Histo1D((get_rand_string(), "Total with shower;Momentum (GeV/c);N entries", 400, 0, 8), "mom" )
    h_k_total_noshower = df_k.Filter("!hasShower").Histo1D((get_rand_string(), "Total no shower;Momentum (GeV/c);N entries", 400, 0, 8), "mom" )

    df_pi = df.Filter("abs(pdg) == 211")
    h_pi_total = df_pi.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);N entries", 400, 0, 8), "mom" )
    h_pi_total_shower = df_pi.Filter("hasShower").Histo1D((get_rand_string(), "Total with shower;Momentum (GeV/c);N entries", 400, 0, 8), "mom" )
    h_pi_total_noshower = df_pi.Filter("!hasShower").Histo1D((get_rand_string(), "Total no shower;Momentum (GeV/c);N entries", 400, 0, 8), "mom" )

    if not tof_only:
        # dE/dx kaon efficiency
        h_k_dedx = apply_efficiency(h_k_total, gr_eff_dedx)
        h_k_dedx.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_k_dedx.SetLineWidth(2)
        h_k_dedx.SetLineStyle(2)

        # dE/dx pion mis-id
        h_pi_dedx = apply_misid(h_pi_total, gr_eff_dedx)
        h_pi_dedx.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_pi_dedx.SetLineWidth(2)
        h_pi_dedx.SetLineStyle(2)

        # dE/dx + TOF 30 kaon efficiency
        h_k_dedx_tof30 = apply_efficiency(h_k_total_shower, gr_eff_dedx_tof30)
        h_k_dedx_track_only = apply_efficiency(h_k_total_noshower, gr_eff_dedx)
        h_k_dedx_tof30.Add(h_k_dedx_track_only)
        h_k_dedx_tof30.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_k_dedx_tof30.SetLineWidth(2)
        h_k_dedx_tof30.SetLineStyle(7)

        # dE/dx + TOF 30 pion mis-id
        h_pi_dedx_tof30 = apply_misid(h_pi_total_shower, gr_eff_dedx_tof30)
        h_pi_dedx_track_only = apply_misid(h_pi_total_noshower, gr_eff_dedx)
        h_pi_dedx_tof30.Add(h_pi_dedx_track_only)
        h_pi_dedx_tof30.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_pi_dedx_tof30.SetLineWidth(2)
        h_pi_dedx_tof30.SetLineStyle(7)

        # dE/dx + TOF 10 kaon efficiency
        h_k_dedx_tof10 = apply_efficiency(h_k_total_shower, gr_eff_dedx_tof10)
        h_k_dedx_tof10.Add(h_k_dedx_track_only)
        h_k_dedx_tof10.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_k_dedx_tof10.SetLineWidth(2)

        # dE/dx + TOF 10 pion mis-id
        h_pi_dedx_tof10 = apply_misid(h_pi_total_shower, gr_eff_dedx_tof10)
        h_pi_dedx_tof10.Add(h_pi_dedx_track_only)
        h_pi_dedx_tof10.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_pi_dedx_tof10.SetLineWidth(2)

        h_k_total.GetYaxis().SetMaxDigits(3)
        canvas1 = create_canvas()
        h_k_total.DrawClone()
        h_k_dedx.DrawClone("same")
        h_k_dedx_tof30.DrawClone("same")
        h_k_dedx_tof10.DrawClone("same")
        h_pi_dedx.DrawClone("same")
        h_pi_dedx_tof30.DrawClone("same")
        h_pi_dedx_tof10.DrawClone("same")

        leg = create_legend(0.2, 0.75, 0.76, 0.91)
        leg.SetNColumns(2)
        leg.SetMargin(0.15)

        h_leg1 = ROOT.TH1F(get_rand_string(), "K total (#varepsilon=100%)", 1, 0, 1)
        leg.AddEntry(h_leg1, h_leg1.GetTitle(), "l")

        h_leg2 = ROOT.TH1F(get_rand_string(), "dE/dx only", 1, 0, 1)
        h_leg2.SetLineStyle(2)
        leg.AddEntry(h_leg2, h_leg2.GetTitle(), "l")

        h_leg3 = ROOT.TH1F(get_rand_string(), "K identified", 1, 0, 1)
        h_leg3.SetLineColor(ROOT.TColor.GetColor("#ff0000"))
        leg.AddEntry(h_leg3, h_leg3.GetTitle(), "l")

        h_leg4 = ROOT.TH1F(get_rand_string(), "dE/dx + TOF 30 ps", 1, 0, 1)
        h_leg4.SetLineStyle(7)
        leg.AddEntry(h_leg4, h_leg4.GetTitle(), "l")

        h_leg5 = ROOT.TH1F(get_rand_string(), "#pi misidentified", 1, 0, 1)
        h_leg5.SetLineColor(ROOT.TColor.GetColor("#0066ff"))
        leg.AddEntry(h_leg5, h_leg5.GetTitle(), "l")

        h_leg6 = ROOT.TH1F(get_rand_string(), "dE/dx + TOF 10 ps", 1, 0, 1)
        h_leg6.SetLineStyle(1)
        leg.AddEntry(h_leg6, h_leg6.GetTitle(), "l")
        leg.DrawClone()
        canvas1.Update()
        return canvas1
    else:
        # TOF 30 kaon efficiency
        h_k_random_track_only = apply_efficiency(h_k_total_noshower, gr_eff_random)
        h_k_tof30 = apply_efficiency(h_k_total_shower, gr_eff_tof30)
        h_k_tof30.Add(h_k_random_track_only)
        h_k_tof30.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_k_tof30.SetLineWidth(2)
        h_k_tof30.SetLineStyle(7)

        # TOF 30 pion mis-id
        h_pi_random_track_only = apply_misid(h_pi_total_noshower, gr_eff_random)
        h_pi_tof30 = apply_misid(h_pi_total_shower, gr_eff_tof30)
        h_pi_tof30.Add(h_pi_random_track_only)
        h_pi_tof30.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_pi_tof30.SetLineWidth(2)
        h_pi_tof30.SetLineStyle(7)

        # TOF 10 kaon efficiency
        h_k_tof10 = apply_efficiency(h_k_total_shower, gr_eff_tof10)
        h_k_tof10.Add(h_k_random_track_only)
        h_k_tof10.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_k_tof10.SetLineWidth(2)

        # TOF 10 pion mis-id
        h_pi_tof10 = apply_misid(h_pi_total_shower, gr_eff_tof10)
        h_pi_tof10.Add(h_pi_random_track_only)
        h_pi_tof10.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_pi_tof10.SetLineWidth(2)

        h_k_total.GetYaxis().SetMaxDigits(3)
        canvas2 = create_canvas()
        h_k_total.DrawClone()
        h_k_tof30.DrawClone("same")
        h_k_tof10.DrawClone("same")
        h_pi_tof30.DrawClone("same")
        h_pi_tof10.DrawClone("same")

        leg = create_legend(0.2, 0.75, 0.76, 0.91)
        leg.SetNColumns(2)
        leg.SetMargin(0.15)

        h_leg1 = ROOT.TH1F(get_rand_string(), "K total (#varepsilon=100%)", 1, 0, 1)
        leg.AddEntry(h_leg1, h_leg1.GetTitle(), "l")

        h_leg2 = ROOT.TH1F(get_rand_string(), "", 1, 0, 1)
        h_leg2.SetLineStyle(2)
        leg.AddEntry(0, "", "")

        h_leg3 = ROOT.TH1F(get_rand_string(), "K identified", 1, 0, 1)
        h_leg3.SetLineColor(ROOT.TColor.GetColor("#ff0000"))
        leg.AddEntry(h_leg3, h_leg3.GetTitle(), "l")

        h_leg4 = ROOT.TH1F(get_rand_string(), "TOF 30 ps", 1, 0, 1)
        h_leg4.SetLineStyle(7)
        leg.AddEntry(h_leg4, h_leg4.GetTitle(), "l")

        h_leg5 = ROOT.TH1F(get_rand_string(), "#pi misidentified", 1, 0, 1)
        h_leg5.SetLineColor(ROOT.TColor.GetColor("#0066ff"))
        leg.AddEntry(h_leg5, h_leg5.GetTitle(), "l")

        h_leg6 = ROOT.TH1F(get_rand_string(), "TOF 10 ps", 1, 0, 1)
        h_leg6.SetLineStyle(1)
        leg.AddEntry(h_leg6, h_leg6.GetTitle(), "l")
        leg.DrawClone()
        canvas2.Update()
        return canvas2


def draw_gen_vs_reco():
    # df analysis
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55").Filter("abs(pdg) == 321")
    df_reco = df_gen.Filter("hasTrack")
    df_reco_shower = df_reco.Filter("hasShower")

    #creation and styling
    h1 = df_gen.Histo1D((get_rand_string(), "Generated;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h2 = df_reco.Histo1D((get_rand_string(), "Has track;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h3 = df_reco_shower.Histo1D((get_rand_string(), "Has track and shower;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    n_gen_low = df_gen.Filter("mom < 3").Count()
    n_reco_low = df_reco.Filter("mom < 3").Count()
    n_reco_shower_low = df_reco_shower.Filter("mom < 3").Count()
    n_gen_total = df_gen.Count()
 
    h2.SetLineColor(ROOT.TColor.GetColor("#20bf55"))
    h3.SetLineColor(ROOT.TColor.GetColor("#01baef"))
    h1.SetLineWidth(4)
    h2.SetLineWidth(4)
    h3.SetLineWidth(4)

    #Calculating
    print(f"Fraction of generated below 3 GeV {100*n_gen_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track below 3 GeV {100*n_reco_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track+shower below 3 GeV {100*n_reco_shower_low.GetValue()/n_gen_total.GetValue()}")


    # plotting
    c1 = create_canvas()
    h1.GetYaxis().SetTitleOffset(1.4)
    h1.DrawClone()
    h2.DrawClone("same")
    h3.DrawClone("same")
    leg = create_legend()
    leg.AddEntry(h1.GetPtr(), h1.GetTitle(), "l")
    leg.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg.Draw()
    c1.Update()

    c2 = create_canvas()
    h2.Divide(h1.GetPtr())
    h3.Divide(h1.GetPtr())
    h2.Draw("histo")
    h2.GetYaxis().SetTitleOffset(1.4)
    h3.Draw("histo same")
    leg2 = create_legend()
    leg2.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg2.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg2.Draw()
    c2.Update()
    #splitline{Kaons from Z#rightarrowq#bar{q}}{E_{cm} = 250 GeV/c^{2}}

    input("wait")


def test():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Define("pt", "sqrt(mcPx*mcPx + mcPy*mcPy)")\
            .Define("costheta", "cos(atan2(pt, mcPz))")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("abs(pdg) == 2212")
    df_reco = df_gen.Filter("hasTrack")

    # h_gen_test = df_gen.Histo2D((get_rand_string(), "Generated;pt;pz", 1000, 0, 4, 1000, -2, 2), "pt", "mcPz" )
    # h_reco_test = df_reco.Histo2D((get_rand_string(), "Has track;pt;pz", 1000, 0, 4, 1000, -2, 2), "pt", "mcPz" )
    ROOT.gStyle.SetPalette(1)
    c1 = create_canvas()
    h_gen_test = df_gen.Histo2D((get_rand_string(), "Generated;cos theta;pt", 200,-1, 1, 200, 0, 4), "costheta", "pt" )
    h_reco_test = df_reco.Histo2D((get_rand_string(), "Has track;cos theta;pt", 200,-1, 1, 200, 0, 4), "costheta", "pt" )
    h_reco_test.Divide(h_gen_test.GetPtr())
    h_reco_test.Draw("colz")
    c1.Update()
    c2 = create_canvas()
    h_gen_test2 = df_gen.Histo2D((get_rand_string(), "Generated;cos theta;p", 200,-1, 1, 200, 0, 4), "costheta", "mom" )
    h_reco_test2 = df_reco.Histo2D((get_rand_string(), "Has track;cos theta;p", 200,-1, 1, 200, 0, 4), "costheta", "mom" )
    h_reco_test2.Divide(h_gen_test2.GetPtr())
    h_reco_test2.Draw("colz")
    c2.Update()

    input("wait")

def draw_hadr_vs_secondary():
    # df analysis
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55").Filter("abs(pdg) == 321")
    df1 = df_gen.Filter("isHadronisationDecay")
    df2 = df_gen.Filter("isBottomQuarkDecay")

    #creation and styling
    h1 = df_gen.Histo1D((get_rand_string(), "Generated;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h2 = df1.Histo1D((get_rand_string(), "Not from b-quark decay;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h3 = df2.Histo1D((get_rand_string(), "From b quark decay;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    n_gen_low = df_gen.Filter("mom < 3").Count()
    n_reco_low = df1.Filter("mom < 3").Count()
    n_reco_shower_low = df2.Filter("mom < 3").Count()
    n_gen_total = df_gen.Count()
 
    h2.SetLineColor(ROOT.TColor.GetColor("#20bf55"))
    h3.SetLineColor(ROOT.TColor.GetColor("#01baef"))
    h1.SetLineWidth(4)
    h2.SetLineWidth(4)
    h3.SetLineWidth(4)

    #Calculating
    print(f"Fraction of generated below 3 GeV {100*n_gen_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track below 3 GeV {100*n_reco_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track+shower below 3 GeV {100*n_reco_shower_low.GetValue()/n_gen_total.GetValue()}")


    # plotting
    c1 = create_canvas()
    h1.GetYaxis().SetTitleOffset(1.3)
    h1.GetYaxis().SetMaxDigits(3)
    h1.GetYaxis().SetNdivisions(512)
    h1.DrawClone()
    h2.DrawClone("same")
    h3.DrawClone("same")
    leg = create_legend()
    leg.AddEntry(h1.GetPtr(), h1.GetTitle(), "l")
    leg.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg.Draw()
    c1.Update()

    c2 = create_canvas()
    h2.Divide(h1.GetPtr())
    h3.Divide(h1.GetPtr())
    h2.Draw("histo")
    h2.GetYaxis().SetTitleOffset(1.3)
    h2.GetYaxis().SetMaxDigits(3)
    h2.GetYaxis().SetNdivisions(512)
    h3.Draw("histo same")
    leg2 = create_legend()
    leg2.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg2.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg2.Draw()
    c2.Update()
    #splitline{Kaons from Z#rightarrowb#bar{b}}{E_{cm} = 250 GeV/c^{2}}

    input("wait")


def draw_prim_vs_secondary():
    # df analysis
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55").Filter("abs(pdg) == 321").Filter("hasTrack")
    df1 = df_gen.Filter("isInRecoPrimaryVertex")
    df2 = df_gen.Filter("isInRecoSecondaryVertex")

    #creation and styling
    h1 = df_gen.Histo1D((get_rand_string(), "Total reconstructed (has track);Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h2 = df1.Histo1D((get_rand_string(), "Reconstructed in primary vertex;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h3 = df2.Histo1D((get_rand_string(), "Reconstructed in secondary vertex;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    n_gen_low = df_gen.Filter("mom < 3").Count()
    n_reco_low = df1.Filter("mom < 3").Count()
    n_reco_shower_low = df2.Filter("mom < 3").Count()
    n_gen_total = df_gen.Count()
 
    h2.SetLineColor(ROOT.TColor.GetColor("#20bf55"))
    h3.SetLineColor(ROOT.TColor.GetColor("#01baef"))
    h1.SetLineWidth(4)
    h2.SetLineWidth(4)
    h3.SetLineWidth(4)

    #Calculating
    print(f"Fraction of generated below 3 GeV {100*n_gen_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track below 3 GeV {100*n_reco_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track+shower below 3 GeV {100*n_reco_shower_low.GetValue()/n_gen_total.GetValue()}")


    # plotting
    c1 = create_canvas()
    h1.GetYaxis().SetTitleOffset(1.3)
    h1.GetYaxis().SetMaxDigits(3)
    h1.GetYaxis().SetNdivisions(512)
    h1.DrawClone()
    h2.DrawClone("same")
    h3.DrawClone("same")
    leg = create_legend()
    leg.AddEntry(h1.GetPtr(), h1.GetTitle(), "l")
    leg.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg.Draw()
    c1.Update()

    c2 = create_canvas()
    h2.Divide(h1.GetPtr())
    h3.Divide(h1.GetPtr())
    h2.Draw("histo")
    h2.GetYaxis().SetTitleOffset(1.3)
    h2.GetYaxis().SetMaxDigits(3)
    h2.GetYaxis().SetNdivisions(512)
    h3.Draw("histo same")
    leg2 = create_legend()
    leg2.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg2.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg2.Draw()
    c2.Update()
    #splitline{Kaons from Z#rightarrowb#bar{b}}{E_{cm} = 250 GeV/c^{2}}

    input("wait")


def get_particles_fractions():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root").Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55")
    leg = create_legend()
    n_total = df.Count()
    n_pi = df.Filter("abs(pdg) == 211").Count()
    n_k = df.Filter("abs(pdg) == 321").Count()
    n_p = df.Filter("abs(pdg) == 2212").Count()

    print("Pi fraction:", 100*n_pi.GetValue()/n_total.GetValue())
    print("K fraction:", 100*n_k.GetValue()/n_total.GetValue())
    print("P fraction:", 100*n_p.GetValue()/n_total.GetValue())
    input("wait")

# get_particles_fractions()
# draw_gen_vs_reco()
# draw_hadr_vs_secondary()
# test()
c1 = get_efficiency_graphs()

# c2 = get_efficiency_graphs(tof_only=True)
input("wait")