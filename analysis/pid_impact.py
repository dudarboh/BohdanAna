#!/usr/bin/env python3

"""Produces momentum distribution plots for chapter 8"""

import numpy as np
import ROOT
from utils import *

ROOT.EnableImplicitMT()


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


def get_efficiency_graphs():
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
    gr_eff_dedx_tof30 = convert_graph_sp_to_eff ( combine_two_graphs(gr_dedx_sp, gr_tof30_sp) )
    gr_eff_dedx_tof10 = convert_graph_sp_to_eff ( combine_two_graphs(gr_dedx_sp, gr_tof10_sp) )
    gr_eff_random = ROOT.TGraph()
    for i in range( gr_eff_dedx.GetN() ):
        x = gr_eff_dedx.GetPointX(i)
        gr_eff_random.SetPoint(i, x, 0.5)

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df_k = df.Filter("!isSimulated && !isOverlay")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Filter("hasTrack")\
            .Filter("abs(pdg) == 321")
            
    h_k_total = df_k.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    h_k_total_shower = df_k.Filter("hasShower").Histo1D((get_rand_string(), "Total with shower;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    h_k_total_noshower = df_k.Filter("!hasShower").Histo1D((get_rand_string(), "Total no shower;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )

    df_pi = df.Filter("!isSimulated && !isOverlay")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Filter("hasTrack")\
            .Filter("abs(pdg) == 211")

    h_pi_total = df_pi.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    h_pi_total_shower = df_pi.Filter("hasShower").Histo1D((get_rand_string(), "Total with shower;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    h_pi_total_noshower = df_pi.Filter("!hasShower").Histo1D((get_rand_string(), "Total no shower;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )


    h_k_dedx_track_only = apply_efficiency(h_k_total_noshower, gr_eff_dedx)
    h_k_random_track_only = apply_efficiency(h_k_total_noshower, gr_eff_random)
    h_k_dedx = apply_efficiency(h_k_total, gr_eff_dedx)
    h_k_dedx.SetLineColor(DEDX_COLOR)
    h_k_dedx.SetLineWidth(2)

    h_k_dedx_tof30 = apply_efficiency(h_k_total_shower, gr_eff_dedx_tof30)
    h_k_dedx_tof30.Add(h_k_dedx_track_only)
    h_k_dedx_tof30.SetLineColor( ROOT.TColor.GetColor("#c712a9") )
    h_k_dedx_tof30.SetLineWidth(2)

    h_k_dedx_tof10 = apply_efficiency(h_k_total_shower, gr_eff_dedx_tof10)
    h_k_dedx_tof10.Add(h_k_dedx_track_only)
    h_k_dedx_tof10.SetLineColor( ROOT.TColor.GetColor("#cd54b9") )
    h_k_dedx_tof10.SetLineWidth(2)

    h_k_tof30 = apply_efficiency(h_k_total_shower, gr_eff_tof30)
    h_k_tof30.Add(h_k_random_track_only)
    h_k_tof30.SetLineColor( ROOT.TColor.GetColor("#001d75") )
    h_k_tof30.SetLineWidth(2)

    h_k_tof10 = apply_efficiency(h_k_total_shower, gr_eff_tof10)
    h_k_tof10.Add(h_k_random_track_only)
    h_k_tof10.SetLineColor( ROOT.TColor.GetColor("#0091ea") )
    h_k_tof10.SetLineWidth(2)

    ############## PIONS
    h_pi_dedx_track_only = apply_misid(h_pi_total_noshower, gr_eff_dedx)
    h_pi_random_track_only = apply_misid(h_pi_total_noshower, gr_eff_random)
    h_pi_dedx = apply_misid(h_pi_total, gr_eff_dedx)
    h_pi_dedx.SetLineColor(DEDX_COLOR)
    h_pi_dedx.SetLineWidth(2)

    h_pi_dedx_tof30 = apply_misid(h_pi_total_shower, gr_eff_dedx_tof30)
    h_pi_dedx_tof30.Add(h_pi_dedx_track_only)
    h_pi_dedx_tof30.SetLineColor( ROOT.TColor.GetColor("#c712a9") )
    h_pi_dedx_tof30.SetLineWidth(2)

    h_pi_dedx_tof10 = apply_misid(h_pi_total_shower, gr_eff_dedx_tof10)
    h_pi_dedx_tof10.Add(h_pi_dedx_track_only)
    h_pi_dedx_tof10.SetLineColor( ROOT.TColor.GetColor("#cd54b9") )
    h_pi_dedx_tof10.SetLineWidth(2)

    h_pi_tof30 = apply_misid(h_pi_total_shower, gr_eff_tof30)
    h_pi_tof30.Add(h_pi_random_track_only)
    h_pi_tof30.SetLineColor( ROOT.TColor.GetColor("#001d75") )
    h_pi_tof30.SetLineWidth(2)

    h_pi_tof10 = apply_misid(h_pi_total_shower, gr_eff_tof10)
    h_pi_tof10.Add(h_pi_random_track_only)
    h_pi_tof10.SetLineColor( ROOT.TColor.GetColor("#0091ea") )
    h_pi_tof10.SetLineWidth(2)


    canvas = create_canvas()
    h_k_total.Draw()
    h_k_dedx.Draw("same")
    h_k_dedx_tof30.Draw("same")
    h_k_dedx_tof10.Draw("same")
    h_k_tof30.Draw("same")
    h_k_tof10.Draw("same")
    h_pi_dedx.Draw("same")
    h_pi_dedx.SetLineColor(3)
    h_pi_dedx_tof30.Draw("same")
    h_pi_dedx_tof30.SetLineColor(ROOT.kGreen+1)

    canvas.Update()

    canvas2 = create_canvas()
    h_pi_total.Draw()
    h_pi_dedx.Draw("same")
    h_pi_dedx_tof30.Draw("same")
    h_pi_dedx_tof10.Draw("same")
    h_pi_tof30.Draw("same")
    h_pi_tof10.Draw("same")

    canvas2.Update()

    input("wait")


def main():
    get_efficiency_graphs()

    # df analysis
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Filter("abs(pdg) == 321")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Define("pt", "sqrt(mcPx*mcPx + mcPy*mcPy)")\
            .Define("costheta", "cos(atan2(pt, mcPz))")
    df_gen = df.Filter("!isSimulated && !isOverlay")
    df_reco = df_gen.Filter("hasTrack")
    df_reco_shower = df_reco.Filter("hasShower")

    #creation and styling
    h1 = df_gen.Histo1D((get_rand_string(), "Generated;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    h2 = df_reco.Histo1D((get_rand_string(), "Has track;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    h3 = df_reco_shower.Histo1D((get_rand_string(), "Has track and shower;Momentum (GeV/c);N entries", 500, 0, 10), "mom" )
    n_gen_low = df_gen.Filter("mom < 3").Count()
    n_reco_low = df_reco.Filter("mom < 3").Count()
    n_reco_shower_low = df_reco_shower.Filter("mom < 3").Count()
    n_gen_total = df_gen.Count()
 
    # h_gen_test = df_gen.Histo2D((get_rand_string(), "Generated;pt;pz", 1000, 0, 4, 1000, -2, 2), "pt", "mcPz" )
    # h_reco_test = df_reco.Histo2D((get_rand_string(), "Has track;pt;pz", 1000, 0, 4, 1000, -2, 2), "pt", "mcPz" )
    ROOT.gStyle.SetPalette(1)
    h_gen_test = df_gen.Histo2D((get_rand_string(), "Generated;pt;pz", 500,-1, 1, 500, 0, 2), "costheta", "pt" )
    h_reco_test = df_reco.Histo2D((get_rand_string(), "Has track;pt;pz", 500,-1, 1, 500, 0, 2), "costheta", "pt" )

    h_reco_test.Divide(h_gen_test.GetPtr())
    h_reco_test.Draw("colz")
    input("wait")
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


    # c3 = create_canvas()
    # h_test.Draw("histo")
    # leg3 = create_legend()
    # leg3.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    # leg3.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    # leg3.Draw()
    # c3.Update()

    input("wait")










    print(f"Calculating for dEdx")
    dedx_graphs = { "pik" : {}, "kp" : {} }
    _, _, _, dedx_graphs["pik"]["TOF"] = calculate_pid_graphs(histos["TOF"][30][pion]["log"], histos["TOF"][30][kaon]["log"])
    dedx_graphs["pik"]["TOF"].SetTitle("TOF 30 ps;Momentum (GeV/c);#pi/K separation power")
    _, _, _, dedx_graphs["pik"]["dEdx"] = calculate_pid_graphs(histos["dEdx"][pion], histos["dEdx"][kaon])
    dedx_graphs["pik"]["dEdx"].SetTitle("dE/dx")
    c3 = draw_dedx_sep_powers(dedx_graphs["pik"]["TOF"], dedx_graphs["pik"]["dEdx"])
    c3.Modified()
    c3.Update()

    _, _, _, dedx_graphs["kp"]["TOF"] = calculate_pid_graphs(histos["TOF"][30][kaon]["log"], histos["TOF"][30][proton]["log"])
    dedx_graphs["kp"]["TOF"].SetTitle("TOF 30 ps;Momentum (GeV/c);K/p separation power")
    _, _, _, dedx_graphs["kp"]["dEdx"] = calculate_pid_graphs(histos["dEdx"][kaon], histos["dEdx"][proton])
    dedx_graphs["kp"]["dEdx"].SetTitle("dE/dx")
    c4 = draw_dedx_sep_powers(dedx_graphs["kp"]["TOF"], dedx_graphs["kp"]["dEdx"])
    c4.Modified()
    c4.Update()

    resolution_graphs = { "pik" : {}, "kp" : {} }
    for RES in RESOLUTIONS:
        print(f"Calculating for {RES} resolution")
        _, _, _, resolution_graphs["pik"][RES] = calculate_pid_graphs(histos["TOF"][RES][pion]["lin"], histos["TOF"][RES][kaon]["lin"])
        _, _, _, resolution_graphs["kp"][RES] = calculate_pid_graphs(histos["TOF"][RES][kaon]["lin"], histos["TOF"][RES][proton]["lin"])

    c1 = draw_resolution_sep_powers(resolution_graphs["pik"])
    resolution_graphs["pik"][0].SetTitle(";Momentum (GeV/c);#pi/K separation power")
    c1.Modified()
    c1.Update()
    c2 = draw_resolution_sep_powers(resolution_graphs["kp"])
    resolution_graphs["kp"][0].SetTitle(";Momentum (GeV/c);K/p separation power")
    c2.Modified()
    c2.Update()

    input("wait")

if __name__ == "__main__":
    main()