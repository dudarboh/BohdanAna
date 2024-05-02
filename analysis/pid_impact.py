#!/usr/bin/env python3

"""Produces separation power plots"""

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

def convert_p_value_to_sep_power(p_value):
    '''Return separation power equivalent to the given p-value (1 - efficiency)'''
    return 2*ROOT.Math.gaussian_quantile_c(p_value, 1)

def convert_sep_power_to_eff(sep_power):
    '''Return efficiency from the given sep power'''
    return ROOT.Math.gaussian_cdf(0.5*sep_power)

def find_optimal_cut(h1, h2, debug=False):
    '''Calcualte the x value between the two histograms where eff=1-misid'''
    cut_mass = []
    efficiencies = []
    mis_ids = []
    diff = []
    n_signal = h1.GetEntries()
    n_background = h2.GetEntries()

    for bin in range(1, h1.GetXaxis().GetNbins()):
        bin_x = h1.GetXaxis().GetBinUpEdge(bin)
        # integral includes first/last bins
        eff = h1.Integral(0, bin)/n_signal if n_signal !=0 else 0
        mis_id = h2.Integral(0, bin)/n_background if n_background !=0 else 0
        if mis_id > eff:
            # h1 is the right histogram and h2 is the left! Flip the effieincy!
            eff, mis_id = 1-eff, 1-mis_id

        cut_mass.append( bin_x )
        efficiencies.append(eff)
        mis_ids.append(mis_id)
        diff.append( abs(1-eff - mis_id) )

    optimal_idx = np.argmin( diff )

    optimal_cut = cut_mass[optimal_idx]
    optimal_eff = efficiencies[optimal_idx]
    optimal_mis_id = mis_ids[optimal_idx]

    # print(f"Cut: {round(optimal_cut, 2)} / eff: {round(optimal_eff, 2)} / mis id: {round(optimal_mis_id, 2)}")
    if debug:
        draw_optimal_cut(h1, h2, optimal_cut)

    return optimal_cut, optimal_eff, optimal_mis_id

def draw_optimal_cut(h1, h2, cut):
    '''Draw two histograms and the calculated cut wehere eff=1-misid. Used for debugging'''
    # my default pi/k/p colours
    h1, h2 = h1.Clone(), h2.Clone()
    h1.Scale(1./h1.GetEntries())
    h2.Scale(1./h2.GetEntries())

    h1_fill = h1.Clone()
    h2_fill = h2.Clone()
    for bin in range(1, h1_fill.GetXaxis().GetNbins() ):
        if h1_fill.GetBinCenter(bin) < cut:
            h1_fill.SetBinContent(bin, 0.)
        else:
            h2_fill.SetBinContent(bin, 0.)

    margin = 0.22
    ROOT.gStyle.SetPadLeftMargin(0.9*margin)
    ROOT.gStyle.SetPadRightMargin(0.1*margin)
    ROOT.gStyle.SetPadTopMargin(0.3*margin)
    ROOT.gStyle.SetPadBottomMargin(0.7*margin)
    canvas = ROOT.TCanvas(get_rand_string(),"", 600, 600)

    h1.Draw("hist")
    h1.GetXaxis().SetTitle("Mass^{2} (GeV^{2}/c^{4})")
    h1.GetYaxis().SetTitle("Normalised N entries")
    # h1.GetYaxis().SetTitleOffset(1.2)
    h1.GetYaxis().SetMaxDigits(3)
    h2.Draw("hist same")
    h1_fill.Draw("hist f same")
    h2_fill.Draw("hist f same")

    ymax = 1.05*max(h1.GetMaximum(), h2.GetMaximum())
    h1.SetMaximum(ymax)

    h1.SetLineWidth(4)
    h1.SetLineColor(pion.color)
    h1_fill.SetFillStyle(3001)
    h1_fill.SetFillColor(pion.color)
    h1_fill.SetLineColor(pion.color)

    h2.SetLineWidth(4)
    h2.SetLineColor(kaon.color)
    h2_fill.SetFillStyle(3001)
    h2_fill.SetFillColor(kaon.color)
    h2_fill.SetLineColor(kaon.color)

    line = ROOT.TLine(cut, 0., cut, ymax)
    line.SetLineColor(15)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw()
    
    canvas.Update()
    input("wait")


def calculate_pid_graphs(h1, h2):
    '''Return graphs of the separation power, efficiency, mis-id, cut values versus momentum based on two 2D histograms'''
    gr_cut = ROOT.TGraph()
    gr_eff = ROOT.TGraph()
    gr_misid = ROOT.TGraph()
    gr_sep_power = ROOT.TGraph()

    for i in range(1, h1.GetXaxis().GetNbins() ):
        x_low, x, x_up = h1.GetXaxis().GetBinLowEdge(i), h1.GetXaxis().GetBinCenter(i), h1.GetXaxis().GetBinUpEdge(i)
        h1_proj = h1.ProjectionY("h1_proj", i, i)
        h2_proj = h2.ProjectionY("h2_proj", i, i)

        # print(f"Momentum: {x_low:.2f} -- {x_up:.2f}")
        cut, eff, misid = find_optimal_cut(h1_proj, h2_proj)
        sep_power = convert_p_value_to_sep_power(1-eff)

        gr_cut.SetPoint(i-1, x, cut)
        gr_eff.SetPoint(i-1, x, eff)
        gr_misid.SetPoint(i-1, x, misid)
        gr_sep_power.SetPoint(i-1, x, sep_power)
    return gr_cut, gr_eff, gr_misid, gr_sep_power

def draw_resolution_sep_powers(graphs):
    '''Draw the separation power graphs for different TOF resolutions'''
    canvas = create_canvas(0.38, 0.4, 0.65)
    canvas.SetTicky(False)

    legend = create_legend(0.4, 0.63, 0.98, 0.94)
    for i, (res, gr) in enumerate(graphs.items()):
        if i == 0:
            # canvas.DrawFrame(0., MIN_SEP_POWER, 19., MAX_SEP_POWER)
            gr.Draw("ALP")
            gr.GetYaxis().SetTitleOffset(0.9)
            gr.GetXaxis().SetRangeUser(0, 19)
            gr.GetXaxis().SetNdivisions(506)
            gr.GetYaxis().SetRangeUser(MIN_SEP_POWER, MAX_SEP_POWER)
            canvas.Modified()
            canvas.Update()
            # draw an axis on the right side
            x_pos = canvas.GetUxmax()
            axis_eff = ROOT.TGaxis(x_pos, canvas.GetUymin(), x_pos-0.001, canvas.GetUymax(), MIN_SEP_POWER, MAX_SEP_POWER)
            axis_eff.SetTitleColor( ROOT.gStyle.GetTitleColor("Y") )
            axis_eff.SetTitleFont( ROOT.gStyle.GetTitleFont("Y") )
            axis_eff.SetTitleSize( ROOT.gStyle.GetTitleSize("Y") )
            axis_eff.CenterTitle(True)
            axis_eff.SetTitleOffset(2.2)
            axis_eff.SetTitle("Efficiency (%)")
            axis_eff.SetLabelColor( ROOT.gStyle.GetLabelColor("Y") )
            axis_eff.SetLabelFont( ROOT.gStyle.GetLabelFont("Y") )
            axis_eff.SetLabelOffset(-0.14)
            axis_eff.SetLabelSize( ROOT.gStyle.GetLabelSize("Y") )
            axis_eff.SetTickLength(0.03)
            for j in range(MIN_SEP_POWER, MAX_SEP_POWER+1):
                axis_eff.ChangeLabel(j+1, -1, -1, -1, ROOT.kBlack, -1, f"{convert_sep_power_to_eff(j)*100:.2f}")
            axis_eff.DrawClone()
        else:
            gr.Draw("LPsame")
        gr.SetLineColor(COLORS_RESOLUTION[i])
        gr.SetMarkerColor(COLORS_RESOLUTION[i])
        gr.SetLineWidth(4)
        gr.SetMarkerStyle(20)
        legend.AddEntry(gr, f"{res} ps","lp")

    legend.DrawClone()

    # latex = ROOT.TLatex()
    # latex.SetTextFont(52)
    # latex.DrawLatex(12, 6.06, "ILD preliminary")

    canvas.Update()
    return canvas

def draw_dedx_sep_powers(gr_tof, gr_dedx):
    '''Draw the separation power graphs for TOF and dEdx'''
    canvas = create_canvas(0.38, 0.4, 0.65)
    canvas.SetLogx()
    canvas.SetTicky(False)

    legend = create_legend(0.35, 0.68, 0.74, 0.84)
    gr_tof.Draw("ALP")
    gr_tof.GetXaxis().SetRangeUser(0.5, 20)
    gr_tof.GetXaxis().SetMoreLogLabels()
    gr_tof.GetXaxis().SetNoExponent()
    gr_tof.GetXaxis().SetNdivisions(506)
    gr_tof.GetYaxis().SetTitleOffset(0.9)
    gr_tof.GetYaxis().SetRangeUser(MIN_SEP_POWER, MAX_SEP_POWER)
    gr_tof.SetLineColor(COLORS_DEDX[0])
    gr_tof.SetMarkerColor(COLORS_DEDX[0])
    gr_tof.SetLineWidth(4)
    gr_tof.SetMarkerStyle(20)
    legend.AddEntry(gr_tof, gr_tof.GetTitle(),"lp")
    canvas.Modified()
    canvas.Update()
    # draw an axis on the right side
    # NOTE: FREAKING ROOT BUG https://root-forum.cern.ch/t/getumin-getumax-show-wrong-results-for-the-canvases-with-the-log-scale/58867
    x_pos = 10**canvas.GetUxmax()
    axis_eff = ROOT.TGaxis(x_pos, canvas.GetUymin(), x_pos-0.001, canvas.GetUymax(), MIN_SEP_POWER, MAX_SEP_POWER)
    axis_eff.SetTitleColor( ROOT.gStyle.GetTitleColor("Y") )
    axis_eff.SetTitleFont( ROOT.gStyle.GetTitleFont("Y") )
    axis_eff.SetTitleSize( ROOT.gStyle.GetTitleSize("Y") )
    axis_eff.CenterTitle(True)
    axis_eff.SetTitleOffset(2.2)
    axis_eff.SetTitle("Efficiency (%)")
    axis_eff.SetLabelColor( ROOT.gStyle.GetLabelColor("Y") )
    axis_eff.SetLabelFont( ROOT.gStyle.GetLabelFont("Y") )
    axis_eff.SetLabelOffset(-0.14)
    axis_eff.SetLabelSize( ROOT.gStyle.GetLabelSize("Y") )
    axis_eff.SetTickLength(0.03)
    for j in range(MIN_SEP_POWER, MAX_SEP_POWER+1):
        axis_eff.ChangeLabel(j+1, -1, -1, -1, ROOT.kBlack, -1, f"{convert_sep_power_to_eff(j)*100:.2f}")
    axis_eff.DrawClone()
    gr_dedx.Draw("LPsame")
    gr_dedx.SetLineColor(COLORS_DEDX[2])
    gr_dedx.SetMarkerColor(COLORS_DEDX[2])
    gr_dedx.SetLineWidth(4)
    gr_dedx.SetMarkerStyle(20)
    legend.AddEntry(gr_dedx, gr_dedx.GetTitle(),"lp")

    legend.DrawClone()

    # latex = ROOT.TLatex()
    # latex.SetTextFont(52)
    # latex.DrawLatex(12, 6.06, "ILD preliminary")

    canvas.Update()
    return canvas


def find_optimal_cut_2d(h1_2d, h2_2d):
    '''Return'''
    # I will try tomorrow to feed two 2D histograms and find optimal line cut
    # return optimal_cut, optimal_eff, optimal_mis_id
    pass

def main():

    df_init = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna_Zqq.root").Filter("quarksToPythia == 44")
    df = df_init.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)").Filter("abs(pdg) == 321")
    h1 = df.Histo1D((get_rand_string(), "Produced;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h2 = df.Filter("isReconstructed").Histo1D((get_rand_string(), "Reconstructed;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h3 = df.Filter("isReconstructed && hasTrack").Histo1D((get_rand_string(), "Reconstructed && Track;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h4 = df.Filter("isReconstructed && hasTrack && hasShower").Histo1D((get_rand_string(), "Reconstructed && Track && Shower;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h5 = df.Filter("isReconstructed && hasTrack && hasShower && !isSimulated").Histo1D((get_rand_string(), "Reconstructed && Track && Shower && not sim;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h6 = df.Filter("isReconstructed && hasTrack && hasShower && !isSimulated && !isOverlay").Histo1D((get_rand_string(), "Reconstructed && Track && Shower && not sim && not overlay;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h7 = df.Filter("isReconstructed && hasTrack && hasShower && !isSimulated && !isOverlay && isInTrueSecondaryVertex").Histo1D((get_rand_string(), "Reconstructed && Track && Shower && not sim && not overlay && isInTrueSecondaryVertex ;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )
    h8 = df.Filter("isReconstructed && hasTrack && hasShower && !isSimulated && !isOverlay && isInTrueSecondaryVertex && isInRecoSecondaryVertex").Histo1D((get_rand_string(), "Reconstructed && Track && Shower && not sim && not overlay && isInTrueSecondaryVertex ;Momentum (GeV/c);N entries", 1000, 0, 20), "mom" )


    h2.SetLineColor(ROOT.kRed+2)
    h3.SetLineColor(3)
    h4.SetLineColor(4)
    h5.SetLineColor(5)
    h6.SetLineColor(6)
    h7.SetLineColor(7)
    h8.SetLineColor(8)

    h1.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h4.Draw("same")
    h5.Draw("same")
    h6.Draw("same")
    h7.Draw("same")
    h8.Draw("same")
    input("wait")




    # Get all histos first to utilise lazy RDataFrame behaviour
    histos = {"TOF" : {}, "dEdx" : {}, "Combined30" : {}}
    for RES in RESOLUTIONS:
        histos["TOF"][RES] = {}
        df = df_init.Define("beta", f"{TRACK_LENGTH_COLUMN}/(tofClosest{RES}*299.792458)")\
                    .Define("mass2", f"{MOMENTUM_COLUMN}*{MOMENTUM_COLUMN}*( 1./(beta*beta) - 1.)")
        for p in particles:
            histos["TOF"][RES][p] = {"lin": {}, "log" : {}}
            df_p = df.Filter(f"abs(pdg) == {p.pdg}")
            h_lin = df_p.Histo2D((get_rand_string(), "", N_MOMENTUM_BINS, MIN_MOMENTUM, MAX_MOMENTUM, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")
            h_log = df_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")
            histos["TOF"][RES][p] = {"lin": h_lin, "log" : h_log}
    for p in particles:
        df_p = df_init.Filter(f"abs(pdg) == {p.pdg}").Filter("dEdx > 0.")
        # h = df_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, MIN_DEDX, MAX_DEDX), MOMENTUM_COLUMN, "dEdx")
        h = df_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), MOMENTUM_COLUMN, "dEdx")
        histos["dEdx"][p] = h


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