import ROOT
import numpy as np

from my_utilities import *

ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
# ROOT.EnableImplicitMT()

def extract_data(df, tof_column="tofClosest0"):
    df = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")
    df = df.Define("inverseBeta", "1/beta")
    df = df.Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")
    return df
    # df_mass = df.Filter("mass2 > 0.").Define("mass", "sqrt(mass2)")
    # return df_mass

def get_pdg_histogram(df, pdg="all", y_name="mass2"):
    n_x_bins, x_min, x_max = 70, 0, 20
    if y_name == "mass2":
        n_y_bins, y_min, y_max = 3000, -10, 10
    elif y_name == "beta":
        n_y_bins, y_min, y_max = 3000, 0, 2

    df_filtered = df.Filter(f"abs(pdg) == {pdg}") if pdg != "all" else df
    h = df_filtered.Histo2D((f"h_{pdg}", "", n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "harmonicMomToEcal_IKF_zedLambda",y_name)
    return h

def analyse_optimal_cut(h1, h2, debug=False):
    '''
        h1, h2 - 1D histograms in a given momentum slice
    '''
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

    print(f"Cut: {round(optimal_cut, 2)} / eff: {round(optimal_eff, 2)} / mis id: {round(optimal_mis_id, 2)}")
    if debug:
        draw_optimal_cut(h1, h2, optimal_cut)

    return optimal_cut, optimal_eff, optimal_mis_id

def draw_optimal_cut(h1, h2, cut):
    # my default pi/k/p colours
    colors = ['#1b9e77', '#d95f02', '#7570b3']
    colors = [ ROOT.TColor.GetColor(c) for c in colors]

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
    h1.SetLineColor(colors[0])
    h1_fill.SetFillStyle(3001)
    h1_fill.SetFillColor(colors[0])
    h1_fill.SetLineColor(colors[0])

    h2.SetLineWidth(4)
    h2.SetLineColor(colors[1])
    h2_fill.SetFillStyle(3001)
    h2_fill.SetFillColor(colors[1])
    h2_fill.SetLineColor(colors[1])

    line = ROOT.TLine(cut, 0., cut, ymax)
    line.SetLineColor(15)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw()
    
    canvas.Update()
    input("wait")

def get_separation_power(p_value):
    '''p-value == 1 - efficiency'''
    return 2*ROOT.Math.gaussian_quantile_c(p_value, 1)


def analyse_pid(h1, h2):
    '''
        h1,h2 - 2D histogram of something vs momentum
    '''
    gr_cut = ROOT.TGraph()
    gr_eff = ROOT.TGraph()
    gr_misid = ROOT.TGraph()
    gr_sep_power = ROOT.TGraph()

    for i in range(1, h1.GetXaxis().GetNbins() ):
        x_low, x, x_up = h1.GetXaxis().GetBinLowEdge(i), h1.GetXaxis().GetBinCenter(i), h1.GetXaxis().GetBinUpEdge(i)
        h1_proj = h1.ProjectionY("h1_proj", i, i)
        h2_proj = h2.ProjectionY("h2_proj", i, i)

        print(f"Momentum: {x_low:.2f} -- {x_up:.2f}")
        cut, eff, misid = analyse_optimal_cut(h1_proj, h2_proj, debug=True)
        sep_power = get_separation_power(1-eff)

        gr_cut.SetPoint(i-1, x, cut)
        gr_eff.SetPoint(i-1, x, eff)
        gr_misid.SetPoint(i-1, x, misid)
        gr_sep_power.SetPoint(i-1, x, sep_power)
    return gr_cut, gr_eff, gr_misid, gr_sep_power

def draw_sep_powers(graphs):
    # colors = ["#03045e","#023e8a","#0077b6","#0096c7","#00b4d8","#48cae4"]
    # colors = ["#690000", "#850e0f", "#a21d19", "#c02b25", "#df3831", "#ff463d"]
    colors = ['#0444b3', '#0068cc', '#008add', '#1caaea', '#61caf4', '#98e8ff']
    colors = [ ROOT.TColor.GetColor(c) for c in colors]

    margin = 0.22
    ROOT.gStyle.SetPadLeftMargin(0.8*margin)
    ROOT.gStyle.SetPadRightMargin(0.2*margin)
    ROOT.gStyle.SetPadTopMargin(0.3*margin)
    ROOT.gStyle.SetPadBottomMargin(0.7*margin)
    canvas = ROOT.TCanvas(get_rand_string(),"", 600, 600)
    canvas.SetGridx(True)
    canvas.SetGridy(True)

    legend = ROOT.TLegend(0.4, 0.63, 0.98, 0.94)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for i, (res, gr) in enumerate(graphs.items()):
        if i == 0:
            gr.Draw("ALP")
            gr.SetTitle(";Momentum (GeV/c); Z")
            gr.GetXaxis().SetRangeUser(0, 19)
            gr.GetYaxis().SetRangeUser(0, 6)
            gr.GetYaxis().SetTitleOffset(1.1)
        else:
            gr.Draw("LPsame")
        gr.SetLineColor(colors[i])
        gr.SetMarkerColor(colors[i])
        gr.SetLineWidth(4)
        gr.SetMarkerStyle(20)
        legend.AddEntry(gr, f"{res} ps","lp")

    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextFont(52)
    latex.DrawLatex(12, 6.06, "ILD preliminary")

    canvas.Update()
    return canvas, legend


def main():
    df_init = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
                  .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
                  .Filter("tofClosest0 > 6.")
    graphs_sep_power_pik = {}
    graphs_sep_power_kp = {}

# [0, 20, 40, 60, 80, 100]
    for res in [30]:
        df = extract_data(df_init, tof_column=f"tofClosest{res}")

        y_name = "mass2"
        h_all = get_pdg_histogram(df, "all", y_name)
        h_pi = get_pdg_histogram(df, "211", y_name)
        h_k = get_pdg_histogram(df, "321", y_name)
        h_p = get_pdg_histogram(df, "2212", y_name)

        # c_all = draw_2d_plot(h_all, 1e6)
        # c_pi=  draw_2d_plot(h_pi, 1e6)
        # c_k = draw_2d_plot(h_k, 1e6)
        # c_p = draw_2d_plot(h_p, 1e6)
        # input("wait")

        gr_cut_pik, gr_eff_pik, gr_misid_pik, gr_sep_power_pik = analyse_pid(h_pi, h_k)
        gr_cut_kp, gr_eff_kp, gr_misid_kp, gr_sep_power_kp = analyse_pid(h_k, h_p)

        graphs_sep_power_pik[res] = gr_sep_power_pik
        graphs_sep_power_kp[res] = gr_sep_power_kp

    c1, leg1 = draw_sep_powers(graphs_sep_power_pik)
    input("wait")
    c2, leg2 = draw_sep_powers(graphs_sep_power_kp)
    input("wait")

main()