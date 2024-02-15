import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)
ROOT.EnableImplicitMT()
from utils import *

def get_pdg_histogram(df, pdg="all", y_name="mass2"):
    # with dE/dx we are better use log binning
    n_x_bins = 30
    x_min = -0.3 # 0.1 GeV
    x_max = 1.3 # ~ 20 GeV
    mom_bins = [ 10**(x_min + (x_max-x_min)*i/n_x_bins ) for i in range(n_x_bins+1) ]

    if y_name == "mass2":
        n_y_bins, y_min, y_max = 3000, -3, 3
        y_bins = [ y_min + (y_max - y_min)*i/n_y_bins for i in range(n_y_bins+1) ]
    elif y_name == "dEdx":
        n_y_bins, y_min, y_max = 3000, 0, 1e-6
        y_bins = [ y_min + (y_max - y_min)*i/n_y_bins for i in range(n_y_bins+1) ]


    df_filtered = df.Filter(f"abs(pdg) == {pdg}") if pdg != "all" else df
    h = df_filtered.Histo2D( (f"h_{y_name}_{pdg}", "",n_x_bins, np.array(mom_bins), n_y_bins, np.array(y_bins) ), "harmonicMomToEcal_IKF_zedLambda",y_name)
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
        cut, eff, misid = analyse_optimal_cut(h1_proj, h2_proj)
        sep_power = get_separation_power(1-eff)

        gr_cut.SetPoint(i-1, x, cut)
        gr_eff.SetPoint(i-1, x, eff)
        gr_misid.SetPoint(i-1, x, misid)
        gr_sep_power.SetPoint(i-1, x, sep_power)
    return gr_cut, gr_eff, gr_misid, gr_sep_power


colors = ['#1b9e77', '#d95f02', '#7570b3']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

resolution = 30 # ps
df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
                .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
                .Filter("tofClosest0 > 6. && dEdx > 0.")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/(tofClosest{resolution}*299.792458)")\
                .Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")

graphs_sep_power_pik = {}
graphs_sep_power_kp = {}


h_all = get_pdg_histogram(df, "all", "mass2")
h_pi = get_pdg_histogram(df, "211", "mass2")
h_k = get_pdg_histogram(df, "321", "mass2")
h_p = get_pdg_histogram(df, "2212", "mass2")

h_all_dedx = get_pdg_histogram(df, "all", "dEdx")
h_pi_dedx = get_pdg_histogram(df, "211", "dEdx")
h_k_dedx = get_pdg_histogram(df, "321", "dEdx")
h_p_dedx = get_pdg_histogram(df, "2212", "dEdx")

_, _, _, gr_sp_tof = analyse_pid(h_pi, h_k)
_, _, _, gr_sp_dedx = analyse_pid(h_pi_dedx, h_k_dedx)
gr_sp_combined = ROOT.TGraph()
gr_sp_combined.SetTitle("; momentum (GeV/c); #pi/K separation power")
for i in range( gr_sp_tof.GetN() ):
    gr_sp_combined.SetPoint(i, gr_sp_tof.GetPointX(i), np.sqrt(gr_sp_tof.GetPointY(i)*gr_sp_tof.GetPointY(i) + gr_sp_dedx.GetPointY(i)*gr_sp_dedx.GetPointY(i) ) )

canvas = create_canvas()
canvas.SetLogx()

gr_sp_combined.Draw("ACP")
gr_sp_combined.GetXaxis().SetRangeUser(0.5, 20)
gr_sp_combined.GetXaxis().SetNoExponent()
gr_sp_combined.GetXaxis().SetMoreLogLabels()
gr_sp_combined.GetYaxis().SetRangeUser(0, 6)

gr_sp_dedx.Draw("CPsame")
gr_sp_tof.Draw("CPsame")
gr_sp_tof.SetLineColor( ROOT.TColor.GetColor("#0077b6") )
gr_sp_dedx.SetLineColor( ROOT.kRed+2 )
gr_sp_combined.SetLineColor(ROOT.kViolet-3)
gr_sp_tof.SetLineWidth(3)
gr_sp_dedx.SetLineWidth(3)
gr_sp_combined.SetLineWidth(3)
gr_sp_tof.SetMarkerColor( ROOT.TColor.GetColor("#0077b6") )
gr_sp_dedx.SetMarkerColor( ROOT.kRed+2 )
gr_sp_combined.SetMarkerColor(ROOT.kViolet-3)
gr_sp_tof.SetMarkerStyle(20)
gr_sp_dedx.SetMarkerStyle(20)
gr_sp_combined.SetMarkerStyle(20)
canvas.Update()


legend = create_legend()
legend.AddEntry(gr_sp_tof, f"TOF {resolution} ps")
legend.AddEntry(gr_sp_dedx, "dE/dx")
legend.AddEntry(gr_sp_combined, "TOF #oplus dE/dx")
legend.Draw()
canvas.Update()
input("wait")