import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

ROOT.gStyle.SetCanvasPreferGL(True)
plt.rcParams.update({'font.size':18, 'text.usetex':True})

particles = {}
particles['pion'] = {'mass' : 0.13957039, 'color': '#1b9e77', 'legend' : "#pi^{#pm}"}
particles['kaon'] = {'mass' : 0.493677, 'color': '#d95f02', 'legend' : "K^{#pm}"}
particles['proton'] = {'mass' : 0.93827208816, 'color': '#7570b3', 'legend' : "p "}
SPEED_OF_LIGHT = 299.792458 # mm / ns

def get_mass(true_mass, momentum, track_length, relative_momentum_resolution, relative_length_resolution, tof_resolution):
    '''Return mass of the particle in momentum bins for the plot'''

    # Deduce true TOF
    tof = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + true_mass*true_mass/(momentum*momentum) ) # ns

    # smear everything
    smeared_momentum = momentum*(1. + relative_momentum_resolution)
    smeared_track_length = track_length*(1. + relative_length_resolution)
    smeared_tof = tof + tof_resolution

    smeared_beta = smeared_track_length/(smeared_tof*SPEED_OF_LIGHT)
    mass = smeared_momentum*np.sqrt( 1/(smeared_beta*smeared_beta) - 1) # GeV

    return np.nan_to_num(mass, nan=0), tof

def get_filled_graph(x, y_min, y_max):
    n_points = len(x)
    graph = ROOT.TGraph(2*n_points)
    for i in range(n_points):
        graph.SetPoint(i, x[i], y_max[i])
        graph.SetPoint(n_points + i, x[n_points-i-1], y_min[n_points-i-1])
    return graph


def get_mass_plot():
    # We always work in momentum bins
    momentum = np.arange(0.1, 8, 0.01)

    # # STAR DATA from the paper https://arxiv.org/abs/nucl-ex/0308022
    track_length = 2200. # mm at eta = 1.
    relative_momentum_resolution = 0.013
    relative_length_resolution = 0.002
    tof_resolution = 0.100 # ns

    text= r''' \textbf{STAR apparatus}
    $ L = 2200 \ \mathrm{mm}$
    $ \Delta p/p = 0.013 $
    $ \Delta L / L = 0.002 $
    $ \Delta T = 100 \ \mathrm{ps} $'''

    # ILD roughly or STAR with modern timing
    # track_length = 2200. # mm
    # relative_momentum_resolution = 0.013
    # relative_length_resolution = 0.002
    # tof_resolution = 0.020 # ns

    # text= r''' \textbf{STAR apparatus}
    # \textbf{modern timing}
    # $ L = 2200 \ \mathrm{mm}$
    # $ \Delta p/p = 0.013 $
    # $ \Delta L / L = 0.002 $
    # $ \Delta T = 20 \ \mathrm{ps} $'''


    #pions
    mass_bands = {}
    graphs = {}
    canvas = ROOT.TCanvas("name", "title")
    for name in particles:
        mass_bands[name, 'true'], _ = get_mass(particles[name]['mass'], momentum, track_length, 0., 0., 0.)
        graphs[name, 'true'] = ROOT.TGraph(len(momentum), momentum, mass_bands[name, 'true'])

        mass_bands[name, 'up_dpdl'], _ = get_mass(particles[name]['mass'], momentum, track_length, relative_momentum_resolution, -relative_length_resolution, 0.)
        mass_bands[name, 'down_dpdl'], _ = get_mass(particles[name]['mass'], momentum, track_length, -relative_momentum_resolution, relative_length_resolution, 0.)
        graphs[name, 'dl_band'] = get_filled_graph(momentum, mass_bands[name, 'down_dpdl'], mass_bands[name, 'up_dpdl'])

        mass_bands[name, 'up_total'], _ = get_mass(particles[name]['mass'], momentum, track_length, relative_momentum_resolution, -relative_length_resolution, tof_resolution)
        mass_bands[name, 'down_total'], _ = get_mass(particles[name]['mass'], momentum, track_length, -relative_momentum_resolution, relative_length_resolution, -tof_resolution)
        graphs[name, 'total_band'] = get_filled_graph(momentum, mass_bands[name, 'down_total'], mass_bands[name, 'up_total'])
        graphs[name, 'total_band'].SetFillColorAlpha(ROOT.kRed, 0.2)

        graphs[name, 'true'].Draw("AL" if name == 'pion' else "Lsame")
        graphs[name, 'true'].SetLineColor( ROOT.TColor.GetColor(particles[name]['color']) )
        graphs[name, 'dl_band'].Draw("fsame")
        graphs[name, 'dl_band'].SetFillColorAlpha(ROOT.TColor.GetColor(particles[name]['color']), 0.2)
        graphs[name, 'total_band'].Draw("fsame")
        graphs[name, 'total_band'].SetFillColorAlpha(ROOT.TColor.GetColor(particles[name]['color']), 0.2)
        if name == 'pion':
            graphs[name, 'true'].GetXaxis().SetTitle("Momentum (GeV/c)")
            graphs[name, 'true'].GetYaxis().SetTitle("Mass (GeV/c^{2})")
            graphs[name, 'true'].GetXaxis().SetRangeUser(0, 8)
            graphs[name, 'true'].GetYaxis().SetRangeUser(0, 3.)
        
    canvas.Update()
    input("wait")

    # canvas.SetTopMargin(0.12)
    # canvas.SetBottomMargin(0.14)
    # canvas.SetLeftMargin(0.18)
    # canvas.SetRightMargin(0.07)


def mass_vs_dtdl(x_axis='dt'):
    momentum = 0.8 # GeV
    track_length = 2000. # mm

    if x_axis == 'dt':
        x_min, x_max, x_title = -100, 100, '#Delta T (ps)'
        legend = ROOT.TLegend(0.2, 0.6, 0.59, 0.76)
        # tof_res is expected to be in ns as an input to get_mass()
        mom_res, len_res, tof_res = 0., 0., np.arange(x_min/1000., x_max/1000. + 0.001, 0.001) # ns
    elif x_axis == 'dl':
        x_min, x_max, x_title = -30, 30, '#Delta L (mm)'
        legend = ROOT.TLegend(0.51, 0.6, 0.9, 0.76)
        # len_res is expected to be relative as an input to get_mass()
        mom_res, len_res, tof_res = 0., np.arange(x_min, x_max + 0.1, 0.1)/track_length, 0.

    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetMargin(0.15)
    legend.SetTextAlign(32)
    legend.SetTextFont(42)
    canvas = ROOT.TCanvas("name", "title")
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.14)
    canvas.SetLeftMargin(0.18)
    canvas.SetRightMargin(0.07)

    graphs = {}
    lines = {}

    tof_true_average = 0.
    for name in particles:
        mass, tof_true = get_mass(particles[name]['mass'], momentum, track_length, mom_res, len_res, tof_res)
        tof_true_average += tof_true/3.
        if x_axis == 'dt':
            graphs[name] = ROOT.TGraph(len(tof_res), tof_res*1000, mass)
        elif x_axis == 'dl':
            graphs[name] = ROOT.TGraph(len(len_res), len_res*track_length, mass)

        graphs[name].SetLineColor(ROOT.TColor.GetColor(particles[name]['color']))
        graphs[name].SetLineWidth(3)
        graphs[name].Draw("AL" if name =='pion' else "Lsame")
        if name=='pion':
            graphs[name].GetYaxis().SetTitle('Mass (GeV/c^{2})')
            graphs[name].GetYaxis().SetTitleOffset(1.2)
            graphs[name].GetYaxis().SetRangeUser(-0.1, 1.65)
            graphs[name].GetXaxis().SetTitle(x_title)
            graphs[name].GetXaxis().SetRangeUser(x_min, x_max)
            

        lines[name] = ROOT.TLine(x_min, particles[name]['mass'], x_max, particles[name]['mass'])
        lines[name].SetLineColor(ROOT.TColor.GetColor(particles[name]['color']))
        lines[name].SetLineWidth(2)
        lines[name].SetLineStyle(9)
        lines[name].Draw()
        align_len = 10-len(f'{tof_true:.2f}')
        legend.AddEntry(graphs[name], "#font[62]{" + f"{particles[name]['legend']}" + "}  T_{true}:" + f"{tof_true:>{align_len}.2f} ns", "l")
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    momentum_text = f"p = {momentum}" + " #frac{GeV}{c}"
    track_length_text = f"L = {track_length:.0f}" + " mm"
    if x_axis == 'dt':
        latex.DrawLatex(-70, 1.45, momentum_text)
        latex.DrawLatex(11, 1.45, track_length_text)
    elif x_axis == 'dl':
        latex.DrawLatex(-20, 1.45, momentum_text)
        latex.DrawLatex(4, 1.45, track_length_text)


    canvas.Update()
    if x_axis == 'dt':
        top_axis = ROOT.TGaxis(canvas.GetUxmin(),canvas.GetUymax(),
                            canvas.GetUxmax(), canvas.GetUymax(),x_min/(tof_true_average*1000)*100,x_max/(tof_true_average*1000)*100,510,"-L")
        top_axis.SetTitle("#delta T (%)")
        top_axis.Draw()
    else:
        top_axis = ROOT.TGaxis(canvas.GetUxmin(),canvas.GetUymax(),
                            canvas.GetUxmax(), canvas.GetUymax(),x_min/track_length*100,x_max/track_length*100,510,"-L")
        top_axis.SetTitle("#delta L (%)")
    top_axis.SetLabelSize(0.06)
    top_axis.SetTitleSize(0.07)
    top_axis.SetLabelFont(42)
    top_axis.SetTitleFont(42)
    top_axis.SetTitleOffset(0.78)
    top_axis.Draw()

    canvas.Update()
    input("wait")


get_mass_plot()
# mass_vs_dtdl('dl')