import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from my_utilities import *

ROOT.gStyle.SetCanvasPreferGL(True)
plt.rcParams.update({'font.size':18, 'text.usetex':True})

class Particle:
    def __init__(self, name, legend, color, mass, momentum, track_length):
        self.name = name
        self.legend = legend
        self.color = color
        self.mass = mass
        self.momentum = momentum
        self.track_length = track_length
        self.tof = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + (mass/momentum)**2 ) if momentum !=0 else -1
        self.beta = track_length/(SPEED_OF_LIGHT * self.tof) if self.tof != -1 else -1
        self.gamma = 1/np.sqrt(1 - self.beta**2)

def get_filled_graph(x, y_min, y_max):
    n_points = len(x)
    graph = ROOT.TGraph(2*n_points)
    for i in range(n_points):
        graph.SetPoint(i, x[i], y_max[i])
        graph.SetPoint(n_points + i, x[n_points-i-1], y_min[n_points-i-1])
    return graph

momentum, track_length = 4., 2000.
pion = Particle('pion', '#pi^{#pm}', ROOT.TColor.GetColor('#1b9e77'), 0.13957039, momentum, track_length)
kaon = Particle('kaon', 'K^{#pm}', ROOT.TColor.GetColor('#d95f02'), 0.493677, momentum, track_length)
proton = Particle('proton', 'p', ROOT.TColor.GetColor('#7570b3'), 0.93827208816, momentum, track_length)

particles = [pion, kaon, proton]

margin = 0.22
ROOT.gStyle.SetPadLeftMargin(0.8*margin)
ROOT.gStyle.SetPadRightMargin(0.2*margin)
ROOT.gStyle.SetPadTopMargin(0.3*margin)
ROOT.gStyle.SetPadBottomMargin(0.7*margin)
canvas = ROOT.TCanvas(get_rand_string(),"", 600, 600)



def get_appox_uncertainty(particle, relative_momentum_resolution, relative_length_resolution, tof_res_in_ns):
    dm = particle.mass * np.sqrt(relative_momentum_resolution**2 + (relative_length_resolution**2 + (tof_res_in_ns/particle.tof)**2)*particle.gamma**4)
    return np.nan_to_num(particle.mass + dm, nan=0), np.nan_to_num(particle.mass - dm, nan=0)


def get_mass(particle, relative_momentum_resolution, relative_length_resolution, tof_res_in_ns):
    '''Return mass of the particle in momentum bins for the plot'''

    # smear everything
    smeared_momentum = particle.momentum*(1. + relative_momentum_resolution)
    smeared_track_length = particle.track_length*(1. + relative_length_resolution)
    smeared_tof = particle.tof + tof_res_in_ns

    smeared_beta = smeared_track_length/(smeared_tof*SPEED_OF_LIGHT)
    mass2 = smeared_momentum*smeared_momentum*(1./(smeared_beta*smeared_beta) - 1.) # GeV^2/c^4
    mass = np.sqrt(mass2) # GeV/c^2
    return np.nan_to_num(mass, nan=0)



def get_mass_plot():
    # We always work in momentum bins
    momentum = np.arange(0.1, 8, 0.01)

    # # STAR DATA from the paper https://arxiv.org/abs/nucl-ex/0308022
    track_length = 2200. # mm at eta = 1.
    relative_momentum_resolution = 0.013
    relative_length_resolution = 0.002
    tof_resolution = 0.100 # ns
    # tof_resolution = 0.020 # ns


    mass_bands = {}
    graphs = {}
    canvas = ROOT.TCanvas("name", "title")
    canvas.SetTopMargin(0.02)
    canvas.SetBottomMargin(0.16)
    canvas.SetLeftMargin(0.18)
    canvas.SetRightMargin(0.02)

    fill_alpha = .3

    legend = ROOT.TLegend(0.2, 0.72, 0.66, 0.96)
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetMargin(0.4)
    # legend.SetColumnSeparation(0.25)
    # legend.SetTextAlign(32)
    legend.SetTextFont(42)

    for name in particles:
        mass_bands[name, 'true'] = get_mass(particles[name]['mass'], momentum, track_length, 0., 0., 0.)["mass"]
        mass_bands[name, 'up_dpdl'] = get_mass(particles[name]['mass'], momentum, track_length, relative_momentum_resolution, -relative_length_resolution, 0.)["mass"]
        mass_bands[name, 'down_dpdl'] = get_mass(particles[name]['mass'], momentum, track_length, -relative_momentum_resolution, relative_length_resolution, 0.)["mass"]
        mass_bands[name, 'up_total'] = get_mass(particles[name]['mass'], momentum, track_length, relative_momentum_resolution, -relative_length_resolution, tof_resolution)["mass"]
        mass_bands[name, 'down_total'] = get_mass(particles[name]['mass'], momentum, track_length, -relative_momentum_resolution, relative_length_resolution, -tof_resolution)["mass"]

        graphs[name, 'true'] = ROOT.TGraph(len(momentum), momentum, mass_bands[name, 'true'])
        graphs[name, 'true'].Draw("AL" if name == 'pion' else "Lsame")
        graphs[name, 'true'].SetLineColor( particles[name]['color'] )
        graphs[name, 'true'].SetLineWidth(3)

        graphs[name, 'dl_band'] = get_filled_graph(momentum, mass_bands[name, 'down_dpdl'], mass_bands[name, 'up_dpdl'])
        graphs[name, 'dl_band'].SetFillColorAlpha(particles[name]['color'], fill_alpha)
        graphs[name, 'dl_band'].Draw("fsame")

        graphs[name, 'total_band'] = get_filled_graph(momentum, mass_bands[name, 'down_total'], mass_bands[name, 'up_total'])
        graphs[name, 'total_band'].SetFillColorAlpha(particles[name]['color'], fill_alpha)
        graphs[name, 'total_band'].Draw("fsame")
        if name == 'pion':
            graphs[name, 'true'].GetXaxis().SetTitle("Momentum (GeV/c)")
            graphs[name, 'true'].GetYaxis().SetTitle("Mass (GeV/c^{2})")
            graphs[name, 'true'].GetYaxis().SetTitleOffset(1.2)
            graphs[name, 'true'].GetXaxis().SetRangeUser(0, 8)
            graphs[name, 'true'].GetYaxis().SetRangeUser(0, 2.3)
        

    # pi
    gr_legend_pi = ROOT.TGraph()
    gr_legend_pi.SetFillColor(particles['pion']['color'])
    legend.AddEntry(gr_legend_pi, "#font[62]{#pi}", "f")

    # true mass
    gr_legend_true = ROOT.TGraph()
    gr_legend_true.SetLineWidth(3)
    gr_legend_true.SetLineColor(ROOT.kBlack)
    legend.AddEntry(gr_legend_true, "true mass", "l")

    # K
    gr_legend_k = ROOT.TGraph()
    gr_legend_k.SetFillColor(particles['kaon']['color'])
    legend.AddEntry(gr_legend_k, "#font[62]{K}", "f")

    # DL
    gr_legend_dl = ROOT.TGraph()
    gr_legend_dl.SetFillColorAlpha(ROOT.kBlack, 2*fill_alpha)
    legend.AddEntry(gr_legend_dl, "#DeltaL", "f")

    # p
    gr_legend_p = ROOT.TGraph()
    gr_legend_p.SetFillColor(particles['proton']['color'])
    legend.AddEntry(gr_legend_p, "#font[62]{p}", "f")

    # DL DT
    gr_legend_dldt = ROOT.TGraph()
    gr_legend_dldt.SetFillColor(ROOT.kBlack)
    gr_legend_dldt.SetFillColorAlpha(ROOT.kBlack, fill_alpha)
    legend.AddEntry(gr_legend_dldt, "#DeltaL \u2295 #DeltaT", "f")

    legend.Draw()

    latex = ROOT.TLatex()
    # latex = ROOT.TMathText()
    latex.SetTextFont(42)
    latex.DrawLatex(4.65, 2.08, f"L = {track_length:.0f}" + " mm")
    latex.DrawLatex(4.65, 1.91, "#Deltap/p = " + f"{relative_momentum_resolution:.3f}")
    latex.DrawLatex(4.65, 1.74, "#DeltaL/L = " + f"{relative_length_resolution:.3f}")
    latex.DrawLatex(4.65, 1.57, "#DeltaT = " + f"{1000*tof_resolution:.0f}" + " ps")
    canvas.Update()
    input("wait")

def mass_vs_dtdl(x_axis='dt'):
    momentum = 4 # GeV
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
        vars = get_mass(particles[name]['mass'], momentum, track_length, mom_res, len_res, tof_res)
        mass, tof_true = vars['mass'], vars['tof_true']
        tof_true_average += tof_true/3.
        if x_axis == 'dt':
            graphs[name] = ROOT.TGraph(len(tof_res), tof_res*1000, mass)
        elif x_axis == 'dl':
            graphs[name] = ROOT.TGraph(len(len_res), len_res*track_length, mass)

        graphs[name].SetLineColor(particles[name]['color'])
        graphs[name].SetLineWidth(3)
        graphs[name].Draw("AL" if name =='pion' else "Lsame")
        if name=='pion':
            graphs[name].GetYaxis().SetTitle('Mass (GeV/c^{2})')
            graphs[name].GetYaxis().SetTitleOffset(1.2)
            graphs[name].GetYaxis().SetRangeUser(-0.1, 1.65)
            graphs[name].GetXaxis().SetTitle(x_title)
            graphs[name].GetXaxis().SetRangeUser(x_min, x_max)
            

        lines[name] = ROOT.TLine(x_min, particles[name]['mass'], x_max, particles[name]['mass'])
        lines[name].SetLineColor(particles[name]['color'])
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


# get_mass_plot()
mass_vs_dtdl('dt')