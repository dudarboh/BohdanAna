import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

plt.rcParams.update({'font.size':18, 'text.usetex':True})

particles = {}
particles['pion'] = {'mass' : 0.13957039, 'color': '#1b9e77', 'legend' : "#pi^{#pm}", 'legend' : "#pi^{#pm}"}
particles['kaon'] = {'mass' : 0.493677, 'color': '#d95f02', 'legend' : "K^{#pm}"}
particles['proton'] = {'mass' : 0.93827208816, 'color': '#7570b3', 'legend' : "p"}
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

    return np.nan_to_num(mass, nan=-0.1)

def get_mass_plot():
    # We always work in momentum bins
    momentum = np.arange(0.1, 8, 0.01)

    band_names = ['true', 'dp', 'dp+dl', 'total']
    band_sides = ['plus', 'minus']

    # # STAR DATA from the paper https://arxiv.org/abs/nucl-ex/0308022
    # track_length = 2200. # mm at eta = 1.
    # relative_momentum_resolution = 0.013
    # relative_length_resolution = 0.002
    # tof_resolution = 0.100 # ns

    # text= r''' \textbf{STAR apparatus}
    # $ L = 2200 \ \mathrm{mm}$
    # $ \Delta p/p = 0.013 $
    # $ \Delta L / L = 0.002 $
    # $ \Delta T = 100 \ \mathrm{ps} $'''

    # ILD roughly or STAR with modern timing
    track_length = 2200. # mm
    relative_momentum_resolution = 0.013
    relative_length_resolution = 0.002
    tof_resolution = 0.020 # ns

    text= r''' \textbf{STAR apparatus}
    \textbf{modern timing}
    $ L = 2200 \ \mathrm{mm}$
    $ \Delta p/p = 0.013 $
    $ \Delta L / L = 0.002 $
    $ \Delta T = 20 \ \mathrm{ps} $'''


    #pions
    mass_bands = {}

    for name in particles:
        mass_bands[name, 'true'] = get_mass(particles[name]['mass'], momentum, track_length, 0., 0., 0.)
        mass_bands[name, 'up_dpdl'] = get_mass(particles[name]['mass'], momentum, track_length, relative_momentum_resolution, -relative_length_resolution, 0.)
        mass_bands[name, 'down_dpdl'] = get_mass(particles[name]['mass'], momentum, track_length, -relative_momentum_resolution, relative_length_resolution, 0.)
        mass_bands[name, 'up_total'] = get_mass(particles[name]['mass'], momentum, track_length, relative_momentum_resolution, -relative_length_resolution, tof_resolution)
        mass_bands[name, 'down_total'] = get_mass(particles[name]['mass'], momentum, track_length, -relative_momentum_resolution, relative_length_resolution, -tof_resolution)


    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(momentum, mass_bands['pion', 'true'], color=particles['pion']['color'], label="$m_{\pi, \mathrm{true}}$")
    ax.fill_between(momentum, mass_bands['pion', 'down_dpdl'], mass_bands['pion', 'up_dpdl'], color=particles['pion']['color'], hatch='////', label="$m_{\pi}(p \pm \Delta p, \ell \pm \Delta \ell, \mathrm{TOF})$", alpha=0.2)
    ax.fill_between(momentum, mass_bands['pion', 'down_total'], mass_bands['pion', 'up_total'], color=particles['pion']['color'], hatch='....', label="$m_{\pi}(p \pm \Delta p, \ell \pm \Delta \ell, \mathrm{TOF} \pm \Delta \mathrm{TOF})$", alpha=0.2)

    ax.plot(momentum, mass_bands['kaon', 'true'], color=particles['kaon']['color'], label="$m_{\kaon, \mathrm{true}}$")
    ax.fill_between(momentum, mass_bands['kaon', 'down_dpdl'], mass_bands['kaon', 'up_dpdl'], color=particles['kaon']['color'], hatch='////', label="$m_{\kaon}(p \pm \Delta p, \ell \pm \Delta \ell, \mathrm{TOF})$", alpha=0.2)
    ax.fill_between(momentum, mass_bands['kaon', 'down_total'], mass_bands['kaon', 'up_total'], color=particles['kaon']['color'], hatch='....', label="$m_{\kaon}(p \pm \Delta p, \ell \pm \Delta \ell, \mathrm{TOF} \pm \Delta \mathrm{TOF})$", alpha=0.2)

    ax.plot(momentum, mass_bands['proton', 'true'], color=particles['proton']['color'], label="$m_{\proton, \mathrm{true}}$")
    ax.fill_between(momentum, mass_bands['proton', 'down_dpdl'], mass_bands['proton', 'up_dpdl'], color=particles['proton']['color'], hatch='////', label="$m_{\proton}(p \pm \Delta p, \ell \pm \Delta \ell, \mathrm{TOF})$", alpha=0.2)
    ax.fill_between(momentum, mass_bands['proton', 'down_total'], mass_bands['proton', 'up_total'], color=particles['proton']['color'], hatch='....', label="$m_{\proton}(p \pm \Delta p, \ell \pm \Delta \ell, \mathrm{TOF} \pm \Delta \mathrm{TOF})$", alpha=0.2)

    ax.set_xlabel('Momentum ($\mathrm{GeV}/c$)')
    ax.set_ylabel('Mass ($\mathrm{GeV}/c^{2}$)')
    ax.set_title('')
    ax.set_xlim(0.0, 5.1)
    ax.xaxis.set_ticks(np.arange(0, 5.1, 1))
    ax.set_ylim(0.0, 1.7)
    ax.yaxis.set_ticks(np.arange(0.0, 1.7, 0.2))
    ax.grid()

    pion_legend = mpatches.Patch(color=particles['pion']['color'], label='$\pi^{\pm}$')
    kaon_legend = mpatches.Patch(color=particles['kaon']['color'], label='$K^{\pm}$')
    proton_legend = mpatches.Patch(color=particles['proton']['color'], label='$p$')

    true_legend = mlines.Line2D([], [], color='black', label='true mass')
    dpdl_legend = mpatches.Patch(hatch='////', label='$\Delta L$', fill=False)
    total_legend = mpatches.Patch(hatch='....', label='$\Delta L \oplus \Delta T$', fill=False)

    ax.legend(handles=[pion_legend, kaon_legend, proton_legend, true_legend, dpdl_legend, total_legend], ncol=2)
    ax.text(0.2, 1.1, text, linespacing=1.5, bbox=dict(facecolor='white', edgecolor='black', pad=.5, boxstyle='round'))
    plt.show()

def mass_vs_dtdl(x_axis='dt'):
    momentum = 4 # GeV
    track_length = 2000. # mm

    if x_axis == 'dt':
        x_min, x_max, x_title = -100, 100, '#Delta T (ps)'
        legend = ROOT.TLegend(0.2, 0.67, 0.76, 0.93)
        # tof_res is expected to be in ns as an input to get_mass()
        mom_res, len_res, tof_res = 0., 0., np.arange(x_min/1000., x_max/1000., 0.001) # ns
    elif x_axis == 'dl':
        x_min, x_max, x_title = -30, 30, '#Delta L (mm)'
        legend = ROOT.TLegend(0.72, 0.67, 0.995, 0.93)
        # len_res is expected to be relative as an input to get_mass()
        mom_res, len_res, tof_res = 0., np.arange(x_min, x_max, 0.1)/track_length, 0.

    legend.SetFillStyle(0)
    legend.SetBorderSize(0)

    canvas = ROOT.TCanvas("name", "title")
    graphs = {}
    lines = {}

    for name in particles:
        mass = get_mass(particles[name]['mass'], momentum, track_length, mom_res, len_res, tof_res)
        mass[mass < 0] = 0
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
            graphs[name].GetYaxis().SetRangeUser(-0.1, 1.5)
            graphs[name].GetXaxis().SetTitle(x_title)
            graphs[name].GetXaxis().SetRangeUser(x_min, x_max)
            

        lines[name] = ROOT.TLine(x_min, particles[name]['mass'], x_max, particles[name]['mass'])
        lines[name].SetLineColor(ROOT.TColor.GetColor(particles[name]['color']))
        lines[name].SetLineWidth(2)
        lines[name].SetLineStyle(9)
        lines[name].Draw()

        legend.AddEntry(graphs[name], particles[name]['legend'], "l")
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    momentum_text = f"p = {momentum}" + " #frac{GeV}{c}"
    track_length_text = f"L = {track_length:.0f}" + " mm"
    if x_axis == 'dt':
        latex.DrawLatex(-0.3, 1.35, momentum_text)
        latex.DrawLatex(-0.3, 1.2, track_length_text)
    elif x_axis == 'dl':
        latex.DrawLatex(-12, 1.35, momentum_text)
        latex.DrawLatex(-12, 1.2, track_length_text)

    canvas.Update()
    input("wait")


# get_mass_plot()
mass_vs_dtdl('dt')