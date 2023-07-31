import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

plt.rcParams.update({'font.size':18, 'text.usetex':True})



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

    particles = {}
    particles['pion'] = {'mass' : 0.13957039, 'color': '#1b9e77'}
    particles['kaon'] = {'mass' : 0.493677, 'color': '#d95f02'}
    particles['proton'] = {'mass' : 0.93827208816, 'color': '#7570b3'}

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

def mass_vs_dt():
    momentum = 4 # GeV
    track_length = 2000. # mm
    tof_resolution = np.arange(-0.100, 0.100, 0.001) # ns

    particles = {}
    particles['pion'] = {'mass' : 0.13957039, 'color': '#1b9e77'}
    particles['kaon'] = {'mass' : 0.493677, 'color': '#d95f02'}
    particles['proton'] = {'mass' : 0.93827208816, 'color': '#7570b3'}

    fig, ax = plt.subplots(figsize=(8, 8))
    for name in particles:
        mass = get_mass(particles[name]['mass'], momentum, track_length, 0., 0., tof_resolution)
        mass[mass < 0] = 0

        ax.plot(tof_resolution*1000, mass, color=particles[name]['color'])

    ax.set_xlabel('$\Delta T (\mathrm{ps})$')
    ax.set_ylabel('Mass ($\mathrm{GeV}/c^{2}$)')
    ax.set_xlim(-100, 100)
    ax.set_ylim(-0.1, 1.2)
    ax.yaxis.set_ticks(np.arange(-0.1, 1.2, 0.1))
    ax.grid()

    text= r'''$ p = 4 \ \mathrm{GeV / c}$    
    $ L = 2000 \ \mathrm{mm}$'''

    pion_legend = mpatches.Patch(color=particles['pion']['color'], label='$\pi^{\pm}$')
    kaon_legend = mpatches.Patch(color=particles['kaon']['color'], label='$K^{\pm}$')
    proton_legend = mpatches.Patch(color=particles['proton']['color'], label='$p$')

    ax.legend(handles=[proton_legend, kaon_legend, pion_legend])
    ax.text(30, 0.1, text, linespacing=1.5, bbox=dict(facecolor='white', edgecolor='black', pad=.5, boxstyle='round'))
    plt.show()


def mass_vs_dl():
    momentum = 4 # GeV
    track_length = 2000. # mm
    len_resolution = np.arange(-30, 30, 0.1) # mm

    particles = {}
    particles['pion'] = {'mass' : 0.13957039, 'color': '#1b9e77'}
    particles['kaon'] = {'mass' : 0.493677, 'color': '#d95f02'}
    particles['proton'] = {'mass' : 0.93827208816, 'color': '#7570b3'}

    fig, ax = plt.subplots(figsize=(8, 8))
    for name in particles:
        mass = get_mass(particles[name]['mass'], momentum, track_length, 0., len_resolution/track_length, 0.)
        mass[mass < 0] = 0

        ax.plot(len_resolution, mass, color=particles[name]['color'])

    ax.set_xlabel('$\Delta L (\mathrm{mm})$')
    ax.set_ylabel('Mass ($\mathrm{GeV}/c^{2}$)')
    ax.set_xlim(-30, 30)
    ax.set_ylim(-0.1, 1.2)
    ax.yaxis.set_ticks(np.arange(-0.1, 1.2, 0.1))
    ax.grid()

    text= r'''$ p = 4 \ \mathrm{GeV / c}$    
    $ L = 2000 \ \mathrm{mm}$'''

    pion_legend = mpatches.Patch(color=particles['pion']['color'], label='$\pi^{\pm}$')
    kaon_legend = mpatches.Patch(color=particles['kaon']['color'], label='$K^{\pm}$')
    proton_legend = mpatches.Patch(color=particles['proton']['color'], label='$p$')

    ax.legend(handles=[proton_legend, kaon_legend, pion_legend])
    ax.text(-25, 0.1, text, linespacing=1.5, bbox=dict(facecolor='white', edgecolor='black', pad=.5, boxstyle='round'))
    plt.show()


get_mass_plot()
# mass_vs_dl()