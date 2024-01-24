import ROOT
import numpy as np
from random import choice
from string import ascii_letters

SPEED_OF_LIGHT = 299.792458 # mm / ns
PION_MASS = 0.13957039 # GeV/ c^{2}
KAON_MASS = 0.493677 # GeV/ c^{2}
PROTON_MASS = 0.93827208816 # GeV/ c^{2}

latex = ROOT.TLatex()
latex.SetTextFont(52)
latex.SetTextSize(0.04)

class Particle:
    def __init__(self, name, legend, color, mass, momentum, track_length):
        self.name = name
        self.legend = legend
        self.color = color
        self.mass = mass
        self.mass2 = mass*mass
        self.momentum = momentum
        self.track_length = track_length
        self.tof = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + (mass/momentum)**2 )
        self.beta = track_length/(SPEED_OF_LIGHT * self.tof)
        self.gamma = 1/np.sqrt(1 - self.beta**2)

def create_list_of_particles(momentum, track_length):
    pion = Particle('pion', ' #pi^{#pm}', ROOT.TColor.GetColor('#1b9e77'), PION_MASS, momentum, track_length)
    kaon = Particle('kaon', ' K^{#pm}', ROOT.TColor.GetColor('#d95f02'), KAON_MASS, momentum, track_length)
    proton = Particle('proton', ' p', ROOT.TColor.GetColor('#7570b3'), PROTON_MASS, momentum, track_length)
    return [pion, kaon, proton]

def get_rand_string():
    return ''.join(choice(ascii_letters) for i in range(16))


def get_filled_graph(x, y_min, y_max):
    '''
    Return a TGraph filled between y_min and y_max.
    '''
    n_points = len(x)
    graph = ROOT.TGraph(2*n_points)
    for i in range(n_points):
        graph.SetPoint(i, x[i], y_max[i])
        graph.SetPoint(n_points + i, x[n_points-i-1], y_min[n_points-i-1])
    return graph


def get_linear_uncertainty(particle, dp, dl, dt, func="m"):
    '''
    dp, dl - relative. dt - absolute
    '''
    common_uncrt = np.sqrt(dp**2 + (dl**2 + (dt/particle.tof)**2)*particle.gamma**4)
    f0, df = 0, 0
    if func == "m":
        f0 = particle.mass
        df = f0 * common_uncrt
    elif func == "m2":
        f0 = particle.mass*particle.mass
        df = 2 * f0 * common_uncrt
    return np.nan_to_num(f0 - df, nan=0), np.nan_to_num(f0 + df, nan=0)


def get_uncertainty(particle, dp, dl, dt, func="m"):
    '''
    dp, dl - relative. dt - absolute
    Return mass of the particle in momentum bins for the plot
    '''
    mom_down, mom_up = particle.momentum*(1. - dp), particle.momentum*(1. + dp)
    track_length_down, track_length_up = particle.track_length*(1. - dl), particle.track_length*(1. + dl)
    tof_down, tof_up = particle.tof - dt, particle.tof + dt

    beta_down, beta_up = track_length_down/(tof_up*SPEED_OF_LIGHT), track_length_up/(tof_down*SPEED_OF_LIGHT) 
    m2_down, m2_up = mom_down*mom_down*(1./(beta_up*beta_up) - 1.), mom_up*mom_up*(1./(beta_down*beta_down) - 1.)
    if func == "m":
        return np.nan_to_num(np.sqrt(m2_down), nan=0), np.nan_to_num(np.sqrt(m2_up), nan=0)
    return m2_down, m2_up



######################## DRAWING #################

def draw_2d_plot(h, maximum=1e4):
    '''
    Draw a 2D histogram in the appropriate styling and pause.
    '''
    margin = 0.33
    ROOT.gStyle.SetPadLeftMargin(0.6*margin)
    ROOT.gStyle.SetPadRightMargin(0.4*margin)
    ROOT.gStyle.SetPadTopMargin(0.35*margin)
    ROOT.gStyle.SetPadBottomMargin(0.65*margin)
    canvas = ROOT.TCanvas(get_rand_string(),"",600,600)

    h.Draw("colz")
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleOffset(1.4)

    print(f"The highest bin is: {h.GetMaximum()}", end=" ")
    h.SetMinimum(1)
    h.SetMaximum(maximum)
    print(f", the maximum is set to: {h.GetMaximum()}")

    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.Update()

    palette = h.GetListOfFunctions().FindObject("palette")
    x_min = h.GetXaxis().GetXmin()
    x_max = h.GetXaxis().GetXmax()
    palette.SetX1(x_min + 1.01*(x_max-x_min))
    palette.SetX2(x_min + 1.05*(x_max-x_min))
    palette.SetLabelOffset(0.001)

    canvas.Modified()
    canvas.Update()
    return canvas

