import ROOT
import numpy as np
from random import choice
from string import ascii_letters

ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)


SPEED_OF_LIGHT = 299.792458 # mm / ns
PION_MASS = 0.13957039 # GeV/ c^{2}
KAON_MASS = 0.493677 # GeV/ c^{2}
PROTON_MASS = 0.93827208816 # GeV/ c^{2}

latex = ROOT.TLatex()
latex.SetTextFont(52)
latex.SetTextSize(0.04)

class Particle:
    def __init__(self, name, legend, color, mass, momentum, track_length, pdg):
        self.name = name
        self.legend = legend
        self.color = color
        self.mass = mass
        self.mass2 = mass*mass
        self.momentum = momentum
        self.track_length = track_length
        self.pdg = pdg
        self.tof = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + (mass/momentum)**2 )
        self.beta = track_length/(SPEED_OF_LIGHT * self.tof)
        self.gamma = 1/np.sqrt(1 - self.beta**2)
        self.legend_graph = ROOT.TGraph()
        self.legend_graph.SetFillColor(self.color)

def create_list_of_particles(momentum, track_length):
    pion = Particle('pion', ' #pi^{#pm}', ROOT.TColor.GetColor('#1b9e77'), PION_MASS, momentum, track_length, 211)
    kaon = Particle('kaon', ' K^{#pm}', ROOT.TColor.GetColor('#d95f02'), KAON_MASS, momentum, track_length, 321)
    proton = Particle('proton', ' p', ROOT.TColor.GetColor('#7570b3'), PROTON_MASS, momentum, track_length, 2212)
    return [pion, kaon, proton]

particles = create_list_of_particles(1., 1) # just for colours! overwrite when studying mass_uncertainty!
(pion, kaon, proton) = particles

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
    '''Draw a 2D histogram in the appropriate styling and pause.'''
    canvas = create_canvas(0.33, 0.58, 0.65)

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
    palette.SetMaxDigits(3)
    palette.SetLabelOffset(0.006)

    canvas.Modified()
    canvas.Update()
    return canvas

def create_legend(x1=0.2, y1=0.75, x2=0.76, y2=0.91):
    legend = ROOT.TLegend(x1, y1, x2, y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(62)
    return legend

def create_canvas(margin=0.22, left_margin_fraction=0.8, bottom_margin_fraction=0.7):
    '''
    Create a square canvas 600x600 with equal horizontal and vertical margins to ensure also a square plot inside the axes!
    '''
    ROOT.gStyle.SetPadLeftMargin( left_margin_fraction * margin )
    ROOT.gStyle.SetPadRightMargin( (1 - left_margin_fraction) * margin )
    ROOT.gStyle.SetPadBottomMargin(bottom_margin_fraction * margin)
    ROOT.gStyle.SetPadTopMargin( (1 - bottom_margin_fraction) * margin)
    canvas = ROOT.TCanvas(get_rand_string(), "", 600, 600)
    return canvas

def draw_vertical_mass_lines(maxy):
    lines = {}
    for p in particles:
        lines[p] = ROOT.TLine(p.mass2, 0., p.mass2, maxy)
        lines[p].SetLineColor(15)
        lines[p].SetLineWidth(2)
        lines[p].SetLineStyle(1)
        lines[p].Draw()
    return lines





##################### RMS90 ###################

def interval_quantile_(x, quant=0.9):
    '''Calculate the shortest interval that contains the desired quantile'''
    # the number of possible starting points
    n_low = int(len(x) * (1. - quant))
    # the number of events contained in the quantil
    n_quant = len(x) - n_low

    # Calculate all distances in one go
    distances = x[-n_low:] - x[:n_low]
    i_start = np.argmin(distances)

    return i_start, i_start + n_quant

def fit90(x): 
    x = np.sort(x)
    n10percent = int(round(len(x)*0.1))
    n90percent = len(x) - n10percent
    
    start, end = interval_quantile_(x, quant=0.9)
    
    rms90 = np.std(x[start:end])
    mean90 = np.mean(x[start:end])
    mean90_err = rms90/np.sqrt(n90percent)
    rms90_err = rms90/np.sqrt(2*n90percent)   # estimator in root
    return mean90, rms90, mean90_err, rms90_err
