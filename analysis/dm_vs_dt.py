import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import my_utils
ROOT.gStyle.SetCanvasPreferGL(True)
plt.rcParams.update({'font.size':18, 'text.usetex':True})

SPEED_OF_LIGHT = 299.792458 # mm / ns

def vary_dt(particle, dt):
    '''
    dt is an np.array indicating absolute uncertainties in ns
    Return an array mass of the particle vs dt uncertainty
    '''
    reco_tof = particle.tof + dt
    reco_beta = particle.track_length/(reco_tof*SPEED_OF_LIGHT)
    reco_mass2 = particle.momentum*particle.momentum*(1/(reco_beta*reco_beta) - 1) # GeV^2/c^4
    return np.nan_to_num( np.sqrt(reco_mass2) , nan=0)

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

legend = ROOT.TLegend(0.2, 0.6, 0.59, 0.76)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextFont(62)
# legend.SetMargin(0.15)
# legend.SetTextAlign(32)




graphs = {}
true_masses = {}
slopes ={}
dt = np.arange(-100.,101.)

y_min_limit, y_max_limit = -0.1, 1.65
for p in particles:
    m = vary_dt(p, dt/1000 )
    graphs[p.name] = ROOT.TGraph(len(dt), dt, m)
    graphs[p.name].SetLineColor(p.color)
    graphs[p.name].SetLineWidth(3)
    if p.name=='pion':
        graphs[p.name].Draw("AL")
        graphs[p.name].GetYaxis().SetTitle('Mass (GeV/c^{2})')
        graphs[p.name].GetYaxis().SetTitleOffset(1.2)
        graphs[p.name].GetYaxis().SetRangeUser(y_min_limit, y_max_limit)
        graphs[p.name].GetXaxis().SetTitle("#Delta T (ps)")
        graphs[p.name].GetXaxis().SetRangeUser(dt[0], dt[-1])
    else:
        graphs[p.name].Draw("Lsame")

    true_masses[p.name] = ROOT.TLine(dt[0], p.mass, dt[-1], p.mass)
    true_masses[p.name].SetLineColor(p.color)
    true_masses[p.name].SetLineWidth(2)
    true_masses[p.name].SetLineStyle(9)
    # true_masses[p.name].Draw()

    mid_idx = int(len(dt)/2)
    slope = (m[mid_idx+1] - m[mid_idx-1]) / (dt[mid_idx+1] - dt[mid_idx-1])
    y_line = p.mass + slope * dt
    slopes[p.name] = ROOT.TLine(dt[y_line > y_min_limit][0], y_line[y_line > y_min_limit][0], dt[y_line < y_max_limit][-1], y_line[y_line < y_max_limit][-1])
    slopes[p.name].SetLineColor(p.color)
    slopes[p.name].SetLineWidth(2)
    slopes[p.name].SetLineStyle(9)
    slopes[p.name].Draw()


    legend.AddEntry(graphs[p.name], p.legend, "l")
legend.Draw()

latex = ROOT.TLatex()
latex.SetTextFont(42)
momentum_text = f"p = {momentum}" + " #frac{GeV}{c}"
track_length_text = f"L = {track_length:.0f}" + " mm"
latex.DrawLatex(-70, 1.45, momentum_text)
latex.DrawLatex(11, 1.45, track_length_text)


canvas.Update()
input("wait")