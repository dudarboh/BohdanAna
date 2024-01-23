import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import my_utils

ROOT.gStyle.SetCanvasPreferGL(True)
plt.rcParams.update({'font.size':18, 'text.usetex':True})

def get_appox_uncertainty(particle, dp, dl, dt):
    '''
    dp, dl - relative. dt - absolute
    '''
    dm = particle.mass * np.sqrt(dp**2 + (dl**2 + (dt/particle.tof)**2)*particle.gamma**4)
    return np.nan_to_num(particle.mass - dm, nan=0), np.nan_to_num(particle.mass + dm, nan=0)


def get_uncertainty(particle, dp, dl, dt):
    '''
    dp, dl - relative. dt - absolute
    Return mass of the particle in momentum bins for the plot
    '''
    mom_down, mom_up = particle.momentum*(1. - dp), particle.momentum*(1. + dp)
    track_length_down, track_length_up = particle.track_length*(1. - dl), particle.track_length*(1. + dl)
    tof_down, tof_up = particle.tof - dt, particle.tof + dt

    beta_down, beta_up = track_length_down/(tof_up*SPEED_OF_LIGHT), track_length_up/(tof_down*SPEED_OF_LIGHT) 
    m2_down, m2_up = mom_down*mom_down*(1./(beta_up*beta_up) - 1.), mom_up*mom_up*(1./(beta_down*beta_down) - 1.)

    m_down, m_up = np.nan_to_num(np.sqrt(m2_down), nan=0), np.nan_to_num(np.sqrt(m2_up), nan=0)
    return m_down, m_up


momentum, track_length = np.arange(0.1, 8, 0.01), 2200.
# # STAR DATA from the paper https://arxiv.org/abs/nucl-ex/0308022
dp = 0.013 # relative
dl = 0.002 # relative
dt = 0.100 # ns
# dt = 0.020 # ns

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
graphs = {}

fill_alpha = .3
legend = ROOT.TLegend(0.2, 0.72, 0.93, 0.92)
legend.SetNColumns(2)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetMargin(0.3)
legend.SetColumnSeparation(0.2)
# legend.SetTextAlign(32)
legend.SetTextFont(62)


# calculating
for p in particles:
    m_down_approx, m_up_approx = get_appox_uncertainty(p, dp, dl, dt)
    m_down, m_up = get_uncertainty(p, dp, dl, dt)
    graphs[p, 'approx'] = get_filled_graph(momentum, m_down_approx, m_up_approx)
    graphs[p, 'total'] = get_filled_graph(momentum, m_down, m_up)
    graphs[p, 'mass'] = ROOT.TGraph(len(momentum), momentum, np.full_like(momentum, p.mass))

    if p.name == 'pion':
        graphs[p, 'mass'].Draw("AL")
        graphs[p, 'mass'].GetXaxis().SetTitle("Momentum (GeV/c)")
        graphs[p, 'mass'].GetYaxis().SetTitle("Mass (GeV/c^{2})")
        graphs[p, 'mass'].GetYaxis().SetTitleOffset(1.2)
        graphs[p, 'mass'].GetXaxis().SetRangeUser(0, 8)
        graphs[p, 'mass'].GetYaxis().SetRangeUser(0, 2.3)
    else:
        graphs[p, 'mass'].Draw("L same")

    graphs[p, 'mass'].SetLineColor( p.color )
    graphs[p, 'approx'].SetLineWidth(3)
    graphs[p, 'approx'].SetLineStyle(9)
    graphs[p, 'approx'].SetLineColor(p.color)
    # graphs[p, 'approx'].SetFillColorAlpha(p.color, 0.25)
    graphs[p, 'total'].SetFillColorAlpha(p.color, fill_alpha)
    graphs[p, 'approx'].Draw("Lsame")
    graphs[p, 'total'].Draw("fsame")


#row1 col1
gr_legend_true = ROOT.TGraph()
gr_legend_true.SetLineColor(ROOT.kBlack)
legend.AddEntry(gr_legend_true, "true mass", "l")

#row1 col2
gr_legend_pi = ROOT.TGraph()
gr_legend_pi.SetFillColor(pion.color)
legend.AddEntry(gr_legend_pi, pion.legend, "f")

#row2 col1
gr_legend_dl = ROOT.TGraph()
gr_legend_dl.SetLineStyle(9)
gr_legend_dl.SetLineWidth(3)
legend.AddEntry(gr_legend_dl, "#frac{#Deltap}{p} \u2295 #gamma^{2}(#frac{#DeltaL}{L} \u2295 #frac{#DeltaT}{T})", "l")

#row2 col2
gr_legend_k = ROOT.TGraph()
gr_legend_k.SetFillColor(kaon.color)
legend.AddEntry(gr_legend_k, kaon.legend, "f")

#row3 col1
gr_legend_dldt = ROOT.TGraph()
gr_legend_dldt.SetFillColor(ROOT.kBlack)
gr_legend_dldt.SetFillColorAlpha(ROOT.kBlack, fill_alpha)
legend.AddEntry(gr_legend_dldt, "m(p#pm#Deltap, L#mp#DeltaL, T#pm#DeltaT)", "f")

#row3 col2
gr_legend_p = ROOT.TGraph()
gr_legend_p.SetFillColor(proton.color)
legend.AddEntry(gr_legend_p, proton.legend, "f")


legend.Draw()

latex = ROOT.TLatex()
# latex = ROOT.TMathText()
latex.SetTextFont(52)
latex.SetTextSize(0.04)
x_text_pos, y0_text_pos, dy = 0.44, 1.53, 0.14
latex.DrawLatex(x_text_pos, y0_text_pos, f"L = {track_length:.0f}" + " mm")
latex.DrawLatex(x_text_pos, y0_text_pos - dy, "#Deltap/p = " + f"{dp:.3f}")
latex.DrawLatex(x_text_pos, y0_text_pos - 2*dy, "#DeltaL/L = " + f"{dl:.3f}")
latex.DrawLatex(x_text_pos, y0_text_pos - 3*dy, "#DeltaT = " + f"{1000*dt:.0f}" + " ps")
canvas.Update()
input("wait")
