import ROOT
import numpy as np
from utils import *

def bias_from_fit():
    m2_rec = np.array([1.95304e-02, 2.45192e-01, 8.85596e-01])
    m2_err_rec = np.array([1.27751e-06, 1.48563e-05, 3.45700e-05])
    m2_true = np.array([p.mass2 for p in particles])

    m2_rec = np.sqrt(m2_rec)
    m2_true = np.sqrt(m2_true)

    diff = (m2_rec - m2_true)*1000

    canvas = create_canvas()
    frame1 = canvas.DrawFrame(0., 0., 1., 1.1*max(diff), "; Mass_{true} (GeV/c^{2}); Mass_{reco} - Mass_{true} (MeV/c^{2})")
    frame1.GetYaxis().SetMaxDigits(3)
    gr = ROOT.TGraph(3, m2_true, diff)
    gr.Draw("PL")
    canvas.Update()
    input("wait")

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
        .Filter("abs(pdg) == 321").Filter("tofSETFront0 > 6.")\
        .Define("beta", "trackLengthToSET_IKF_zedLambda/(tofSETFront0*299.792458)")\
        .Define("mom2", "harmonicMomToSET_IKF_zedLambda*harmonicMomToSET_IKF_zedLambda")\
        .Define("mom", "sqrt(mom2)")\
        .Define("mass2", "mom2*(1./(beta*beta) - 1.)")\
        .Define("mass", "sqrt(mass2)*1000.")

h2d = df.Histo2D( (get_rand_string(), "", 100, 0, 5, 100, KAON_MASS*1000 - 10 , KAON_MASS*1000 + 10 ), "mom", "mass")
h2d.GetXaxis().SetTitle("Momentum (GeV/c)")
h2d.GetYaxis().SetTitle("Mass (GeV/c^{2})")
canvas = draw_2d_plot(h2d, h2d.GetMaximum())
canvas.SetLogz(False)
line = ROOT.TLine(0., 493.677, 5., 493.677)
line.SetLineStyle(9)
line.SetLineWidth(3)
line.SetLineColor(ROOT.kRed+1)
line.Draw()
canvas.Update()
input("wait")
