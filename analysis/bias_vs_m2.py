import ROOT
import numpy as np
from utils import *

m2_rec = np.array([1.95304e-02, 2.45192e-01, 8.85596e-01])
m2_err_rec = np.array([1.27751e-06, 1.48563e-05, 3.45700e-05])
m2_true = np.array([p.mass2 for p in particles])

m2_rec = np.sqrt(m2_rec)
m2_true = np.sqrt(m2_true)

diff = m2_rec - m2_true

canvas = create_canvas()
frame1 = canvas.DrawFrame(0., 0., 1., 1.1*max(diff), "Mass^{2} bias vs particles mass^{2}; Mass^{2}_{true} (GeV^{2}/c^{4}); Mass^{2}_{reco} - Mass^{2}_{true} (GeV^{2}/c^{4})")

gr = ROOT.TGraph(3, m2_true, diff)
gr.Draw("PL")
canvas.Update()
input("wait")
