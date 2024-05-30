import ROOT
import numpy as np
from utils import *
ROOT.EnableImplicitMT()

dark24 = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']

colors = [ROOT.TColor.GetColor(c) for c in dark24]

df_z = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root").Filter("!isSimulated && !isOverlay")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
df_h = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_dst.root").Filter("!isSimulated && !isOverlay")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")



h_z_bb_pi_tot = df_z.Filter("quarksToPythia == 11")\
                .Filter("abs(pdg) == 321")\
                .Histo1D((get_rand_string(), "Z #rightarrow b#bar{b} total;Momentum (GeV/c);N entries", 200, 0, 10), "mom")

h_z_bb_pi_hadr = df_z.Filter("quarksToPythia == 11")\
                .Filter("abs(pdg) == 321")\
                .Filter("isHadronisationDecay")\
                .Histo1D((get_rand_string(), "From hadronisation;Momentum (GeV/c);N entries", 200, 0, 10), "mom")

h_z_bb_pi_nohadr = df_z.Filter("quarksToPythia == 11")\
                .Filter("abs(pdg) == 321")\
                .Filter("isCharmQuarkDecay")\
                .Histo1D((get_rand_string(), "Not from hadronisation;Momentum (GeV/c);N entries", 200, 0, 10), "mom")

h_z_bb_pi_nohadr2 = df_z.Filter("quarksToPythia == 11")\
                .Filter("abs(pdg) == 321")\
                .Filter("isBottomQuarkDecay")\
                .Histo1D((get_rand_string(), "Not from hadronisation;Momentum (GeV/c);N entries", 200, 0, 10), "mom")


h_z_bb_pi_tot.SetLineColor(colors[0])
h_z_bb_pi_hadr.SetLineColor(colors[1])
h_z_bb_pi_nohadr.SetLineColor(colors[2])
h_z_bb_pi_nohadr2.SetLineColor(colors[3])
c1 = create_canvas()
leg = create_legend()
h_z_bb_pi_tot.Draw("")
leg.AddEntry(h_z_bb_pi_tot.GetPtr(), h_z_bb_pi_tot.GetTitle(), "l")

h_z_bb_pi_hadr.Draw("same")
leg.AddEntry(h_z_bb_pi_hadr.GetPtr(), h_z_bb_pi_hadr.GetTitle(), "l")

h_z_bb_pi_nohadr.Draw("same")
leg.AddEntry(h_z_bb_pi_nohadr.GetPtr(), h_z_bb_pi_nohadr.GetTitle(), "l")
h_z_bb_pi_nohadr2.Draw("same")
leg.AddEntry(h_z_bb_pi_nohadr2.GetPtr(), h_z_bb_pi_nohadr2.GetTitle(), "l")

input("wait")

df_h_pi = df_h.Filter("abs(pdg) == 211")

df_bb = df_z_pi.Filter("quarksToPythia == 55")
df_cc = df_z_pi.Filter("quarksToPythia == 44")
df_ss = df_z_pi.Filter("quarksToPythia == 33")
df_uu = df_z_pi.Filter("quarksToPythia == 22")
df_dd = df_z_pi.Filter("quarksToPythia == 11")
df_h_bb = df_h_pi.Filter("higgsDaughters == 505")
df_h_ww = df_h_pi.Filter("higgsDaughters == 2424")
df_h_gg = df_h_pi.Filter("higgsDaughters == 2121")
df_h_cc = df_h_pi.Filter("higgsDaughters == 404")
df_h_zz = df_h_pi.Filter("higgsDaughters == 2323")


h_bb = df_bb.Histo1D((get_rand_string(), "Z #rightarrow b#bar{b};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_cc = df_cc.Histo1D((get_rand_string(), "Z #rightarrow c#bar{c};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_ss = df_ss.Histo1D((get_rand_string(), "Z #rightarrow s#bar{s};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_uu = df_uu.Histo1D((get_rand_string(), "Z #rightarrow u#bar{u};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_dd = df_dd.Histo1D((get_rand_string(), "Z #rightarrow d#bar{d};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_h_bb = df_h_bb.Histo1D((get_rand_string(), "ZH #rightarrow #nu_{#mu,#tau}#bar{#nu_{#mu,#tau}}b#bar{b};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_h_ww = df_h_ww.Histo1D((get_rand_string(), "ZH #rightarrow #nu_{#mu,#tau}#bar{#nu_{#mu,#tau}}WW;Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_h_gg = df_h_gg.Histo1D((get_rand_string(), "ZH #rightarrow #nu_{#mu,#tau}#bar{#nu_{#mu,#tau}}gg;Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_h_cc = df_h_cc.Histo1D((get_rand_string(), "ZH #rightarrow #nu_{#mu,#tau}#bar{#nu_{#mu,#tau}}c#bar{c};Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")
h_h_zz = df_h_zz.Histo1D((get_rand_string(), "ZH #rightarrow #nu_{#mu,#tau}#bar{#nu_{#mu,#tau}}ZZ;Momentum (GeV/c);Fraction of particles", 500, 0, 10), "mom")

h_bb.SetLineColor(colors[0])
h_cc.SetLineColor(colors[1])
h_ss.SetLineColor(colors[2])
h_uu.SetLineColor(colors[3])
h_dd.SetLineColor(colors[4])
h_h_bb.SetLineColor(colors[5])
h_h_ww.SetLineColor(colors[6])
h_h_gg.SetLineColor(colors[7])
h_h_cc.SetLineColor(colors[8])
h_h_zz.SetLineColor(colors[9])

c1 = create_canvas()
leg = create_legend()
h_bb.Scale(1./h_bb.GetEntries())
h_bb.Draw("histo")
leg.AddEntry(h_bb.GetPtr(), h_bb.GetTitle(), "l")

h_cc.Scale(1./h_cc.GetEntries())
h_cc.Draw("histo same")
leg.AddEntry(h_cc.GetPtr(), h_cc.GetTitle(), "l")

h_ss.Scale(1./h_ss.GetEntries())
h_ss.Draw("histo same")
leg.AddEntry(h_ss.GetPtr(), h_ss.GetTitle(), "l")

h_uu.Scale(1./h_uu.GetEntries())
h_uu.Draw("histo same")
leg.AddEntry(h_uu.GetPtr(), h_uu.GetTitle(), "l")

h_dd.Scale(1./h_dd.GetEntries())
h_dd.Draw("histo same")
leg.AddEntry(h_dd.GetPtr(), h_dd.GetTitle(), "l")

h_h_bb.Scale(1./h_h_bb.GetEntries())
h_h_bb.Draw("histo same")
leg.AddEntry(h_h_bb.GetPtr(), h_h_bb.GetTitle(), "l")

h_h_ww.Scale(1./h_h_ww.GetEntries())
h_h_ww.Draw("histo same")
leg.AddEntry(h_h_ww.GetPtr(), h_h_ww.GetTitle(), "l")

h_h_gg.Scale(1./h_h_gg.GetEntries())
h_h_gg.Draw("histo same")
leg.AddEntry(h_h_gg.GetPtr(), h_h_gg.GetTitle(), "l")

h_h_cc.Scale(1./h_h_cc.GetEntries())
h_h_cc.Draw("histo same")
leg.AddEntry(h_h_cc.GetPtr(), h_h_cc.GetTitle(), "l")

h_h_zz.Scale(1./h_h_zz.GetEntries())
h_h_zz.Draw("histo same")
leg.AddEntry(h_h_zz.GetPtr(), h_h_zz.GetTitle(), "l")

leg.Draw()
c1.Update()
input("wait")



