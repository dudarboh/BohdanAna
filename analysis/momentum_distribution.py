import ROOT
import numpy as np
from utils import *
ROOT.EnableImplicitMT()

dark24 = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']

colors = [ROOT.TColor.GetColor(c) for c in dark24]

# 2f_Z_hadronic_eLpR_dst.root
# 2f_Z_hadronic_eRpL_dst.root
# 4f_WW_semileptonic_eLpR_dst.root
# 4f_WW_semileptonic_eRpL_dst.root
# BohdanAna.root
# higgs_n23n23h_eLpR_dst.root
# higgs_n23n23h_eRpL_dst.root


def general():
    histos = {}

    n_events_z = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root").Count()
    n_events_ww = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root").Count()
    n_events_higgs = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h.root").Count()

    df_z = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_ww = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_higgs = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    histos["Zqq"] = df_z.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWqq"] = df_ww.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["H"] = df_higgs.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    of_z = df_z.Filter("mom > 10").Count()
    of_ww = df_ww.Filter("mom > 10").Count()
    of_h = df_higgs.Filter("mom > 10").Count()

    z_overflow = 100.*of_z.GetValue()/histos["Zqq"].GetEntries()
    ww_overflow = 100.*of_ww.GetValue()/histos["WWqq"].GetEntries()
    h_overflow = 100.*of_h.GetValue()/histos["H"].GetEntries()

    histos["Zqq"].Scale(1./n_events_z.GetValue())
    histos["WWqq"].Scale(1./n_events_ww.GetValue())
    histos["H"].Scale(1./n_events_higgs.GetValue())

    histos["Zqq"].SetTitle("Z #rightarrow q#bar{q}" + f" (overflow: {z_overflow:.1f} %)" + ";Momentum (GeV/c);N kaons per event")
    histos["WWqq"].SetTitle("WW #rightarrow l#nu_{l} q#bar{q}" + f" (overflow: {ww_overflow:.1f} %)")
    histos["H"].SetTitle("ZH #rightarrow #nu_{#mu,#tau} #bar{#nu}_{#mu,#tau} H" + f" (overflow: {h_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")

    histos["Zqq"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["Zqq"].GetYaxis().SetMaxDigits(3)
    histos["Zqq"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetFillStyle(0)
    leg.SetMargin(0.15)
    c.Update()
    input("wait")


# general()

def zqq():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_bb = df_events.Filter("quarksToPythia == 55").Count()
    n_events_cc = df_events.Filter("quarksToPythia == 44").Count()
    n_events_ss = df_events.Filter("quarksToPythia == 33").Count()
    n_events_uu = df_events.Filter("quarksToPythia == 22").Count()
    n_events_dd = df_events.Filter("quarksToPythia == 11").Count()

    histos["Zbb"] = df.Filter("quarksToPythia == 55").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zcc"] = df.Filter("quarksToPythia == 44").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zss"] = df.Filter("quarksToPythia == 33").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zuu"] = df.Filter("quarksToPythia == 22").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zdd"] = df.Filter("quarksToPythia == 11").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    of_z_bb = df.Filter("quarksToPythia == 55").Filter("mom > 10").Count()
    of_z_cc = df.Filter("quarksToPythia == 44").Filter("mom > 10").Count()
    of_z_ss = df.Filter("quarksToPythia == 33").Filter("mom > 10").Count()
    of_z_uu = df.Filter("quarksToPythia == 22").Filter("mom > 10").Count()
    of_z_dd = df.Filter("quarksToPythia == 11").Filter("mom > 10").Count()

    bb_overflow = 100.*of_z_bb.GetValue()/histos["Zbb"].GetEntries()
    cc_overflow = 100.*of_z_cc.GetValue()/histos["Zcc"].GetEntries()
    ss_overflow = 100.*of_z_ss.GetValue()/histos["Zss"].GetEntries()
    uu_overflow = 100.*of_z_uu.GetValue()/histos["Zuu"].GetEntries()
    dd_overflow = 100.*of_z_dd.GetValue()/histos["Zdd"].GetEntries()

    histos["Zbb"].Scale(1./n_events_bb.GetValue())
    histos["Zcc"].Scale(1./n_events_cc.GetValue())
    histos["Zss"].Scale(1./n_events_ss.GetValue())
    histos["Zuu"].Scale(1./n_events_uu.GetValue())
    histos["Zdd"].Scale(1./n_events_dd.GetValue())

    histos["Zbb"].SetTitle("Z #rightarrow b#bar{b}" + f" (overflow: {bb_overflow:.1f} %)" + ";Momentum (GeV/c);N kaons per event")
    histos["Zcc"].SetTitle("Z #rightarrow c#bar{c}" + f" (overflow: {cc_overflow:.1f} %)")
    histos["Zss"].SetTitle("Z #rightarrow s#bar{s}" + f" (overflow: {ss_overflow:.1f} %)")
    histos["Zuu"].SetTitle("Z #rightarrow u#bar{u}" + f" (overflow: {uu_overflow:.1f} %)")
    histos["Zdd"].SetTitle("Z #rightarrow d#bar{d}" + f" (overflow: {dd_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["Zbb"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["Zbb"].GetYaxis().SetMaxDigits(3)
    histos["Zbb"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")

# zqq()
    
def ww():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 2212")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_ud = df_events.Filter("quarksToPythia == 12 || quarksToPythia == 21").Count()
    n_events_us = df_events.Filter("quarksToPythia == 32 || quarksToPythia == 23").Count()
    n_events_cd = df_events.Filter("quarksToPythia == 14 || quarksToPythia == 41").Count()
    n_events_cs = df_events.Filter("quarksToPythia == 34 || quarksToPythia == 43").Count()

    histos["WWud"] = df.Filter("quarksToPythia == 12 || quarksToPythia == 21").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWus"] = df.Filter("quarksToPythia == 32 || quarksToPythia == 23").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWcd"] = df.Filter("quarksToPythia == 14 || quarksToPythia == 41").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWcs"] = df.Filter("quarksToPythia == 34 || quarksToPythia == 43").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    of_ww_ud = df.Filter("quarksToPythia == 12 || quarksToPythia == 21").Filter("mom > 10").Count()
    of_ww_us = df.Filter("quarksToPythia == 32 || quarksToPythia == 23").Filter("mom > 10").Count()
    of_ww_cd = df.Filter("quarksToPythia == 14 || quarksToPythia == 41").Filter("mom > 10").Count()
    of_ww_cs = df.Filter("quarksToPythia == 34 || quarksToPythia == 43").Filter("mom > 10").Count()


    ud_overflow = 100.*of_ww_ud.GetValue()/histos["WWud"].GetEntries()
    us_overflow = 100.*of_ww_us.GetValue()/histos["WWus"].GetEntries()
    cd_overflow = 100.*of_ww_cd.GetValue()/histos["WWcd"].GetEntries()
    cs_overflow = 100.*of_ww_cs.GetValue()/histos["WWcs"].GetEntries()

    histos["WWud"].Scale(1./n_events_ud.GetValue())
    histos["WWus"].Scale(1./n_events_us.GetValue())
    histos["WWcd"].Scale(1./n_events_cd.GetValue())
    histos["WWcs"].Scale(1./n_events_cs.GetValue())

    histos["WWud"].SetTitle("WW #rightarrow l#nu_{l} ud" + f" (overflow: {ud_overflow:.1f} %)" + ";Momentum (GeV/c);N protons per event")
    histos["WWus"].SetTitle("WW #rightarrow l#nu_{l} us" + f" (overflow: {us_overflow:.1f} %)")
    histos["WWcd"].SetTitle("WW #rightarrow l#nu_{l} cd" + f" (overflow: {cd_overflow:.1f} %)")
    histos["WWcs"].SetTitle("WW #rightarrow l#nu_{l} cs" + f" (overflow: {cs_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["WWud"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["WWud"].GetYaxis().SetMaxDigits(3)
    histos["WWud"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")
# ww()ll 


def higgs():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_bb = df_events.Filter("higgsDaughters == 505").Count()
    n_events_ww = df_events.Filter("higgsDaughters == 2424").Count()
    n_events_gg = df_events.Filter("higgsDaughters == 2121 || higgsDaughters == 909").Count()
    n_events_cc = df_events.Filter("higgsDaughters == 404").Count()
    n_events_zz = df_events.Filter("higgsDaughters == 2323").Count()
    # n_events_ss = df_events.Filter("higgsDaughters == 303").Count()

    histos["Hbb"] = df.Filter("higgsDaughters == 505").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hww"] = df.Filter("higgsDaughters == 2424").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hgg"] = df.Filter("higgsDaughters == 2121 || higgsDaughters == 909").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hcc"] = df.Filter("higgsDaughters == 404").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hzz"] = df.Filter("higgsDaughters == 2323").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    # histos["Hss"] = df.Filter("higgsDaughters == 303").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    #overflows
    of_hbb = df.Filter("higgsDaughters == 505").Filter("mom > 10").Count()
    of_hww = df.Filter("higgsDaughters == 2424").Filter("mom > 10").Count()
    of_hgg = df.Filter("higgsDaughters == 2121 || higgsDaughters == 909").Filter("mom > 10").Count()
    of_hcc = df.Filter("higgsDaughters == 404").Filter("mom > 10").Count()
    of_hzz = df.Filter("higgsDaughters == 2323").Filter("mom > 10").Count()
    # of_hss = df.Filter("higgsDaughters == 303").Filter("mom > 10").Count()

    hbb_overflow = 100.*of_hbb.GetValue()/histos["Hbb"].GetEntries()
    hww_overflow = 100.*of_hww.GetValue()/histos["Hww"].GetEntries()
    hgg_overflow = 100.*of_hgg.GetValue()/histos["Hgg"].GetEntries()
    hcc_overflow = 100.*of_hcc.GetValue()/histos["Hcc"].GetEntries()
    hzz_overflow = 100.*of_hzz.GetValue()/histos["Hzz"].GetEntries()
    # hss_overflow = 100.*of_hss.GetValue()/histos["Hss"].GetEntries()

    histos["Hbb"].Scale(1./n_events_bb.GetValue())
    histos["Hww"].Scale(1./n_events_ww.GetValue())
    histos["Hgg"].Scale(1./n_events_gg.GetValue())
    histos["Hcc"].Scale(1./n_events_cc.GetValue())
    histos["Hzz"].Scale(1./n_events_zz.GetValue())
    # histos["Hss"].Scale(1./n_events_ss.GetValue())

    histos["Hbb"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} b#bar{b}" + f" (overflow: {hbb_overflow:.1f} %)" + ";Momentum (GeV/c);N kaons per event")
    histos["Hww"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} WW" + f" (overflow: {hww_overflow:.1f} %)")
    histos["Hgg"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} gg" + f" (overflow: {hgg_overflow:.1f} %)")
    histos["Hcc"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} c#bar{c}" + f" (overflow: {hcc_overflow:.1f} %)")
    histos["Hzz"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} ZZ" + f" (overflow: {hzz_overflow:.1f} %)")
    # histos["Hss"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} s#bar{s}" + f" (overflow: {hss_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["Hbb"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["Hbb"].GetYaxis().SetMaxDigits(3)
    histos["Hbb"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")
higgs()