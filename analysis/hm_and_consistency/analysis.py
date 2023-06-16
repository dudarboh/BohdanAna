import ROOT
import numpy as np
# ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetLegendFillColor(0)


def plot_2d(n_seg1, n_seg2):
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/hm/trk_len_hm*.root")
    df = df.Filter("tof > 0")#.Filter("abs(pdg) == 2212") # can't be 0 or negative
    df = df.Define("beta", "trackLength/(tof*299.792458)").Filter("beta >= 0 && beta <= 1")
    df = df.Define("mass", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000") # in MeV

    # not so default binning anymore, but for old plots: 30, 0, 15, 200, -100, 1300.

    df = df.Filter(f"nSegments > {n_seg1} && nSegments < {n_seg2}")
    h = df.Histo2D(("h", f"{n_seg1} < N segments < {n_seg2}; momentum [GeV]; Mass [MeV]", 500, 0, 15, 500, -100, 1300.), "momentumHM","mass")
    print(h.GetMaximum())

    def draw_lines():
        lines = {}
        pdgs = [211, 321, 2212]
        m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
        for pdg in pdgs:
            lines[pdg] = ROOT.TLine(0., m_pdg[pdg], 15., m_pdg[pdg])
            lines[pdg].SetLineColor(2)
            lines[pdg].SetLineWidth(1)
            lines[pdg].SetLineStyle(9)
            lines[pdg].Draw()
        return lines


    canvas = ROOT.TCanvas(f"{n_seg1}_n_seg_{n_seg2}")
    h.Draw("colz")
    canvas.Update()
    palette = h.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.93)
    palette.SetX2NDC(0.95)
    lines = draw_lines()
    h.SetStats(0)
    h.SetMinimum(0.1)
    h.SetMaximum(3000)
    canvas.SetLogz()
    canvas.SetGridy(0)
    canvas.Update()
    input("wait")
    canvas.Print(f"{canvas.GetName()}.png")

def plot_1d():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/hm/trk_len_hm*.root")
    df = df.Filter("tof > 0") # can't be 0 or negative
    df = df.Define("beta", "trackLength/(tof*299.792458)").Filter("beta >= 0 && beta <= 1")
    df = df.Define("massIP", "momentumIP*sqrt( 1./(beta*beta) - 1.)*1000") # in MeV
    df = df.Define("massCalo", "momentumCalo*sqrt( 1./(beta*beta) - 1.)*1000") # in MeV
    df = df.Define("massHM", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000") # in MeV

    df = df.Filter("massIP != 0 && massCalo != 0 && massHM != 0")
    # Binning for some plots, i made: (1000, -100, 1300.)

    h_ip = df.Histo1D(("h_ip", "At IP; mass [MeV]; N entries", 1000, -100, 1300.), "massIP")
    h_calo = df.Histo1D(("h_calo", "At Calo; mass [MeV]; N entries", 1000, -100, 1300.), "massCalo")
    h_hm = df.Histo1D(("h_hm", "HM; mass [MeV]; N entries", 1000, -100, 1300.), "massHM")


    max_z = 1.05*max(h_ip.GetMaximum(), h_calo.GetMaximum(), h_hm.GetMaximum())
    # colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kGreen+2]

    def draw_lines():
        lines = {}
        pdgs = [211, 321, 2212]
        m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
        for pdg in pdgs:
            lines[pdg] = ROOT.TLine(m_pdg[pdg], 0., m_pdg[pdg], max_z)
            lines[pdg].SetLineColor(8)
            lines[pdg].SetLineWidth(2)
            lines[pdg].SetLineStyle(9)
            lines[pdg].Draw()
        return lines


    canvas = ROOT.TCanvas()

    h_ip.Draw()
    lines1 = draw_lines()
    h_ip.SetLineColor(1)
    h_ip.SetStats(0)
    h_ip.SetMinimum(0.1)
    h_ip.SetMaximum(max_z)

    h_calo.Draw("same")
    h_calo.SetLineColor(2)
    h_calo.SetStats(0)

    h_hm.Draw("same")
    h_hm.SetLineColor(4)
    h_hm.SetStats(0)

    legend = canvas.BuildLegend()
    legend.SetFillColor(0)
    canvas.SetGridy(0)
    canvas.Update()
    input("wait")


def plot_var():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/hm/trk_len_hm*.root")
    df = df.Filter("tof > 0").Define("beta", "trackLength/(tof*299.792458)").Filter("beta >= 0 && beta <= 1")
    h = df.Histo2D(("h", "N Bad Z;N Total segments; N bad z segments", 1000, 0, 1000, 300, 0, 300), "nSegments", "nBadZ")
    canvas = ROOT.TCanvas()
    h.Draw("colz")
    palette = h.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.93)
    palette.SetX2NDC(0.95)
    canvas.Update()
    input("wait")


def check_outliers():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/hm/trk_len_hm*.root")
    df = df.Filter("tof > 0").Filter("abs(pdg) == 211") # can't be 0 or negative
    df = df.Define("beta", "trackLength/(tof*299.792458)").Filter("beta >= 0 && beta <= 1")
    df = df.Define("mass", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000") # in MeV
    df_bad = df.Filter(" mass >= 200.")

    # Binning for some plots, i made: (1000, -100, 1300.)

    h_total = df.Histo1D(("h_total", "All pions; N Bad Z segments; N entries", 500, 0, 40), "tof")
    h_bad = df_bad.Histo1D(("h_bad", "Bad pions; TOF [ns]; Fraction of bad pions", 500, 0, 40), "tof")

    h_bad.Divide(h_total.GetPtr())

    canvas = ROOT.TCanvas()

    h_bad.Draw()
    canvas.Update()
    input("wait")


def plot_ptpz():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/trk_len.root")
    df = df.Filter("massDefault != 0 && massTanL != 0 && massZ != 0 && abs(pdg) == 321")
    df = df.Define("pt", "hypot(mom_ip._x,mom_ip._y)")
    # Make binning consistent with old plots
    # df.AsNumpy(["pt", "mom_ip._z", "m"])

    h_default = df.Define("m", "massDefault*1000 - 0.493677*1000").Profile2D(("h_default", "Default; pt [GeV]; pz [GeV]", 600, 0, 15, 600, -15, 15, -300, 300), "pt","mom_ip._z", "m")
    h_tanl = df.Define("m", "massTanL*1000 - 0.493677*1000").Profile2D(("h_tanL", "TanL; pt [GeV]; pz [GeV]", 600, 0, 15, 600, -15, 15, -300, 300), "pt","mom_ip._z", "m")
    h_z = df.Define("m", "massZ*1000 - 0.493677*1000").Profile2D(("h_z", "Z; pt [GeV]; pz [GeV]", 600, 0, 15, 600, -15, 15, -300, 300), "pt","mom_ip._z", "m")


    canvas = ROOT.TCanvas()

    h_default.Draw("colz")
    h_default.SetStats(0)
    h_default.SetMinimum(-200)
    h_default.SetMaximum(200)
    canvas.Update()
    input("wait")

    canvas = ROOT.TCanvas()
    h_tanl.Draw("colz")
    h_tanl.SetStats(0)
    h_tanl.SetMinimum(-200)
    h_tanl.SetMaximum(200)
    canvas.Update()
    input("wait")

    canvas = ROOT.TCanvas()
    h_z.Draw("colz")
    h_z.SetStats(0)
    h_z.SetMinimum(-200)
    h_z.SetMaximum(200)
    canvas.Update()
    input("wait")


# algorithms = ["momentumIP" ,"momentumCalo" ,"momentumHM"]
# for alg in algorithms:
    # plot_2d(alg)

# for i in range(100):
#     plot_2d(n_seg1=i*10, n_seg2=(i+1)*10)

# plot_2d(n_seg1=0, n_seg2=240)

# plot_var()
check_outliers()