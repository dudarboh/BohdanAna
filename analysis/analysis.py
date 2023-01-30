import ROOT
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)

def plot_2d():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/trk_len.root")
    df = df.Filter("massDefault != 0 && massTanL != 0 && massZ != 0")
    # Make binning consistent with old plots

    h_default = df.Define("m", "massDefault*1000").Histo2D(("h_default", "Default; momentum [GeV]; Mass [GeV]", 30, 0, 15, 200, -100, 1300.), "momentum","m")
    h_tanl = df.Define("m", "massTanL*1000").Histo2D(("h_tanL", "TanL; momentum [GeV]; Mass [GeV]", 30, 0, 15, 200, -100, 1300.), "momentum","m")
    h_z = df.Define("m", "massZ*1000").Histo2D(("h_z", "Z; momentum [GeV]; Mass [GeV]", 30, 0, 15, 200, -100, 1300.), "momentum","m")


    max_z = 1.05*max(h_default.GetMaximum(), h_tanl.GetMaximum(), h_z.GetMaximum())
    # colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kGreen+2]

    def draw_lines():
        lines = {}
        pdgs = [211, 321, 2212]
        m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
        for pdg in pdgs:
            lines[pdg] = ROOT.TLine(0., m_pdg[pdg], 15., m_pdg[pdg])
            lines[pdg].SetLineColor(8)
            lines[pdg].SetLineWidth(3)
            lines[pdg].SetLineStyle(9)
            lines[pdg].Draw()
        return lines


    canvas = ROOT.TCanvas()
    canvas.Divide(2, 2)

    pad = canvas.cd(1)
    pad.SetLogz()
    h_default.Draw("colz")
    lines1 = draw_lines()
    h_default.SetLineColor(1)
    h_default.SetStats(0)
    h_default.SetMinimum(0.1)
    h_default.SetMaximum(max_z)

    pad = canvas.cd(2)
    pad.SetLogz()
    h_tanl.Draw("colz")
    lines2 = draw_lines()
    h_tanl.SetLineColor(2)
    h_tanl.SetStats(0)
    h_tanl.SetMinimum(0.1)
    h_tanl.SetMaximum(max_z)

    pad = canvas.cd(3)
    pad.SetLogz()
    h_z.Draw("colz")
    lines3 = draw_lines()
    h_z.SetLineColor(4)
    h_z.SetStats(0)
    h_z.SetMinimum(0.1)
    h_z.SetMaximum(max_z)

    canvas.BuildLegend()
    canvas.Update()
    input("wait")

def plot_1d():

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/track_length/trk_len.root")
    df = df.Filter("massDefault != 0 && massTanL != 0 && massZ != 0")
    # Make binning consistent with old plots

    h_default = df.Define("m", "massDefault*1000").Histo1D(("h_default", "Default; mass [MeV]; N entries", 1000, -100, 1300.), "m")
    h_tanl = df.Define("m", "massTanL*1000").Histo1D(("h_tanL", "TanL; mass [MeV]; N entries", 1000, -100, 1300.), "m")
    h_z = df.Define("m", "massZ*1000").Histo1D(("h_z", "Z; mass [MeV]; N entries", 1000, -100, 1300.), "m")


    max_z = 1.05*max(h_default.GetMaximum(), h_tanl.GetMaximum(), h_z.GetMaximum())
    # colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kGreen+2]

    def draw_lines():
        lines = {}
        pdgs = [211, 321, 2212]
        m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
        for pdg in pdgs:
            lines[pdg] = ROOT.TLine(m_pdg[pdg], 0., m_pdg[pdg], max_z)
            lines[pdg].SetLineColor(8)
            lines[pdg].SetLineWidth(3)
            lines[pdg].SetLineStyle(9)
            lines[pdg].Draw()
        return lines


    canvas = ROOT.TCanvas()

    h_default.Draw()
    # lines1 = draw_lines()
    h_default.SetLineColor(1)
    h_default.SetStats(0)
    h_default.SetMinimum(0.1)
    h_default.SetMaximum(max_z)

    h_tanl.Draw("same")
    h_tanl.SetLineColor(2)
    h_tanl.SetStats(0)

    h_z.Draw("same")
    h_z.SetLineColor(4)
    h_z.SetStats(0)

    canvas.BuildLegend()
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


plot_ptpz()
