import ROOT
import numpy as np
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetLegendFillColor(0)
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

    legend = canvas.BuildLegend()
    legend.SetFillColor(0)
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

    legend = canvas.BuildLegend()
    legend.SetFillColor(0)
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


def plot_track_parameters(filename="evt_1_pfo_6"):
    file = ROOT.TFile(f"../build/{filename}.root", "READ")
    gr = {}
    for key in file.GetListOfKeys():
        if(key.GetName() == "omega"):
            gr_temp = key.ReadObj()
            x = np.array(gr_temp.GetX())
            y = abs( np.array(gr_temp.GetY()) )
            gr[key.GetName()] = ROOT.TGraph(len(x), x, y)
            gr[key.GetName()].SetTitle("#Omega")
        else:
            gr[key.GetName()] = key.ReadObj()            
        if "mc" in key.GetName():
            gr[key.GetName()].SetLineColor(2)

    mg = {}
    mg["phi"] = ROOT.TMultiGraph()
    mg["phi"].SetTitle("#varphi Reco vs MC;Hit number;#varphi")
    mg["phi"].Add(gr["phi"])
    mg["phi"].Add(gr["mc_phi"])
    mg["omega"] = ROOT.TMultiGraph()
    mg["omega"].SetTitle("#Omega Reco vs MC;Hit number;#Omega (1/mm)")
    mg["omega"].Add(gr["omega"])
    mg["omega"].Add(gr["mc_omega"])
    mg["pt"] = ROOT.TMultiGraph()
    mg["pt"].SetTitle("p_{T} Reco vs MC;Hit number;pt (GeV)")
    mg["pt"].Add(gr["pt"])
    mg["pt"].Add(gr["mc_pt"])
    mg["tanl"] = ROOT.TMultiGraph()
    mg["tanl"].SetTitle("tan #lambda Reco vs MC;Hit number;tan #lambda")
    mg["tanl"].Add(gr["tanl"])
    mg["tanl"].Add(gr["mc_tanl"])
    mg["pz"] = ROOT.TMultiGraph()
    mg["pz"].SetTitle("p_{z} Reco vs MC;Hit number;pz (GeV)")
    mg["pz"].Add(gr["pz"])
    mg["pz"].Add(gr["mc_pz"])
    mg["z"] = ROOT.TMultiGraph()
    mg["z"].SetTitle("z Reco vs MC;Hit number;z (mm)")
    mg["z"].Add(gr["z"])
    mg["z"].Add(gr["mc_z"])
    mg["trk_len"] = ROOT.TMultiGraph()
    mg["trk_len"].SetTitle("Cumulative track length TanL(blue) vs Default(red) vs \"Z\";Hit number;trk_len (mm)")
    mg["trk_len"].Add(gr["trk_len_tanl"])
    gr["trk_len_tanl"].SetLineColor(4)
    mg["trk_len"].Add(gr["trk_len_default"])
    gr["trk_len_default"].SetLineColor(2)
    mg["trk_len"].Add(gr["trk_len_z"])
    mg["trk_len_diff"] = ROOT.TMultiGraph()

    y ={}
    x = np.array(gr["trk_len_z"].GetX())
    y["z"] = np.array(gr["z"].GetY())
    y["mc_z"] = np.array(gr["mc_z"].GetY())
    y["phi"] = np.array(gr["z"].GetY())
    y["mc_phi"] = np.array(gr["mc_z"].GetY())
    y["trk_len_z"] = np.array(gr["trk_len_z"].GetY())
    y["trk_len_default"] = np.array(gr["trk_len_default"].GetY())
    y["trk_len_tanl"] = np.array(gr["trk_len_tanl"].GetY())

    gr["z_diff"] = ROOT.TGraph(len(x), x, y["z"] - y["mc_z"])
    gr["z_diff"].SetTitle("Z difference Reco - MC ;Hit number; #Delta Z (mm)")
    gr["phi_diff"] = ROOT.TGraph(len(x), x, y["phi"] - y["mc_phi"])
    gr["phi_diff"].SetTitle("#varphi difference Reco - MC ;Hit number; #Delta #varphi")

    gr["trk_len_diff_default"] = ROOT.TGraph(len(x), x, y["trk_len_default"] - y["trk_len_z"])
    gr["trk_len_diff_default"].SetLineColor(2)
    gr["trk_len_diff_default"].SetTitle("Diff track length DEFAULT")

    gr["trk_len_diff_tanl"] = ROOT.TGraph(len(x), x, y["trk_len_tanl"] - y["trk_len_z"])
    gr["trk_len_diff_tanl"].SetLineColor(4)
    gr["trk_len_diff_tanl"].SetTitle("Diff track length TANL")

    mg["trk_len_diff"].SetTitle("Track length diff TanL(Blue) and Default(Red) vs Z option;Hit number;#Delta trk_len (mm)")
    mg["trk_len_diff"].Add(gr["trk_len_diff_default"])
    mg["trk_len_diff"].Add(gr["trk_len_diff_tanl"])

    canvas = ROOT.TCanvas(f"{filename}", f"{filename}", 1280, 720)
    canvas.Divide(1, 3)

    pad = canvas.cd(1)
    mg["phi"].Draw("AL")

    pad = canvas.cd(2)
    mg["omega"].Draw("AL")

    pad = canvas.cd(3)
    mg["tanl"].Draw("AL")
    canvas.Update()
    input("wait")
    canvas.Print(f"{filename}_omega_tanl.png")

    pad = canvas.cd(2)
    mg["pt"].Draw("AL")

    pad = canvas.cd(3)
    mg["pz"].Draw("AL")
    canvas.Update()
    input("wait")
    canvas.Print(f"{filename}_mom.png")

    pad = canvas.cd(2)
    mg["phi"].Draw("AL")

    pad = canvas.cd(3)
    mg["z"].Draw("AL")
    canvas.Update()
    input("wait")
    canvas.Print(f"{filename}_phi_z.png")

    pad = canvas.cd(2)
    gr["phi_diff"].Draw("AL")

    pad = canvas.cd(3)
    gr["z_diff"].Draw("AL")
    canvas.Update()
    input("wait")
    canvas.Print(f"{filename}_phi_z_diff.png")

    pad = canvas.cd(2)
    mg["trk_len"].Draw("AL")

    pad = canvas.cd(3)
    mg["trk_len_diff"].Draw("AL")
    canvas.Update()
    input("wait")
    canvas.Print(f"{filename}_trk_len.png")    


filenames = ["evt_1_pfo_6", "evt_1_pfo_12", "evt_2_pfo_1", "evt_3_pfo_11", "evt_3_pfo_19", "evt_6_pfo_10"]
for filename in filenames:
    plot_track_parameters(filename)
