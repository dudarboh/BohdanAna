import ROOT
ROOT.EnableImplicitMT()

class Algorithm:
    def __init__(self, name, momentum, track_length, tof):
        self.name = name
        self.mom = momentum
        self.len = track_length
        self.tof = tof

colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]

algoClosest30 = Algorithm("Single ECAL hit (30 ps)", "harmonicMomToEcal", "trackLengthToEcal", "tofClosest30")
algoSET50 = Algorithm("Two SET strips (50 ps)", "harmonicMomToSET", "trackLengthToSET", "tofSET50")
algoAverage100 = Algorithm("Ten ECAL hits (100 ps)", "harmonicMomToEcal", "trackLengthToEcal", "tofAverage100")
algorithms = [algoClosest30, algoSET50, algoAverage100]

def draw_lines_1d(maxy):
    lines = {}
    pdgs = [211, 321, 2212]
    m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
    for pdg in pdgs:
        lines[pdg] = ROOT.TLine(m_pdg[pdg], 0., m_pdg[pdg], maxy)
        lines[pdg].SetLineColor(15)
        lines[pdg].SetLineWidth(2)
        lines[pdg].SetLineStyle(9)
        lines[pdg].Draw()
    return lines


df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6. && tofSET0 > 0.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")

histos = []
for alg in algorithms:
    # NOTE beta quality cut significantly distords TOF gaussian distributions! (makes sense...)
    # df_beta = df.Define("beta", f"{alg.len}/({alg.tof}*299.792458)").Filter("beta >= 0 && beta <= 1")
    tof_true = "tofClosest0" if "SET" not in alg.tof else "tofSET0"
    df_res = df.Define("res", f"({alg.tof} - {tof_true})*1000")
    h = df_res.Histo1D((f"{alg.name}", f"{alg.name}; #Delta TOF (ps); N entries", 1500, -300, 300), "res")
    histos.append(h)


canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.18)
canvas.SetRightMargin(0.05)
canvas.SetGridx(True)
maxy = 0
for i, h in enumerate(histos):
    h.GetYaxis().SetTitleOffset(1.6)
    h.SetLineColor(colors[i])
    h.SetMarkerColor(colors[i])
    h.SetMarkerStyle(0)
    h.SetLineWidth(3)
    maxy = max( maxy, 1.25*h.GetMaximum() )


histos[0].SetMaximum(maxy)
histos[0].Draw() # draw axis for lines, god bless...

legend = ROOT.TLegend(0.21, 0.81, 0.6, 0.94)

for i, h in enumerate(histos):
    h.Draw("same")
    legend.AddEntry(h.GetPtr())

legend.Draw()

fits = []
legend_fits = ROOT.TLegend(0.66, 0.81, 0.85, 0.94)
for i, h in enumerate(histos):
    f = ROOT.TF1(f"f_{i}", "gaus", -50, 50)
    f.SetLineColor(colors[i])
    f.SetLineWidth(2)
    h.Fit(f, "RN")
    fits.append(f)
    f.SetNpx(1000)
    f.Draw("same")
    legend_fits.AddEntry(f, "#sigma_{fit} = " + f"{round(f.GetParameter(2), 1)}" + " ps", "l")


legend_fits.Draw("same")
canvas.Update()
input("Finish")
