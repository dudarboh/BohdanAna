import ROOT
ROOT.EnableImplicitMT()

class Algorithm:
    def __init__(self, name, momentum, track_length, tof):
        self.name = name
        self.mom = momentum
        self.len = track_length
        self.tof = tof

colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]

algoClosest0 = Algorithm("Closest 0 ps", "harmonicMomToEcal", "trackLengthToEcal", "tofClosest0")
algoClosest30 = Algorithm("Closest 30 ps", "harmonicMomToEcal", "trackLengthToEcal", "tofClosest30")
algoSET50 = Algorithm("SET 50 ps", "harmonicMomToSET", "trackLengthToSET", "tofSET50")
algoAverage100 = Algorithm("Average 100 ps", "harmonicMomToEcal", "trackLengthToEcal", "tofAverage100")
# algorithms = [algoClosest0, algoClosest30, algoSET50, algoAverage100]
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
df = df.Filter("tofClosest0 > 6. && tofSET0 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")

histos = []
for alg in algorithms:
    df_beta = df.Define("beta", f"{alg.len}/({alg.tof}*299.792458)").Filter("beta >= 0 && beta <= 1")
    df_mass = df_beta.Define("mass", f"{alg.mom}*sqrt( 1./(beta*beta) - 1.)*1000")

    h = df_mass.Histo1D((f"{alg.name}", f"{alg.name}; mass [MeV]; N entries", 2000, 0, 1100), "mass")
    histos.append(h)


canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.18)
maxy = 0
for i, h in enumerate(histos):
    h.GetYaxis().SetTitleOffset(1.4)
    h.SetLineColor(colors[i])
    h.SetMarkerColor(colors[i])
    h.SetMarkerStyle(0)
    h.SetLineWidth(3)
    maxy = max( maxy, 1.05*h.GetMaximum() )


histos[0].SetMaximum(1.05*maxy)
histos[0].Draw() # draw axis for lines, god bless...
lines = draw_lines_1d(1.05*maxy)

legend = ROOT.TLegend(0.6, 0.7, 0.93, 0.93)
legend.AddEntry(lines[211], "True #pi/K/p mass", "l")

for i, h in enumerate(histos):
    h.Draw("same")
    legend.AddEntry(h.GetPtr())

legend.Draw()
canvas.Update()
input("Finish")