import ROOT
import numpy as np

canvas = ROOT.TCanvas()
canvas.SetGridx(True)
canvas.SetGridy(True)

data = -np.load("classical_algorithm_10_layers.npy")

n_bins = 200
x_min = -120
x_max = 120

# n_bins = 400
# x_min = -300
# x_max = 300


h_average = ROOT.TH1F("Analytic average", "Analytic average; reconstructed TOF - true TOF [ps]; N entries", n_bins, x_min, x_max)
h_average.FillN( len(data), data, np.ones_like(data) )
h_average.SetLineColor(ROOT.kRed)
h_average.SetMarkerColor(ROOT.kRed)

gaus_average = ROOT.TF1("g_average", "gaus", x_min, x_max)
gaus_average.SetNpx(500)
gaus_average.SetTitle("Gaussian fit")
gaus_average.SetParNames("Fit Const", "Fit #mu", "Fit #sigma")
gaus_average.SetLineColor(ROOT.kRed+1)

h_average.Fit(gaus_average, "0", "", -20, 20)
h_average.Draw()
gaus_average.Draw("same")

data = -np.load("epic_regression_10_layers_training_10_layers.npy")
h_ml = ROOT.TH1F("ML approach", "ML approach; reconstructed TOF - true TOF [ps]; N entries", n_bins, x_min, x_max)
h_ml.FillN(len(data), data, np.ones_like(data))
h_ml.SetLineColor(ROOT.kBlue)
h_ml.SetMarkerColor(ROOT.kBlue)


gaus_ml = ROOT.TF1("g_ml", "gaus", x_min, x_max)
gaus_ml.SetTitle("Gaussian fit")
gaus_ml.SetNpx(500)
gaus_ml.SetParNames("Fit Const", "Fit #mu", "Fit #sigma")
gaus_ml.SetLineColor(ROOT.kBlue +1)

h_ml.Fit(gaus_ml, "0", "", -15, 5)
h_ml.Draw("same")
gaus_ml.Draw("same")

h_average.GetYaxis().SetRangeUser(0., 1.3*max(h_average.GetMaximum(), h_ml.GetMaximum()))

leg_average = ROOT.TLegend(0.19, 0.65, 0.58, 0.94)
leg_average.SetBorderSize(0)
leg_average.SetFillStyle(0)
leg_average.SetMargin(0.1)
entry1 = leg_average.AddEntry(h_average, h_average.GetTitle(), "l")
entry2 = leg_average.AddEntry("", f"Histo mean: {h_average.GetMean():.2f}", "")
entry3 = leg_average.AddEntry("", f"Histo std. dev.: {h_average.GetStdDev():.2f}", "")
entry4 = leg_average.AddEntry(gaus_average, gaus_average.GetTitle(), "l")
entry5 = leg_average.AddEntry("", f"Fit #mu: {gaus_average.GetParameter(1):.2f}", "")
entry6 = leg_average.AddEntry("", f"Fit #sigma: {gaus_average.GetParameter(2):.2f}", "")
entry1.SetTextColor(ROOT.kRed)
entry2.SetTextColor(ROOT.kRed)
entry3.SetTextColor(ROOT.kRed)
entry4.SetTextColor(ROOT.kRed+1)
entry5.SetTextColor(ROOT.kRed+1)
entry6.SetTextColor(ROOT.kRed+1)
leg_average.Draw()

leg_ml = ROOT.TLegend(0.595, 0.65, 0.985, 0.94)
leg_ml.SetBorderSize(0)
leg_ml.SetFillStyle(0)
leg_ml.SetMargin(0.1)
entry1 = leg_ml.AddEntry(h_ml, h_ml.GetTitle(), "l")
entry2 = leg_ml.AddEntry("", f"Histo mean: {h_ml.GetMean():.2f}", "")
entry3 = leg_ml.AddEntry("", f"Histo std. dev.: {h_ml.GetStdDev():.2f}", "")
entry4 = leg_ml.AddEntry(gaus_ml, gaus_ml.GetTitle(), "l")
entry5 = leg_ml.AddEntry("", f"Fit #mu: {gaus_ml.GetParameter(1):.2f}", "")
entry6 = leg_ml.AddEntry("", f"Fit #sigma: {gaus_ml.GetParameter(2):.2f}", "")
entry1.SetTextColor(ROOT.kBlue)
entry2.SetTextColor(ROOT.kBlue)
entry3.SetTextColor(ROOT.kBlue)
entry4.SetTextColor(ROOT.kBlue+1)
entry5.SetTextColor(ROOT.kBlue+1)
entry6.SetTextColor(ROOT.kBlue+1)
leg_ml.Draw()

latex = ROOT.TLatex()
latex.SetTextFont(52)
latex.DrawLatex(-4.34, 1640., "ILD work in progress")

canvas.Update()
input("Finish")
