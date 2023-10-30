import ROOT
import numpy as np

def interval_quantile_(x, quant=0.9):
    '''Calculate the shortest interval that contains the desired quantile'''
    # the number of possible starting points
    n_low = int(len(x) * (1. - quant))
    # the number of events contained in the quantil
    n_quant = len(x) - n_low

    # Calculate all distances in one go
    distances = x[-n_low:] - x[:n_low]
    i_start = np.argmin(distances)

    return i_start, i_start + n_quant


def fit90(x): 
    x = np.sort(x)
    n10percent = int(round(len(x)*0.1))
    n90percent = len(x) - n10percent
    
    start, end = interval_quantile_(x, quant=0.9)
    
    rms90 = np.std(x[start:end])
    mean90 = np.mean(x[start:end])
    mean90_err = rms90/np.sqrt(n90percent)
    rms90_err = rms90/np.sqrt(2*n90percent)   # estimator in root
    return mean90, rms90, mean90_err, rms90_err


canvas = ROOT.TCanvas()
canvas.SetGridx(True)
canvas.SetGridy(True)

data = -np.load("classical_algorithm_10_layers.npy")
mean90_average, rms90_average, _, _ = fit90(data)
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
h_average.Draw()

data = -np.load("epic_regression_10_layers_training_10_layers.npy")
mean90_ml, rms90_ml, _, _ = fit90(data)
h_ml = ROOT.TH1F("ML approach", "ML approach; reconstructed TOF - true TOF [ps]; N entries", n_bins, x_min, x_max)
h_ml.FillN(len(data), data, np.ones_like(data))
h_ml.SetLineColor(ROOT.kBlue)
h_ml.SetMarkerColor(ROOT.kBlue)
h_ml.Draw("same")

h_average.GetYaxis().SetRangeUser(0., 1.3*max(h_average.GetMaximum(), h_ml.GetMaximum()))

leg_average = ROOT.TLegend(0.17, 0.65, 0.56, 0.94)
leg_average.SetBorderSize(0)
leg_average.SetFillStyle(0)
leg_average.SetMargin(0.1)
entry1 = leg_average.AddEntry(h_average, h_average.GetTitle(), "l")
entry2 = leg_average.AddEntry("", "Mean_{90}: " + f"{mean90_average:.2f}", "")
entry3 = leg_average.AddEntry("", "RMS_{90}: " + f"{rms90_average:.2f}", "")
entry1.SetTextColor(ROOT.kRed)
entry2.SetTextColor(ROOT.kRed)
entry3.SetTextColor(ROOT.kRed)
leg_average.Draw()

leg_ml = ROOT.TLegend(0.57, 0.65, 0.96, 0.94)
leg_ml.SetBorderSize(0)
leg_ml.SetFillStyle(0)
leg_ml.SetMargin(0.1)
entry1 = leg_ml.AddEntry(h_ml, h_ml.GetTitle(), "l")
entry2 = leg_ml.AddEntry("", "Mean_{90}: " + f"{mean90_ml:.2f}", "")
entry3 = leg_ml.AddEntry("", "RMS_{90}: " + f"{rms90_ml:.2f}", "")
entry1.SetTextColor(ROOT.kBlue)
entry2.SetTextColor(ROOT.kBlue)
entry3.SetTextColor(ROOT.kBlue)
leg_ml.Draw()

latex = ROOT.TLatex()
latex.SetTextFont(52)
latex.DrawLatex(-4.34, 1640., "ILD work in progress")

canvas.Update()
input("Finish")
