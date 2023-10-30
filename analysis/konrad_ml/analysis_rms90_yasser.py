import ROOT
import numpy as np

def rms90(hist):
    hist0 = hist.Clone()
    nBins = hist.GetNbinsX()

    integralHist = -10.0
    IntMean90 = 0.9
    IntOneSigma = 1.0 - np.exp( -1.0 )
    MeanFirstBin = 1
    MeanLastBin = 2
    localIntegral = 0.0
    localMean = 0.0
    localMean2 = 0.0
    localRMS = 0.0
    RMS90 = 1000000.0
    RMSOneSigma = 1000000.0

    for firstBin in range(1, nBins - 1):
        localMean = 0.0
        localMean2 = 0.0
        MeanFirstBin = firstBin
        MeanLastBin = firstBin
        for lastBin in range(firstBin, nBins):
            if ( hist.Integral( firstBin , lastBin ) <= IntMean90 ) and ( hist.Integral( firstBin , lastBin + 1 ) > IntMean90 ):
                MeanLastBin = lastBin

        if not ( ( hist.Integral( MeanFirstBin , MeanLastBin ) <= IntMean90 ) and ( hist.Integral( MeanFirstBin , MeanLastBin + 1 ) > IntMean90 ) ):
            continue

        localIntegral = hist.Integral( MeanFirstBin , MeanLastBin )
        hist0.GetXaxis().SetRangeUser( hist0.GetBinLowEdge( MeanFirstBin ) , hist0.GetBinLowEdge( MeanLastBin ) + hist0.GetBinWidth( MeanLastBin ) )
        localMean = hist0.GetMean()
        localRMS = hist0.GetStdDev()
        if ( hist0.GetBinLowEdge( MeanLastBin ) + hist0.GetBinWidth( MeanLastBin ) - hist0.GetBinLowEdge( MeanFirstBin ) ) < lowestRange :
        {
            lowestRange = hist0.GetBinLowEdge( MeanLastBin ) + hist0.GetBinWidth( MeanLastBin ) - hist0.GetBinLowEdge( MeanFirstBin );
            Mean90 = hist0.GetMean();
            RMS90 = hist0.GetStdDev();
            Mean90Err = hist0.GetMeanError();
            RMS90Err = hist0.GetStdDevError();
            lowerLimit90 = hist.GetXaxis().GetBinLowEdge( MeanFirstBin );
            upperLimit90 = hist.GetXaxis().GetBinLowEdge( MeanLastBin ) + hist.GetBinWidth( MeanLastBin );
            lowerBin90 = MeanFirstBin;
            upperBin90 = MeanLastBin + 1;
        }
    }
    MeanFirstBin = 1;
    MeanLastBin = 2;
    for ( int firstBin = 1 ; firstBin < nBins - 1 ; ++firstBin )
    {
        localMean = 0.0;
        localMean2 = 0.0;
        for ( int lastBin = firstBin ; lastBin < nBins ; ++lastBin )
        {
        
            if ( ( hist.Integral( firstBin , lastBin ) <= IntOneSigma ) && ( hist.Integral( firstBin , lastBin + 1 ) > IntOneSigma ) )
            {            
                MeanFirstBin = firstBin;
                MeanLastBin = lastBin;
            }
        }
        localIntegral = hist.Integral( MeanFirstBin , MeanLastBin );
        for ( int iBin = MeanFirstBin ; iBin <= MeanLastBin ; ++iBin )
        {
            localMean += ( hist.GetBinContent( iBin ) ) * ( hist.GetXaxis().GetBinLowEdge( iBin ) + ( hist.GetXaxis().GetBinWidth( iBin ) / 2.0 ) ) / localIntegral;
            localMean2 += ( hist.GetBinContent( iBin ) ) * pow( hist.GetXaxis().GetBinLowEdge( iBin ) + ( hist.GetXaxis().GetBinWidth( iBin ) / 2.0 ) , 2 ) / localIntegral;
        }
        localRMS = sqrt( localMean2 - pow( localMean , 2 ) );
    
        if ( localRMS < RMSOneSigma )
        {
            MeanOneSigma = localMean;
            RMSOneSigma = localRMS;
            MeanOneSigmaErr = fabs( MeanOneSigma ) / sqrt( hist.GetEntries() );
            RMSOneSigmaErr = RMSOneSigma / sqrt( hist.GetEntries() );
            lowerLimitOneSigma = hist.GetXaxis().GetBinLowEdge( MeanFirstBin );
            upperLimitOneSigma = hist.GetXaxis().GetBinLowEdge( MeanLastBin ) + hist.GetXaxis().GetBinWidth( MeanLastBin );        
        }
    }
}


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
