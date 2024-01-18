import ROOT
from random import choice
from string import ascii_letters

def get_rand_string():
    return ''.join(choice(ascii_letters) for i in range(16))

def draw_2d_plot(h, maximum=1e4):
    '''
    Draw a 2D histogram in the appropriate styling and pause.
    '''
    margin = 0.33
    ROOT.gStyle.SetPadLeftMargin(0.6*margin)
    ROOT.gStyle.SetPadRightMargin(0.4*margin)
    ROOT.gStyle.SetPadTopMargin(0.35*margin)
    ROOT.gStyle.SetPadBottomMargin(0.65*margin)
    canvas = ROOT.TCanvas(get_rand_string(),"",600,600)

    h.Draw("colz")
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleOffset(1.4)

    print(f"The highest bin is: {h.GetMaximum()}", end=" ")
    h.SetMinimum(1)
    h.SetMaximum(maximum)
    print(f", the maximum is set to: {h.GetMaximum()}")

    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.Update()

    palette = h.GetListOfFunctions().FindObject("palette")
    x_min = h.GetXaxis().GetXmin()
    x_max = h.GetXaxis().GetXmax()
    palette.SetX1(x_min + 1.01*(x_max-x_min))
    palette.SetX2(x_min + 1.05*(x_max-x_min))
    palette.SetLabelOffset(0.001)

    canvas.Modified()
    canvas.Update()
    return canvas
