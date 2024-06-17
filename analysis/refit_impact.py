from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

# refittedOmegaIP
# refittedTanLambdaIP
# refittedD0IP
# refittedZ0IP
# refittedPhiIP
# refittedOmegaErrIP
# refittedTanLambdaErrIP
# refittedD0ErrIP
# refittedZ0ErrIP
# refittedPhiErrIP
# refittedRecoIpPx
# refittedRecoIpPy
# refittedRecoIpPz


def refit_impact(variable="omega", particle=kaon):
    ROOT.gStyle.SetOptTitle(True)
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
            .Filter("tofClosest0 > 6.")

    if particle == kaon:
        df = df.Filter("abs(pdg) == 321")
    elif particle == proton:
        df = df.Filter("abs(pdg) == 2212")
    else:
        return 0

    df = df.Define("pTrue", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    if variable == "omega":
        title_default = ";Momentum (GeV/c);#Omega_{default} - #Omega_{true} (1/mm)"
        title_refitted = ";Momentum (GeV/c);#Omega_{refitted} - #Omega_{true} (1/mm)"
        title_default_pull = ";(#Omega_{default} - #Omega_{true}) / #sigma^{default}_{#Omega};N entries"
        title_refitted_pull = ";(#Omega_{refitted} - #Omega_{true}) / #sigma^{refitted}_{#Omega};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.0001, 0.0001

        y_column_default = "omegaIP - omegaTrue"
        y_column_refitted = "refittedOmegaIP - omegaTrue"
        y_column_default_pull = "(omegaIP - omegaTrue)/omegaErrIP"
        y_column_refitted_pull = "(refittedOmegaIP - omegaTrue)/refittedOmegaErrIP"
    elif variable == "tanLambda":
        title_default = ";Momentum (GeV/c);tan#lambda_{default} - tan#lambda_{true}"
        title_refitted = ";Momentum (GeV/c);tan#lambda_{refitted} - tan#lambda_{true}"
        title_default_pull = ";(tan#lambda_{default} - tan#lambda_{true}) / #sigma^{default}_{tan#lambda};N entries"
        title_refitted_pull = ";(tan#lambda_{refitted} - tan#lambda_{true}) / #sigma^{refitted}_{tan#lambda};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.01, 0.01

        y_column_default = "tanLambdaIP - tanLambdaTrue"
        y_column_refitted = "refittedTanLambdaIP - tanLambdaTrue"
        y_column_default_pull = "(tanLambdaIP - tanLambdaTrue)/tanLambdaErrIP"
        y_column_refitted_pull = "(refittedTanLambdaIP - tanLambdaTrue)/refittedTanLambdaErrIP"
    elif variable == "d0":
        title_default = ";Momentum (GeV/c);d_{0}_{default} - d_{0}_{true} (mm)"
        title_refitted = ";Momentum (GeV/c);d_{0}_{refitted} - d_{0}_{true} (mm)"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.01, 0.01
        y_column_default = "d0IP - d0True"
        y_column_refitted = "refittedD0IP - d0True"


    h_default = df.Define("diff", y_column_default).Histo2D((get_rand_string(), title_default, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_refitted = df.Define("diff", y_column_refitted).Histo2D((get_rand_string(), title_refitted, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_default_pull = df.Define("diff", y_column_default_pull).Histo1D((get_rand_string(), title_default_pull, 200, -5, 5), "diff")
    h_refitted_pull = df.Define("diff", y_column_refitted_pull).Histo1D((get_rand_string(), title_refitted_pull, 200, -5, 5), "diff")

    h_default.GetYaxis().SetTitleOffset(1.8)
    h_refitted.GetYaxis().SetTitleOffset(1.8)
    h_refitted_pull.SetLineColor(ROOT.kRed)

    h_default.GetYaxis().SetTitleOffset(1.8)
    h_refitted.GetYaxis().SetTitleOffset(1.8)

    c1 = draw_2d_plot(h_default, 0.4, 0.65, 0.65)
    c2 = draw_2d_plot(h_refitted, 0.4, 0.65, 0.65)
    c3 = create_canvas()
    h_default_pull.Draw()
    h_refitted_pull.Draw("same")
    h_default_pull.GetYaxis().SetMaxDigits(3)
    h_default_pull.SetMinimum( 0. )
    h_default_pull.SetMaximum( 1.05*max(h_default_pull.GetMaximum(), h_refitted_pull.GetMaximum()) )
    latex.DrawLatexNDC(0.2, 0.8, "Std Dev: "+f"{h_default_pull.GetStdDev():.3f}")
    latex.DrawLatexNDC(0.2, 0.75, "#color[2]{Std Dev: " + f"{h_refitted_pull.GetStdDev():.3f}" + "}")
    c3.Update()

    input("wait")

refit_impact()