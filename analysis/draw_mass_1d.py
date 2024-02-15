import ROOT
ROOT.EnableImplicitMT()
from utils import *

def draw_lines_1d(maxy):
    lines = {}
    for p in particles:
        lines[p] = ROOT.TLine(p.mass2, 0., p.mass2, maxy)
        lines[p].SetLineColor(15)
        lines[p].SetLineWidth(2)
        lines[p].SetLineStyle(9)
        lines[p].Draw()
    return lines


# total range
n_bins, x_min, x_max = 500, -0.3, 1.2

#pion peak
n_bins, x_min, x_max = 500, -10e-3, 50e-3
#kaon peak
# n_bins, x_min, x_max = 500, 0.2, 0.28
# #proton peak
# n_bins, x_min, x_max = 500, 0.8, 0.95


def main():
    df_init = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
                  .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
                  .Filter("tofClosest0 > 6.")

    h_phi_lambda = df_init.Define("beta", "trackLengthToEcal_SHA_phiLambda_IP/(tofClosest0*299.792458)")\
                .Define("mass2", "(recoIpPx*recoIpPx+recoIpPy*recoIpPy+recoIpPz*recoIpPz)*( 1./(beta*beta) - 1.)")\
                .Histo1D((get_rand_string(), "", n_bins, x_min, x_max),"mass2")
                # .Define("mass", "sqrt(mass2)*1000")\

    h_phi_zed = df_init.Define("beta", "trackLengthToEcal_SHA_phiZed_IP/(tofClosest0*299.792458)")\
                .Define("mass2", "(recoIpPx*recoIpPx+recoIpPy*recoIpPy+recoIpPz*recoIpPz)*( 1./(beta*beta) - 1.)")\
                .Histo1D((get_rand_string(), "", n_bins, x_min, x_max),"mass2")
                # .Define("mass", "sqrt(mass2)*1000")\

    h_zed_lambda = df_init.Define("beta", "trackLengthToEcal_SHA_zedLambda_IP/(tofClosest0*299.792458)")\
                .Define("mass2", "(recoIpPx*recoIpPx+recoIpPy*recoIpPy+recoIpPz*recoIpPz)*( 1./(beta*beta) - 1.)")\
                .Histo1D((get_rand_string(), "", n_bins, x_min, x_max),"mass2")
                # .Define("mass", "sqrt(mass2)*1000")\

    canvas = create_canvas(margin=0.33, left_margin_fraction=0.7, bottom_margin_fraction=0.7)
    legend = create_legend()
    legend.SetTextFont(42)
    h_phi_lambda.Draw("L")
    h_phi_lambda.SetLineColor(ROOT.TColor.GetColor('#000000'))
    h_phi_lambda.SetLineWidth(3)
    h_phi_lambda.SetTitle(";Mass^{2} (GeV^{2}/c^{4}); N entries")
    h_phi_lambda.GetXaxis().SetMaxDigits(3)
    h_phi_lambda.GetXaxis().SetTitleOffset(1.1)
    h_phi_lambda.GetYaxis().SetMaxDigits(3)
    h_phi_lambda.SetMaximum( 1.05*max(h_phi_lambda.GetMaximum(), h_phi_zed.GetMaximum(), h_zed_lambda.GetMaximum()) )
    legend.AddEntry(h_phi_lambda.GetPtr(), "#frac{#Delta#varphi}{#Omega}#sqrt{1+tan^{2}#lambda}", "l")

    h_phi_zed.Draw("L same")
    h_phi_zed.SetLineWidth(3)
    h_phi_zed.SetLineColor(ROOT.TColor.GetColor('#688E26'))
    legend.AddEntry(h_phi_zed.GetPtr(), "#sqrt{#frac{#Delta#varphi}{#Omega}^{2} + #Deltaz^{2}}", "l")

    h_zed_lambda.Draw("L same")
    h_zed_lambda.SetLineWidth(3)
    h_zed_lambda.SetLineColor(ROOT.TColor.GetColor('#FAA613'))
    legend.AddEntry(h_zed_lambda.GetPtr(), "#frac{#Deltaz}{tan#lambda}#sqrt{1+tan^{2}#lambda}", "l")

    legend.Draw()
    lines = draw_lines_1d(h_phi_lambda.GetMaximum())
    legend.AddEntry(lines[pion], "true #pi^{#pm} mass", "l")

    canvas.Update()
    input("wait")


main()