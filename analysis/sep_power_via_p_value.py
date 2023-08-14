import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)
ROOT.EnableImplicitMT()

colors = ['#1b9e77', '#d95f02', '#7570b3']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

def get_p_value(h1, h2):
    '''
        h1 should be to the left, (let's assume signal)
        h2 should be to the right (let's assume background)
    '''
    # as bins in mass: 500, -100, 1300.
    # i changed it to 3000 -0 6000 , because with 100 ps ALL masses move to higher values
    efficiencies = []
    mis_ids = []
    cut_mass = []
    diff = []
    n_signal = h1.GetEntries()
    n_background = h2.GetEntries()

    for cut in range(1, 3000):
        eff = h1.Integral(0, cut)/n_signal
        mis_id = h2.Integral(0, cut)/n_background

        cut_mass.append( h1.GetXaxis().GetBinLowEdge(cut+1)  )
        efficiencies.append(eff)
        mis_ids.append(mis_id)
        diff.append( abs(1-eff - mis_id) )

    cut_mass = np.array(cut_mass)
    efficiencies = np.array(efficiencies)
    mis_ids = np.array(mis_ids)
    diff = np.array(diff)

    p_value = mis_ids[ np.argmin(diff) ]
    sep_power = 2*ROOT.Math.gaussian_quantile_c(p_value, 1)
    print(f"Found cut at mass: {round(cut_mass[np.argmin(diff)], 2)}\
            with (efficiency): {round(1-efficiencies[np.argmin(diff)], 2)} and mis id: {round(p_value, 2)} and sep power {round(sep_power, 2)}")
    return sep_power


# df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/dEdx/final.root")
df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna2.root")
df = df.Filter("tofClosest0 > 6.")
n_mom_bins = 70
n_mass_bins = 3000

def get_sep_power_graph(df, tof_column="tofClosest0"):

    df_total = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_total = df_total.Histo2D(("h_total", "All; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    df_pion = df.Filter("abs(pdg) == 211")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_pion = df_pion.Histo2D(("h_pion", "#pi; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    df_kaon = df.Filter("abs(pdg) == 321")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_kaon = df_kaon.Histo2D(("h_kaon", "K; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    df_proton = df.Filter("abs(pdg) == 2212")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_proton = df_proton.Histo2D(("h_proton", "proton; momentum [GeV]; Mass [MeV]", n_mom_bins, 0, 20, n_mass_bins, 0, 6000.), "harmonicMomToEcal_IKF_zedLambda","mass")

    # DRAW FANCY 2D plot
    # ROOT.gStyle.SetPadRightMargin(0.12)
    # canvas = ROOT.TCanvas("c_2d_m_vs_p_total",
    #                         "",
    #                         int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
    #                         int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    # h_2d_total.Draw("colz")

    # h_2d_total.SetMinimum(1)
    # h_2d_total.SetMaximum(100000)
    # canvas.SetLogz()
    # canvas.SetGridx(0)
    # canvas.SetGridy(0)
    # canvas.Update()
    # palette = h_2d_total.GetListOfFunctions().FindObject("palette")
    # palette.SetX1NDC(0.89)
    # palette.SetX2NDC(0.91)
    # canvas.Modified()
    # canvas.Update()

    # canvas2 = ROOT.TCanvas("c_proj_m_vs_p",
    #                         "",
    #                         int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
    #                         int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )

    gr_sp_pik = ROOT.TGraph()
    gr_sp_kp = ROOT.TGraph()
    gr_sp_pik.SetTitle(f"{tof_column}; momentum [GeV]; #pi/K separation power")
    gr_sp_kp.SetTitle("title; momentum [GeV]; K/p separation power")
    # gr_sp_pik.SetPoint(0, 0, 0)
    # gr_sp_kp.SetPoint(0, 0, 0)

    for i in range(1, n_mom_bins):
        h_proj_pion = h_2d_pion.ProjectionY("h_pion_projection", i, i)
        h_proj_kaon = h_2d_kaon.ProjectionY("h_kaon_projection", i, i)
        h_proj_proton = h_2d_proton.ProjectionY("h_proton_projection", i, i)

        print(f"Computing for momentum range: {round(h_2d_total.GetXaxis().GetBinLowEdge(i), 2)} to {round(h_2d_total.GetXaxis().GetBinLowEdge(i+1), 2)}")
        gr_sp_pik.SetPoint(i, h_2d_total.GetXaxis().GetBinCenter(i), get_p_value(h_proj_pion, h_proj_kaon) )
        gr_sp_kp.SetPoint(i, h_2d_total.GetXaxis().GetBinCenter(i), get_p_value(h_proj_kaon, h_proj_proton))

        # h_proj_pion.Draw()
        # h_proj_pion.SetLineColor(colors[0])
        # h_proj_kaon.Draw("Lsame")
        # h_proj_kaon.SetLineColor(colors[1])
        # h_proj_proton.Draw("Lsame")
        # h_proj_proton.SetLineColor(colors[2])
        # canvas2.SetLogy()
        # canvas2.Update()

    # canvas3 = ROOT.TCanvas("c_sep_power",
    #                         "",
    #                         int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
    #                         int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    return gr_sp_kp
    # gr_sp_pik.Draw("APL")
    # gr_sp_kp.Draw("PLsame")
    # gr_sp_kp.SetLineColor(4)
    # canvas3.Update()
    # input("wait")


# get_sep_power_graph(df, tof_column="tofClosest90")
# input("EEEEEEEEEND")

canvas = ROOT.TCanvas("c",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )

canvas.SetGridx(True)
canvas.SetGridy(True)

gr_sp={}
colors = ["#03045e","#023e8a","#0077b6","#0096c7","#00b4d8","#48cae4"]
# colors = ["#690000", "#850e0f", "#a21d19", "#c02b25", "#df3831", "#ff463d"]
colors = [ ROOT.TColor.GetColor(c) for c in colors[::-1]]

legend = ROOT.TLegend()
for i, res in enumerate( [0, 10, 30, 50, 70, 90] ):
    gr_sp[res] = get_sep_power_graph(df, tof_column=f"tofClosest{int(res)}")
    gr_sp[res].Draw("AL" if i == 0 else "Lsame")
    gr_sp[res].SetLineColor(colors[i])
    gr_sp[res].SetLineWidth(4)
    legend.AddEntry(gr_sp[res], f"{res} ps","l")

gr_sp[0].GetXaxis().SetRangeUser(0, 19)
gr_sp[0].GetYaxis().SetRangeUser(0, 6)
legend.Draw()
canvas.Update()
input("wait")
