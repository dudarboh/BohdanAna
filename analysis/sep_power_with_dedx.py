import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)
ROOT.EnableImplicitMT()

def get_p_value(h1, h2):
    '''
        h1 should be to the left, (let's assume signal)
        h2 should be to the right (let's assume background)
    '''
    # as bins in mass: 500, -100, 1300.
    # i changed it to 3000 -0 6000 , because with 100 ps ALL masses move to higher values and truncation happens which gives wrong results at hig momentum
    efficiencies = []
    mis_ids = []
    cut_x = []
    diff = []
    n_signal = h1.GetEntries()
    n_background = h2.GetEntries()
    print()

    for cut in range(1, h1.GetXaxis().GetNbins()):
        eff = h1.Integral(0, cut)/n_signal
        mis_id = h2.Integral(0, cut)/n_background
        if eff < mis_id:
            eff, mis_id = mis_id, eff

        cut_x.append( h1.GetXaxis().GetBinLowEdge(cut+1)  )
        efficiencies.append(eff)
        mis_ids.append(mis_id)
        diff.append( abs(1-eff - mis_id) )

    cut_x = np.array(cut_x)
    efficiencies = np.array(efficiencies)
    mis_ids = np.array(mis_ids)
    diff = np.array(diff)

    p_value = mis_ids[ np.argmin(diff) ]
    sep_power = 2*ROOT.Math.gaussian_quantile_c(p_value, 1)
    print(f"Computed p-value for {h1.GetName()} and {h2.GetName()}")
    print(f"Found cut at x: {cut_x[np.argmin(diff)]}\
            with (efficiency): {round(1-efficiencies[np.argmin(diff)], 2)} and mis id: {round(p_value, 2)} and sep power {round(sep_power, 2)}")
    return sep_power

colors = ['#1b9e77', '#d95f02', '#7570b3']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6. && dEdx > 0.")
n_mass_bins = 3000

# with dE/dx we are better use log binning
n_mom_bins = 30
min_x = -0.3 # 0.1 GeV
max_x = 1.3 # ~ 20 GeV
mom_bins = [ 10**(min_x + (max_x-min_x)*i/n_mom_bins ) for i in range(n_mom_bins+1) ]

n_dedx_bins = 3000
min_y = 0
max_y = 1e-6
dedx_bins = [ min_y + (max_y - min_y)*i/n_dedx_bins for i in range(n_dedx_bins+1) ]

def get_sep_power_graph(df, tof_column="tofClosest30"):
    df = df.Define("mom", "harmonicMomToEcal_IKF_zedLambda")\
           .Define("trk_len", "trackLengthToEcal_IKF_zedLambda")\
           .Define("beta", f"trk_len/({tof_column}*299.792458)")\
           .Filter("beta >= 0 && beta <= 1")\
           .Define("mass", "mom*sqrt( 1./(beta*beta) - 1.)*1000")\

    ### TOF 2D HISTOGRAMS ###
    h_tof_total = df.Histo2D(("h_total", "All; momentum [GeV]; Mass [MeV]", n_mom_bins, np.array(mom_bins), 3000, 0, 6000.), "mom","mass")
    h_tof_pion = df.Filter("abs(pdg) == 211").Histo2D(("h_pion", "#pi; momentum [GeV]; Mass [MeV]", n_mom_bins, np.array(mom_bins), 3000, 0, 6000.), "mom","mass")
    h_tof_kaon = df.Filter("abs(pdg) == 321").Histo2D(("h_kaon", "K; momentum [GeV]; Mass [MeV]", n_mom_bins, np.array(mom_bins), 3000, 0, 6000.), "mom","mass")
    h_tof_proton = df.Filter("abs(pdg) == 2212").Histo2D(("h_proton", "p; momentum [GeV]; Mass [MeV]", n_mom_bins, np.array(mom_bins), 3000, 0, 6000.), "mom","mass")
    ### END ###

    ### dEdx 2D HISTOGRAMS ###
    h_dedx_pion = df.Filter("abs(pdg) == 211").Histo2D(("h_dedx_pion", "title; momentum [GeV]; dEdx [GeV/mm]", n_mom_bins, np.array(mom_bins), n_dedx_bins, np.array(dedx_bins)), "mom", "dEdx")
    h_dedx_kaon = df.Filter("abs(pdg) == 321").Histo2D(("h_dedx_kaon", "title; momentum [GeV]; dEdx [GeV/mm]", n_mom_bins, np.array(mom_bins), n_dedx_bins, np.array(dedx_bins)), "mom", "dEdx")
    h_dedx_proton = df.Filter("abs(pdg) == 2212").Histo2D(("h_dedx_proton", "title; momentum [GeV]; dEdx [GeV/mm]", n_mom_bins, np.array(mom_bins), n_dedx_bins, np.array(dedx_bins)), "mom", "dEdx")
    ### END ###

    # canvas2 = ROOT.TCanvas("c_proj_vs_p",
    #                         "",
    #                         int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
    #                         int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )


    ### Derive pik separation power for TOF, dE/dx and combine both ###
    gr_sp_pik = ROOT.TGraph()
    gr_sp_pik.SetTitle("TOF 30 ps; momentum [GeV]; #pi/K separation power")

    gr_sp_pik_dedx = ROOT.TGraph()
    gr_sp_pik_dedx.SetTitle("dE/dx; momentum [GeV]; #pi/K separation power")

    gr_sp_pik_combined = ROOT.TGraph()
    gr_sp_pik_combined.SetTitle("Combined TOF #oplus dE/dx; momentum [GeV]; #pi/K separation power")

    gr_sp_kp = ROOT.TGraph()
    gr_sp_kp.SetTitle("TOF 30 ps; momentum [GeV]; K/p separation power")

    gr_sp_kp_dedx = ROOT.TGraph()
    gr_sp_kp_dedx.SetTitle("dE/dx; momentum [GeV]; K/p separation power")

    gr_sp_kp_combined = ROOT.TGraph()
    gr_sp_kp_combined.SetTitle("Combined TOF #oplus dE/dx; momentum [GeV]; K/p separation power")



    for i in range(1, n_mom_bins):
        ### PI/K ###
        print(f"Computing for momentum range: {round(h_tof_total.GetXaxis().GetBinLowEdge(i), 2)} to {round(h_tof_total.GetXaxis().GetBinLowEdge(i+1), 2)}")
        x = h_tof_total.GetXaxis().GetBinCenter(i)
        h_proj_tof_pion = h_tof_pion.ProjectionY("h_pion_proj_tof", i, i)
        h_proj_tof_kaon = h_tof_kaon.ProjectionY("h_kaon_proj_tof", i, i)
        sp_tof = get_p_value(h_proj_tof_pion, h_proj_tof_kaon)
        gr_sp_pik.SetPoint(i-1, x, sp_tof )

        h_proj_dedx_pion = h_dedx_pion.ProjectionY("h_pion_proj_dedx", i, i)
        h_proj_dedx_kaon = h_dedx_kaon.ProjectionY("h_kaon_proj_dedx", i, i)
        sp_dedx = get_p_value(h_proj_dedx_pion, h_proj_dedx_kaon)
        gr_sp_pik_dedx.SetPoint(i-1, x, sp_dedx )

        gr_sp_pik_combined.SetPoint(i-1, x, np.sqrt( sp_tof*sp_tof + sp_dedx*sp_dedx) )
        ### END ###

        ### K/P ###
        h_proj_tof_kaon = h_tof_kaon.ProjectionY("h_kaon_proj_tof", i, i)
        h_proj_tof_proton = h_tof_proton.ProjectionY("h_proton_proj_tof", i, i)
        sp_tof = get_p_value(h_proj_tof_kaon, h_proj_tof_proton)
        gr_sp_kp.SetPoint(i-1, x, sp_tof )

        h_proj_dedx_kaon = h_dedx_kaon.ProjectionY("h_kaon_proj_dedx", i, i)
        h_proj_dedx_proton = h_dedx_proton.ProjectionY("h_proton_proj_dedx", i, i)
        sp_dedx = get_p_value(h_proj_dedx_kaon, h_proj_dedx_proton)
        gr_sp_kp_dedx.SetPoint(i-1, x, sp_dedx )

        gr_sp_kp_combined.SetPoint(i-1, x, np.sqrt( sp_tof*sp_tof + sp_dedx*sp_dedx) )
        ### END ###

        # h_proj_tof_pion.Draw()
        # h_proj_tof_pion.SetLineColor(colors[0])
        # h_proj_tof_kaon.Draw("Lsame")
        # h_proj_tof_kaon.SetLineColor(colors[1])
        # canvas2.SetLogy()
        # canvas2.Update()
        # input("wait")

        # h_proj_dedx_pion.Draw()
        # h_proj_dedx_pion.SetLineColor(colors[0])
        # h_proj_dedx_kaon.Draw("Lsame")
        # h_proj_dedx_kaon.SetLineColor(colors[1])
        # canvas2.SetLogy()
        # canvas2.Update()
        # input("wait")



    ROOT.gStyle.SetOptTitle(0)
    canvas3 = ROOT.TCanvas("c_sep_power",
                            "",
                            int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                            int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    canvas3.SetGridx()
    canvas3.SetGridy()
    canvas3.SetLogx()

    gr_sp_pik.Draw("ACP")
    gr_sp_pik_dedx.Draw("CPsame")
    gr_sp_pik_combined.Draw("CPsame")
    gr_sp_pik.SetLineColor( ROOT.TColor.GetColor("#0077b6") )
    gr_sp_pik_dedx.SetLineColor( ROOT.kRed+2 )
    gr_sp_pik_combined.SetLineColor(ROOT.kViolet-3)
    gr_sp_pik.SetLineWidth(3)
    gr_sp_pik_dedx.SetLineWidth(3)
    gr_sp_pik_combined.SetLineWidth(3)
    gr_sp_pik.SetMarkerColor( ROOT.TColor.GetColor("#0077b6") )
    gr_sp_pik_dedx.SetMarkerColor( ROOT.kRed+2 )
    gr_sp_pik_combined.SetMarkerColor(ROOT.kViolet-3)
    gr_sp_pik.SetMarkerStyle(20)
    gr_sp_pik_dedx.SetMarkerStyle(20)
    gr_sp_pik_combined.SetMarkerStyle(20)

    gr_sp_pik.GetXaxis().SetRangeUser(0.5, 20)
    gr_sp_pik.GetXaxis().SetNoExponent()
    gr_sp_pik.GetXaxis().SetMoreLogLabels()
    gr_sp_pik.GetYaxis().SetRangeUser(0, 6)

    canvas3.BuildLegend()
    canvas3.Update()
    input("wait")

    gr_sp_kp.Draw("ACP")
    gr_sp_kp_dedx.Draw("CPsame")
    gr_sp_kp_combined.Draw("CPsame")
    gr_sp_kp.SetLineColor( ROOT.TColor.GetColor("#0077b6") )
    gr_sp_kp.SetLineWidth(3)
    gr_sp_kp_dedx.SetLineColor( ROOT.kRed+2 )
    gr_sp_kp_combined.SetLineColor(ROOT.kViolet-3)
    gr_sp_kp_dedx.SetLineWidth(3)
    gr_sp_kp_combined.SetLineWidth(3)
    gr_sp_kp.SetMarkerColor( ROOT.TColor.GetColor("#0077b6") )
    gr_sp_kp_dedx.SetMarkerColor( ROOT.kRed+2 )
    gr_sp_kp_combined.SetMarkerColor(ROOT.kViolet-3)
    gr_sp_kp.SetMarkerStyle(20)
    gr_sp_kp_dedx.SetMarkerStyle(20)
    gr_sp_kp_combined.SetMarkerStyle(20)

    gr_sp_kp.GetXaxis().SetRangeUser(0.5, 20)
    gr_sp_kp.GetXaxis().SetNoExponent()
    gr_sp_kp.GetXaxis().SetMoreLogLabels()
    gr_sp_kp.GetYaxis().SetRangeUser(0, 6)
    canvas3.BuildLegend()
    canvas3.Update()
    input("wait")


get_sep_power_graph(df, tof_column="tofClosest30")
input("EEEEEEEEEND")