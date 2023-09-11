import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
# ROOT.EnableImplicitMT()


code='''

RVec <double> dToImpact(const RVec<XYZVector>& hits, const XYZVector& rImpact){
    auto getDistance = [&](const XYZVector& hit) { return (hit - rImpact).R(); };
    return Map(hits, getDistance);
}


RVec <double> dToLine(const RVec <XYZVector>& hits, const XYZVector& vtx, const XYZVector& p){
    RVec <double> distance;
    for (const auto& hit:hits){
        double d = (hit - vtx).Cross(p.Unit()).R();
        distance.push_back(d);
    }
    return distance;
}

RVec <bool> selectHits(const RVec<double>& dToLine, const RVec<int>& layer_hit, bool only_closest=true, int n_layers=10, double cyl_cut=9999.){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    if(!only_closest){
        for (int i=0; i < nHits; ++i)
            if (dToLine[i] < cyl_cut && layer_hit[i] < n_layers)
                selected[i] = true;
        return selected;
    }

    for (int layer=0; layer<n_layers; ++layer){
        map<int, double> layer_hits;
        for (int i=0; i < nHits; ++i)
            if(dToLine[i] < cyl_cut && layer_hit[i] == layer)
                layer_hits[i] = dToLine[i];

        int min_idx =(*min_element(layer_hits.begin(), layer_hits.end(),
                        [](const auto& l, const auto& r) { return l.second < r.second; }) ).first ;
        selected[min_idx] = true;
    }
    return selected;
}

double tofAvg(const RVec <double>& tofHit, const RVec<double>& dToImpact){
    int nHits = tofHit.size();
    double tofSum = 0.;
    for(int i=0; i < nHits; ++i) tofSum += tofHit[i] - dToImpact[i]/SPEED_OF_LIGHT;
    return tofSum/nHits;
}

'''
ROOT.gInterpreter.Declare(code)


colors = [ROOT.TColor.GetColor('#1b9e77'), ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#e7298a')]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Define("")


#PLOT 2D
def plot_2d(track_length_column="trackLengthToEcal_IKF_zedLambda", tof_column="tofClosest0"):
    print("1")
    df_beta = df.Define("mom", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)").\
                Define("beta", f"{track_length_column}/({tof_column}*299.792458)").Filter("beta >= 0 && beta <= 1")
    df_mass = df_beta.Define("mass", "mom*sqrt( 1./(beta*beta) - 1.)*1000")
    print("2")
    h = df_mass.Histo2D((f"h_{track_length_column}", f"{track_length_column}; momentum [GeV]; Mass [MeV]", 500, 0, 15, 500, -100, 1300.), "mom","mass")
    print("3")

    canvas = ROOT.TCanvas(f"{track_length_column}",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    print("4")

    h.Draw("colz")
    print("5")
    h.SetMinimum(1)
    h.SetMaximum(10000)
    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.SetRightMargin(0.12)
    canvas.Update()
    print("6")
    palette = h.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.89)
    palette.SetX2NDC(0.91)
    canvas.Modified()
    canvas.Update()
    print("7")
    input("wait")
    canvas.Close()
    print("8")


#PLOT 1D
def plot_1d():
    algorithms = ["trackLengthToEcal_IKF_zedLambda",
                "trackLengthToEcal_SHA_zedLambda_ECAL",
                "trackLengthToEcal_SHA_phiLambda_IP"]

    histos = []
    for name in algorithms:
        df_beta = df.Define("mom", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)").\
                    Define("beta", f"{name}/(tofClosest0*299.792458)").Filter("beta >= 0 && beta <= 1")
        df_mass = df_beta.Define("mass", "mom*sqrt( 1./(beta*beta) - 1.)*1000")
        h = df_mass.Histo1D((f"h_{name}", f"{name}; mass [MeV]; N entries", 2000, 0, 1100), "mass")
        histos.append(h)

    ROOT.gStyle.SetPadLeftMargin(0.18)
    hs = ROOT.THStack()
    canvas = ROOT.TCanvas("mass_1d_track_lenghts",
                            "",
                            int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                            int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    legend = ROOT.TLegend()
    for i, (name, h) in enumerate(zip(algorithms, histos)):
        h.Draw("" if i == 0 else "Lsame")
        h.SetLineColor(colors[i])
        h.SetLineWidth(3)
        legend.AddEntry(h.GetPtr(),name,"l")

    legend.Draw()
    canvas.Modified()
    canvas.Update()
    input("wait")


plot_1d()
# plot_2d("trackLengthToEcal_SHA_phiLambda_IP", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_phiZed_IP", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_zedLambda_IP", "tofClosest0")

# plot_2d("trackLengthToEcal_SHA_phiLambda_ECAL", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_phiZed_ECAL", "tofClosest0")
# plot_2d("trackLengthToEcal_SHA_zedLambda_ECAL", "tofClosest0")

# plot_2d("trackLengthToEcal_IKF_phiLambda", "tofClosest0")
# plot_2d("trackLengthToEcal_IKF_phiZed", "tofClosest0")
# plot_2d("trackLengthToEcal_IKF_zedLambda", "tofClosest0")