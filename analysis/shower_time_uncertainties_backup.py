import ROOT
import numpy as np
from sklearn import linear_model
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
# ROOT.EnableImplicitMT()

ROOT.gInterpreter.Declare('#include "tof.hpp"')
df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/tof_studies.root")\
    .Filter("if (rdfentry_ % 1000000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")\
    .Filter("sqrt(recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy + recoCaloPz*recoCaloPz) < 10.")\
    .Range(100000)

# Only n_layers layers!
n_layers = 10
time_resolution = 50
time_res_cut = 3.6*time_resolution
df = df.Redefine("dl", f"dl[layerHit<{n_layers}]")\
        .Redefine("dToImpact", f"dToImpact[layerHit<{n_layers}]")\
        .Redefine("dToLine", f"dToLine[layerHit<{n_layers}]")\
        .Redefine("tHit", f"tHit[layerHit<{n_layers}]")\
        .Redefine("tSurface", f"tSurface[layerHit<{n_layers}]")\
        .Redefine("layerHit", f"layerHit[layerHit<{n_layers}]")

h = df.Histo1D(("h", "title", 200, 0, 5.), "dToLine")
h.Draw()
input("wait")


#smear
df = df.Define("tHit50ps", "smear(rdfslot_, tHit, gaus50)")\
        .Define("tHit100ps", "smear(rdfslot_, tHit, gaus100)")\
        .Define("tSurface50ps", "getTimeAtSurface(tHit50ps, dToImpact)")\
        .Define("tSurface100ps", "getTimeAtSurface(tHit100ps, dToImpact)")

df = df.Define("CylMask", "selectCylinderHits(dToLine, 6.)")\
        .Define("CyltSurface", f"tSurface{time_resolution}ps[CylMask]")\
        .Define("dtCylinder", "( Mean( CyltSurface ) - tofClosest0)*1000")\

df = df.Define("cylMedianMask", f"selectMedianHits(CyltSurface, {time_res_cut})")\
        .Define("dtCylMedian", "( Mean( CyltSurface[cylMedianMask] ) - tofClosest0)*1000")\

df = df.Define("medianMask", f"selectMedianHits(tSurface{time_resolution}ps, {time_res_cut})")\
        .Define("dToLineMedianMask", "dToLine[medianMask]")\
        .Define("tSurfaceMedianMask", f"tSurface{time_resolution}ps[medianMask]")\
        .Define("medianCylMask", "selectCylinderHits(dToLineMedianMask, 6.)")\
        .Define("medianCyltSurface", "tSurfaceMedianMask[medianCylMask]")\
        .Define("dtMedianCylinder", "( Mean( medianCyltSurface ) - tofClosest0)*1000")\

df = df.Define("frankMask", "selectFrankHits(dToLine, layerHit)")\
        .Define("dtFrank", f"( Mean( tSurface{time_resolution}ps[frankMask] ) - tofClosest0)*1000")

df = df.Define("CylMaskKonrad", "selectCylinderHits(dToLine, 7.74)")\
        .Define("CyltSurfaceKonrad", f"tSurface50ps[CylMaskKonrad]")\
        .Define("cylMedianMaskKonrad", f"selectMedianHits(CyltSurfaceKonrad, 292.0)")\
        .Define("dtCylMedianKonrad", "( Mean( CyltSurfaceKonrad[cylMedianMaskKonrad] ) - tofClosest0)*1000")\

df = df.Define("CylMaskKonrad2", "selectCylinderHits(dToLine, 9.36)")\
        .Define("CyltSurfaceKonrad2", f"tSurface50ps[CylMaskKonrad2]")\
        .Define("cylMedianMaskKonrad2", f"selectMedianHits(CyltSurfaceKonrad2, 179)")\
        .Define("dtCylMedianKonrad2", "( Mean( CyltSurfaceKonrad2[cylMedianMaskKonrad2] ) - tofClosest0)*1000")\


data = df.AsNumpy(["dtCylinder", "dtCylMedian", "dtMedianCylinder", "dtFrank", "dtCylMedianKonrad", "dtCylMedianKonrad2"])

mu90, rms90, _, _ = fit90( data["dtCylinder"] )
print("dtCylinder:", mu90," / ", rms90)
mu90, rms90, _, _ = fit90( data["dtCylMedian"] )
print("dtCylMedian:", mu90," / ", rms90)

mu90, rms90, _, _ = fit90( data["dtCylMedianKonrad"] )
print("dtCylMedianKonrad:", mu90," / ", rms90)

mu90, rms90, _, _ = fit90( data["dtCylMedianKonrad2"] )
print("dtCylMedianKonrad2:", mu90," / ", rms90)


mu90, rms90, _, _ = fit90( data["dtMedianCylinder"] )
print("dtMedianCylinder:", mu90," / ", rms90)

mu90, rms90, _, _ = fit90( data["dtFrank"] )
print("dtFrank:", mu90," / ", rms90)





input("wait")






def selectFrankHits(layer, dPerp):
    '''take only closest hit per layer and return index mask'''
    mask = np.full(dPerp.shape, False)
    for i in range(0, 10):
        indicesInLayer = np.where(layer == i)[0]
        if indicesInLayer.size == 0:
            continue
        idx = indicesInLayer[ np.argmin( dPerp[indicesInLayer] ) ]
        mask[idx] = True
    return mask

def selectCylinderHits(dPerp):
    '''take only hits within some radii and return index mask'''
    r = 6.
    mask = dPerp < r
    while not mask.any():
        r += 0.05
        mask = dPerp < r
    return mask

def selectRANSACHits(dTotal, tHit, threshold=None):
    if tHit.size < 4:
        return np.array([True]*tHit.size)
    ransac = linear_model.RANSACRegressor(residual_threshold=threshold)
    ransac.fit(dTotal[:, None], tHit[:, None])
    return ransac.inlier_mask_

def selectCustomRANSACHits(tHit, threshold=200.):
    nHits = tHit.size
    if nHits % 2 == 0:
        # A trick to pick the left element if N elements is even.
        median = np.median( np.insert(tHit, 0, -999999) )
    else:
        median = np.median(tHit)
    mask = abs(tHit - median) < threshold
    return mask


def plot_shower(df, i, check_algorithm=None):
    #extract the data
    shower = df.Filter(f"rdfentry_ == {i}")

    true_tof = shower["tofClosest0"][0]
    layer = np.array([ x for x in shower["layerHit"][0] ])
    dTotal = np.array([ x for x in shower["dToImpact"][0] ])
    dAlong = np.array([ x for x in shower["dl"][0] ])
    dPerp = np.array([ x for x in shower["dToLine"][0] ])
    tHit = np.array([ x for x in shower["tHit"][0] ])
    tSurface = np.array([ x for x in shower["tSurface"][0] ])

    dt = (tSurface - true_tof)*1000

    # default plotting settings
    max_plot_distance = 1.1*max(max(dAlong), max(dPerp))
    max_plot_time = 250

    margin=0.25
    left_margin_fraction=0.8
    bottom_margin_fraction=0.7
    ROOT.gStyle.SetPadLeftMargin( left_margin_fraction * margin )
    ROOT.gStyle.SetPadRightMargin( (1 - left_margin_fraction) * margin )
    ROOT.gStyle.SetPadBottomMargin(bottom_margin_fraction * margin)
    ROOT.gStyle.SetPadTopMargin( (1 - bottom_margin_fraction) * margin)
    canvas = ROOT.TCanvas("c", "title", 1800, 600)
    canvas.Divide(3, 1)

    if check_algorithm == "Frank":
        good_hits_mask = selectFrankHits(layer, dPerp)
        bad_hits_mask = np.invert(good_hits_mask)
    elif check_algorithm == "RANSAC":
        good_hits_mask = selectRANSACHits(dTotal, tSurface, 0.13)
        bad_hits_mask = np.invert(good_hits_mask)
    elif check_algorithm == "Cylinder":
        good_hits_mask = selectCylinderHits(dPerp)
        bad_hits_mask = np.invert(good_hits_mask)


    if check_algorithm:
        n_good_points = len(tSurface[good_hits_mask])
        n_bad_points = len(tSurface[bad_hits_mask])
    # Plot three distributions
    pad1 = canvas.cd(1)
    frame1 = pad1.DrawFrame(-10, -max_plot_time, max_plot_distance, max_plot_time, "Time vs longitudinal distance; d_{along} (mm); TOF_{hit} - TOF_{true} (ps)")
    frame1.GetYaxis().SetTitleOffset(1.4)
    frame1.GetXaxis().SetTitleOffset(1.1)
    if not check_algorithm:
        grAlong = ROOT.TGraph(len(tSurface), dAlong, dt)
        grAlong.Draw("P")
    else:
        if n_good_points > 0:
            grAlongGood = ROOT.TGraph(n_good_points, dAlong[good_hits_mask], dt[good_hits_mask])
            grAlongGood.SetMarkerColor(ROOT.kGreen+2)
            grAlongGood.Draw("P")
        if n_bad_points > 0:
            grAlongBad = ROOT.TGraph(n_bad_points, dAlong[bad_hits_mask], dt[bad_hits_mask])
            grAlongBad.SetMarkerColor(ROOT.kRed+2)
            grAlongBad.Draw("Psame")



    pad2 = canvas.cd(2)
    frame2 = pad2.DrawFrame(-10, -max_plot_time, max_plot_distance, max_plot_time, "Time vs transverse distance; d_{perp} (mm); TOF_{hit} - TOF_{true} (ps)")
    frame2.GetYaxis().SetTitleOffset(1.4)
    frame2.GetXaxis().SetTitleOffset(1.1)
    if not check_algorithm:
        grPerp = ROOT.TGraph(len(tSurface), dPerp, dt)
        grPerp.Draw("P")
    else:
        if n_good_points > 0:
            grPerpGood = ROOT.TGraph(len(tSurface[good_hits_mask]), dPerp[good_hits_mask], dt[good_hits_mask])
            grPerpGood.SetMarkerColor(ROOT.kGreen+2)
            grPerpGood.Draw("P")
        if n_bad_points > 0:
            grPerpBad = ROOT.TGraph(len(tSurface[bad_hits_mask]), dPerp[bad_hits_mask], dt[bad_hits_mask])
            grPerpBad.SetMarkerColor(ROOT.kRed+2)
            grPerpBad.Draw("Psame")

    pad3 = canvas.cd(3)
    frame3 = pad3.DrawFrame(-10, -10, max_plot_distance, max_plot_distance, "Longit. vs transv. distance; d_{perp} (mm); d_{along} (mm)")
    frame3.GetYaxis().SetTitleOffset(1.4)
    frame3.GetXaxis().SetTitleOffset(1.1)
    if not check_algorithm:
        grXY = ROOT.TGraph(len(tSurface), dPerp, dAlong)
        grXY.Draw("P")
    else:
        if n_good_points > 0:
            grXYGood = ROOT.TGraph(len(tSurface[good_hits_mask]), dPerp[good_hits_mask], dAlong[good_hits_mask])
            grXYGood.SetMarkerColor(ROOT.kGreen+2)
            grXYGood.Draw("P")
        if n_bad_points > 0:
            grXYBad = ROOT.TGraph(len(tSurface[bad_hits_mask]), dPerp[bad_hits_mask], dAlong[bad_hits_mask])
            grXYBad.SetMarkerColor(ROOT.kRed+2)
            grXYBad.Draw("Psame")

    canvas.Update()
    canvas.Print(f"./shower_examples/shower_example_{i}_{check_algorithm}.png")
    input("wait")

def doRANSACScan(showers):

    canvas = create_canvas()
    gr = ROOT.TGraph()
    gr.SetTitle("RANSAC optimal cut scan; cut value (ps); RMS_{90} of the residual")

    cuts = [100 + i*10 for i in range(20)] # in ps
    residuals = {cut: [] for cut in cuts}

    for i, (dTotal, dPerp, t, true_tof) in enumerate(zip(showers["dToImpact"], showers["dToLine"], showers["tSurface50ps"], showers["tofClosest0"])):
        # preselect cylinder
        dPerp = np.array( [x for x in dPerp] )
        mask_cyl = selectCylinderHits(dPerp)
        dTotal = np.array( [x for x in dTotal] )[mask_cyl]
        t = np.array( [x for x in t] )[mask_cyl]*1e3

        for cut in cuts:
            # good_hits_mask = selectRANSACHits(dTotal, t, cut)
            good_hits_mask = selectCustomRANSACHits(t, cut)
            residual = np.mean( t[good_hits_mask] ) - true_tof*1e3
            residuals[cut].append(residual)
        if (i+1) % 20 == 0:
            print(i)
            # print(rms90)
            for j, (cut, res) in enumerate(residuals.items()):
                _, rms90, _, _ = fit90( np.array(res) )
                gr.SetPoint(j, cut, rms90)
            gr.Draw("APL")
            gr.GetXaxis().SetNdivisions(506)
            gr.SetMinimum(15.)
            gr.SetMaximum(20.)
            canvas.Update()
        if i % 1000 == 0:
            canvas.Print(f"RANSAC{i}_scan_50ps_precylinder.png")

def compare_methods(showers):

    RANSAC_CUT = 200.
    TIME_RESOLUTION = 50 # only 50 and 100
    canvas = create_canvas(0.4, 0.7, 0.85)
    methods = ["Cylinder only", "cyl + RANSAC", "cyl + Median cut"]
    h = ROOT.TH1F("h", "RANSAC optimal cut scan;; RMS_{90} of the residual", len(methods), 0, len(methods))
    if TIME_RESOLUTION == 50:
        # h.SetMinimum(14.)
        # h.SetMaximum(30.)
        h.SetMinimum(15.)
        h.SetMaximum(20.)

    elif TIME_RESOLUTION == 100:
        h.SetMinimum(25.)
        h.SetMaximum(40.)

    for i, m in enumerate(methods):
        h.GetXaxis().SetBinLabel(i+1, m)

    residuals = { m:[] for m in methods}

    for i, (dTotal, dPerp, layer, t, true_tof) in enumerate( zip( showers["dToImpact"], showers["dToLine"], showers["layerHit"], showers[f"tSurface{TIME_RESOLUTION}ps"], showers["tofClosest0"] ) ):
        dTotal = np.array( [x for x in dTotal] )
        dPerp = np.array( [x for x in dPerp] )
        layer = np.array( [x for x in layer] )
        t = np.array( [x for x in t] )*1e3

        # mask_ransac = selectRANSACHits(dTotal, t, RANSAC_CUT)
        # residual = np.mean( t[mask_ransac] ) - true_tof*1e3
        # residuals["RANSAC only"].append(residual)

        # mask_default = selectFrankHits(layer, dPerp)
        # residual = np.mean( t[mask_default] ) - true_tof*1e3
        # residuals["Default only"].append(residual)

        # t_ransac = t[mask_ransac]
        # mask_ransac_default = selectFrankHits(layer[mask_ransac], dPerp[mask_ransac])
        # residual = np.mean( t_ransac[mask_ransac_default] ) - true_tof*1e3
        # residuals["RANSAC + default"].append(residual)

        # t_default = t[mask_default]
        # mask_default_ransac = selectRANSACHits(dTotal[mask_default], t_default, RANSAC_CUT)
        # residual = np.mean( t_default[mask_default_ransac]) - true_tof*1e3
        # residuals["default + RANSAC"].append(residual)

        # t_ransac = t[mask_ransac]
        # mask_ransac_cyl = selectCylinderHits(dPerp[mask_ransac])
        # residual = np.mean( t_ransac[mask_ransac_cyl] ) - true_tof*1e3
        # residuals["RANSAC + cyl"].append(residual)

        mask_cyl = selectCylinderHits(dPerp)
        t_cyl = t[mask_cyl]
        residual = np.mean( t_cyl ) - true_tof*1e3
        residuals["Cylinder only"].append(residual)

        mask_cyl_ransac1 = selectRANSACHits(dTotal[mask_cyl], t_cyl, RANSAC_CUT)
        residual = np.mean( t_cyl[mask_cyl_ransac1]) - true_tof*1e3
        residuals["cyl + RANSAC"].append(residual)

        t_cyl = t[mask_cyl]
        mask_cyl_ransac2 = selectCustomRANSACHits(t_cyl, RANSAC_CUT)
        residual = np.mean( t_cyl[mask_cyl_ransac2]) - true_tof*1e3
        residuals["cyl + Median cut"].append(residual)


        if i > 20:
            print(i)
            for j, (method, res) in enumerate(residuals.items()):
                _, rms90, _, _ = fit90( np.array(res) )
                h.SetBinContent(j+1, rms90)
                if i % 500 == 0:
                    print(method, rms90)
            h.Draw()
            h.GetXaxis().LabelsOption("v")
            canvas.Update()
        if i % 1000 == 0:
            
            canvas.Print(f"compare_methods2_{i}_50ps_fixed.png")


def find_optimal_cuts():
    '''Scans perpendicular distance and dt to median cuts'''
    canvas = create_canvas(0.4, 0.7, 0.85)

    file = ROOT.TFile("/nfs/dust/ilc/user/dudarboh/tof/tof_studies.root")
    tree = file.treename

    for event in tree:
        layer = np.array([ x for x in shower["layerHit"][0] ])
        dTotal = np.array([ x for x in shower["dToImpact"][0] ])
        dAlong = np.array([ x for x in shower["dl"][0] ])
        dPerp = np.array([ x for x in shower["dToLine"][0] ])
        tHit = np.array([ x for x in shower["tHit"][0] ])
        tSurface = np.array([ x for x in shower["tSurface"][0] ])


    methods = ["Cylinder only", "cyl + RANSAC", "cyl + Median cut"]
    h = ROOT.TH1F("h", "RANSAC optimal cut scan;; RMS_{90} of the residual", len(methods), 0, len(methods))
    if TIME_RESOLUTION == 50:
        # h.SetMinimum(14.)
        # h.SetMaximum(30.)
        h.SetMinimum(15.)
        h.SetMaximum(20.)

    elif TIME_RESOLUTION == 100:
        h.SetMinimum(25.)
        h.SetMaximum(40.)

    for i, m in enumerate(methods):
        h.GetXaxis().SetBinLabel(i+1, m)

    residuals = { m:[] for m in methods}

    for i, (dTotal, dPerp, layer, t, true_tof) in enumerate( zip( showers["dToImpact"], showers["dToLine"], showers["layerHit"], showers[f"tSurface{TIME_RESOLUTION}ps"], showers["tofClosest0"] ) ):
        dTotal = np.array( [x for x in dTotal] )
        dPerp = np.array( [x for x in dPerp] )
        layer = np.array( [x for x in layer] )
        t = np.array( [x for x in t] )*1e3

        # mask_ransac = selectRANSACHits(dTotal, t, RANSAC_CUT)
        # residual = np.mean( t[mask_ransac] ) - true_tof*1e3
        # residuals["RANSAC only"].append(residual)

        # mask_default = selectFrankHits(layer, dPerp)
        # residual = np.mean( t[mask_default] ) - true_tof*1e3
        # residuals["Default only"].append(residual)

        # t_ransac = t[mask_ransac]
        # mask_ransac_default = selectFrankHits(layer[mask_ransac], dPerp[mask_ransac])
        # residual = np.mean( t_ransac[mask_ransac_default] ) - true_tof*1e3
        # residuals["RANSAC + default"].append(residual)

        # t_default = t[mask_default]
        # mask_default_ransac = selectRANSACHits(dTotal[mask_default], t_default, RANSAC_CUT)
        # residual = np.mean( t_default[mask_default_ransac]) - true_tof*1e3
        # residuals["default + RANSAC"].append(residual)

        # t_ransac = t[mask_ransac]
        # mask_ransac_cyl = selectCylinderHits(dPerp[mask_ransac])
        # residual = np.mean( t_ransac[mask_ransac_cyl] ) - true_tof*1e3
        # residuals["RANSAC + cyl"].append(residual)

        mask_cyl = selectCylinderHits(dPerp)
        t_cyl = t[mask_cyl]
        residual = np.mean( t_cyl ) - true_tof*1e3
        residuals["Cylinder only"].append(residual)

        mask_cyl_ransac1 = selectRANSACHits(dTotal[mask_cyl], t_cyl, RANSAC_CUT)
        residual = np.mean( t_cyl[mask_cyl_ransac1]) - true_tof*1e3
        residuals["cyl + RANSAC"].append(residual)

        t_cyl = t[mask_cyl]
        mask_cyl_ransac2 = selectCustomRANSACHits(t_cyl, RANSAC_CUT)
        residual = np.mean( t_cyl[mask_cyl_ransac2]) - true_tof*1e3
        residuals["cyl + Median cut"].append(residual)


        if i > 20:
            print(i)
            for j, (method, res) in enumerate(residuals.items()):
                _, rms90, _, _ = fit90( np.array(res) )
                h.SetBinContent(j+1, rms90)
                if i % 500 == 0:
                    print(method, rms90)
            h.Draw()
            h.GetXaxis().LabelsOption("v")
            canvas.Update()
        if i % 1000 == 0:
            
            canvas.Print(f"compare_methods2_{i}_50ps_fixed.png")


# doRANSACScan(df)
# compare_methods(df)
# for i in range(58, 78):
    # plot_shower(df, i)
    # plot_shower(df, i, "Cylinder")
