import numpy as np
import ROOT

def get_mass(momentum, tof, track_length):
    speed_of_light = 299.792458 # mm / ns
    beta = track_length/(tof*speed_of_light)
    return momentum*np.sqrt( 1/(beta*beta) - 1) # GeV


def old_python_plot():
    p = 1.7
    l0 = 2000. # mm
    t0 = 7.12057 # ns
    m0 = 493.677

    t_lim = l0/c
    t = np.arange(t_lim, 10., 0.00001)

    l_lim = t0*c
    l = np.arange(2000., l_lim, 0.1)

    # vs t
    # beta = l0/(t*c)
    # m = p/beta * np.sqrt(1. - beta**2) * 1000.
    # plt.plot(t*1000 - t0*1000, m - m0)
    # plt.xlabel(r'$t - t_{true}$, [ps]')

    # vs l
    beta = l/(t0*c)
    m = p/beta * np.sqrt(1. - beta**2) * 1000.
    plt.plot(l - l0, m - m0)
    plt.xlabel(r'$l - l_{true}$, [mm]')


    plt.ylabel(r"$m - m_{true}$, [MeV]")
    plt.grid(True)
    plt.show()


def dm_vs_dl():
    c = ROOT.TCanvas()

    colors = [ROOT.TColor.GetColor('#fc9272'), ROOT.TColor.GetColor('#de2d26')]

    # Assume pion/kaon particle with momentum (1 -- 10 GeV) and track_length 2 m.
    speed_of_light = 299.792458 # mm / ns
    mass = 0.13957039 # GeV
    # mass = 0.493677 # GeV
    momentum = np.arange(0.1, 8, 0.01) # GeV
    track_length = 2000 # mm

    # Deduce TOF
    tof = track_length/speed_of_light * np.sqrt( 1 + mass*mass/(momentum*momentum) ) # ns


    # Reconstruct mass back, using +- 5,10% of the track_length
    vary_track_length = np.array([mod*track_length for mod in [0.99, 0.995, 1., 1.005, 1.01]]) # mm

    reco_mass = {}
    gr = {'-10' : ROOT.TGraph(), '-5' : ROOT.TGraph(), '0' : ROOT.TGraph(), '+5' : ROOT.TGraph(), '+10' : ROOT.TGraph()}

    m = {}
    m['-10'] = get_mass(momentum, tof, vary_track_length[0])
    m['-5'] = get_mass(momentum, tof, vary_track_length[1])
    m['0'] = get_mass(momentum, tof, vary_track_length[2])
    m['+5'] = get_mass(momentum, tof, vary_track_length[3])
    m['+10'] = get_mass(momentum, tof, vary_track_length[4])

    n_points = len(momentum)
    for key, mass in m.items():
        counter = 0
        for i, (x, y) in enumerate(zip(momentum, mass)):
            if np.isnan(y):
                continue
            gr[key].SetPoint(counter, x, y*1000)
            counter += 1

    gr['0'].Draw('AL')
    gr['0'].SetLineWidth(5)

    gr['-5'].Draw('L same')
    gr['-5'].SetLineWidth(5)
    gr['-5'].SetLineColor(colors[0])
    gr['+5'].Draw('L same')
    gr['+5'].SetLineWidth(5)
    gr['+5'].SetLineColor(colors[0])

    gr['-10'].Draw('L same')
    gr['-10'].SetLineWidth(5)
    gr['-10'].SetLineColor(colors[1])
    gr['+10'].Draw('L same')
    gr['+10'].SetLineWidth(5)
    gr['+10'].SetLineColor(colors[1])

    gr['0'].GetYaxis().SetRangeUser(0, 600)
    gr['0'].SetTitle("True track length (2 m); momentum (GeV); mass (MeV)")
    gr['+5'].SetTitle("#pm 0.5 % (#pm 10 mm)")
    gr['+10'].SetTitle("#pm 1 % (#pm 20 mm)")

    c.BuildLegend()
    c.Update()
    input("wait")


def dm_vs_dt():
    c = ROOT.TCanvas()

    colors = [ROOT.TColor.GetColor('#fc9272'), ROOT.TColor.GetColor('#de2d26')]

    # Assume pion/kaon particle with momentum (1 -- 10 GeV) and track_length 2 m.
    speed_of_light = 299.792458 # mm / ns
    mass = 0.13957039 # GeV
    # mass = 0.493677 # GeV
    momentum = np.arange(0.1, 8, 0.01) # GeV

    tof = 6 # ns

    # Deduce TOF
    track_length = momentum*speed_of_light*tof * np.sqrt( 1 /(mass*mass + momentum*momentum) ) # mm

    # Reconstruct mass back, using +- 5,10% of the track_length
    vary_tof = np.array([mod*tof for mod in [0.99, 0.995, 1., 1.005, 1.01]]) # ns

    reco_mass = {}
    gr = {'-10' : ROOT.TGraph(), '-5' : ROOT.TGraph(), '0' : ROOT.TGraph(), '+5' : ROOT.TGraph(), '+10' : ROOT.TGraph()}

    m = {}
    m['-10'] = get_mass(momentum, vary_tof[0], track_length)
    m['-5'] = get_mass(momentum, vary_tof[1], track_length)
    m['0'] = get_mass(momentum, vary_tof[2], track_length)
    m['+5'] = get_mass(momentum, vary_tof[3], track_length)
    m['+10'] = get_mass(momentum, vary_tof[4], track_length)

    n_points = len(momentum)
    for key, mass in m.items():
        counter = 0
        for i, (x, y) in enumerate(zip(momentum, mass)):
            if np.isnan(y):
                continue
            gr[key].SetPoint(counter, x, y*1000)
            counter += 1

    gr['0'].Draw('AL')
    gr['0'].SetLineWidth(5)

    gr['-5'].Draw('L same')
    gr['-5'].SetLineWidth(5)
    gr['-5'].SetLineColor(colors[0])
    gr['+5'].Draw('L same')
    gr['+5'].SetLineWidth(5)
    gr['+5'].SetLineColor(colors[0])

    gr['-10'].Draw('L same')
    gr['-10'].SetLineWidth(5)
    gr['-10'].SetLineColor(colors[1])
    gr['+10'].Draw('L same')
    gr['+10'].SetLineWidth(5)
    gr['+10'].SetLineColor(colors[1])

    gr['0'].GetYaxis().SetRangeUser(0, 600)
    gr['0'].SetTitle(f"True TOF ({tof} ns); momentum (GeV); mass (MeV)")
    text = f"#pm 0.5 % (#pm {5*tof} ps)"
    gr['+5'].SetTitle(text)
    text = f"#pm 1 % (#pm {10*tof} ps)"
    gr['+10'].SetTitle(text)

    c.BuildLegend()
    c.Update()
    input("wait")