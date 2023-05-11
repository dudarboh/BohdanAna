import numpy as np
import ROOT

def get_mass(momentum, tof, track_length):
    speed_of_light = 299.792458 # mm / ns
    beta = track_length/(tof*speed_of_light)
    return momentum*np.sqrt( 1/(beta*beta) - 1) # GeV


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

# gr[0].Draw('AL')
# grshade[0].SetFillStyle(3013)
# grshade[0].SetFillColor(16)
# grshade[0].Draw("af")
# input("wait")

# for i,trk_len in enumerate(vary_track_length):
#     n_points = 0
#     for j, (x, y) in enumerate(zip(momentum, reco_mass[i])):
#         if y == np.nan :
#             continue
#         gr[i].SetPoint(n_points, x, y)

#         n_points += 1

# gr[2].Draw('AL')
# gr[2].SetTitle("True track length (2 m); momentum (GeV); mass (MeV)")
# gr[2].SetLineColor(colors[0])

# gr[1].Draw('Lsame')
# gr[1].SetLineColor(colors[0])
# gr[3].Draw('Lsame')
# gr[3].SetLineColor(colors[0])

# gr[0].Draw('Lsame')
# gr[0].SetLineColor(colors[1])
# gr[4].Draw('Lsame')
# gr[4].SetLineColor(colors[1])


# input("wait")




#  TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

#    c1->SetGrid();
#    c1->DrawFrame(0,0,2.2,12);
   
#    const Int_t n = 20;
#    Double_t x[n], y[n],ymin[n], ymax[n];
#    Int_t i;
#    for (i=0;i<n;i++) {
#      x[i] = 0.1+i*0.1;
#      ymax[i] = 10*sin(x[i]+0.2);
#      ymin[i] = 8*sin(x[i]+0.1);
#      y[i] = 9*sin(x[i]+0.15);
#    }
#    TGraph *grmin = new TGraph(n,x,ymin);
#    TGraph *grmax = new TGraph(n,x,ymax);
#    TGraph *gr    = new TGraph(n,x,y);
#    TGraph *grshade = new TGraph(2*n);
#    for (i=0;i<n;i++) {
#       grshade->SetPoint(i,x[i],ymax[i]);
#       grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
#    }
#    grshade->SetFillStyle(3013);
#    grshade->SetFillColor(16);
#    grshade->Draw("f");
#    grmin->Draw("l");
#    grmax->Draw("l");
#    gr->SetLineWidth(4);
#    gr->SetMarkerColor(4);
#    gr->SetMarkerStyle(21);
#    gr->Draw("CP");







# p = 1.7
# l0 = 2000. # mm
# t0 = 7.12057 # ns
# m0 = 493.677

# t_lim = l0/c
# t = np.arange(t_lim, 10., 0.00001)

# l_lim = t0*c
# l = np.arange(2000., l_lim, 0.1)

# # vs t
# # beta = l0/(t*c)
# # m = p/beta * np.sqrt(1. - beta**2) * 1000.
# # plt.plot(t*1000 - t0*1000, m - m0)
# # plt.xlabel(r'$t - t_{true}$, [ps]')

# # vs l
# beta = l/(t0*c)
# m = p/beta * np.sqrt(1. - beta**2) * 1000.
# plt.plot(l - l0, m - m0)
# plt.xlabel(r'$l - l_{true}$, [mm]')


# plt.ylabel(r"$m - m_{true}$, [MeV]")
# plt.grid(True)
# plt.show()
