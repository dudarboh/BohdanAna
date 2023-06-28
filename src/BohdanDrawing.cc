#include "BohdanDrawing.h"
#include "TOF.h"
#include "EVENT/TrackState.h"
#include "UTIL/Operators.h"
#include "UTIL/TrackTools.h"
#include "marlinutil/DDMarlinCED.h"
#include "marlinutil/CalorimeterHitType.h"

#include "TF1.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <numeric>
#include <cstdlib>
	
using CLHEP::RandGauss;
using dd4hep::rec::Vector3D;

static int picture_counter = 1;

TStyle* getMyStyle(){
    //Mostly CMS style with some ILD style comments
    TStyle* myStyle = new TStyle("myStyle", "My Style");

    // For the canvas:
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetCanvasColor(kWhite);  // ILD 10
    myStyle->SetCanvasDefH(600); //Height of canvas
    myStyle->SetCanvasDefW(600); //Width of canvas
    myStyle->SetCanvasDefX(0);   //Position on screen
    myStyle->SetCanvasDefY(0);

    // For the Pad:
    myStyle->SetPadBorderMode(0);
    myStyle->SetPadGridX(false);
    myStyle->SetPadGridY(false);
    myStyle->SetPadColor(kWhite); // ILD 10
    myStyle->SetGridColor(0);
    myStyle->SetGridStyle(3);
    myStyle->SetGridWidth(1);

    // For the frame:
    myStyle->SetFrameBorderMode(0);
    myStyle->SetFrameFillColor(0); // ILD 10
    myStyle->SetFrameLineWidth(1); // ILD 2
    myStyle->SetFrameBorderSize(1);
    myStyle->SetFrameFillStyle(0);
    myStyle->SetFrameLineColor(1);
    myStyle->SetFrameLineStyle(1);

    // For the histo:
    // myStyle->SetHistFillColor(1);
    // myStyle->SetHistFillStyle(0);
    myStyle->SetHistLineColor(1);
    myStyle->SetHistLineStyle(0);
    myStyle->SetHistLineWidth(1); // ILD 2
    // myStyle->SetLegoInnerR(Float_t rad = 0.5);
    myStyle->SetNumberContours(256); //default 20

    myStyle->SetEndErrorSize(2);
    myStyle->SetErrorX(0.);

    myStyle->SetMarkerStyle(20);

    //For the fit/function:
    myStyle->SetOptFit(1); // ILD 0
    myStyle->SetFitFormat("5.4g");
    myStyle->SetFuncColor(2);
    myStyle->SetFuncStyle(1);
    myStyle->SetFuncWidth(1); // ILD 2

    //For the date:
    myStyle->SetOptDate(0);
    // myStyle->SetDateX(Float_t x = 0.01);
    // myStyle->SetDateY(Float_t y = 0.01);

    // For the statistics box:
    myStyle->SetOptFile(0);
    myStyle->SetOptStat(0); // mean and RMS: SetOptStat("mr");
    myStyle->SetStatColor(kWhite); // ILD 10
    myStyle->SetStatFont(42);
    myStyle->SetStatFontSize(0.025); // ILD 0.07
    myStyle->SetStatTextColor(1);
    myStyle->SetStatFormat("6.4g");
    myStyle->SetStatBorderSize(1); // ILD 0
    myStyle->SetStatH(0.1);
    myStyle->SetStatW(0.15);
    // myStyle->SetStatStyle(Style_t style = 1001);
    // myStyle->SetStatX(Float_t x = 0);
    // myStyle->SetStatY(Float_t y = 0);

    // Margins:
    myStyle->SetPadTopMargin(0.05); // ILD 0.08
    myStyle->SetPadBottomMargin(0.13); // ILD 0.18
    myStyle->SetPadLeftMargin(0.16); // ILD 0.17
    myStyle->SetPadRightMargin(0.02); // ILD 0.08

    // For the Global title:
    myStyle->SetOptTitle(0);    // 0=No Title
    myStyle->SetTitleFont(42);
    myStyle->SetTitleColor(1);
    myStyle->SetTitleTextColor(1);
    myStyle->SetTitleFillColor(10); // ILD 0
    myStyle->SetTitleFontSize(0.05);
    // myStyle->SetTitleH(0); // Set the height of the title box
    // myStyle->SetTitleW(0); // Set the width of the title box
    // myStyle->SetTitleX(0); // Set the position of the title box
    // myStyle->SetTitleY(0.985); // Set the position of the title box
    // myStyle->SetTitleStyle(Style_t style = 1001);
    // myStyle->SetTitleBorderSize(2); // ILD 0

    // For the axis titles:
    myStyle->SetTitleColor(1, "XYZ");
    myStyle->SetTitleFont(42, "XYZ");
    myStyle->SetTitleSize(0.06, "XYZ"); // ILD 0.07
    // myStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // myStyle->SetTitleYSize(Float_t size = 0.02);
    myStyle->SetTitleXOffset(0.9); // ILD 1
    myStyle->SetTitleYOffset(1.25); // ILD 1.1
    // myStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    myStyle->SetLabelColor(1, "XYZ");
    myStyle->SetLabelFont(42, "XYZ");
    myStyle->SetLabelOffset(0.007, "XYZ"); // ILD 0.015
    myStyle->SetLabelSize(0.05, "XYZ"); // ILD 0.06

    // For the axis:
    myStyle->SetAxisColor(1, "XYZ");
    myStyle->SetStripDecimals(kTRUE);
    myStyle->SetTickLength(0.03, "XYZ");
    myStyle->SetNdivisions(510, "XYZ"); // ILD 506
    myStyle->SetPadTickX(0);  // 0=Text labels (and tics) only on bottom, 1=Text labels on top and bottom. ILD 1
    myStyle->SetPadTickY(1);

    // Change for log plots:
    myStyle->SetOptLogx(0);
    myStyle->SetOptLogy(0);
    myStyle->SetOptLogz(0);

    // Postscript options:
    myStyle->SetPaperSize(20.,20.);
    // myStyle->SetLineScalePS(Float_t scale = 3);
    // myStyle->SetLineStyleString(Int_t i, const char* text);
    // myStyle->SetHeaderPS(const char* header);
    // myStyle->SetTitlePS(const char* pstitle);

    // myStyle->SetBarOffset(Float_t baroff = 0.5);
    // myStyle->SetBarWidth(Float_t barwidth = 0.5);
    // myStyle->SetPaintTextFormat(const char* format = "g");
    // myStyle->SetTimeOffset(Double_t toffset);
    // myStyle->SetHistMinimumZero(kTRUE);

    // myStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    myStyle->SetPalette(1,0); // from ILD

    return myStyle;
}

void displayPFO(EVENT::ReconstructedParticle* pfo, IMPL::TrackStateImpl tsStdReco, IMPL::TrackStateImpl tsEasy){
    std::vector<EVENT::Track*> tracks;
    if ( not pfo->getTracks().empty() )
        tracks = getSubTracks(pfo->getTracks()[0]);
    for(auto* track: tracks){
        auto hits = track->getTrackerHits();
        for (auto* hit : hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int size = 4; // larger point
            unsigned long color = 0x000000;
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }
    //plot ecal state for fun
    const TrackState* tsCalo = tracks.back()->getTrackState( TrackState::AtCalorimeter );
    auto posCalo = tsCalo->getReferencePoint();
    ced_hit(posCalo[0], posCalo[1], posCalo[2] + tsCalo->getZ0(), 0, 4, 0x000000); // track state at calorimeterer


    std::vector<EVENT::Cluster*> clusters = pfo->getClusters();
    for(auto* cluster: clusters){
        auto hits = cluster->getCalorimeterHits();
        for (auto* hit: hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int size = 4; // larger point
            unsigned long color = 0x000000;
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }

    auto plotHelixFromTrackState = [](TrackStateImpl ts, unsigned long color){
        double omega = ts.getOmega();
        // double tanL = ts.getTanLambda();
        double phi = ts.getPhi();
        double d0 = ts.getD0();
        double z0 = ts.getZ0();
        auto ref = ts.getReferencePoint();
        std::array<double, 3> mom = UTIL::getTrackMomentum(&ts, 3.5);

        //
        int charge = omega > 0 ? 1 : -1;
        double x = -d0*std::sin(phi) + ref[0];
        double y = d0*std::cos(phi) + ref[1];
        double z = z0 + ref[2];


        int type = 0; // point
        int size = 8; // larger point
        ced_hit(x, y, z, type, size, color);
        size = 2;
        DDMarlinCED::drawHelix( 3.5, charge, x, y, z, mom[0], mom[1], mom[2], type, size, color, 0, 2000, 2600, 0);
        DDMarlinCED::drawHelix( -3.5, -charge, x, y, z, mom[0], mom[1], mom[2], type, size, color, 0, 2000, 2600, 0);
    };

    plotHelixFromTrackState(tsStdReco, 0xfc0505);
    plotHelixFromTrackState(tsEasy, 0x0905fc);
}


void displayFTDSimHits(EVENT::LCEvent* evt){
    EVENT::LCCollection* hits = evt->getCollection("FTDCollection");

    for (int i=0; i < hits->getNumberOfElements(); i++) {
        SimTrackerHit* hit = static_cast<EVENT::SimTrackerHit*>(hits->getElementAt(i));
        auto pos = hit->getPosition();
        int type = 0; // point
        int size = 4; // larger point
        unsigned long color = 0x000000;
        ced_hit(pos[0], pos[1], pos[2], type, size, color);
    }
}

void plotCanvas(EVENT::Cluster* cluster, Vector3D posAtEcal, Vector3D momAtEcal, MCParticle* mc){
    std::cout<<"PDG: "<<mc->getPDG()<<std::endl;
    Vector3D mom(mc->getMomentum());
    std::cout<<"Momentum: "<<mom.r()<<" ( "<<mom.rho()<<", "<<mom.z()<<")"<<std::endl;

    auto hits = cluster->getCalorimeterHits();
    auto frankHits = selectFrankEcalHits(cluster, posAtEcal, momAtEcal, 10);

    //fill maps
    std::map <CalorimeterHit*, double> hit2dToTrack;
    std::map <CalorimeterHit*, double> hit2time;
    std::map <CalorimeterHit*, double> hit2timeSmeared;
    for(auto hit : hits){
        bool isECALHit = ( CHT( hit->getType() ).caloID() == CHT::ecal );
        if ( (!isECALHit) ) continue;
        hit2dToTrack[hit] = (Vector3D(hit->getPosition()) - posAtEcal).r();
        hit2time[hit] = hit->getTime();
        hit2timeSmeared[hit] = RandGauss::shoot(hit->getTime(), 0.1);
    }
    if ( hit2time.empty() ) return;

    TStyle* myStyle = getMyStyle();
    myStyle->cd();
    gROOT->ForceStyle();
    TCanvas c = TCanvas("c", "All ECAL hits");

    std::vector<double> x_all, y_all, y_all_smeared;
    for(auto hit : hits){
        bool isECALHit = ( CHT( hit->getType() ).caloID() == CHT::ecal );
        if ( (!isECALHit) ) continue;
        x_all.push_back( hit2dToTrack[hit] );
        y_all.push_back( hit2time[hit] );
        y_all_smeared.push_back( hit2timeSmeared[hit] );
    }
    std::vector<double> y_err(x_all.size(), 0.1);
    TGraphErrors gr( x_all.size(), x_all.data(), y_all.data(), nullptr, y_err.data() );
    gr.GetXaxis()->SetTitle("d to track (mm)");
    gr.GetYaxis()->SetTitle("Hit time (ns)");
    gr.Draw("AP");
    TLatex t;
    t.SetNDC();
    t.DrawLatex(0.2,  0.87, Form("#splitline{PDG: %d}{p: %.1f GeV}", mc->getPDG(), Vector3D(mc->getMomentum()).r() ));
    TLegend* legend = new TLegend(0.57, 0.84, 0.91, 0.93);
    legend->AddEntry(&gr, "ECAL hits (true time)", "pe");
    legend->Draw();
    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step1.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }

    gr.SetMarkerColor(15);
    gr.SetLineColor(15);
    TGraph gr2( x_all.size(), x_all.data(), y_all_smeared.data());
    gr2.SetTitle("100 ps smeared hit time;d to track (mm); Hit time (ns)");
    gr2.Draw("Psame");

    legend->Clear();
    legend->AddEntry(&gr, "Previous step", "pe");
    legend->AddEntry(&gr2, "Smear time 100 ps", "p");
    legend->Draw();

    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step2.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }

    std::vector<double> x_frank, y_frank_smeared;
    for(auto hit : frankHits){
        bool isECALHit = ( CHT( hit->getType() ).caloID() == CHT::ecal );
        if ( (!isECALHit) ) continue;
        x_frank.push_back( hit2dToTrack[hit] );
        y_frank_smeared.push_back( hit2timeSmeared[hit] );
    }
 
    TGraph gr3( x_frank.size(), x_frank.data(), y_frank_smeared.data() );
    TF1 f("f", Form("%f + (x-%f)/%f", y_all_smeared.front(), x_all.front(), CLHEP::c_light), 0., *std::max_element(x_all.begin(), x_all.end()) );
    gr3.SetTitle("Selected hits;d to track (mm); Hit time (ns)");
    gr2.SetMarkerColor(15);
    gr2.SetLineColor(15);
    gr2.Draw("P");
    gr3.Draw("Psame");

    f.Draw("same");

    legend->Clear();
    legend->AddEntry(&gr2, "Previous step", "p");
    legend->AddEntry(&gr3, "Selected hits", "p");
    legend->AddEntry(&f, "Speed of light", "l");
    legend->Draw();

    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step3.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }



    std::vector<double> y_frank_smeared_corr;
    for(int i=0; i < x_frank.size(); i++){
        y_frank_smeared_corr.push_back( y_frank_smeared[i] - x_frank[i]/CLHEP::c_light );
    }

    TGraph gr4( x_frank.size(), x_frank.data(), y_frank_smeared_corr.data() );
    double average = std::reduce(y_frank_smeared_corr.begin(), y_frank_smeared_corr.end()) / double(y_frank_smeared_corr.size());
    TF1 f2("f", Form("%f", average), 0., *std::max_element(x_all.begin(), x_all.end()) );

    gr4.SetTitle("Corrected hit time;d to track (mm); Corrected hit time (ns)");
    gr.Draw("AP");
    gr.SetLineStyle(0);
    gr.SetLineColor(0);
    gr.SetMarkerStyle(0);
    gr.SetLineColor(0);
    gr4.Draw("Psame");

    f2.Draw("same");

    legend->Clear();
    legend->AddEntry(&gr4, "Time corrected", "p");
    legend->AddEntry(&f2, Form("Average (TOF) = %.1f ns", average), "l");
    legend->Draw();

    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step4.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }
    picture_counter++;

}