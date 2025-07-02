#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCrystalBall.h"
#include "TCanvas.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TChain.h"
#include "ROOT/RDataFrame.hxx"
#include "TH1D.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "RooPlot.h"

// THE AIM OF THIS CODE IS TO FIT SIMULATED DATA IN ORDER TO SAVE THE 
// PARAMETERS OF THE DISTRIBUTIONS AND SUBSEQUENTLY USE THEM IN THE DATA FITTING

// HERE BOTH THE SIGNAL AND THE BACKGROUND ARE FITTED
// BREM 1 CATEGORY

using namespace RooFit;

// declaration of functions
ROOT::RDF::RNode ApplyFiltersSim(ROOT::RDF::RNode df);
void PlotFitPlots (RooRealVar D0_M, RooDataSet data, RooAbsPdf& model, TCanvas* canvas, int iBrem, std::string name = "fit" );


void fit(){

    ////////////////////////
    ////////// SIGNAL //////
    ////////////////////////

    // Opening the file and defining the RDataFrame
    std::cout << "Defining RDataFrame ..." << '\n';
    std::string pathToFile =  "../../Simulation/D02Kpiee/Kpiee_tree.root";
    ROOT::RDataFrame df_s("DecayTree", pathToFile);
    
    // Defining new columns to have variables in MeV
    // {Recall: Rapidsim output is in GeV}
    auto df_sig = df_s    .Define("D0_M_MeV",     [](double x) { return x * 1000; }, {"D0_M"})
                            // Jpsi is a strange name to call the "dilepton couple"
                          .Define("Jpsi_M_MeV",   [](double x) { return x * 1000; }, {"Jpsi_M"});

    // Filtering on dilepton invariant mass to observe just 1 region                     
    auto df_S = ApplyFiltersSim(df_sig);

    TCanvas* canvas = new TCanvas();

    // PLOTTING HISTOGRAMS OF DATA WITH FIT
    RooRealVar D0_M("D0_M_MeV", "D0 Mass", 1600, 2150, "MeV");

    df_S.Snapshot("DecayTree", "temp_unbinned.root", {"D0_M_MeV"});
    TFile* fileIn = TFile::Open("temp_unbinned.root");
    TTree* tree = (TTree*) fileIn->Get("DecayTree");
    RooDataSet data("data", "data", tree, RooArgSet(D0_M));

    std::cout << "Defining Fit model..." << '\n';

    // Defining the fit model: Johnson function
    RooRealVar mean =  RooRealVar("mean", "mean", 1846.5, 1800, 1900);
    RooRealVar sigma =  RooRealVar("sigma", "sigma", 48.7, 0, 100);
    RooRealVar gamma =  RooRealVar("gamma", "gamma", 0.19336, -100, 100);
    RooRealVar delta =  RooRealVar("delta", "delta", 1.4173, -1000, 1000);
    RooJohnson signal =  RooJohnson("signal","signal", D0_M, mean, sigma, gamma, delta);

    // Plot fit and data
    PlotFitPlots(D0_M, data, signal, canvas, 1, "fitSignal");

    // Saving parameters of the fit
    RooArgSet* params = signal.getParameters(D0_M);
    params->writeToFile("paramsSignal.txt");


    ////////////////////////
    ////// BACKGROUND //////
    ////////////////////////

    // Opening the file and defining the RDataFrame
    std::cout << "Defining RDataFrame ..." << '\n';
    RooRealVar D0_M_bkg("D0_M_MeV", "D0 Mass", 1400, 2000, "MeV");

    std::string pathToFile_bkg =  "../../Simulation/D02Kpipienu/Kpipienu_tree.root";
    ROOT::RDataFrame df_b ("DecayTree", pathToFile_bkg);
    
    // Defining new columns to have variables in MeV
    auto df_bkg = df_b  .Define("D0_M_MeV",     [](double x) { return x * 1000; }, {"D0_M_lm2em"})
                        .Define("Jpsi_M_MeV",   [](double x) { return x * 1000; }, {"Jpsi_M_lm2em"});

    // Filtering on dilepton invariant mass to observe just 1 bin
    auto df_B = ApplyFiltersSim(df_bkg);

    df_B.Snapshot("DecayTree", "temp_unbinned.root", {"D0_M_MeV"});
    TFile* fileIn_bkg = TFile::Open("temp_unbinned.root");
    TTree* tree_bkg = (TTree*) fileIn_bkg ->Get("DecayTree");
    RooDataSet data_bkg("data_bkg", "data", tree_bkg, RooArgSet(D0_M_bkg));

    std::cout << "Defining Fit model for bkg..." << '\n';

    // Defining the fit model: Bukin function
    RooRealVar mean_bkg("mean_bkg", "Mean Mass", 1700, 1500, 1800);
    RooRealVar width ("width", "Width",  50, 1, 200);
    RooRealVar asym("asym", "asymm", -0.2, -1, 1);
    RooRealVar rhoL("rhoL", "rhoL", -0.1, -10, 10);
    RooRealVar rhoR("rhoR", "rhoR", 0.1, -10, 10);

    RooBukinPdf bkg("bukin", "Bukin PDF", D0_M_bkg, mean_bkg, width, asym, rhoL, rhoR);    
    
    // Plot fit and data
    PlotFitPlots(D0_M_bkg, data_bkg, bkg, canvas, 1, "fitBkg");

    // Saving parameters of the fit
    RooArgSet* params_bkg = bkg.getParameters(D0_M_bkg);
    params_bkg->writeToFile("paramsBkg.txt");

    std::cout << "Finish!" << '\n';

}

// ######################################################### //
// #################### Functions body ##################### //
// ######################################################### //

//APPLY FILTERS to SIM: don't need to filter on brem1/0 cathegory bc it's already done in the simulated data with the setting
ROOT::RDF::RNode ApplyFiltersSim(ROOT::RDF::RNode df){

    std::cout << "Hello, I'm filtering bin 3 events on RapidSim data..." << '\n';
    return df.Filter("Jpsi_M_MeV > 675 && Jpsi_M_MeV < 875", "Jpsi_M_MeV");
}


// FUNCTION TO FIT AND PLOT SIM DATA 
void PlotFitPlots (RooRealVar D0_M, RooDataSet data, RooAbsPdf& model, TCanvas* canvas, int iBrem, std::string name = "fitSignal") {

    std::cout << "Plotting fit plots ..." << '\n';

    // dividing the canvas in two pads: one for the fit and one for the pull distribution
    TPad upperPad("upperPad", "upperPad", .02, .2525, .998, .995);
    TPad lowerPad("lowerPad", "lowerPad", .02, .005, .998, .250);

    upperPad.SetBottomMargin(0.1);
    upperPad.SetBorderMode(0);
    lowerPad.SetTopMargin(0.001);
    lowerPad.SetBottomMargin(0.45);
    lowerPad.SetBorderMode(0);
    lowerPad.SetFillStyle(4000);

    lowerPad.Draw();
    upperPad.Draw();

    upperPad.cd();
    // performing the fit
    RooNLLVar * nll = new RooNLLVar("nll","nll", model ,data);
    RooMinimizer m1(*nll) ;
    m1.setVerbose(kFALSE);
    m1.setPrintLevel(3);
    m1.setStrategy(1);
    m1.migrad();
    m1.hesse();

    // Saving the fit result
    RooFitResult* result = m1.save();
    result->Print("v");

    unsigned int i_minimization = 0;

    // plot
    RooPlot* frame = D0_M.frame();
    data.plotOn(frame);
    model.plotOn(frame);

    frame->GetXaxis()->SetTitleSize(0);

    frame->Draw();

    // legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);

    int entries1 = data.sumEntries();
    leg->AddEntry(frame->getObject(0), "RapidSim Data", "lep");
    leg->AddEntry((TObject*)0, Form("\tEntries = %d", entries1), "");
    leg->AddEntry(frame->getObject(1), "Fit model", "l");
    leg->Draw();


    // lower pad for the pull distribution
    lowerPad.cd();
    RooHist* hpull0 = frame->pullHist();
    RooPlot* frameP = D0_M.frame();
    frameP -> SetTitle("Pull Distribution");
    frameP->addPlotable(hpull0, "P");
    hpull0-> SetMarkerSize(0.7);
    frameP->GetXaxis()->SetLabelSize(0.15);
    frameP->GetYaxis()->SetLabelSize(0.12);
    frameP->GetXaxis()->SetTitleSize(0.2);
    frameP->GetXaxis()->SetTitle("m(D^{0}) [MeV]");
    frameP->GetYaxis()->SetTitleSize(0.2);
    frameP->Draw();

    canvas->SaveAs((name + ".png").c_str());
    canvas->Clear();

}


