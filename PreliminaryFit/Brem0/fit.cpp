#include "ROOT/RDataFrame.hxx"
#include "RooCrystalBall.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TTree.h"

// THE AIM OF THIS CODE IS TO FIT SIMULATED DATA IN ORDER TO SAVE THE
// PARAMETERS OF THE DISTRIBUTIONS AND SUBSEQUENTLY USE THEM IN THE DATA FITTING

// HERE BOTH THE SIGNAL AND THE BACKGROUND ARE FITTED
// BREM 0 CATEGORY

using namespace RooFit;

// declaration of functions
ROOT::RDF::RNode ApplyFiltersSim(ROOT::RDF::RNode df);
void PlotFitPlots(RooRealVar D0_M, RooDataSet data, RooAbsPdf &model,
                  TCanvas *canvas, int iBrem, std::string name = "fit");
void PlotFitPlots(RooRealVar D0_M, RooDataHist data, RooAbsPdf &model,
                  TCanvas *canvas, int iBrem, RooFitResult *result = nullptr);

void fit() {

  ////////////////////////
  ////////// SIGNAL //////
  ////////////////////////

  // Opening the file and defining the RDataFrame
  std::cout << "Defining RDataFrame ..." << '\n';
  std::string pathToFile = "../../Simulation/D02Kpiee/Kpiee_tree.root";
  ROOT::RDataFrame df_s("DecayTree", pathToFile);

  // Defining new columns to have variables in MeV
  // {Recall: Rapidsim output is in GeV}
  // Jpsi is a strange name to call the "dilepton couple"
  auto df_sig =
      df_s.Define("D0_M_MeV", [](double x) { return x * 1000; }, {"D0_M"})
          .Define("Jpsi_M_MeV", [](double x) { return x * 1000; }, {"Jpsi_M"});

  // Filtering on dilepton invariant mass to observe just 1 region
  auto df_S = ApplyFiltersSim(df_sig);

  TCanvas *canvas = new TCanvas();

  // PLOTTING HISTOGRAMS OF DATA WITH FIT
  RooRealVar D0_M("D0_M_MeV", "D0 Mass", 1600, 2150, "MeV");

  df_S.Snapshot("DecayTree", "temp_unbinned.root", {"D0_M_MeV"});
  TFile *fileIn = TFile::Open("temp_unbinned.root");
  TTree *tree = (TTree *)fileIn->Get("DecayTree");
  RooDataSet data("data", "data", tree, RooArgSet(D0_M));

  std::cout << "Defining Fit model..." << '\n';

  // Defining the fit model: Crystal Ball function
  RooRealVar mean("mean", "Mean Mass", 1865, 1800, 1900);
  RooRealVar sigma("sigma", "Width", 40, 10, 100);
  RooRealVar sigma1("sigma1", "Width1", 20, 1, 100);
  RooRealVar alpha("alpha", "Alpha", 0.7, 0.001, 10);
  RooRealVar n("n", "n", 3, 0, 200);
  RooRealVar alpha1("alpha1", "Alpha1", 0.7, 0.001, 10);
  RooRealVar n1("n1", "n1", 0, 0, 200);
  RooCrystalBall signal("signal", "signal", D0_M, mean, sigma, sigma1,
                        alpha, n, alpha1, n1);

  n1.setConstant(true);
  n.setConstant(true);

  // Plot fit and data
  PlotFitPlots(D0_M, data, signal, canvas, 0, "fitSignal");

  // Saving parameters of the fit
  RooArgSet *params = signal.getParameters(D0_M);
  params->writeToFile("paramsSignal.txt");

  ////////////////////////
  ////// BACKGROUND //////
  ////////////////////////

  // Opening the file and defining the RDataFrame
  std::cout << "Defining RDataFrame ..." << '\n';
  RooRealVar D0_M_bkg("D0_M_MeV", "D0 Mass", 1700, 1900, "MeV");

  std::string pathToFile_bkg = "../../Simulation/D02Kpipipi/Kpipipi_tree.root";
  ROOT::RDataFrame df_b("DecayTree", pathToFile_bkg);

  // Defining new columns to have variables in MeV
  // The strange names of the variables are due to the rapidsim output.
  // This generation was for the MISID BKG: i.e. D0 \to K pi pi pi
  // It was performed considering the mis-identification of 2 pions as
  // electrons lp2ep means that!

  auto df_bkg = df_b.Define("D0_M_MeV", [](double x) { return x * 1000; },
                            {"D0_M_lp2ep_lm2em"})
                    .Define("Jpsi_M_MeV", [](double x) { return x * 1000; },
                            {"Jpsi_M_lp2ep_lm2em"})
                    .Define("D0_P_MeV", [](double x) { return x * 1000; },
                            {"D0_P_lp2ep_lm2em"})
                    .Define("D0_eta__", [](double x) { return x; },
                            {"D0_eta_lp2ep_lm2em"})
                    .Define("lp_P_MeV", [](double x) { return x * 1000; },
                            {"lp_P_lp2ep_lm2em"})
                    .Define("lp_eta__", [](double x) { return x; },
                            {"lp_eta_lp2ep_lm2em"})
                    .Define("lm_P_MeV", [](double x) { return x * 1000; },
                            {"lm_P_lp2ep_lm2em"})
                    .Define("lm_eta__", [](double x) { return x; },
                            {"lm_eta_lp2ep_lm2em"});

  // Filtering on dilepton invariant mass to observe just 1 bin
  auto df_B = ApplyFiltersSim(df_bkg);

  // In this case it is necessary to apply the efficiency weights
  // since the PID distorts the shape so to rely on the output of the
  // simulation, we apply efficiency of reconstruction

  // Opening the file with the efficiencies
  TFile *fileEff = TFile::Open("efficiencies_correct.root", "READ");
  if (!fileEff || fileEff->IsZombie()) {
    std::cerr << "Errore nell'aprire il file ROOT!" << std::endl;
    return;
  }
  // Getting the TH2D histogram for efficiency
  TH2D *histEff = (TH2D *)fileEff->Get("eff_DLLe>3");
  if (!histEff) {
    std::cerr << "Istogramma TH2D non trovato!" << std::endl;
    fileEff->Close();
    return;
  }

  // Normalizing the histogram to have a maximum of 1
  double maxEff = histEff->GetMaximum();
  histEff->Scale(1.0 / maxEff);

  // Defining a lambda function to apply the efficiency weights according to the
  // momentum and eta of the leptons
  auto addWeights = [histEff](double P, double eta) {
    int binX = histEff->GetXaxis()->FindBin(P);
    int binY = histEff->GetYaxis()->FindBin(eta);
    return histEff->GetBinContent(binX, binY);
  };

  // Evaluating the efficiency for the leptons and the total efficiency
  auto df_weighted =
      df_B.Define("eff_lp", addWeights, {"lp_P_MeV", "lp_eta__"})
          .Define("eff_lm", addWeights, {"lm_P_MeV", "lm_eta__"})
          .Define("eff",
                  [](double eff_lp, double eff_lm) { return eff_lp * eff_lm; },
                  {"eff_lp", "eff_lm"});

  // Saving the weighted data in a file
  df_weighted.Snapshot("DecayTree", "temp_unbinned.root", {"D0_M_MeV", "eff"});
  TFile *fileIn_bkg = TFile::Open("temp_unbinned.root");
  TTree *tree_bkg = (TTree *)fileIn_bkg->Get("DecayTree");

  // Creating the RooDataHist from the TTree
  auto histdata =
      df_weighted.Histo1D({"data", "data", 75, 1700, 1900}, "D0_M_MeV", "eff");
  RooDataHist data_bkg("data_bkg", "Weighted histogram", D0_M_bkg,
                       Import(*histdata));

  std::cout << "Defining Fit model for bkg..." << '\n';

  // Defining the fit model: Bukin function
  RooRealVar mean_bkg("mean_bkg", "Peak position", 1860, 1800, 1900);
  RooRealVar width("width", "Width", 40, 10, 100);
  RooRealVar asym("asym", "Asymmetry", 0.0, -1, 1);
  RooRealVar rhoL("rhoL", "Left tail", 0.1, -10, 10);
  RooRealVar rhoR("rhoR", "Right tail", 0.1, -10, 10);

  RooBukinPdf bkg("bukin", "Bukin PDF", D0_M_bkg, mean_bkg, width, asym, rhoL,
                  rhoR);

  // Plot fit and data
  PlotFitPlots(D0_M_bkg, data_bkg, bkg, canvas, 0, res);

  // Saving parameters of the fit
  RooArgSet *params_bkg = bkg.getParameters(D0_M_bkg);
  params_bkg->writeToFile("paramsBkg.txt");

  std::cout << "Finish!" << '\n';
}

// ######################################################### //
// #################### Functions body ##################### //
// ######################################################### //

// APPLY FILTERS to SIM: don't need to filter on brem1/0 cathegory bc it's
// already done in the data
ROOT::RDF::RNode ApplyFiltersSim(ROOT::RDF::RNode df) {

  std::cout << "Hello, I'm filtering bin 3 events on RapidSim data..." << '\n';
  return df.Filter("Jpsi_M_MeV > 675 && Jpsi_M_MeV < 875", "Jpsi_M_MeV");
}

// FUNCTION TO FIT AND PLOT SIM DATA WITH CRYSTAL BALL
void PlotFitPlots(RooRealVar D0_M, RooDataSet data, RooAbsPdf &model,
                  TCanvas *canvas, int iBrem, std::string name = "fitSignal") {

  std::cout << "Plotting fit plots ..." << '\n';

  // Create pads for upper and lower plots
  // The upper pad will show the fit and data, while the lower pad will show the
  // pull distribution
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

  // Perform the fit
  RooNLLVar *nll = new RooNLLVar("nll", "nll", model, data);

  RooMinimizer m1(*nll);
  m1.setVerbose(kFALSE);
  m1.setPrintLevel(3);
  m1.setStrategy(1);
  m1.migrad();
  m1.hesse();

  RooFitResult *result = m1.save();
  result->Print("v");

  unsigned int i_minimization = 0;

  // Plot the fit and data
  RooPlot *frame = D0_M.frame();
  data.plotOn(frame);
  model.plotOn(frame);

  frame->GetXaxis()->SetTitleSize(0);

  frame->Draw();

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  int entries1 = data.sumEntries();
  leg->AddEntry(frame->getObject(0), "RapidSim Data", "lep");
  leg->AddEntry((TObject *)0, Form("\tEntries = %d", entries1), "");
  leg->AddEntry(frame->getObject(1), "Fit model", "l");
  leg->Draw();

  // lower pad for the pull distribution
  lowerPad.cd();
  RooHist *hpull0 = frame->pullHist();
  RooPlot *frameP = D0_M.frame();
  frameP->SetTitle("Pull Distribution");
  frameP->addPlotable(hpull0, "P");
  hpull0->SetMarkerSize(0.7);
  frameP->GetXaxis()->SetLabelSize(0.15);
  frameP->GetYaxis()->SetLabelSize(0.12);
  frameP->GetXaxis()->SetTitleSize(0.2);
  frameP->GetXaxis()->SetTitle("m(D^{0}) [MeV]");
  frameP->GetYaxis()->SetTitleSize(0.2);
  frameP->Draw();

  canvas->SaveAs((name + ".png").c_str());
  canvas->Clear();
  
}

void PlotFitPlots(RooRealVar D0_M, RooDataHist data, RooAbsPdf &model,
                  TCanvas *canvas, int iBrem, RooFitResult *result) {

  std::cout << "Plotting fit plots ..." << '\n';

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

  RooPlot *frame = D0_M.frame();
  data.plotOn(frame);
  model.plotOn(frame);

  frame->GetXaxis()->SetTitleSize(0);

  frame->Draw();

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  int entries1 = data.sumEntries();
  leg->AddEntry(frame->getObject(0), "RapidSim Data", "lep");
  leg->AddEntry(frame->getObject(1), "Bukin pdf", "l");
  leg->Draw();

  lowerPad.cd();
  RooHist *hpull0 = frame->pullHist();
  RooPlot *frameP = D0_M.frame();
  frameP->SetTitle("Pull Distribution");
  frameP->addPlotable(hpull0, "P");
  hpull0->SetMarkerSize(0.7);
  frameP->GetXaxis()->SetLabelSize(0.15);
  frameP->GetYaxis()->SetLabelSize(0.12);
  frameP->GetXaxis()->SetTitleSize(0.2);
  frameP->GetXaxis()->SetTitle("m(D^{0}) [MeV]");
  frameP->GetYaxis()->SetTitleSize(0.2);
  frameP->Draw();

  std::string name = ("fitRSBkgEff_B" + std::to_string(iBrem) + "bin3").c_str();
  canvas->SaveAs((name + ".png").c_str());
  canvas->Clear();
  canvas->Close();

  std::string filename = name + "_fit_results.txt";
  std::ofstream outFile(filename);

  if (outFile.is_open()) {

    outFile << "Fit parameters:\n";
    result->printMultiline(outFile, 1, true);

    outFile << "\nStato della convergenza:\n";
    outFile << "Status: " << result->status() << " (0 = successo)\n";
    outFile << "Qualità della matrice di covarianza: " << result->covQual()
            << " (3 = ottima, 2 = buona, 1 = accettabile, 0 = pessima o non "
               "calcolata)\n";

    outFile << "\nMatrice di covarianza:\n";
    for (int k = 0; k < result->covarianceMatrix().GetNrows(); k++) {
      for (int j = 0; j < result->covarianceMatrix().GetNcols(); j++) {
        outFile << result->covarianceMatrix()(k, j) << "\t";
      }
      outFile << "\n";
    }

    outFile.close();
    std::cout << "Risultati salvati in fit_results.txt" << std::endl;
  } else {
    std::cerr << "Errore: impossibile aprire il file!" << std::endl;
  }
}
