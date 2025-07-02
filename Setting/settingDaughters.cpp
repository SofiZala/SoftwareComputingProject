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
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// THE AIM OF THIS CODE IS TO CREATE 2 ROOTFILES.

// FIRST:  .ROOT FILE CONTAINING A TGRAPH OF THE SIGMA AS A FUNCTION OF RANGE OF
// MOMENTUM FOR EACH DAUGHTER PARTICLE THAT IS A HADRON IN THE DECAY.
// THE ROOTFILE WILL BE USED AS INPUT FOR RAPIDSIM TO SMEAR THE MOMENTUM
// DISTRIBUTIONS OF SUCH DAUGHTER PARTICLES. FOR EACH PARTICLE:
//         1. THE DATASET IS DIVIDED IN MOMENTUM RANGE (5 BINS) WITH THE SAME
//         STATISTICS
//         2. EACH BIN IS FITTED WITH A GAUSSIAN PDF: SIGMA VALUES AND ERRORS
//         ARE EXTRACTED
//         3. CONSTRUCTION OF THE TGRAPH (MEAN VALUE IN THE RANGE OF P - SIGMA
//         FROM THE FIT)

// SECOND: .ROOT FILE CONTAINING THE DISTRIBUTIONS OF (P-PTRUE)/PTRUE FOR THE
// ELECTRONS.

// functions declaration
ROOT::RDF::RNode ApplyFiltersMC(ROOT::RDF::RNode df, int iBrem);
std::vector<double> EvaluatePerc(ROOT::RDF::RNode df_mc, std::string variableP);
std::vector<ROOT::RDF::RResultPtr<TH1D>>
ApplyFiltersPrange(ROOT::RDF::RNode df_mc, std::vector<double> quantiles,
                   std::string variableP, std::string variableSigma);
TGraphErrors *CreateGraph(std::string varName,
                          std::vector<ROOT::RDF::RResultPtr<TH1D>> varInBins,
                          std::vector<double> quantiles, TCanvas *canvas);
TGraphErrors *ProcessParticle(ROOT::RDF::RNode df, const std::string &name,
                              TFile *outFile, TCanvas *canvas);
std::vector<ROOT::RDF::RResultPtr<TH1D>> MakeElectronHisto(ROOT::RDF::RNode df);
void PlotElectronHisto(TFile *outFileElectrons,
                       std::vector<ROOT::RDF::RResultPtr<TH1D>> MCHistos,
                       TCanvas *canvas);

void CreateFolder(const char *folder);
void ChangeDirectory(const char *folder);

void settingDaughters(int const range = -1) {

  // Opening the rootfile containing the true and reconstructed momentum
  std::cout << "Defining RDataFrame for simulated events..." << '\n';
  std::string pathToFileMC = "inputFile.root";
  ROOT::RDataFrame df_mc("DecayTree_preselection", pathToFileMC);

  // Definition of the resolution sigma as: (P-Ptrue) / Ptrue
  auto df_mc_r =
      df_mc.Range(0, range)
          .Define("spip_sigmaP", [](float x, float y) { return (x - y) / y; },
                  {"spip_P", "spip_TRUEP"})
          .Define("Km_sigmaP", [](float x, float y) { return (x - y) / y; },
                  {"Km_P", "Km_TRUEP"})
          .Define("pip_sigmaP", [](float x, float y) { return (x - y) / y; },
                  {"pip_P", "pip_TRUEP"})
          .Define("lm_sigmaP", [](float x, float y) { return (x - y) / y; },
                  {"lm_P", "lm_TRUEP"})
          .Define("lp_sigmaP", [](float x, float y) { return (x - y) / y; },
                  {"lp_P", "lp_TRUEP"});

  TCanvas *canvas = new TCanvas();

  ////////////////////////////////////////////////////////////////
  /////////////////////////////// BREM 1 /////////////////////////
  ////////////////////////////////////////////////////////////////

  // Filtering the data for brem cathegory: B1
  auto df_mc1 = ApplyFiltersMC(df_mc_r, 1);

  // Creating the output directory
  CreateFolder("HadronBrem1");
  ChangeDirectory("HadronBrem1");

  // Creating the output ROOT file for hadrons
  TFile *outFileB1Hadrons = new TFile("Run3_HadronSmearing.root", "RECREATE");

  // List of particles to process
  std::vector<std::string> particleNames = {"spip", "Km", "pip"};
  std::vector<TGraphErrors *> graphs;

  // Process each particle and create TGraphs
  for (const auto &name : particleNames) {
    TGraphErrors *g = ProcessParticle(df_mc1, name, outFileB1Hadrons, canvas);
    graphs.push_back(g);
  }

  outFileB1Hadrons->Close();

  ChangeDirectory("..");

  // Creating the output directory for electrons
  CreateFolder("ElectronBrem1");
  ChangeDirectory("ElectronBrem1");

  // Creating the output ROOT file for electrons
  TFile *outFileB1Electrons =
      new TFile("Run3_ElectronSmearing.root", "RECREATE");
  // Creating histograms for electrons
  std::vector<ROOT::RDF::RResultPtr<TH1D>> MCHistos1 =
      MakeElectronHisto(df_mc1);
  // Plotting the histograms for electrons
  PlotElectronHisto(outFileB1Electrons, MCHistos1, canvas);
  outFileB1Electrons->Close();

  ChangeDirectory("..");

  ////////////////////////////////////////////////////////////////
  /////////////////////////////// BREM 0 /////////////////////////
  ////////////////////////////////////////////////////////////////

  // Filtering the data for brem cathegory: B0
  auto df_mc0 = ApplyFiltersMC(df_mc_r, 0);

  // Creating the output directory
  CreateFolder("HadronBrem0");
  ChangeDirectory("HadronBrem0");

  TFile *outFileB0Hadrons = new TFile("Run3_HadronSmearing.root", "RECREATE");

  // Process each particle and create TGraphs
  for (const auto &name : particleNames) {

    // ProcessParticle() calls all the functions needed to process the data and
    // obtain the final TGraphs
    ProcessParticle(df_mc0, name, outFileB0Hadrons, canvas);
    TGraphErrors *g =
        (TGraphErrors *)outFileB0Hadrons->Get((name + "_TGraph").c_str());
    graphs.push_back(g);
  }

  outFileB0Hadrons->Close();

  ChangeDirectory("..");

  // Creating the output directory for electrons
  CreateFolder("ElectronBrem0");
  ChangeDirectory("ElectronBrem0");

  // Creating the output ROOT file for electrons
  TFile *outFileB0Electrons =
      new TFile("Run3_ElectronSmearing.root", "RECREATE");

  // Creating histograms for electrons
  std::vector<ROOT::RDF::RResultPtr<TH1D>> MCHistos0 =
      MakeElectronHisto(df_mc0);
  // Plotting the histograms for electrons
  PlotElectronHisto(outFileB0Electrons, MCHistos0, canvas);
  outFileB0Electrons->Close();

  canvas->Close();
}

// ######################################################### //
// #################### Functions body ##################### //
// ######################################################### //

// Apply filters to data: filters needed to distinguish brem1/0 cathegory
ROOT::RDF::RNode ApplyFiltersMC(ROOT::RDF::RNode df, int iBrem) {

  std::cout << "Applying brem filter " << iBrem << " on MC data ..." << '\n';

  if (iBrem == 0) {
    std::cout << "Brem 0 filter applied!" << '\n';
    return df.Filter("lp_HASBREM == 0 && lm_HASBREM == 0", "brem0");
  } else if (iBrem == 1) {
    std::cout << "Brem 1 filter applied" << '\n';
    return df.Filter("lp_HASBREM == 1 || lm_HASBREM == 1", "brem1");
  } else {
    std::cout << "No Brem filter applied!" << '\n';
    return df;
  }
}

// This function is used only to call all functions needed to process the data
// and create the TGraphs
TGraphErrors *ProcessParticle(ROOT::RDF::RNode df, const std::string &name,
                              TFile *outFile, TCanvas *canvas) {
  std::cout << name << ": creating histograms..." << '\n';

  std::string variableP = name + "_P";
  std::string variableSigma = name + "_sigmaP";

  // evaluate the quantiles of the momentum distribution to divide the dataset
  // in 5 bins
  auto quantiles = EvaluatePerc(df, variableP);

  // Divide the dataset in 5 bins according to the range of momentum and create
  // histograms
  auto histsInBins =
      ApplyFiltersPrange(df, quantiles, variableP, variableSigma);

  // Create a TGraph from the histograms in the bins
  TGraphErrors *graph = CreateGraph(name, histsInBins, quantiles, canvas);

  // Saving
  outFile->cd();
  graph->SetName((name + "_TGraph").c_str());
  graph->Write();
  return graph;
}

// EvaluatePerc calculates the quantiles of the momentum distribution
std::vector<double> EvaluatePerc(ROOT::RDF::RNode df_mc,
                                 std::string variableP) {

  auto P_values =
      df_mc.Define("dummy", "0").Snapshot("DecayTree_full", "temp.root");

  TFile ftemp("temp.root");
  TTree *t = (TTree *)ftemp.Get("DecayTree_full");

  std::vector<float> all_P;

  float P_value;
  t->SetBranchAddress(variableP.c_str(), &P_value);

  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    all_P.push_back(P_value);
  }
  ftemp.Close();
  gSystem->Unlink("temp.root");

  // Sort according to P values
  std::sort(all_P.begin(), all_P.end());

  // evaluate quantiles
  size_t n = all_P.size();

  double q20 = all_P[n * 0.2];
  double q40 = all_P[n * 0.4];
  double q60 = all_P[n * 0.6];
  double q80 = all_P[n * 0.8];

  std::vector<double> quantiles{q20, q40, q60, q80

  };

  return quantiles;
}

// ApplyFiltersPrange divides the dataset in 5 subsamples, each containing the
// same number of events thanks to the previous evaluation of quantiles. Each
// subsamples form a new dataset which is returned as a TH1D
std::vector<ROOT::RDF::RResultPtr<TH1D>>
ApplyFiltersPrange(ROOT::RDF::RNode df_mc, std::vector<double> quantiles,
                   std::string variableP, std::string variableSigma) {

  std::cout << "Filtering the dataset and dividiving in different bins "
               "according to the range of momentum decided!"
            << '\n';

  auto h_df_bin0 =
      df_mc
          .Filter([quantiles](float p) { return p <= quantiles[0]; },
                  {variableP})
          .Histo1D({"sigmaP",
                    ("sigmaP with P in [ 0 - " +
                     std::to_string(static_cast<int>(quantiles[0])) + " ] ")
                        .c_str(),
                    200, -0.03, 0.03},
                   variableSigma);

  auto h_df_bin1 =
      df_mc
          .Filter(
              [quantiles](float p) {
                return p > quantiles[0] && p <= quantiles[1];
              },
              {variableP})
          .Histo1D({"sigmaP",
                    ("sigmaP with P in [ " +
                     std::to_string(static_cast<int>(quantiles[0])) + "-" +
                     std::to_string(static_cast<int>(quantiles[1])) + " ] ")
                        .c_str(),
                    200, -0.03, 0.03},
                   variableSigma);

  auto h_df_bin2 =
      df_mc
          .Filter(
              [quantiles](float p) {
                return p > quantiles[1] && p <= quantiles[2];
              },
              {variableP})
          .Histo1D({"sigmaP",
                    ("sigmaP with P in [ " +
                     std::to_string(static_cast<int>(quantiles[1])) + "-" +
                     std::to_string(static_cast<int>(quantiles[2])) + " ] ")
                        .c_str(),
                    200, -0.03, 0.03},
                   variableSigma);

  auto h_df_bin3 =
      df_mc
          .Filter(
              [quantiles](float p) {
                return p > quantiles[2] && p <= quantiles[3];
              },
              {variableP})
          .Histo1D({"sigmaP",
                    ("sigmaP with P in [ " +
                     std::to_string(static_cast<int>(quantiles[2])) + "-" +
                     std::to_string(static_cast<int>(quantiles[3])) + " ] ")
                        .c_str(),
                    200, -0.03, 0.03},
                   variableSigma);

  auto h_df_bin4 =
      df_mc
          .Filter([quantiles](float p) { return p > quantiles[3]; },
                  {variableP})
          .Histo1D(
              {"sigmaP",
               ("sigmaP with P in [ " +
                std::to_string(static_cast<int>(quantiles[3])) + " - inf ] ")
                   .c_str(),
               200, -0.03, 0.03},
              variableSigma);

  std::vector<ROOT::RDF::RResultPtr<TH1D>> histo_vector{
      h_df_bin0, h_df_bin1, h_df_bin2, h_df_bin3, h_df_bin4};

  return histo_vector;
}

// CreateGraph creates a TGraph from the histograms in the bins
TGraphErrors *CreateGraph(std::string varName,
                          std::vector<ROOT::RDF::RResultPtr<TH1D>> varInBins,
                          std::vector<double> quantiles, TCanvas *canvas) {

  std::cout << varName + ": Creating TGraph for ..." << '\n';

  if (!varInBins.empty()) {

    CreateFolder(varName.c_str());
    ChangeDirectory(varName.c_str());

    std::vector<double> mean_vec;
    std::vector<double> sigma_vec;
    std::vector<double> sigma_error_vec;

    // For each dataset in the bins, fit a Gaussian and extract the mean and
    // sigma
    for (int i = 0; i < varInBins.size(); i++) {

      std::string histName = varInBins[i]->GetName();
      canvas->cd();

      // variable definition
      RooRealVar sigmaP("sigmaP", "sigmaP", -0.03, 0.03);

      // parameters definition
      RooRealVar mean("mean", "Mean Mass", 0, -0.03, 0.03);
      RooRealVar sigma("sigma", "Width", 0.01, 0, 2);
      RooGaussian gauss("gauss", "gauss", sigmaP, mean, sigma);

      // Create a frame for plotting
      RooPlot *frame = sigmaP.frame();

      // Getting the histogram
      TH1D *h_var = varInBins[i].GetPtr();
      RooDataHist data("data", "data", RooArgList(sigmaP), h_var);

      // Fit and Plot data
      data.plotOn(frame);
      gauss.fitTo(data);
      gauss.plotOn(frame);

      // Extracting fitted parameters and errors
      double fittedMean = mean.getVal();
      double fittedSigma = sigma.getVal();
      double errorMean = mean.getError();
      double errorSigma = sigma.getError();

      // save results
      if (i < 5) {
        mean_vec.push_back(fittedMean);
        sigma_vec.push_back(fittedSigma);
        sigma_error_vec.push_back(sigma.getError());
      }

      std::string histTitle = varInBins[i]->GetTitle();
      frame->GetXaxis()->SetTitle((varName + "_" + histTitle).c_str());
      frame->SetTitle("");
      frame->Draw();

      TLegend *leg = new TLegend(0.55, 0.7, 0.9, 0.9);
      leg->AddEntry(frame->getObject(0), "Data", "lep");
      leg->AddEntry(frame->getObject(1), "Fit Gauss", "l");

      int entries = h_var->GetEntries();
      leg->AddEntry((TObject *)0, Form("\tEnt. = %d", entries), "");
      leg->AddEntry((TObject *)0,
                    Form("\tMean = %.6f #pm %.6f", fittedMean, errorMean), "");
      leg->AddEntry((TObject *)0,
                    Form("\tStd dev= %.6f #pm %.6f", fittedSigma, errorSigma),
                    "");
      leg->SetTextSize(0.03);

      leg->Draw();

      canvas->SaveAs(
          (varName + "_" + histName + "_B1_bin" + std::to_string(i) + ".png")
              .c_str());
      canvas->Clear();
    }

    std::ofstream out_result((varName + "_results.txt").c_str());
    if (out_result.is_open()) {
      out_result << varName + " results:\n";
      out_result << "Range: \t \t \t Mean values: \t \t \t Sigma values: \n"
                 << "Bin 0: [0, " << quantiles[0] << "] MeV \t" << mean_vec[0]
                 << "\t \t" << sigma_vec[0] << "\n"
                 << "Bin 1: [" << quantiles[0] << ", " << quantiles[1]
                 << "] MeV \t" << mean_vec[1] << "\t" << sigma_vec[1] << "\n"
                 << "Bin 2: [" << quantiles[1] << ", " << quantiles[2]
                 << "] MeV \t" << mean_vec[2] << "\t" << sigma_vec[2] << "\n"
                 << "Bin 3: [" << quantiles[2] << ", " << quantiles[3]
                 << "] MeV \t" << mean_vec[3] << "\t" << sigma_vec[3] << "\n"
                 << "Bin 4: [" << quantiles[3] << ", inf] MeV \t" << mean_vec[4]
                 << "\t \t" << sigma_vec[4] << "\n";

    } else {
      std::cerr << "Error opening file for writing." << std::endl;
    }

    double sigma_vals[5];
    double sigma_error_vals[5];
    for (int i = 0; i < 5; i++) {
      sigma_vals[i] = sigma_vec[i];
      sigma_error_vals[i] = sigma_error_vec[i];
    }

    double bin0 = (0 + quantiles[0]) / 2;
    double bin1 = (quantiles[0] + quantiles[1]) / 2;
    double bin2 = (quantiles[1] + quantiles[2]) / 2;
    double bin3 = (quantiles[2] + quantiles[3]) / 2;
    double bin4 = (quantiles[3] + 100000) / 2;

    double bin_cent[5] = {bin0, bin1, bin2, bin3, bin4};

    TGraphErrors *graph =
        new TGraphErrors(5, bin_cent, sigma_vals, 0, sigma_error_vals);

    graph->SetTitle("Sigma vs P; P (MeV); Sigma");

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    TCanvas *c_graph = new TCanvas("c_graph", "Sigma vs P", 800, 600);
    graph->Draw("APL");
    c_graph->SaveAs((varName + "_graphB1.png").c_str());

    ChangeDirectory("..");

    return graph;
  }
  return nullptr;
}

// save histo for electrons
std::vector<ROOT::RDF::RResultPtr<TH1D>>
MakeElectronHisto(ROOT::RDF::RNode df) {

  auto h_lm_sigmaP =
      df.Histo1D({"lm_sigmaP", "lm_sigmaP", 150, -0.8, 0.3}, "lm_sigmaP");
  auto h_lp_sigmaP =
      df.Histo1D({"lp_sigmaP", "lp_sigmaP", 150, -0.8, 0.3}, "lp_sigmaP");

  std::vector<ROOT::RDF::RResultPtr<TH1D>> electronHistos{
      h_lm_sigmaP,
      h_lp_sigmaP,
  };

  return electronHistos;
}

// plot histograms for electrons
void PlotElectronHisto(TFile *outFileElectrons,
                       std::vector<ROOT::RDF::RResultPtr<TH1D>> MCHistos,
                       TCanvas *canvas) {

  if (!MCHistos.empty()) {

    outFileElectrons->cd();
    canvas->cd();
    std::cout << "Plotting histos of simulated data..." << '\n';

    for (int i = 0; i < MCHistos.size(); i++) {

      std::string histTitle = MCHistos[i]->GetTitle();

      MCHistos[i]->SetName((histTitle).c_str());
      MCHistos[i]->Write("");

      // DECOMMENTA SE LI VUOI VEDERE SU CANVAS MENTRE LI SALVA NEL .ROOT
      MCHistos[i]->SetStats(0);
      MCHistos[i]->GetXaxis()->SetTitle((histTitle + "[GeV]").c_str());
      MCHistos[i]->GetXaxis()->SetTitleSize(0.04);
      MCHistos[i]->GetXaxis()->SetLabelSize(0.04);
      MCHistos[i]->SetLineColor(kTeal - 1);
      MCHistos[i]->SetFillColor(kGreen - 10);
      MCHistos[i]->SetFillStyle(3003);
      MCHistos[i]->SetLineWidth(1);
      MCHistos[i]->Draw("hist");
      MCHistos[i]->SetTitle("");

      // LEGEND
      TLegend *leg = new TLegend(0.70, 0.65, 0.92, 0.88);
      // Extracting manually values for mean, std dev and entries
      double mean2 = MCHistos[i]->GetMean();
      double std2 = MCHistos[i]->GetRMS();
      int entries2 = MCHistos[i]->GetEntries();
      TH1D *MCHisto = MCHistos[i].GetPtr();
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.032);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->AddEntry(MCHisto, "MC", "l");
      leg->AddEntry((TObject *)0, Form("\tEnt. = %d", entries2), "");
      leg->AddEntry((TObject *)0, Form("\tMean = %.2f", mean2), "");
      leg->AddEntry((TObject *)0, Form("\tStd dev= %.2f", std2), "");
      leg->Draw("SAME");
      canvas->SaveAs((histTitle + ".png").c_str());
      canvas->Clear();
    }

    std::cout << "MC histograms saved in .root file!" << '\n';
  }
}

void CreateFolder(const char *folder) {
  struct stat info;
  if (stat(folder, &info) != 0) {
    mkdir(folder, 0777);
  }
}

void ChangeDirectory(const char *folder) {
  if (chdir(folder) == 0) {
    std::cout << "Changed directory to: " << folder << '\n';
  } else {
    perror("Error changing directory");
  }
}
