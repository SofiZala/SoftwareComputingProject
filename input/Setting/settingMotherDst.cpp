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

// AIM OF THIS CODE:
// CREATE A .ROOT FILE CONTAINING:
//                          - 1D HISTOGRAM OF PT,
//                          - 1D HISTOGRAM OF ETA,
//                          - 1D HISTOGRAM OF PHI
// OF THE MOTHER PARTICLE Dst WHICH WILL BE USED AS INPUT FOR THE RAPID SIM
// SIMULATION. THE DISTRIBUTIONS ARE SAMPLED ACCORDING TO PARAMETERS OF THE PDFs
// GIVEN IN INPUT FOR PT AND ETA, WHILE FOR PHI THE SAMPLE IS MADE FROM A
// DISTRIBUTION. RAPID SIM WILL USE THIS FILE TO SIMULATE THE Dst PROPERTIES

ROOT::RDF::RNode ApplyFiltersMC(ROOT::RDF::RNode df, int iBrem);
std::vector<ROOT::RDF::RResultPtr<TH1D>> MakeVarsHistos(ROOT::RDF::RNode df);
void ProcessDstInformation(ROOT::RDF::RNode df_mc, int brem = -1);

void CreateFolder(const char *folder);
void ChangeDirectory(const char *folder);

void settingMotherDst(int const range = -1) {

  // Creation and opening of the output file: .root file containing the required
  // histogram

  //////////////////////////////////////////////////////////////////
  //////////////   1D-histogram of PT , ETA AND PHI ////////////////
  /////////////////////////////////////////////////////////////////

  // Open the rootfile ontaining the phi info
  std::string pathToFileMC = "inputFile.root";
  ROOT::RDataFrame df_mc("DecayTree_preselection", pathToFileMC);

  // Filtering the data for brem cathegory: B1 and B0
  auto df_mc1 = ApplyFiltersMC(df_mc, 1);
  auto df_mc0 = ApplyFiltersMC(df_mc, 0);

  // Creation of the output of this code for both categories
  ProcessDstInformation(df_mc1, 1);
  ProcessDstInformation(df_mc0, 0);
}

// ######################################################### //
// #################### Functions body ##################### //
// ######################################################### //

// Apply filters to MC: filters needed for brem1/0 cathegory
ROOT::RDF::RNode ApplyFiltersMC(ROOT::RDF::RNode df, int iBrem) {

  std::cout << "Hello, I'm applying brem filter " << iBrem << " on MC data..."
            << '\n';

  auto df_brem_bin = df;
  if (iBrem == 0) {
    std::cout << "Brem 0 filter applied" << '\n';
    // lp and lm are the names of the electrons (e+ e-)
    // HASBREM is a variable that indicates if the electron has Bremstrahlung or
    // not
    return df_brem_bin =
               df.Filter("lp_HASBREM == 0 && lm_HASBREM == 0", "brem0");
  } else if (iBrem == 1) {
    std::cout << "Brem 1 filter applied" << '\n';
    return df_brem_bin =
               df.Filter("lp_HASBREM == 1 || lm_HASBREM == 1", "brem1");
  } else {
    std::cout << "No Brem filter applied" << '\n';
    return df_brem_bin;
  }
}

// Function to create histograms for the variables of interest
void ProcessDstInformation(ROOT::RDF::RNode df_mc, int brem = -1) {

  // Entering the correct directory for the brem value
  CreateFolder(("DstBrem" + std::to_string(brem)).c_str());
  ChangeDirectory(("DstBrem" + std::to_string(brem)).c_str());

  // Creating the output file for the histograms
  TFile *outFile = new TFile("LHCc14.root", "RECREATE");

  // phi: the phi distribution is sampled from the Dst_TRUEphi variable of the
  // mc data in in the input rootfile
  auto h_phi =
      df_mc.Histo1D({"phi", "Dst_TRUEphi", 25, 0, 2 * 3.14}, "Dst_TRUEphi");

  // For PT and ETA, the distributions are sampled from the Johnson pdf and
  // Gauss pdf respectively according to the parameters read from the file
  // params_Dst_bremX.txt (X=0 or X=1)
  RooRealVar x("x", "pT", 0, 25);
  RooRealVar y("y", "eta", 1.5, 5);

  // Declaration and reading the parameters of the Johnson pdf and Gauss pdf
  // from file and definition of the pdf
  RooRealVar mean_PT("mean_PT", "mean", 0, 10);
  RooRealVar sigma_PT("sigma_PT", "sigma", 0, 10);
  RooRealVar gamma_PT("gamma_PT", "gamma", 0, 1);
  RooRealVar delta_PT("delta_PT", "delta", 0, 1);
  RooRealVar mean_eta("mean_eta", "mean", 1.5, 5);
  RooRealVar sigma_eta("sigma_eta", "sigma", 0, 1);

  RooArgSet params(mean_PT, sigma_PT, gamma_PT, delta_PT, mean_eta, sigma_eta);
  if (params.readFromFile(
          ("../params_Dst_brem" + std::to_string(brem) + ".txt").c_str())) {
    std::cerr << "Error! File not opened.\n";
    return;
  }

  mean_PT.setConstant(true);
  sigma_PT.setConstant(true);
  gamma_PT.setConstant(true);
  delta_PT.setConstant(true);
  mean_eta.setConstant(true);
  sigma_eta.setConstant(true);

  RooJohnson johnson("johnson", "JohnsonSU PDF", x, mean_PT, sigma_PT, gamma_PT,
                     delta_PT);

  // Sampling of events according to the Johnson distribution
  RooDataSet *sample = johnson.generate(x, 10000);

  // Creation of the 1D histogram for PT
  TH1D *h_gen = new TH1D("pT", "Sampled pT; pT; Events", 36, 0, 25);

  for (int i = 0; i < sample->numEntries(); ++i) {
    const RooArgSet *row = sample->get(i);
    h_gen->Fill(row->getRealValue("x"));
  }

  // --- Plot ---
  TCanvas *c = new TCanvas("c", "Sampled Histogram", 800, 600);

  h_gen->Draw();
  c->SaveAs("pT.png");
  c->Clear();

  // Creation of the 1D histogram for eta
  RooGaussian gauss("gauss", "gauss PDF", y, mean_eta, sigma_eta);

  // Sampling of events according to the Gaussian distribution
  RooDataSet *sample1 = gauss.generate(y, 10000);

  // Creation of the 1D histogram for ETA
  TH1D *h_gen1 = new TH1D("eta", "Sampled eta; eta; Events", 36, 0, 5);

  for (int i = 0; i < sample1->numEntries(); ++i) {
    const RooArgSet *row1 = sample1->get(i);
    h_gen1->Fill(row1->getRealValue("y"));
  }

  h_gen1->Draw();
  c->SaveAs("eta.png");
  c->Clear();

  h_phi->Draw();
  c->SaveAs("phi.png");
  c->Close();

  // Writing the histograms to the output file

  outFile->cd();
  h_gen->Write();
  h_gen1->Write();
  h_phi->Write();
  outFile->Close();
  ChangeDirectory("..");
}

// Function to create a folder if it does not exist
void CreateFolder(const char *folder) {
  struct stat info;
  if (stat(folder, &info) != 0) {
    mkdir(folder, 0777);
  }
}

// Function to change the current working directory
void ChangeDirectory(const char *folder) {
  if (chdir(folder) == 0) {
    std::cout << "Changed directory to: " << folder << '\n';
  } else {
    perror("Error changing directory");
  }
}