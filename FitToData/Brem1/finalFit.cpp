#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCrystalBall.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TSystem.h"
#include <sys/stat.h>
#include <unistd.h>

// THE AIM OF THIS CODE IS TO FIT DATA
// USING PARAMETERS FROM PRELIMINARY FITS TO SIGNAL AND BACKGROUND SIMULATIONS

using namespace RooFit;

void Plot(RooRealVar D0_M, RooDataSet data, RooAddPdf model, double nSig,
          double nBkg, double lower, double upper);

void finalFit() {

  // Open the input data file and get the tree
  std::string inputDataFile = "data.root";
  std::cout << "Defining RDataFrame ..." << '\n';
  TFile *inputFile = TFile::Open(inputDataFile.c_str());
  TTree *treeIn = (TTree *)inputFile->Get("DecayTree_BDTcut");

  // Define the observable
  RooRealVar D0_M("D0_M", "m(D^{0}) [MeV]", 1600, 2150, "MeV");
  RooDataSet data("data", "data", D0_M, Import(*treeIn));

  // Define parameters for the fit
  // Signal -> Johnson PDF
  RooRealVar mean = RooRealVar("mean", "mean", 1846.5, 1800, 1900);
  RooRealVar sigma = RooRealVar("sigma", "sigma", 48.7, 0, 100);
  RooRealVar gamma = RooRealVar("gamma", "gamma", 0.19336, -100, 100);
  RooRealVar delta = RooRealVar("delta", "delta", 1.4173, -1000, 1000);
  RooJohnson signal =
      RooJohnson("signal", "signal", D0_M, mean, sigma, gamma, delta);

  // Bkg partially reconstructed -> Bukin PDF
  RooRealVar mean_bkg("mean_bkg", "Mean Mass", 1700, 1500, 1800);
  RooRealVar width("width", "Width", 50, 1, 200);
  RooRealVar asym("asym", "asymm", -0.2, -1, 1);
  RooRealVar rhoL("rhoL", "rhoL", -0.1, -10, 10);
  RooRealVar rhoR("rhoR", "rhoR", 0.1, -10, 10);
  RooBukinPdf partReco("partReco", "partReco", D0_M, mean_bkg, width, asym,
                       rhoL, rhoR);

  // Comb bkg -> Chebychev polynomial
  RooRealVar c0("c0", "coefficient #0", -0.78, -100, 100);
  RooChebychev cheb("cheb", "Chebychev PDF", D0_M, RooArgList(c0));

  // extended model
  RooRealVar nSig("nSig", "Number of signal events", 2074, 1, 1e6);
  RooRealVar nCheb("nCheb", "Number of background events", 1212, 1, 1e6);
  RooRealVar nPartReco("nPartReco", "nPartReco", 671, 1, 1e6);
  RooAddPdf model{"model", "model", RooArgList(signal, cheb, partReco),
                  RooArgList(nSig, nCheb, nPartReco)};

  // Read and set parameters from preliminary fit
  RooArgSet *params = model.getParameters(D0_M);

  if (!gSystem->AccessPathName("../../PreliminaryFit/Brem1/paramsSignal.txt")) {
    params->readFromFile("../../PreliminaryFit/Brem1/paramsSignal.txt");
  } else {
    std::cerr << "File paramsSignal.txt not found!" << std::endl;
  }

  if (!gSystem->AccessPathName("../../PreliminaryFit/Brem1/paramsBkg.txt")) {
    params->readFromFile("../../PreliminaryFit/Brem1/paramsBkg.txt");
  } else {
    std::cerr << "File paramsBkg.txt not found!" << std::endl;
  }

  // Fixing constant parameters
  delta.setConstant(kTRUE);
  gamma.setConstant(kTRUE);
  mean_bkg.setConstant(kTRUE);
  width.setConstant(kTRUE);
  asym.setConstant(kTRUE);
  rhoL.setConstant(kTRUE);
  rhoR.setConstant(kTRUE);

  // Performing the fit until the quality of the fit is good
  RooFitResult *res;
  for (int i = 0; i < 15; ++i) {
    res = model.fitTo(data, RooFit::Save(true), Extended(true));

    // estimated distance to minimum
    float edm = res->edm();

    // quality of covariance matrix
    int covQual = res->covQual();

    // status of the fit
    float status = res->status();

    if (edm < 0.0001 && covQual == 3 && status == 0)
      break;
    std::cout << "------------- RETRYING THE FIT --------------- " << '\n';
  }

  // Print the fit results
  const TMatrixDSym &cor = res->correlationMatrix();
  const TMatrixDSym &cov = res->covarianceMatrix();

  std::cout << "Covariance Matrix: \n";
  cov.Print();
  std::cout << "Correlation Matrix: \n";
  cor.Print();

  // Evaluate the number of sig and bkg events in the signal region
  // Define the signal region as mean +/- 3 sigma
  double lowerBoundSigRegion = mean.getVal() - 3 * sigma.getVal();
  double upperBoundSigRegion = mean.getVal() + 3 * sigma.getVal();
  D0_M.setRange("signalRegion", lowerBoundSigRegion, upperBoundSigRegion);
  // Evaluate the numbers in range
  RooAbsReal *chebIntegral =
      cheb.createIntegral(D0_M, NormSet(D0_M), Range("signalRegion"));
  RooAbsReal *partIntegral =
      partReco.createIntegral(D0_M, NormSet(D0_M), Range("signalRegion"));
  RooAbsReal *sigIntegral =
      signal.createIntegral(D0_M, NormSet(D0_M), Range("signalRegion"));

  double nChebInSigRegion = chebIntegral->getVal() * nCheb.getVal();
  double nPartRecoInSigRegion = partIntegral->getVal() * nPartReco.getVal();
  double nSigInSigRegion = sigIntegral->getVal() * nSig.getVal();

  double totbkg = nChebInSigRegion + nPartRecoInSigRegion;

  // Evaluate errors
  double err_nSigInRegion = sigIntegral->getVal() * nSig.getError();
  double err_nChebInRegion = chebIntegral->getVal() * nCheb.getError();
  double err_nPartRecoInRegion = partIntegral->getVal() * nPartReco.getError();
  double cov_nCheb_nPartReco = cov(2, 3);
  double err_bkgtot =
      sqrt(pow(err_nChebInRegion, 2) + pow(err_nPartRecoInRegion, 2) +
           2 * chebIntegral->getVal() * partIntegral->getVal() *
               cov_nCheb_nPartReco);

  // Plot data, fit and results
  Plot(D0_M, data, model, nSigInSigRegion,
       nChebInSigRegion + nPartRecoInSigRegion, lowerBoundSigRegion,
       upperBoundSigRegion);

  res->Print("v");

  // Save the results to a file
  std::string filename = "FitData_results.txt";
  std::ofstream outFile(filename);

  if (outFile.is_open()) {

    outFile << "Fit parameters:\n";
    res->printMultiline(outFile, 1, true);

    outFile << "Status: " << res->status() << " (0 = success)\n";
    outFile << "Quality of the covariance matrix: " << res->covQual();

    outFile << "\nCovariance matrix:\n";
    for (int k = 0; k < res->covarianceMatrix().GetNrows(); k++) {
      for (int j = 0; j < res->covarianceMatrix().GetNcols(); j++) {
        outFile << res->covarianceMatrix()(k, j) << "\t";
      }
      outFile << "\n";
    }

    outFile.close();
    std::cout << "Fit results saved!" << std::endl;
  } else {
    std::cerr << "Errore: file not opened!" << std::endl;
  }

  // Save the parameters to a file
  params->writeToFile("paramsfit.txt");
}

// ############ FUNCTIONS ############# //

void Plot(RooRealVar D0_M, RooDataSet data, RooAddPdf model, double nSig,
          double nBkg, double lower, double upper) {
  std::cout << "Plotting plots ..." << '\n';

  // Create a canvas with two pads: one for the fit and one for the pull
  // distribution
  TCanvas *c = new TCanvas();
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

  // Upper: fit
  upperPad.cd();
  float binWidth = 17.;
  float min = 1600;
  float max = 2150;
  int nBins = (max - min) / binWidth;
  D0_M.setBins(nBins);
  RooPlot *frame = D0_M.frame();
  data.plotOn(frame, MarkerSize(0.7));

  model.plotOn(frame, Components("signal"), LineColor(kRed),
               LineStyle(kDashed));
  model.plotOn(frame, Components("cheb"), LineColor(kViolet),
               LineStyle(kDashed));
  model.plotOn(frame, Components("partReco"), LineColor(kGreen),
               LineStyle(kDashed));
  model.plotOn(frame);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->SetTitle("");
  frame->Draw();

  // legend for the objects in the frame
  TLegend *leg = new TLegend(0.55, 0.6, 0.9, 0.88);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(frame->getObject(0), "Data", "l");
  leg->AddEntry(frame->getObject(1),
                "Signal D^{0} #rightarrow K^{-}#pi^{+}e^{-}e^{+}", "l");
  leg->AddEntry(frame->getObject(2), "Comb. bkg.", "l");
  leg->AddEntry(frame->getObject(3), "Part. reco.", "l");
  leg->AddEntry(frame->getObject(4), "Fit", "l");
  leg->Draw();

  // legend to print numeric results
  TLegend *leg2 = new TLegend(0.05, 0.65, 0.35, 0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.05);

  TString nRegion = Form("Signal region: [%.0f, %.0f]", lower, upper);
  leg2->AddEntry((TObject *)0, nRegion, "");

  TString nSigStr = Form("N_{sig} = %.1f", nSig);
  leg2->AddEntry((TObject *)0, nSigStr, "");

  TString nBkgStr = Form("N_{bkg} = %.1f", nBkg);
  leg2->AddEntry((TObject *)0, nBkgStr, "");

  double frac = nSig / nBkg;
  TString nSigBkg = Form("N_{sig/bkg} = %.3f", frac);
  leg2->AddEntry((TObject *)0, nSigBkg, "");
  leg2->Draw();

  // lower : pull distribution
  lowerPad.cd();
  RooHist *hPull = frame->pullHist();
  RooPlot *frameP = D0_M.frame(Title("Pull Distribution"));
  frameP->addPlotable(hPull, "P");
  hPull->SetMarkerSize(0.7);
  frameP->GetXaxis()->SetLabelSize(0.2);
  frameP->GetYaxis()->SetLabelSize(0.15);
  frameP->GetXaxis()->SetTitleSize(0.2);
  frameP->GetYaxis()->SetNdivisions(302);
  frameP->GetXaxis()->SetTitle("m(D^{0}) [MeV]");
  frameP->SetMinimum(-5);
  frameP->SetMaximum(5);
  frameP->Draw();
  TLine *line_plus3 = new TLine(min, 3, max, 3);
  TLine *line_minus3 = new TLine(min, -3, max, -3);
  line_plus3->SetLineStyle(2);
  line_minus3->SetLineStyle(2);
  line_plus3->SetLineColor(kGray + 2);
  line_minus3->SetLineColor(kGray + 2);
  line_plus3->Draw("same");
  line_minus3->Draw("same");

  c->SaveAs("FitData.png");
}
