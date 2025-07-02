// THE AIM OF THIS CODE IS TO SET THE INPUT ROOT FILES FOR THE SIMULATION
// According to the chosen brem value, it will overwrite/recreate the files
// containing the smearing options for daughters and the D*+ distribution for
// the mother, in the rapidsim directory used for the simulation.

// function declarations
void Hadron(int brem = -1);
void Electron(int brem = -1);
void Dst(int brem = -1);
void CopyFunction(TFile *input, TFile *output);

void SetInput(int brem = -1) {

  // Check on the input value
  if (brem == 0 || brem == 1) {

    // call the functions to set the input files for the simulation
    Hadron(brem);
    Electron(brem);
    Dst(brem);

  } else {
    std::cerr << "Invalid brem value. Please use 0 or 1." << '\n';
  }
}

// Hadron function
// This function is used to set the input file for the hadron smearing
// configuration
void Hadron(int brem = -1) {

  // Open the input file for hadrons in the correct brem directory
  TFile *inputFileHad = TFile::Open(
      ("HadronBrem" + std::to_string(brem) + "/Run3_HadronSmearing.root")
          .c_str(),
      "READ");

  // check if the file was opened successfully
  if (!inputFileHad || inputFileHad->IsZombie()) {
    std::cerr << "Error! \n";
  }
  std::cout << "File Input: HadronBrem" << brem
            << "/Run3_HadronSmearing.root opened!" << '\n';

  // Create the output file for hadrons in the rapidsim directory
  TFile *outFileHad = new TFile(
      "/opt/RapidSim/rootfiles/smear/Run3_HadronSmearing.root", "RECREATE");

  // Copy the contents of the input file to the output file with the sama name
  CopyFunction(inputFileHad, outFileHad);

  std::cout << "File Output: Run3_HadronSmearing.root saved in smearing "
               "configuration!"
            << '\n';
  std::cout << '\n';
}

// Electron function
// This function is used to set the input file for the electron smearing
// configuration
void Electron(int brem = -1) {

  // Open the input file for electrons in the correct brem directory
  TFile *inputFileEle = TFile::Open(
      ("ElectronBrem" + std::to_string(brem) + "/Run3_ElectronSmearing.root")
          .c_str(),
      "READ");

  if (!inputFileEle || inputFileEle->IsZombie()) {
    std::cerr << "Error! \n";
  }

  std::cout << "File Input:"
            << "ElectronBrem" << brem << "/Run3_ElectronSmearing.root opened!"
            << '\n';

  // Create the output file for electrons in the rapidsim directory
  TFile *outFileEle = new TFile(
      "/opt/RapidSim/rootfiles/smear/Run3_ElectronSmearing.root", "RECREATE");

  // Copy
  CopyFunction(inputFileEle, outFileEle);

  std::cout << "File Output: Run3_ElectronSmearing.root saved in smearing "
               "configuration!"
            << '\n';
  std::cout << '\n';
}

// Dst function
// This function is used to set the input file for the D*+ distribution
// configuration
void Dst(int brem = -1) {

  // Open the input file for D*+ distribution in the correct brem directory
  TFile *inputFileDst = TFile::Open(
      ("DstBrem" + std::to_string(brem) + "/LHCc14.root").c_str(), "READ");

  if (!inputFileDst || inputFileDst->IsZombie()) {
    std::cerr << "Error! \n";
  }

  std::cout << "File Input: DstBrem" << brem << "/LHCc14.root opened!" << '\n';

  // Create the output file for D*+ distribution in the rapidsim directory
  TFile *outFileDst =
      new TFile("/opt/RapidSim/rootfiles/fonll/LHCc14.root", "RECREATE");

  // Copy
  CopyFunction(inputFileDst, outFileDst);

  std::cout << "File Output: LHCc14.root saved in fonll configuration!" << '\n';
  std::cout << '\n';
}

// Copy all objects from input ROOT file to output ROOT file
void CopyFunction(TFile *input, TFile *output) {

  TIter next(input->GetListOfKeys());
  TKey *key;

  while ((key = (TKey *)next())) {
    TObject *obj = key->ReadObj();
    output->cd();
    std::cout << "Copying " << obj->GetName() << "..." << '\n';
    obj->Write(obj->GetName());
  }

  input->Close();
  output->Close();
}