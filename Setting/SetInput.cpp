// THE AIM OF THIS CODE IS TO SET THE INPUT ROOT FILES FOR THE SIMULATION
// According to the chosen brem value, it will overwrite/recreate the files
// containing the smearing options for daughters and the D*+ distribution for
// the mother, in the rapidsim directory used for the simulation.

// function declarations
void SetInputFile(int brem = -1, TString typeOfParticle = "Electron");
void CopyFunction(TFile *input, TFile *output);

void SetInput(int brem = -1)
{

    // Check on the input value
    if (brem == 0 || brem == 1)
    {

        // call the functions to set the input files for the simulation
        SetInputFile(brem, "Hadron");
        SetInputFile(brem, "Electron");
        SetInputFile(brem, "Dst");
    }
    else
    {
        std::cerr << "Invalid brem value. Please use 0 or 1." << '\n';
    }
}

void SetInputFile(int brem = - 1, TString typeOfParticle = "Electron")
{
    TString name, fileRoot, motherDir;
    // Set paths
    if (typeOfParticle == "Hadron")
    {
        name = "HadronBrem";
        fileRoot = "Run3_HadronSmearing.root";
        motherDir = "smear/";
    }
    else if (typeOfParticle == "Electron")
    {
        name = "ElectronBrem";
        fileRoot = "Run3_ElectronSmearing.root";
        motherDir = "smear/";
    }
    else if (typeOfParticle == "Dst")
    {
        name = "DstBrem";
        fileRoot = "LHCc14.root";
        motherDir = "fonll/";
    }else{
        std::cerr<<"Not viable particle\n";
        return;
    }
    // Open the input file in the correct brem directory
    TString path = name + brem + "/" + fileRoot;
    std::cout <<"Path: "<<path<<"\n";
    TFile *inputFile = TFile::Open(path, "READ");

    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error! \n";
    }

    std::cout << "File Input: "<< name << brem <<"/"<<fileRoot<< " opened!" << '\n';

    // Create the output file in the rapidsim directory
    TFile *outFile =
        new TFile("/opt/RapidSim/rootfiles/" + motherDir + fileRoot, "RECREATE");

    // Copy
    CopyFunction(inputFile, outFile);

    std::cout << "File Output: "<< fileRoot <<" saved in " << motherDir <<"configuration!" << '\n';
    std::cout << '\n';
}

// Copy all objects from input ROOT file to output ROOT file
void CopyFunction(TFile *input, TFile *output)
{

    TIter next(input->GetListOfKeys());
    TKey *key;

    while ((key = (TKey *)next()))
    {
        TObject *obj = key->ReadObj();
        output->cd();
        std::cout << "Copying " << obj->GetName() << "..." << '\n';
        obj->Write(obj->GetName());
    }

    input->Close();
    output->Close();
}

