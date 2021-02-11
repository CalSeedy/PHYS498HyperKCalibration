#include <iostream>
#include <dirent.h>
#include "./include/WCSimRootEvent.hh"
#include "./include/WCSimRootGeom.hh"
#include "./include/WCSimRootOptions.hh"

int test() {
	bool verbose = true;
	bool hybrid = false;
    std::vector<std::string> rootfiles;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir("./RootFiles/")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            std::string name(ent->d_name);
            printf("%s\n", name.c_str());
            if (name.find("_flat") == std::string::npos) {
                rootfiles.push_back("./RootFiles/" + name);
            }
        }
        closedir(dir);
    }
    else {
        /* could not open directory */
        perror("");
        return EXIT_FAILURE;
    }
    
    rootfiles.pop_back();
    rootfiles.pop_back();

    for (int idx = 0; idx < rootfiles.size(); idx++) {
        std::string filename = rootfiles.at(idx);
        std::cout << "Analysing Root file: " << filename << std::endl;
        std::string fn(filename);
        fn = fn.substr(0, fn.find(".root"));

		TFile* file = new TFile(filename.c_str(), "read");
		if (!file->IsOpen()) {
			std::cout << "Error, could not open input file: " << filename.c_str() << std::endl;
			return -1;
		}

		// Get the a pointer to the tree from the file
		TTree* tree = (TTree*)file->Get("wcsimT");
		tree->GetCurrentFile();
		// Get the number of events
		int nevent = tree->GetEntries();
		if (verbose) printf("nevent %d\n", nevent);

		// Create a WCSimRootEvent to put stuff from the tree in
		WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
		WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();
		// Set the branch address for reading from the tree
		TBranch* branch = tree->GetBranch("wcsimrootevent");
		branch->SetAddress(&wcsimrootsuperevent);
		// Force deletion to prevent memory leak 
		tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

		TBranch* branch2;
		if (hybrid) {
			branch2 = tree->GetBranch("wcsimrootevent2");
			branch2->SetAddress(&wcsimrootsuperevent2);
			// Force deletion to prevent memory leak 
			tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
		}

		// Options tree - only need 1 "event"
		TTree* opttree = (TTree*)file->Get("wcsimRootOptionsT");
		WCSimRootOptions* opt = 0;
		opttree->SetBranchAddress("wcsimrootoptions", &opt);
		if (verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
		if (opttree->GetEntries() == 0) {
			exit(9);
		}
		opttree->GetEntry(0);
		opt->Print();

		// Geometry tree - only need 1 "event"
		TTree* geotree = (TTree*)file->Get("wcsimGeoT");
		WCSimRootGeom* geo = 0;
		geotree->SetBranchAddress("wcsimrootgeom", &geo);
		if (verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
		if (geotree->GetEntries() == 0) {
			exit(9);
		}
		geotree->GetEntry(0);
		// start with the main "subevent", as it contains most of the info
		// and always exists.
		int num_PMTs = 0;
		int num_mPMTs = 0;
		num_PMTs = geo->GetWCNumPMT(0);
		if (hybrid) num_mPMTs = geo->GetWCNumPMT(1);
		else num_mPMTs = 0;
		int numPMTs = num_mPMTs + num_PMTs;
		std::cout << "Detected" << numPMTs << " (" << num_mPMTs << ", " << num_PMTs << ") PMTs!" << std::endl;



		/*
        std::string out1 = fn.substr(12,fn.length()) + (std::string)"-position.png";
        std::string out2 = fn.substr(12,fn.length()) + (std::string)"-verts.png";
        std::string out3 = fn.substr(12,fn.length()) + (std::string)"-main.png";


        std::cout << out1 << std::endl;
        std::cout << out2 << std::endl;
        std::cout << out3 << std::endl;
		*/
    }


	return 0;
}