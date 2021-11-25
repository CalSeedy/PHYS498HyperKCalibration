#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <THStack.h>
#include <TMath.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

WCSimRootGeom* geo = 0;

bool eventDisplay_closeWall = false;
int nPMTpermPMT = 19;

// Simple example of reading a generated Root file
const int nPMTtypes = 2;
double PMTradius[nPMTtypes];

bool is_hybrid(std::string name) {
	int found = -1;
	int found2 = -1;

	found = name.find("mPMT");
	found2 = name.find("hybrid");
	if ((found != std::string::npos) || (found2 != std::string::npos)) return true;
	else return false;
}


struct Pos
{
	std::vector<double> data;

	Pos() {
		for (int i = 0; i < 3; i++) data.push_back(0.0);
	}

	Pos(std::vector<double> values) {
		for (int i = 0; i < values.size(); i++) {
			data.push_back(values.at(i));
		}
	}

	Pos(const Pos& other) {
		for (int i = 0; i < other.data.size(); i++) {
			data.push_back(other.data.at(i));
		}
	}


	double Get(int ind) {
		if ((ind < data.size()) && (ind >= 0)) {
			return data.at(ind);
		}
		else {
			return std::nan("1");
		}
	}

	void Print() {
		std::cout << "Pos: {";
		for (int i = 0; i < data.size(); i++) {
			std::cout << data.at(i) << ", ";
		}
		std::cout << "}" << std::endl;
	}

	double distanceTo(Pos other) {
		return sqrt((data[0] - other.Get(0))* (data[0] - other.Get(0)) + (data[1] - other.Get(1))*(data[1] - other.Get(1)) + (data[2] - other.Get(2))*(data[2] - other.Get(2)));
	}


	void add(double x) {
		for (int i = 0; i < data.size(); i++) data.at(i) += x;
	}

	void divide(double x) {
		for (int i = 0; i < data.size(); i++) data.at(i) /= x;
	}

	void divide(int x) {
		for (int i = 0; i < data.size(); i++) data.at(i) /= x;
	}

	void add(Pos p) {
		for (int i = 0; i < data.size(); i++) data.at(i) += p.Get(i); 
	}

	friend Pos operator+(const Pos& p1, const Pos& p2) {
		Pos out(p1);
		out.add(p2);
		return out;
	}

	Pos& operator+=(const Pos& p) {
		this->add(p);
		return *this;
	}

	friend Pos operator/(const Pos& p1, const double& x) {
		Pos out(p1);
		out.divide(x);
		return out;
	}

	friend Pos operator/(const Pos& p1, const int& x) {
		Pos out(p1);
		out.divide(x);
		return out;
	}

	Pos& operator/=(const double& c) {
		this->divide(c);
		return *this;
	}
};

std::vector<std::string> getFiles() {
	std::vector<std::string> out;
	// this causes TFile File does not exist error?
	/*
	std::ifstream infile("./rootfiles.txt");
	std::string line;
	while (std::getline(infile, line)) {
		if (line == "") continue;

		out.push_back(line);
	}
	for (auto f : out) std::cout << f << std::endl;
	*/

	////////////////////////////////////////////////////////////////////////////
	// WCSim Hybrid Test Files
	out.push_back("RootFiles/wcsim_hybrid\(10MeVe-\).root");
	out.push_back("RootFiles/wcsim_hybrid\(500MeVmu\).root");
	out.push_back("RootFiles/wcsim_hybrid\(500MeVe-\).root");
	out.push_back("RootFiles/wcsim_hybrid\(10GeVe-\).root");
	out.push_back("RootFiles/wcsim_hybrid\(10GeVmu\).root");
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// atmospheric muons and electrons 
	//out.push_back("/hepstore/pmenendez/HK/atm_emu/wcsim_atm_emu_20kBL_43.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_emu/wcsim_atm_emu_40kBL_1011.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_emu/wcsim_atm_emu_3kmPMT_1000.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_emu/wcsim_atm_emu_5kmPMT_666.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_emu/wcsim_atm_emu_10kmPMT_888.root");
	
	// atmospheric taus
	//out.push_back("/hepstore/pmenendez/HK/atm_tau/wcsim_atm_tau_20kBL_1000.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_tau/wcsim_atm_tau_40kBL_1000.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_tau/wcsim_atm_tau_3kmPMT_1000.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_tau/wcsim_atm_tau_5kmPMT_1000.root");
	//out.push_back("/hepstore/pmenendez/HK/atm_tau/wcsim_atm_tau_10kmPMT_1000.root");
	
	//Renium Source, simulate solar neutrinos
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingA/wcsim_Rn_scalingA_20kBL_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingA/wcsim_Rn_scalingA_40kBL_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingA/wcsim_Rn_scalingA_5kmPMT_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingA/wcsim_Rn_scalingA_10kmPMT_0.root");
	//
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingB/wcsim_Rn_scalingA_20kBL_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingB/wcsim_Rn_scalingA_40kBL_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingB/wcsim_Rn_scalingA_3kmPMT_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingB/wcsim_Rn_scalingA_5kmPMT_0.root");
	//out.push_back("/hepstore/pmenendez/HK/Rn/scalingB/wcsim_Rn_scalingA_10kmPMT_0.root");
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// HK Monte Carlo
	//out.push_back("HK/WCSim_hk_20bl_4200hz_r14374_0hz_435nmCol_middle.root");
	//out.push_back("HK/WCSim_hkhybridmpmt10pc_20bl_3kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle.root");
	//out.push_back("HK/WCSim_hkhybridmpmt10pc_20bl_5kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle.root");
	//out.push_back("HK/WCSim_hkhybridmpmt10pc_20bl_10kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle.root");
	////////////////////////////////////////////////////////////////////////////

	return out;
}


double bias_function(double charge) {
	if (charge < 1000) {
		return 0.990274 * pow(charge, 0) + 0.0010733 * pow(charge, 1) - 4.43408e-06 * pow(charge, 2) + 4.63877e-09 * pow(charge, 3) + 1.50912e-12 * pow(charge, 4) - 5.58521e-15 * pow(charge, 5) + 2.61881e-18 * pow(charge, 6);
	}
	else if (charge < 15000) {
		return 1.11385 * pow(charge, 0) - 0.000388516 * pow(charge, 1) + 9.9336e-08 * pow(charge, 2) - 1.54218e-11 * pow(charge, 3) + 1.37418e-15 * pow(charge, 4) - 6.425e-20 * pow(charge, 5) + 1.21684e-24 * pow(charge, 6);
	}
	else return 0.22682;
}


int main(int argc, char** argv) {
	char* infilename = NULL;
	char* outfilename = "";
	bool verbose = false;//false;
	bool hybrid = true;
	bool bias = false;
	bool ratios = false;
	int nPMTtypes = 2;
	int c = -1;
	int evLimit = 0;
	while ((c = getopt(argc, argv, ":f:o:n:vbr")) != -1) {//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
		switch (c) {
		case 'f':
			infilename = optarg;
			break;
		case 'o':
			outfilename = optarg;
			break;
		case 'v':
			verbose = true;
			break;
		case 'b':
			bias = true;
			break;
		case 'n':
			evLimit = std::stoi(optarg);
			break;
		case 'r':
			ratios = true;
			break;
		}
	}
	auto start = std::chrono::system_clock::now();
	std::vector<std::string> files;
	// Open the file
	if (infilename != NULL) {
		files.push_back((std::string)infilename);
	}
	else {
		files = getFiles();
	}

	TCanvas *c1, *c2, *c3, *c4, *c5, *c6, *c7, *c8, *c9;
	TFile *file, *outfile;
	TH1D *h1, *pmtQ, *dists, *pmtPosX, *pmtPosY, *pmtPosZ, *vertX, *vertY, *vertZ, *mh1, *mpmtQ, *mpmtPosX, *mpmtPosY, *mpmtPosZ, *mvertX, *mvertY, *mvertZ;

	std::vector<std::string> consoletxt;
	for (std::string filename : files) {
		std::cout << "Analysing " << filename << std::endl;
		hybrid = is_hybrid(filename);
		if (verbose) std::cout << "Is hybrid? " << hybrid << std::endl;
		file = new TFile(filename.c_str(), "READ");

		if (!file->IsOpen()) {
			std::cout << "Error, could not open input file: " << filename << std::endl;
			return -1;
		}

		// Get the a pointer to the tree from the file
		TTree* tree = (TTree*)file->Get("wcsimT");
		if (!tree) { 
			std::cerr << "TTree wcsimT missing\n";
			exit(1); 
		}
		// Get the number of events
		int nevent = ((int)tree->GetEntries());//std::min(((int)tree->GetEntries()),100000);
		if (verbose) printf("nevent %d\n", nevent);
		if (evLimit == 0) evLimit = nevent;
		//bool test = filename.find("HK/");
		//if (test != std::string::npos) evLimit = 100;
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

		// Geometry tree - only need 1 "event"
		TTree* geotree = (TTree*)file->Get("wcsimGeoT");
		geotree->SetBranchAddress("wcsimrootgeom", &geo);
		if (verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
		if (geotree->GetEntries() == 0) {
			exit(9);
		}
		geotree->GetEntry(0);
		PMTradius[0] = geo->GetWCPMTRadius();
		PMTradius[1] = geo->GetWCPMTRadius(true);
		int nPMT = geo->GetWCNumPMT();
		int nmPMT = geo->GetWCNumPMT(true);
		std::cout << "Number of PMTs of 1st type = " << nPMT << ", radius = " << PMTradius[0] << std::endl;
		std::cout << "Number of PMTs of 2nd type = " << nmPMT << ", radius = " << PMTradius[1] << std::endl;

		// Options tree - only need 1 "event"
		TTree* opttree = (TTree*)file->Get("wcsimRootOptionsT");
		WCSimRootOptions* opt = 0;
		opttree->SetBranchAddress("wcsimrootoptions", &opt);
		if (verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
		if (opttree->GetEntries() == 0) {
			exit(9);
		}
		opt->Print();
		//opttree->GetEntry(0);
		if (verbose) std::cout << "Starting analysis..." << std::endl;
		// start with the main "subevent", as it contains most of the info
		// and always exists.
		WCSimRootTrigger* wcsimrootevent;
		WCSimRootTrigger* wcsimrootevent2;

		int num_trig = 0;

		//TFile * fOutput = new TFile(OutputFile,"recreate");

		h1 = new TH1D("PMT Hits", "Number of PMT Hits", 30, 0, 30);
		pmtQ = new TH1D("PMT Charge", "Hit PMT Charge", 50, 0, 20);
		dists = new TH1D("PMT Distances (<= 10m)", "Hit PMT Distances (<= 10m)", 100, 0, 1000);
		pmtPosX = new TH1D("X Pos", "X Pos for hit PMTs", 200, -4000, 4000);
		pmtPosY = new TH1D("Y Pos", "Y Pos for hit PMTs", 200, -4000, 4000);
		pmtPosZ = new TH1D("Z Pos", "Z Pos for hit PMTs", 200, -4000, 4000);
		vertX = new TH1D("Vert X", "Event vertex X positions", 100, -5000, 5000);
		vertY = new TH1D("Vert Y", "Event vertex Y positions", 100, -5000, 5000);
		vertZ = new TH1D("Vert Z", "Event vertex Z positions", 100, -5000, 5000);

		mh1 = new TH1D("mPMT Hits", "Number of mPMT Hits", 30, 0, 30);
		mpmtQ = new TH1D("mPMT Charge", "Hit mPMT Charge", 50, 0, 20);
		mpmtPosX = new TH1D("mPMT X Pos", "X Pos for hit mPMTs", 200, -4000, 4000);
		mpmtPosY = new TH1D("mPMT Y Pos", "Y Pos for hit mPMTs", 200, -4000, 4000);
		mpmtPosZ = new TH1D("mPMT Z Pos", "Z Pos for hit mPMTs", 200, -4000, 4000);
		mvertX = new TH1D("mVert X", "Event vertex X positions", 100, -5000, 5000);
		mvertY = new TH1D("mVert Y", "Event vertex Y positions", 100, -5000, 5000);
		mvertZ = new TH1D("mVert Z", "Event vertex Z positions", 100, -5000, 5000);

		std::vector<Pos> positions;
		std::vector<Pos> mPositions;
		std::vector<int> IDs;
		std::vector<int> mIDs;
		std::vector<double> Qdata(nPMT);
		if ((nmPMT % 19 != 0) && (hybrid)) { std::cout << "Leftover PMTs in mPMTs, should be a factor of 19." << std::endl; exit(-1); }
		//std::vector<Pos> mPMTAvgPos;
		if (nmPMT == 0) nmPMT = 1;
		std::vector<double> mQdata(nmPMT);
		if (hybrid) {
			mQdata.resize(nmPMT / 19);
		}
		//assume PMTs in mPMTs are successively made in batches of 19, so 19 continuous indexes make up the mPMT.
		/*
			for (int N = 0; N < nmPMT / 19; N++) {
				std::vector<Pos> mPMTPositions;
				for (int i = 0; i < 19; i++) {
					WCSimRootPMT pmt = geo->GetPMT(19 * N + i, true);
					std::vector<double> pos = {
					pmt.GetPosition(0), pmt.GetPosition(1), pmt.GetPosition(2)
					};
					Pos p(pos);
					mPMTPositions.push_back(p);
				}

				Pos avg;
				for (Pos p : mPMTPositions) {
					avg += p;
				}
				avg /= 19;
				if (verbose) {
					std::cout << "***mPMT#" << N << "Pos Avg.***" << std::endl;
					avg.Print();
				}
				mPMTAvgPos.push_back(avg);
			}
		}
		*/

		std::cout << filename << std::endl;
		int a = filename.rfind("/"); // path/to/folder/file.ext    -- find last slash (/)
		int b = filename.substr(0, a - 1).rfind("/"); // find second to last slash (/)
		std::string fn1 = filename.substr(b + 1, filename.rfind(".") - b - 1); // attempt to grab "file"
		if (fn1.substr(0, 9) == "RootFiles") fn1 = (std::string)"test" + filename.substr(a, filename.rfind(".") - a);
		std::cout << fn1 << std::endl;
		int slash = fn1.rfind("/");
		std::string name = fn1.substr(slash + 1, fn1.size() - slash);
		if (bias) fn1.append("_bias");
		if (bias) name.append("_bias");
		double q1, q2, minX, minY, maxX, maxY, ratio, minR, maxR;
		TFile* dataOut = new TFile((name + (std::string)"-chargeData.root").c_str(), "RECREATE");
		TTree* dataStore = new TTree("charge data", "all charge data");
		dataStore->Branch("q1", &q1, "q1/D");
		dataStore->Branch("q2", &q2, "q2/D");
		if (ratios) dataStore->Branch("ratio", &ratio, "ratio/D");
		dataStore->Branch("hybrid", &hybrid, "hybrid/O");

		double tLimit = 9000.0; // limit lifetime to 9000 ns == 9 um
		std::vector<int> events;
		while (events.size() < evLimit) {
			int ev = std::rand() % nevent;
			while (std::find(events.begin(), events.end(), ev) != events.end()) {
				ev = std::rand() % nevent;
			}
			events.push_back(ev);
			events.shrink_to_fit();
		}
		std::vector<double> filterQ1, filterQ2;
		// Now loop over events
		for (int ev : events)
		{		
			// Read the event from the tree into the WCSimRootEvent instance
			tree->GetEntry(ev);
			if (verbose) std::cout << "Got tree entry" << std::endl;


			wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
			if (hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
			if (verbose) {
				printf("********************************************************");
				printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
					wcsimrootevent->GetHeader()->GetDate());
				printf("Mode %d\n", wcsimrootevent->GetMode());
				printf("Number of subevents %d\n",
					wcsimrootsuperevent->GetNumberOfSubEvents());

				printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
				printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
					wcsimrootevent->GetVtx(1), wcsimrootevent->GetVtx(2));
			}

			if (verbose) {
				printf("Jmu %d\n", wcsimrootevent->GetJmu());
				printf("Npar %d\n", wcsimrootevent->GetNpar());
				printf("Ntrack %d\n", wcsimrootevent->GetNtrack());

			}

			vertX->Fill(wcsimrootevent->GetVtx(0));
			vertY->Fill(wcsimrootevent->GetVtx(1));
			vertZ->Fill(wcsimrootevent->GetVtx(2));

			if (hybrid) {
				mvertX->Fill(wcsimrootevent2->GetVtx(0));
				mvertY->Fill(wcsimrootevent2->GetVtx(1));
				mvertZ->Fill(wcsimrootevent2->GetVtx(2));
			}

			std::vector<float> triggerInfo;
			triggerInfo.clear();
			triggerInfo = wcsimrootevent->GetTriggerInfo();

			std::vector<float> triggerInfo2;
			triggerInfo2.clear();
			if (hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();

			if (verbose) {
				for (int v = 0; v < triggerInfo.size(); v++) {
					std::cout << "Trigger entry #" << v << ", info = " << triggerInfo[v] << std::endl;
				}
				if (hybrid) {
					for (int v = 0; v < triggerInfo2.size(); v++) {
						std::cout << "Trigger2 entry #" << v << ", info = " << triggerInfo2[v] << std::endl;
					}
				}
			}

			int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
			int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
			int ncherenkovhits2 = 0; if (hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();
			int ncherenkovdigihits2 = 0; if (hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits();

			if (verbose) {
				printf("node id: %i\n", ev);
				printf("Ncherenkovhits %d\n", ncherenkovhits);
				printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
				printf("Ncherenkovhits2 %d\n", ncherenkovhits2);
				printf("Ncherenkovdigihits2 %d\n", ncherenkovdigihits2);
				std::cout << "RAW HITS:" << std::endl;
			}

			// Loop through elements in the TClonesArray of WCSimRootCherenkovHits
			for (int pmtType = 0; pmtType < nPMTtypes; pmtType++) {
				if (pmtType == 1 && !hybrid) continue;
				if (verbose) std::cout << "PMT Type = " << pmtType << std::endl;

				TClonesArray* timeArray;//An array of pointers on CherenkovHitsTimes.
				if (pmtType == 0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
				else timeArray = wcsimrootevent2->GetCherenkovHitTimes();

				double totalPe = 0;
				int totalHit = 0;

				int nhits;
				if (pmtType == 0) nhits = ncherenkovhits;
				else nhits = ncherenkovhits2;

				WCSimRootCherenkovHit* wcsimrootcherenkovhit;
				for (int i = 0; i < nhits; i++)
				{
					if (verbose) std::cout << "Hit #" << i << std::endl;

					if (pmtType == 0) wcsimrootcherenkovhit = (wcsimrootevent->GetCherenkovHits())->At(i);
					else wcsimrootcherenkovhit = (wcsimrootevent2->GetCherenkovHits())->At(i);

					int tubeNumber = wcsimrootcherenkovhit->GetTubeID();
					int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
					int peForTube = wcsimrootcherenkovhit->GetTotalPe(1);
					if (verbose) std::cout << "(" << pmtType << ", " << i << ") Got, tube: #" << tubeNumber << "| pe: " << peForTube << std::endl;
					if (pmtType == 0) h1->Fill(peForTube);
					else mh1->Fill(peForTube);

					WCSimRootPMT pmt;
					if (pmtType == 0) pmt = geo->GetPMT(tubeNumber - 1, false);
					else pmt = geo->GetPMT(tubeNumber - 1, true);

					std::vector<double> PMTpos{ 0.0, 0.0, 0.0 };
					for (int j = 0; j < 3; j++) {
						PMTpos.at(j) = pmt.GetPosition(j);
					}

					Pos pos = Pos(PMTpos);
					if (pmtType == 0) {
						positions.push_back(pos);
						IDs.push_back(tubeNumber - 1);
					} else {
						int idx = tubeNumber - 1;
						mPositions.push_back(pos);
						mIDs.push_back(idx);
					}


					if (verbose) std::cout << "Added Pos and ID" << std::endl;

					if (pmtType == 0) {
						pmtPosX->Fill(PMTpos.at(0)); pmtPosY->Fill(PMTpos.at(1)); pmtPosZ->Fill(PMTpos.at(2));
					}
					else {
						mpmtPosX->Fill(PMTpos.at(0)); mpmtPosY->Fill(PMTpos.at(1)); mpmtPosZ->Fill(PMTpos.at(2));
					}

					if (verbose) {
						printf("Total pe: %d times (", peForTube);

						for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
						{
							WCSimRootCherenkovHitTime* HitTime = (WCSimRootCherenkovHitTime*)timeArray->At(j);
							if (verbose) printf("%6.2f ", HitTime->GetTruetime());
						}
						std::cout << ")" << std::endl;
					}

					if (verbose) std::cout << "Position of PMT = " << PMTpos.at(0) << ", " << PMTpos.at(1) << ", " << PMTpos.at(2) << std::endl;

					totalPe += peForTube;
					totalHit++;				
				}
				if (verbose) std::cout << "Total Pe : " << totalPe << std::endl;
			}// End of loop over Cherenkov hits

			// Get the number of digitized hits
			// Loop over sub events
			if (verbose) std::cout << "DIGITIZED HITS:" << std::endl;
			double totalQ = 0.0;
			for (int pmtType = 0; pmtType < nPMTtypes; pmtType++) {
				if (verbose) std::cout << "PMT Type = " << pmtType << std::endl;
				// Grab the big arrays of times and parent IDs
				TClonesArray* timeArray;
				if (pmtType == 0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
				else timeArray = wcsimrootevent2->GetCherenkovHitTimes();

				double totalPe = 0;
				int totalHit = 0;

				int nhits = 0;
				if (pmtType == 0) nhits = ncherenkovdigihits;
				else nhits = ncherenkovdigihits2;

				for (int i = 0; i < nhits; i++)
				{
					TObject* Hit;
					if (pmtType == 0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
					else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

					WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

					int tubeNumber = wcsimrootcherenkovdigihit->GetTubeId();
					double peForTube = wcsimrootcherenkovdigihit->GetQ();
					double t = wcsimrootcherenkovdigihit->GetT();

					if (t > tLimit) continue;
					WCSimRootPMT pmt;
					if (pmtType == 0) {
						pmt = geo->GetPMT(tubeNumber - 1, false);
						pmtQ->Fill(peForTube);
						Qdata[tubeNumber - 1] += peForTube;
					}
					else {
						pmt = geo->GetPMT(tubeNumber - 1, true);
						mpmtQ->Fill(peForTube);
						int idx = (int)std::floor((tubeNumber - 1) / 19);
						//if (((tubeNumber - 1) % 19) != 0) std::cout << "Tube ID (" << tubeNumber - 1 << ") isn't a multiple of 19!" << std::endl;
						mQdata[idx] += peForTube;
					}


					totalQ += peForTube;
					double PMTpos[3] = { pmt.GetPosition(0), pmt.GetPosition(1), pmt.GetPosition(2)};
					if (verbose) {
						printf("Total pe: %1.1f times( ", peForTube);
						std::cout << ")" << std::endl;
						std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << std::endl;
					}
				} // End of loop over Cherenkov digitised hits
				if (verbose) std::cout << "Total Pe : " << totalPe << ", total hit : " << totalHit << std::endl;
			}

			for (int i = 0; i < positions.size(); i++) { // for every standard PMT
				int id1 = IDs.at(i);
				if (hybrid) {
					for (int j = 0; j < mPositions.size(); j++) { //  for every mPMT
						try {
							int id2 = (int)std::floor(mIDs.at(j) / 19);
							if ((Qdata.at(id1) == 0.0) || (mQdata.at(id2) == 0.0)) continue;
							Pos pos1 = positions.at(i);

							Pos pos2 = mPositions.at(j);

							if (verbose) {
								pos1.Print();
								pos2.Print();
							}

							double dist = pos1.distanceTo(pos2);

							if (dist <= 1000) { // within 10 m
								dists->Fill(dist);

								if (dist <= 150) { // within 1.5 m
									q1 = Qdata.at(id1);
									q2 = mQdata.at(id2);
									if (bias) q1 *= bias_function(q1);
									if (ratios) ratio = q1 / q2;
									//filterQ1.push_back(q1 * mult);
									//filterQ2.push_back(q2);
									dataStore->Fill();
								}
							}

							if (q1 < minX) minX = q1;
							if (q1 > maxX) maxX = q1;
							if (q2 < minY) minY = q2;
							if (q2 > maxY) maxY = q2;
							if (ratios) {
								if (ratio < minR) minR = ratio;
								if (ratio > maxR) maxR = ratio;
							}

							if ((i % 10 == 0) && (j == mPositions.size() - 1 ) && (j != 0)) {
								std::system("clear");
								auto now = std::chrono::system_clock::now();
								for (auto st : consoletxt) std::cout << st << std::endl;
								std::cout << "Event: " << ev << std::endl;
								std::cout << filename << "| (" << i << "/ " << positions.size() << ") | Time Running: " << std::chrono::duration<double>(now - start).count() << " s" << std::endl;
							}
						}
						catch (std::out_of_range& err) {
							std::cout << i << ", " << j << " (" << IDs.at(i) << ", " << mIDs.at(j) << ")" << std::endl;
							return 1;
						}
					}
				} else {
					for (int j = i; j < positions.size(); j++) { //  for every other PMT
						try {
							int id2 = IDs.at(j);
							if (id1 != id2) {
								if ((Qdata.at(id1) == 0.0) || (Qdata.at(id2) == 0.0)) continue;
								Pos pos1 = positions.at(i);

								Pos pos2 = positions.at(j);

								if (verbose) {
									pos1.Print();
									pos2.Print();
								}

								double dist = pos1.distanceTo(pos2);

								if (dist <= 1000) { // within 10 m
									dists->Fill(dist);

									if (dist <= 150) { // within 1.5 m
										q1 = Qdata.at(id1);
										q2 = Qdata.at(id2);
										if (bias) q1 *= bias_function(q1);
										//filterQ1.push_back(q1 * mult);
										if (bias) q2 *= bias_function(q2);
										//filterQ2.push_back(q2 * mult);
										if (ratios) ratio = q1 / q2;
										dataStore->Fill();
									}
								}
								if (q1 < minX) minX = q1;
								if (q1 > maxX) maxX = q1;
								if (q2 < minY) minY = q2;
								if (q2 > maxY) maxY = q2;
								if (ratios) {
									if (ratio < minR) minR = ratio;
									if (ratio > maxR) maxR = ratio;
								}
								if ((i % 10 == 0) && (j == positions.size() - 1) && (j != 0)) {
									std::system("clear");
									auto now = std::chrono::system_clock::now();
									for (auto st : consoletxt) std::cout << st << std::endl;
									std::cout << "Event: " << ev << std::endl;
									std::cout << (std::string)filename << "| (" << i << "/ " << positions.size() << ") | Time Running: " << std::chrono::duration<double>(now - start).count() << " s" << std::endl;
								}
							}
						}
						catch (std::out_of_range& err) {
							std::cout << i << ", " << j << " (" << IDs.at(i) << ", " << IDs.at(j) << ")" << std::endl;
							return 1;
						}
					}
				}
			}
			IDs.clear(); positions.clear(); for (double x : Qdata) x = 0.0;
			if (hybrid) { mIDs.clear(); mPositions.clear(); for (double x : mQdata) x = 0.0;}
			dataStore->Write();
			// renitialize super event between loops.
			wcsimrootsuperevent->ReInitialize();
			if (hybrid) wcsimrootsuperevent2->ReInitialize();
		} // End of loop over events
		dataOut->Write();
		dataOut->Close();
		std::cout << "num_trig " << num_trig << "\n";

		std::cout << "Pos size: " << positions.size() << std::endl;
		auto now = std::chrono::system_clock::now();
		std::string s = (std::string)filename + "| Finished at: " + std::to_string(std::chrono::duration<double>(now - start).count()) + " s";
		consoletxt.push_back(s);
		//std::system("clear");
		for (auto st : consoletxt) std::cout << st << std::endl;

		int nHigh = 3;
		int nWide = 1;
		//std::cout << "Zeros (q1, q2): " << count(filterQ1.begin(), filterQ1.end(), 0) << ", " << count(filterQ2.begin(), filterQ2.end(), 0) << std::endl;
		c1 = new TCanvas("c1", "First canvas", nWide * 3840, nHigh * 2160);
		c1->Divide(nWide, nHigh);
		THStack* posXStack = new THStack("posXStack", "PMT X-position; PMT Position, x [cm]; # of PMTs");
		THStack* posYStack = new THStack("posYStack", "PMT Y-position; PMT Position, y [cm]; # of PMTs");
		THStack* posZStack = new THStack("posZStack", "PMT Z-position; PMT Position, z [cm]; # of PMTs");
		pmtPosX->SetFillColor(kBlue);
		pmtPosY->SetFillColor(kBlue);
		pmtPosZ->SetFillColor(kBlue);
		posXStack->Add(pmtPosX); 
		posYStack->Add(pmtPosY);
		posZStack->Add(pmtPosZ);
		if (hybrid) {
			mpmtPosX->SetFillColor(kRed); posXStack->Add(mpmtPosX);
			mpmtPosY->SetFillColor(kRed); posYStack->Add(mpmtPosY);;
			mpmtPosZ->SetFillColor(kRed); posZStack->Add(mpmtPosZ);;
		}
		c1->cd(1); posXStack->Draw();
		c1->cd(2); posYStack->Draw();
		c1->cd(3); posZStack->Draw();
		c1->Draw();

		THStack* vertXStack = new THStack("vertXStack", "Event vertex X-position; Vertex Position, x [cm]; Count");
		THStack* vertYStack = new THStack("vertYStack", "Event vertex Y-position; Vertex Position, y [cm]; Count");
		THStack* vertZStack = new THStack("vertZStack", "Event vertex Z-position; Vertex Position, z [cm]; Count");
		c2 = new TCanvas("c2", "Second canvas", nWide * 3840, nHigh * 2160);
		c2->Divide(nWide, nHigh);
		vertX->SetFillColor(kBlue); vertXStack->Add(vertX);
		if (hybrid) mvertX->SetFillColor(kRed); vertXStack->Add(mvertX);
		vertY->SetFillColor(kBlue); vertYStack->Add(vertY);
		if (hybrid) mvertY->SetFillColor(kRed); vertYStack->Add(mvertY);
		vertZ->SetFillColor(kBlue); vertZStack->Add(vertZ);
		if (hybrid) mvertZ->SetFillColor(kRed); vertZStack->Add(mvertZ);
		c2->cd(1); vertXStack->Draw();
		c2->cd(2); vertYStack->Draw();
		c2->cd(3); vertZStack->Draw();
		c2->Draw(); 

		THStack* hitStack = new THStack("hitStack", "Number of hits on (m)PMTs; # of Hits on (m)PMT; Count");
		THStack* chargeStack = new THStack("chargeStack", "Charge on a Hit (m)PMT; Charge on hit (m)PMT [pe]; Count");
		c3 = new TCanvas("c3", "Third canvas", nWide * 3840, nHigh * 2160);
		c3->Divide(nWide, nHigh);
		c3->cd(1); h1->SetFillColor(kBlue); hitStack->Add(h1);
		if (hybrid) mh1->SetFillColor(kRed); hitStack->Add(mh1);
		hitStack->Draw();
		c3->cd(2); pmtQ->SetFillColor(kBlue); chargeStack->Add(pmtQ);
		if (hybrid) mpmtQ->SetFillColor(kRed); chargeStack->Add(mpmtQ);
		chargeStack->Draw();
		c3->cd(3); dists->GetXaxis()->SetTitle("Distance between hit (m)PMTs [cm]"); dists->GetYaxis()->SetTitle("Count"); dists->SetFillColor(kBlue); dists->Draw();
		c3->Draw();

		//std::cout << "a: " << minX << ", " << maxX << "| b: " << minY << ", " << maxY << "| c: " << minR << ", " << maxR << std::endl;
		if (verbose) std::cout << "making charge plots" << std::endl;
		if (((maxX / minX > 1e100) && (minX >= 0.0)) || ((maxX / minX < -1e100) && (minX <= 0.0)) || minX < 1e-10 && minX >= 0.0) minX = 0.001; 
		if (((maxY / minY > 1e100) && (minY >= 0.0)) || ((maxY / minY < -1e100) && (minY <= 0.0)) || minY < 1e-10 && minY >= 0.0) minY = 0.001; 
		if (((maxR / minR > 1e100) && (minR >= 0.0)) || ((maxR / minR < -1e100) && (minR <= 0.0)) || minR < 1e-10 && minR >= 0.0) minR = 0.001; 
		//std::cout << "Max (q1, q2): " << maxX << ", " << maxY << "| Min (q1, q2): " << minX << ", " << minY << std::endl;
		/*
		pmtQvQ = new TH2D("pmtQvQ", ((std::string)"Q2 vs Q1 for hit (m)PMTs (" + name + (std::string)")").c_str(), nbins, minX, maxX, nbins, minY, maxY);
		nbins = 100;
		QvQProfile = new TProfile("QvQ Profile", ((std::string)"Q2 vs Q1 for hit (m)PMTs (" + name + (std::string)")").c_str(), nbins, minX, maxX, minY, maxY);

		for (int i = 0; i < filterQ1.size(); i++) {
			pmtQvQ->Fill(filterQ1.at(i), filterQ2.at(i));
			QvQProfile->Fill(filterQ1.at(i), filterQ2.at(i));
		}
		*/

		TFile* data = new TFile((name + (std::string)"-chargeData.root").c_str(), "READ");
		dataStore = (TTree*)data->Get("charge data");
		if (outfilename == "") outfilename = (fn1 + (std::string)"-analysis.root").c_str();
		std::cout << "File " << outfilename << " is open for writing" << std::endl;
		outfile = new TFile(outfilename, "RECREATE");
		int nbins = 300;
		double logMinX = log(minX);
		double logMaxX = log(maxX);
		double logMinY = log(minY);
		double logMaxY = log(maxY);
		double logMinR = log(minR);
		double logMaxR = log(maxR);
		if (logMaxX / logMinX > 1e100) logMinX = 0.0;
		if (logMaxY / logMinY > 1e100) logMinY = 0.0;
		if (logMaxR / logMinR > 1e100) logMinR = 0.0;
		c4 = new TCanvas("c4", "Fourth canvas", 3840 * 2, 2160);
		c4->Divide(2, 1);
		char* yaxis; 
		if (hybrid) yaxis = "Charge on hit mPMT, Q2";
		else yaxis = "Charge on other hit PMT, Q2";
		char* args;
		asprintf(&args, "q2:q1>>pmtQvQ(%i, 0.0, %g, %i, %g, %g)", nbins, std::min(maxX / 3.0 , 2000.0), nbins, minY, std::min(maxY / 3.0, 1500.0));
		c4->cd(1); 
		dataStore->Draw(args, "", "COLZ");
		TH2D* pmtQvQ = (TH2D*)gPad->GetPrimitive("pmtQvQ");
		pmtQvQ->SetContour(100);
		pmtQvQ->GetXaxis()->SetTitle("Charge on hit PMT, Q1");
		pmtQvQ->GetYaxis()->SetTitle(yaxis);
		c4->cd(2); 
		nbins = 100;
		delete args;
		asprintf(&args, "q2:q1>>QvQProfile(%i, 0.0, %g, %i, %g, %g)", nbins, std::min(maxX / 3.0, 2000.0), nbins, minY, std::min(maxY / 3.0, 1500.0));
		dataStore->Draw(args, "", "PROF");
		TProfile* QvQProfile = (TProfile*)gPad->GetPrimitive("QvQProfile");
		QvQProfile->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); QvQProfile->GetYaxis()->SetTitle(yaxis);
		c4->Draw();

		//std::cout << "Max (q1, q2): " << maxX << ", " << maxY << "| Min (q1, q2): " << minX << ", " << minY << std::endl;
		if (hybrid) yaxis = "Charge for hit mPMT, Q2";
		else yaxis = "Charge for other hit PMT, Q2";
		c5 = new TCanvas("c5", "Fifth canvas", 3840 * 2, 2160);
		c5->Divide(2, 1);
		c5->cd(1);
		nbins = 300;
		delete args;
		asprintf(&args, "log(q2):log(q1)>>pmtLogQvQ(%i, %g, %g, %i, %g, %g)", nbins, logMinX, logMaxX, nbins, logMinY, logMaxY);
		dataStore->Draw(args, "", "COLZ");
		TH2D* pmtLogQvQ = (TH2D*)gPad->GetPrimitive("pmtLogQvQ");
		pmtLogQvQ->SetContour(100); 
		pmtLogQvQ->GetXaxis()->SetTitle("Log Q for hit PMT, ln(Q1)"); 
		pmtLogQvQ->GetYaxis()->SetTitle(yaxis);

		c5->cd(2);
		nbins = 100;
		delete args;
		asprintf(&args, "log(q2):log(q1)>>LogQvQProfile(%i, %g, %g, %i, %g, %g)", nbins, logMinX, logMaxX, nbins, logMinY, logMaxY);
		dataStore->Draw(args, "", "PROF");
		TProfile* LogQvQProfile = (TProfile*)gPad->GetPrimitive("LogQvQProfile");
		LogQvQProfile->GetXaxis()->SetTitle("Log Q for hit PMT, ln(Q1)"); 
		LogQvQProfile->GetYaxis()->SetTitle(yaxis);
		c5->Draw();

		pmtQvQ->Write();
		QvQProfile->Write();
		pmtLogQvQ->Write();
		LogQvQProfile->Write();

		if (hybrid) {
			mh1->Write();
			mpmtQ->Write();
			mpmtPosX->Write();
			mpmtPosY->Write();
			mpmtPosZ->Write();
			mvertX->Write();
			mvertY->Write();
			mvertZ->Write();
		}

		TH2D* ratioQvQ, * ratioLogQvQ, * altRatioQvQ, * altRatioLogQvQ;
		TProfile* ratioQvQProfile, * ratioLogQvQProfile, * altRatioQvQProfile, * altRatioLogQvQProfile;
		if (ratios) {

			//ratioQvQ = new TH2D("ratioQvQ", ((std::string)"Ratio of (m)PMT charges (" + name + (std::string)")").c_str(), nbins, minX, maxX, nbins, minY, maxY);
			//ratioQvQProfile = new TProfile("QvQ Ratio Profile", ((std::string)"Ratio of (m)PMT charges (" + name + (std::string)")").c_str(), nbins, minX, maxX, minY, maxY);
			
			c6 = new TCanvas("c6", "Sixth canvas", 3840 * 2, 2160);
			c6->Divide(2, 1);
			yaxis = "Charge Ratio for hit (m)PMTs, Q1 / Q2";
			c6->cd(1);
			nbins = 100;
			delete args;
			asprintf(&args, "ratio:q2>>ratioQvQ(%i, %g, %g, %i, %g, %g)", nbins, 0.0, std::min(maxY / 3.0, 2000.0), nbins, minR, maxR);
			dataStore->Draw(args, "", "COLZ");
			ratioQvQ = (TH2D*)gPad->GetPrimitive("ratioQvQ");
			ratioQvQ->SetContour(100); 
			ratioQvQ->GetXaxis()->SetTitle("Charge on hit mPMT, Q2"); 
			ratioQvQ->GetYaxis()->SetTitle(yaxis);
			c6->cd(2);
			nbins = 50;
			delete args;
			asprintf(&args, "ratio:q2>>ratioQvQProfile(%i, %g, %g, %i, %g, %g)", nbins, 0.0, std::min(maxY / 3.0, 2000.0), nbins, minR, maxR);
			dataStore->Draw(args, "", "PROF");
			ratioQvQProfile = (TProfile*)gPad->GetPrimitive("ratioQvQProfile");
			ratioQvQProfile->GetXaxis()->SetTitle("Charge on hit mPMT, Q2"); 
			ratioQvQProfile->GetYaxis()->SetTitle(yaxis);
			c6->Draw();

			//std::cout << "Max (q1, q2): " << maxX << ", " << maxY << "| Min (q1, q2): " << minX << ", " << minY << std::endl;
			c7 = new TCanvas("c7", "Seventh canvas", 3840 * 2, 2160);
			c7->Divide(2, 1);
			yaxis = "Charge Ratio, Q1 / Q2";
			c7->cd(1);
			nbins = 100;
			delete args;
			asprintf(&args, "ratio:log(q2)>>ratioLogQvQ(%i, %g, %g, %i, %g, %g)", nbins, logMinY, logMaxY, nbins, minR, maxR);
			dataStore->Draw(args, "", "COLZ");
			ratioLogQvQ = (TH2D*)gPad->GetPrimitive("ratioLogQvQ");
			ratioLogQvQ->SetContour(100); 
			ratioLogQvQ->GetXaxis()->SetTitle("Charge on hit mPMT, ln(Q2)"); 
			ratioLogQvQ->GetYaxis()->SetTitle(yaxis);

			c7->cd(2);
			nbins = 50;
			delete args;
			asprintf(&args, "ratio:log(q2)>>ratioLogQvQProfile(%i, %g, %g, %i, %g, %g)", nbins, logMinY, logMaxY, nbins, minR, maxR);
			dataStore->Draw(args, "", "PROF");
			ratioLogQvQProfile = (TProfile*)gPad->GetPrimitive("ratioLogQvQProfile");
			ratioLogQvQProfile->GetXaxis()->SetTitle("Charge on hit mPMT, ln(Q2)"); 
			ratioLogQvQProfile->GetYaxis()->SetTitle(yaxis);
			c7->Draw();


			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Alternate X-axis

			c8 = new TCanvas("c8", "Eigth canvas", 3840 * 2, 2160);
			c8->Divide(2, 1);
			yaxis = "Charge Ratio for hit (m)PMTs, Q1 / Q2";
			c8->cd(1);
			nbins = 100;
			delete args;
			asprintf(&args, "ratio:q1>>altRatioQvQ(%i, %g, %g, %i, %g, %g)", nbins, 0.0, std::min(maxX / 3.0, 2000.0), nbins, minR, maxR);
			dataStore->Draw(args, "", "COLZ");
			altRatioQvQ = (TH2D*)gPad->GetPrimitive("altRatioQvQ");
			altRatioQvQ->SetContour(100); 
			altRatioQvQ->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); 
			altRatioQvQ->GetYaxis()->SetTitle(yaxis);

			c8->cd(2);
			nbins = 50;
			delete args;
			asprintf(&args, "ratio:q1>>altRatioQvQProfile(%i, %g, %g, %i, %g, %g)", nbins, 0.0, std::min(maxX / 3.0, 2000.0), nbins, minR, maxR);
			dataStore->Draw(args, "", "PROF");
			altRatioQvQProfile = (TProfile*)gPad->GetPrimitive("altRatioQvQProfile");
			altRatioQvQProfile->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); 
			altRatioQvQProfile->GetYaxis()->SetTitle(yaxis);
			c8->Draw();

			//std::cout << "Max (q1, q2): " << maxX << ", " << maxY << "| Min (q1, q2): " << minX << ", " << minY << std::endl;

			c9 = new TCanvas("c9", "Nineth canvas", 3840 * 2, 2160);
			c9->Divide(2, 1);
			yaxis = "Charge Ratio, Q1 / Q2";
			c9->cd(1);
			nbins = 100;
			delete args;
			asprintf(&args, "ratio:log(q1)>>altRatioLogQvQ(%i, %g, %g, %i, %g, %g)", nbins, logMinX, logMaxX, nbins, minR, maxR);
			dataStore->Draw(args, "", "COLZ");
			altRatioLogQvQ = (TH2D*)gPad->GetPrimitive("altRatioLogQvQ");
			altRatioLogQvQ->SetContour(100); 
			altRatioLogQvQ->GetXaxis()->SetTitle("Charge on hit PMT, ln(Q1)"); 
			altRatioLogQvQ->GetYaxis()->SetTitle(yaxis);

			c9->cd(2);
			nbins = 50;
			delete args;
			asprintf(&args, "ratio:log(q1)>>altRatioLogQvQProfile(%i, %g, %g, %i, %g, %g)", nbins, logMinX, logMaxX, nbins, minR, maxR);
			dataStore->Draw(args, "", "PROF");
			altRatioLogQvQProfile = (TProfile*)gPad->GetPrimitive("altRatioLogQvQProfile");
			altRatioLogQvQProfile->GetXaxis()->SetTitle("Charge on hit PMT, ln(Q1)"); 
			altRatioLogQvQProfile->GetYaxis()->SetTitle(yaxis);
			c9->Draw();
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
		data->Close();
		std::string out1 = fn1 + (std::string)"-positions.png";
		std::string out2 = fn1 + (std::string)"-verts.png";
		std::string out3 = fn1 + (std::string)"-main.png";
		std::string out4 = fn1 + (std::string)"-charges.png";
		std::string out5 = fn1 + (std::string)"-logCharges.png";

		c1->Print(out1.c_str());
		c2->Print(out2.c_str());
		c3->Print(out3.c_str());
		c4->Print(out4.c_str());
		c5->Print(out5.c_str());
		
		if (ratios) {
			std::string out6 = fn1 + (std::string)"-chargeRatio.png";
			std::string out7 = fn1 + (std::string)"-logChargeRatio.png";
			c6->Print(out6.c_str());
			c7->Print(out7.c_str());

			std::string out8 = fn1 + (std::string)"-chargeRatio_alternate.png";
			std::string out9 = fn1 + (std::string)"-logChargeRatio_alternate.png";
			c8->Print(out8.c_str());
			c9->Print(out9.c_str());
		}
		

		h1->Write();
		pmtQ->Write();
		dists->Write();
		pmtPosX->Write();
		pmtPosY->Write();
		pmtPosZ->Write();
		vertX->Write();
		vertY->Write();
		vertZ->Write();

		if (ratios) {
			ratioQvQ->Write();
			ratioQvQProfile->Write();
			ratioLogQvQ->Write();
			ratioLogQvQProfile->Write();
			altRatioQvQ->Write();
			altRatioLogQvQ->Write();
			altRatioQvQProfile->Write();
			altRatioLogQvQProfile->Write();
		}

		outfile->Close();
		outfilename = "";
		delete dataOut;
		delete c1, c2, c3, c4, c5;
		if (ratios) delete c6, c7, c8, c9;
	}
	
	//delete file, outfile;
	//delete h1, pmtQ, dists, pmtPosX, pmtPosY, pmtPosZ, vertX, vertY, vertZ;
	//delete pmtQvQ, pmtLogQvQ, QvQProfile, LogQvQProfile;
	//if (hybrid) delete mh1, mpmtQ, mpmtPosX, mpmtPosY, mpmtPosZ, mvertX, mvertY, mvertZ;
	return 0;
}
