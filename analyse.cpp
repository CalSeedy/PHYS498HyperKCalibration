#include <ctime>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
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
	if ((found != std::string::npos) && (found2 != std::string::npos)) return true;
	else return false;
}


struct Pos
{
	std::vector<double> data;

	Pos(std::vector<double> values) {
		for (int i = 0; i < values.size(); i++) {
			data.push_back(values.at(i));
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
};



int main(int argc, char** argv) {

	char* filename = NULL;
	char* outfilename = NULL;
	bool verbose = false;//false;
	bool hybrid = true;
	bool Gamma = false;//Is the mother particle a gamma or another particle?
	float cvacuum = 3e8 / 1e9;//speed of light, in meter per ns.
	float nindex = 1.373;//refraction index of water
	int nPMTtypes = 2;
	int c = -1;
	while ((c = getopt(argc, argv, "f:o:v")) != -1) {//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
		switch (c) {
		case 'f':
			filename = optarg;
			break;
		case 'o':
			outfilename = optarg;
			break;
		case 'v':
			verbose = true;
			break;
		}
	}
	auto start = std::chrono::system_clock::now();
	TFile* file;
	// Open the file
	if (filename == NULL) {
		filename = "./RootFiles/wcsim_hybrid\(10MeVe-\).root";
	}

	file = new TFile(filename, "read");

	if (!file->IsOpen()) {
		std::cout << "Error, could not open input file: " << filename << std::endl;
		return -1;
	}
	std::cout << "Analysing " << filename << std::endl;
	// Get the a pointer to the tree from the file
	TTree* tree = (TTree*)file->Get("wcsimT");

	// Get the number of events
	int nevent = ((int)tree->GetEntries());//std::min(((int)tree->GetEntries()),100000);
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
	opttree->GetEntry(0);
	opt->Print();

	// start with the main "subevent", as it contains most of the info
	// and always exists.
	WCSimRootTrigger* wcsimrootevent;
	WCSimRootTrigger* wcsimrootevent2;

	int num_trig = 0;

	//TFile * fOutput = new TFile(OutputFile,"recreate");

	TH1D* h1 = new TH1D("PMT Hits", "Number of PMT Hits", 30, 0, 30);
	TH1D* pmtQ = new TH1D("PMT Charge", "Hit PMT Charge", 50, 0, 20);
	TH1D* dists = new TH1D("PMT Distances (<= 10m)", "Hit PMT Distances (<= 10m)", 100, 0, 1000);
	TH1D* pmtPosX = new TH1D("X Pos", "X Pos for hit PMTs", 200, -4000, 4000);
	TH1D* pmtPosY = new TH1D("Y Pos", "Y Pos for hit PMTs", 200, -4000, 4000);
	TH1D* pmtPosZ = new TH1D("Z Pos", "Z Pos for hit PMTs", 200, -4000, 4000);
	TH1D* vertX = new TH1D("Vert X", "Event vertex X positions", 1000, -5000, 5000);
	TH1D* vertY = new TH1D("Vert Y", "Event vertex Y positions", 1000, -5000, 5000);
	TH1D* vertZ = new TH1D("Vert Z", "Event vertex Z positions", 1000, -5000, 5000);

	TH1D* mh1 = new TH1D("mPMT Hits", "Number of mPMT Hits", 30, 0, 30);
	TH1D* mpmtQ = new TH1D("mPMT Charge", "Hit mPMT Charge", 50, 0, 20);
	TH1D* mpmtPosX = new TH1D("mPMT X Pos", "X Pos for hit mPMTs", 200, -4000, 4000);
	TH1D* mpmtPosY = new TH1D("mPMT Y Pos", "Y Pos for hit mPMTs", 200, -4000, 4000);
	TH1D* mpmtPosZ = new TH1D("mPMT Z Pos", "Z Pos for hit mPMTs", 200, -4000, 4000);
	TH1D* mvertX = new TH1D("Vert X", "Event vertex X positions", 1000, -5000, 5000);
	TH1D* mvertY = new TH1D("Vert Y", "Event vertex Y positions", 1000, -5000, 5000);
	TH1D* mvertZ = new TH1D("Vert Z", "Event vertex Z positions", 1000, -5000, 5000);

	std::vector<Pos> positions;
	std::vector<Pos> mPositions;
	std::vector<int> IDs;
	std::vector<int> mIDs;
	std::vector<double> Qdata(nPMT);
	std::vector<double> mQdata(nmPMT);

	std::vector<std::string> consoletxt;
	// Now loop over events
	for (int ev = 0; ev < nevent; ev++)
	{

		// Read the event from the tree into the WCSimRootEvent instance
		tree->GetEntry(ev);

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
		
		//event is the same??
		/*
		if (hybrid) {
			mvertX->Fill(wcsimrootevent2->GetVtx(0));
			mvertY->Fill(wcsimrootevent2->GetVtx(1));
			mvertZ->Fill(wcsimrootevent2->GetVtx(2));
		}*/

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

		// Now read the tracks in the event
		// Get the number of tracks
		int ntrack = wcsimrootevent->GetNtrack();
		if (verbose) printf("ntracks=%d\n", ntrack);

		double particleStart[3];
		double particleStop[3];
		double particleDir[3];

		// Loop through elements in the TClonesArray of WCSimTracks
		for (int i = 0; i < ntrack; i++)
		{
			TObject* element = (wcsimrootevent->GetTracks())->At(i);
			WCSimRootTrack* wcsimroottrack = (WCSimRootTrack*)(element);

			if (i == (ntrack - 1)) {//Mother particle
				for (int j = 0; j < 3; j++) {
					particleStart[j] = wcsimroottrack->GetStart(j);
					particleStop[j] = wcsimroottrack->GetStop(j);
					particleDir[j] = wcsimroottrack->GetDir(j);
				}
			}
			if (verbose) {
				printf("Track ipnu: %d\n", wcsimroottrack->GetIpnu());
				printf("Track parent ID: %d\n", wcsimroottrack->GetParenttype());
				printf("Track energy: %f\n", wcsimroottrack->GetE());
				printf("Track momentum: %f\n", wcsimroottrack->GetP());
				printf("Track mass: %f\n", wcsimroottrack->GetM());

				for (int j = 0; j < 3; j++) {
					printf("Track start: %d %f\n", j, wcsimroottrack->GetStart(j));
					printf("Track dir: %d %f\n", j, wcsimroottrack->GetDir(j));
				}
			}
		}  // End of loop over tracks

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
					positions.push_back(PMTpos);
					IDs.push_back(tubeNumber - 1);
				}
				else {
					mPositions.push_back(PMTpos);
					mIDs.push_back(tubeNumber - 1);
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

			int nhits;
			if (pmtType == 0) nhits = ncherenkovdigihits;
			else nhits = ncherenkovdigihits2;

			// Loop through elements in the TClonesArray of WCSimRootCherenkovHits
			for (int i = 0; i < nhits; i++)
			{
				TObject* Hit;
				if (pmtType == 0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
				else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

				WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit =
					dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

				int tubeNumber = wcsimrootcherenkovdigihit->GetTubeId();
				double peForTube = wcsimrootcherenkovdigihit->GetQ();

				WCSimRootPMT pmt;
				if (pmtType == 0) {
					pmt = geo->GetPMT(tubeNumber - 1, false);
					pmtQ->Fill(peForTube);
					Qdata[tubeNumber - 1] += peForTube;
				}
				else {
					pmt = geo->GetPMT(tubeNumber - 1, true);
					mpmtQ->Fill(peForTube);
					mQdata[tubeNumber - 1] += peForTube;
				}


				totalQ += peForTube;
				double PMTpos[3] = { pmt.GetPosition(0), pmt.GetPosition(1), pmt.GetPosition(2), };
				if (verbose) {
					printf("Total pe: %1.1f times( ", peForTube);
					std::cout << ")" << std::endl;
					std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << std::endl;
				}

				double distance_pmt_vertex = 0;
				int mPMT_number = pmt.GetmPMTNo();
				int pmt_number_in_mPMT = pmt.GetmPMT_PMTNo();
			} // End of loop over Cherenkov digitised hits
			if (verbose) std::cout << "Total Pe : " << totalPe << ", total hit : " << totalHit << std::endl;
		}

		// reinitialize super event between loops.
		wcsimrootsuperevent->ReInitialize();
		if (hybrid) wcsimrootsuperevent2->ReInitialize();
	} // End of loop over events
	std::cout << "num_trig " << num_trig << "\n";

	std::vector<double> filterQ1, filterQ2;

	for (int i = 0; i < positions.size(); i++) { // for every standard PMT
		if (hybrid) {
			for (int j = 0; j < mPositions.size(); j++) { //  for every mPMT
				try {
					Pos pos1 = positions.at(i);
					double x1 = pos1.Get(0);
					double y1 = pos1.Get(1);
					double z1 = pos1.Get(2);

					Pos pos2 = mPositions.at(j);
					double x2 = pos2.Get(0);
					double y2 = pos2.Get(1);
					double z2 = pos2.Get(2);

					if (verbose) {
						pos1.Print();
						pos2.Print();
					}

					double dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

					if (dist <= 1000) {
						dists->Fill(dist);

						if (dist <= 150) {
							filterQ1.push_back(Qdata.at(IDs.at(i)));
							filterQ2.push_back(mQdata.at(mIDs.at(j)));
						}
					}

					if ((i % 5 == 0) && (j % (mPositions.size() - 1) == 0)) {
						//std::system("clear");
						auto now = std::chrono::system_clock::now();
						for (auto st : consoletxt) std::cout << st << std::endl;
						std::cout << filename << "| " << i << ", " << j << "| Time Running: " << std::chrono::duration<double>(now - start).count() << " s" << std::endl;
					}

				}
				catch (std::out_of_range& err) {
					std::cout << i << ", " << j << " (" << IDs.at(i) << ", " << mIDs.at(j) << ")" << std::endl;
					return 1;
				}
			}
		}
		else {
			for (int j = i; j < positions.size(); j++) { //  for every other PMT
				try {
					if (IDs.at(i) != IDs.at(j)) {
						Pos pos1 = positions.at(i);
						double x1 = pos1.Get(0);
						double y1 = pos1.Get(1);
						double z1 = pos1.Get(2);

						Pos pos2 = positions.at(j);
						double x2 = pos2.Get(0);
						double y2 = pos2.Get(1);
						double z2 = pos2.Get(2);

						if (verbose) {
							pos1.Print();
							pos2.Print();
						}

						double dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

						if (dist <= 1000) {
							dists->Fill(dist);

							if (dist <= 150) {
								filterQ1.push_back(Qdata.at(IDs.at(i)));
								filterQ2.push_back(Qdata.at(IDs.at(j)));
							}
						}

						if ((i % 5 == 0) && (j % (positions.size() - 1) == 0)) {
							//std::system("clear");
							auto now = std::chrono::system_clock::now();
							for (auto st : consoletxt) std::cout << st << std::endl;
							std::cout << (std::string)filename << "| " << i << ", " << j << "| Time Running: " << std::chrono::duration<double>(now - start).count() << " s" << std::endl;
						}
					}
				}
				catch (std::out_of_range& err) {
					std::cout << i << ", " << j << " (" << IDs.at(i) << ", " << IDs.at(j) << ")" << std::endl;
					return 1;
				}


			}
		}
		if (i == (positions.size() - 1)) {
			auto now = std::chrono::system_clock::now();
			std::string s = (std::string)filename + "| Finished at: " + std::to_string(std::chrono::duration<double>(now - start).count()) + " s";
			consoletxt.push_back(s);
		}
	}
	//std::system("clear");
	for (auto st : consoletxt) std::cout << st << std::endl;

	int nHigh = 3;
	int nWide = 1;

	TCanvas* c1 = new TCanvas("c1", "First canvas", nWide * 1920, nHigh * 1080);
	c1->Divide(nWide, nHigh);
	c1->cd(1); pmtPosX->GetXaxis()->SetTitle("PMT Position, x [cm]"); pmtPosX->GetYaxis()->SetTitle("# of PMTs"); pmtPosX->SetFillColor(kBlue); pmtPosX->Draw();
	c1->cd(2); pmtPosY->GetXaxis()->SetTitle("PMT Position, y [cm]"); pmtPosY->GetYaxis()->SetTitle("# of PMTs"); pmtPosY->SetFillColor(kBlue); pmtPosY->Draw();
	c1->cd(3); pmtPosZ->GetXaxis()->SetTitle("PMT Position, z [cm]"); pmtPosZ->GetYaxis()->SetTitle("# of PMTs"); pmtPosZ->SetFillColor(kBlue); pmtPosZ->Draw();
	if (hybrid) {
		c1->cd(1); mpmtPosX->SetFillColor(kRed); mpmtPosX->Draw("SAME");
		c1->cd(2); mpmtPosY->SetFillColor(kRed); mpmtPosY->Draw("SAME");
		c1->cd(3); mpmtPosZ->SetFillColor(kRed); mpmtPosZ->Draw("SAME");
	}
	c1->Draw();


	TCanvas* c2 = new TCanvas("c2", "Second canvas", nWide * 1920, nHigh * 1080);
	c2->Divide(nWide, nHigh);
	c2->cd(1); vertX->GetXaxis()->SetTitle("Event Vertex, x [cm]"); vertX->GetYaxis()->SetTitle("Count"); vertX->SetFillColor(kBlue); vertX->Draw();
	if (hybrid) mvertX->SetFillColor(kRed); mvertX->Draw("SAME");
	c2->cd(2); vertY->GetXaxis()->SetTitle("Event Vertex, y [cm]"); vertY->GetYaxis()->SetTitle("Count"); vertY->SetFillColor(kBlue); vertY->Draw();
	if (hybrid) mvertY->SetFillColor(kRed); mvertY->Draw("SAME");
	c2->cd(3); vertZ->GetXaxis()->SetTitle("Event Vertex, z [cm]"); vertZ->GetYaxis()->SetTitle("Count"); vertZ->SetFillColor(kBlue); vertZ->Draw();
	if (hybrid) mvertZ->SetFillColor(kRed); mvertZ->Draw("SAME");
	c2->Draw();


	TCanvas* c3 = new TCanvas("c3", "Third canvas", nWide * 1920, nHigh * 1080);
	c3->Divide(nWide, nHigh);
	c3->cd(1); h1->GetXaxis()->SetTitle("Hits on a PMT"); h1->GetYaxis()->SetTitle("Count"); h1->SetFillColor(kBlue); h1->Draw();
	if (hybrid) mh1->SetFillColor(kRed); mh1->Draw("SAME");
	c3->cd(2); pmtQ->GetXaxis()->SetTitle("Charge on hit PMT [eV]"); pmtQ->GetYaxis()->SetTitle("Count"); pmtQ->SetFillColor(kBlue); pmtQ->Draw();
	if (hybrid) mpmtQ->SetFillColor(kRed); mpmtQ->Draw("SAME");
	c3->cd(3); dists->GetXaxis()->SetTitle("Distance between hit PMTs [cm]"); dists->GetYaxis()->SetTitle("Count"); dists->SetFillColor(kBlue); dists->Draw();
	c3->Draw();


	double maxQ1 = double(filterQ1.at(std::distance(filterQ1.begin(), std::max_element(std::begin(filterQ1), std::end(filterQ1)))));
	double maxQ2 = double(filterQ2.at(std::distance(filterQ2.begin(), std::max_element(std::begin(filterQ2), std::end(filterQ2)))));
	double minQ1 = double(filterQ1.at(std::distance(filterQ1.begin(), std::min_element(std::begin(filterQ1), std::end(filterQ1)))));
	double minQ2 = double(filterQ2.at(std::distance(filterQ2.begin(), std::min_element(std::begin(filterQ2), std::end(filterQ2)))));
	maxQ1 *= 1.1;
	maxQ2 *= 1.1;
	if (minQ1 > 0.0) minQ1 = 0.0;
	else minQ1 *= 1.1;

	if (minQ2 > 0.0) minQ2 = 0.0;
	else minQ2 *= 1.1;

	std::string fn(filename);
	fn = fn.substr(0, fn.find(".root"));

	int nbins = 200;
	TH2D* pmtQvQ = new TH2D("QvQ", ((std::string)"Q2 vs Q1 for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1, nbins, minQ2, maxQ2);
	nbins = 50;
	TProfile* QvQProfile = new TProfile("QvQ Profile", ((std::string)"Q2 vs Q1 for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1, minQ2, maxQ2);

	for (int i = 0; i < filterQ1.size(); i++) {
		pmtQvQ->Fill(filterQ1.at(i), filterQ2.at(i));
		QvQProfile->Fill(filterQ1.at(i), filterQ2.at(i));
	}

	TCanvas* c4 = new TCanvas("c4", "Fourth canvas", 1920 * 2, 1080);
	c4->Divide(2, 1);
	c4->cd(1); pmtQvQ->SetContour(100); pmtQvQ->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); pmtQvQ->GetYaxis()->SetTitle("Charge on other hit PMT, Q2"); pmtQvQ->Draw("COLZ");
	c4->cd(2); QvQProfile->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); QvQProfile->GetYaxis()->SetTitle("Charge on other hit PMT, Q2"); QvQProfile->Draw();
	c4->Draw();


	for (int i = 0; i < filterQ1.size(); i++) {
		filterQ1.at(i) = log(filterQ1.at(i));
	}
	for (int i = 0; i < filterQ2.size(); i++) {
		filterQ2.at(i) = log(filterQ2.at(i));
	}

	maxQ1 = double(filterQ1.at(std::distance(filterQ1.begin(), std::max_element(std::begin(filterQ1), std::end(filterQ1)))));
	maxQ2 = double(filterQ2.at(std::distance(filterQ2.begin(), std::max_element(std::begin(filterQ2), std::end(filterQ2)))));
	minQ1 = double(filterQ1.at(std::distance(filterQ1.begin(), std::min_element(std::begin(filterQ1), std::end(filterQ1)))));
	minQ2 = double(filterQ2.at(std::distance(filterQ2.begin(), std::min_element(std::begin(filterQ2), std::end(filterQ2)))));
	maxQ1 *= 1.1;
	maxQ2 *= 1.1;
	if (minQ1 > 0.0) minQ1 = 0.0;
	else minQ1 *= 1.1;

	if (minQ2 > 0.0) minQ2 = 0.0;
	else minQ2 *= 1.1;

	nbins = 200;
	TH2D* pmtLogQvQ = new TH2D("LogQvQ", ((std::string)"ln(Q2) vs ln(Q1) for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1, nbins, minQ2, maxQ2);
	nbins = 100;
	TProfile* LogQvQProfile = new TProfile("Log QvQ Profile", ((std::string)"ln(Q2) vs ln(Q1) for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1, minQ2, maxQ2);

	for (int i = 0; i < filterQ1.size(); i++) {
		pmtLogQvQ->Fill(filterQ1.at(i), filterQ2.at(i));
		//LogQvQProfile->Fill(filterQ1.at(i), filterQ2.at(i));

	}


	TCanvas* c5 = new TCanvas("c5", "Fifth canvas", 1920 * 2, 1080);
	c5->Divide(2, 1);
	c5->cd(1); pmtLogQvQ->SetContour(100); pmtLogQvQ->GetXaxis()->SetTitle("Log Q for hit PMT, ln(Q1)"); pmtLogQvQ->GetYaxis()->SetTitle("Log Q for hit mPMT, ln(Q2)"); pmtLogQvQ->Draw("COLZ");
	c5->cd(2); LogQvQProfile->GetXaxis()->SetTitle("Log Q for hit PMT, ln(Q1)"); LogQvQProfile->GetYaxis()->SetTitle("Log Q for hit mPMT, ln(Q2)"); LogQvQProfile->Draw();
	c5->Draw();


	std::string out1 = fn.substr(12, fn.length()) + (std::string)"-positions.png";
	std::string out2 = fn.substr(12, fn.length()) + (std::string)"-verts.png";
	std::string out3 = fn.substr(12, fn.length()) + (std::string)"-main.png";
	std::string out4 = fn.substr(12, fn.length()) + (std::string)"-charges.png";
	std::string out5 = fn.substr(12, fn.length()) + (std::string)"-logCharges.png";


	c1->Print(out1.c_str());
	c2->Print(out2.c_str());
	c3->Print(out3.c_str());
	c4->Print(out4.c_str());
	c5->Print(out5.c_str());

	//if (outfilename == NULL) sprintf(outfilename, "out.root");
	/*
	TFile* outfile = new TFile(outfilename, "RECREATE");
	std::cout << "File " << outfilename << " is open for writing" << std::endl;


	outfile->Close();
	*/
	return 0;
}
