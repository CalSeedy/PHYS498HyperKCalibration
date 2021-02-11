#include <algorithm>
#include <ctime>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdio.h>     
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "./include/WCSimRootEvent.hh"
#include "./include/WCSimRootGeom.hh"
#include "./include/WCSimRootOptions.hh"


struct Pos
{
	std::vector<double> data;

	Pos(std::vector<double> values) {
		for (int i = 0; i < values.size(); i++) {
			data.push_back(values.at(i));
		}
	}

	Pos() {
		for (int i = 0; i < 3; i++) {
			data.push_back(0.0);
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

bool is_hybrid(std::string filename) {
	if (filename.find("mPMT") != std::string::npos)	return true;
	else return false;
}



// Simple example of reading a generated Root file
int readPMT(bool verbose = false)//const char* filename = "wcsim.root", bool verbose = false)
{
	std::vector<std::string> rootfiles;
	DIR* dir;
	struct dirent* ent;
	if ((dir = opendir("./RootFiles/")) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir(dir)) != NULL) {
			std::string name(ent->d_name);
			//printf("%s\n", name.c_str());
			if ((name.find("_flat") == std::string::npos) || (name.find("/") == std::string::npos)) {
				rootfiles.push_back("./RootFiles/" + name);
			}
		}
		closedir(dir);
	}
	else {
		/* could not open directory */
		perror("");
		std::cout << "Failed to open file!" << std::endl;
		return EXIT_FAILURE;
	}
	rootfiles.pop_back();
	rootfiles.pop_back();

	for (auto rf : rootfiles) {
		std::cout << rf << std::endl;
	}

	std::ofstream txtfile;
	txtfile.open("output.txt");

	auto start = std::chrono::system_clock::now();
	//std::stringstream consoletxt;
	std::vector<std::string> consoletxt;
	for (int idx = 0; idx < rootfiles.size(); idx++) {
		std::string filename = rootfiles.at(idx);
		std::cout << "Analysing Root file: " << filename.c_str() << std::endl;
		bool hybrid = is_hybrid(filename);
		double totalQ = 0.0;
		txtfile << "----" << filename.c_str() << "----" << std::endl;
		// Clear global scope
		//gROOT->Reset();

		// Open the file
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
		file->Close();

		TH1D* h1 = new TH1D("PMT Hits", "Number of PMT Hits", 30, 0, 30);
		TH1D* pmtQ = new TH1D("PMT Charge", "Hit PMT Charge", 50, 0, 20);
		TH1D* dists = new TH1D("PMT Distances (<= 10m)", "Hit PMT Distances (<= 10m)", 100, 0, 1000);
		TH1D* pmtPosX = new TH1D("X Pos", "X Pos for hit PMTs", 200, -4500, 4500);
		TH1D* pmtPosY = new TH1D("Y Pos", "Y Pos for hit PMTs", 200, -4500, 4500);
		TH1D* pmtPosZ = new TH1D("Z Pos", "Z Pos for hit PMTs", 200, -4500, 4500);
		TH1D* vertX = new TH1D("Vert X", "Event vertex X positions", 200, -100, 100);
		TH1D* vertY = new TH1D("Vert Y", "Event vertex Y positions", 200, -100, 100);
		TH1D* vertZ = new TH1D("Vert Z", "Event vertex Z positions", 200, -100, 100);


		TH1D* mpmtQ = new TH1D("mPMT Charge", "Hit mPMT Charge", 50, 0, 20);
		TH1D* mpmtDists = new TH1D("mPMT Distances (<= 10m)", "Hit mPMT Distances (<= 10m)", 100, 0, 1000);
		TH1D* mpmtPosX = new TH1D("mPMT X Pos", "X Pos for hit mPMTs", 200, -4500, 4500);
		TH1D* mpmtPosY = new TH1D("mPMT Y Pos", "Y Pos for hit mPMTs", 200, -4500, 4500);
		TH1D* mpmtPosZ = new TH1D("mPMT Z Pos", "Z Pos for hit mPMTs", 200, -4500, 4500);

		if (!hybrid) delete mpmtQ, mpmtDists, mpmtPosX, mpmtPosY, mpmtPosZ;

		std::vector<Pos> positions;
		std::vector<int> IDs, IDs_mPMT;
		std::vector<double> Qdata(num_PMTs + num_mPMTs); //TArrayD* Qdata  = new TArrayD(off);
		std::vector<bool> is_mPMT;
		for (double& v : Qdata) { v = 0.0; }
		int num_trig = 0;
		//std::cout << "Detected " << numPMTs << " PMTs!" << std::endl;
		double maxQ1(0.0), maxQ2(0.0), minQ1(0.0), minQ2(0.0);

		WCSimRootTrigger* wcsimrootevent;
		WCSimRootTrigger* wcsimrootevent2;
		// Now loop over events
		for (int ev = 0; ev < nevent; ev++)
		{
			// Read the event from the tree into the WCSimRootEvent instance
			tree->GetEntry(ev);
			wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
			if (hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);

			std::vector<float> triggerInfo;
			triggerInfo.clear();
			triggerInfo = wcsimrootevent->GetTriggerInfo();

			std::vector<float> triggerInfo2;
			triggerInfo2.clear();
			if (hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();


			vertX->Fill(wcsimrootevent->GetVtx(0));
			vertY->Fill(wcsimrootevent->GetVtx(1));
			vertZ->Fill(wcsimrootevent->GetVtx(2));


			txtfile << "Event: " << ev + 1 << "----" << std::endl;
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
			// Now read the tracks in the event

			// Get the number of tracks
			int ntrack = wcsimrootevent->GetNtrack();
			if (verbose) printf("ntracks=%d\n", ntrack);

			int i;
			// Loop through elements in the TClonesArray of WCSimTracks
			for (i = 0; i < ntrack; i++)
			{
				TObject* element = (wcsimrootevent->GetTracks())->At(i);

				//WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
				WCSimRootTrack* wcsimroottrack = (WCSimRootTrack*)(element);

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
			txtfile << "Tracks: " << i << std::endl;


			int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
			int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
			int ncherenkovhits2 = 0; if (hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();
			int ncherenkovdigihits2 = 0; if (hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits();

			if (verbose) {
				printf("node id: %i\n", ev);
				printf("Ncherenkovhits (normal) %d\n", ncherenkovhits);
				printf("Ncherenkovdigihits (normal) %d\n", ncherenkovdigihits);
				printf("Ncherenkovhits (mPMT) %d\n", ncherenkovhits2);
				printf("Ncherenkovdigihits (mPMT) %d\n", ncherenkovdigihits2);
				std::cout << "RAW HITS:" << std::endl;
			}

			for (int pmtType = 0; pmtType < 2; pmtType++) {
				// Grab the big arrays of times and parent IDs
				TClonesArray* timeArray;//An array of pointers on CherenkovHitsTimes.
				if (pmtType == 0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
				else timeArray = wcsimrootevent2->GetCherenkovHitTimes();

				int totalPe = 0;
				// Loop through elements in the TClonesArray of WCSimRootCherenkovHits
				int nhits;
				if (pmtType == 0) nhits = ncherenkovhits;
				else nhits = ncherenkovhits2;

				for (i = 0; i < nhits; i++)
				{
					WCSimRootCherenkovHit* wcsimrootcherenkovhit;
					if (pmtType == 0) wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>((wcsimrootevent->GetCherenkovHits())->At(i));
					else wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>((wcsimrootevent2->GetCherenkovHits())->At(i));
					TObject* Hit = (wcsimrootevent->GetCherenkovHits())->At(i);

					int tubeNumber = wcsimrootcherenkovhit->GetTubeID();

					int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
					int peForTube = wcsimrootcherenkovhit->GetTotalPe(1);

					WCSimRootPMT pmt = geo->GetPMT(tubeNumber - 1, (bool)pmtType);
					totalPe += peForTube;
					h1->Fill(peForTube);

					if (verbose) printf("Total pe: %d times( ", peForTube);

					for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
					{
						WCSimRootCherenkovHitTime* HitTime =
							dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));

						if (verbose) printf("%6.2f ", HitTime->GetTruetime());
					}
					if (verbose) std::cout << ")" << std::endl;

				} // End of loop over Cherenkov hits
				if (verbose) std::cout << "Total Pe : " << totalPe << std::endl;
				txtfile << "Total Pe : " << totalPe << std::endl;
			}

			// Look at digitized hit info
			// Get the number of digitized hits
			// Loop over sub events

			if (verbose) std::cout << "DIGITIZED HITS:" << std::endl;
			for (int pmtType = 0; pmtType < 2; pmtType++) {
				int nhits;
				if (pmtType == 0) nhits = ncherenkovdigihits;
				else nhits = ncherenkovdigihits2;
				if (verbose) std::cout << "PMT type = " << pmtType << std::endl;

				TClonesArray* timeArray;
				if (pmtType == 0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
				else timeArray = wcsimrootevent2->GetCherenkovHitTimes();

				double totalPe = 0;
				int totalHit = 0;


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
					is_mPMT.push_back(pmtType);
					IDs.push_back(tubeNumber - 1);
					double peForTube = wcsimrootcherenkovdigihit->GetQ();
					if (pmtType == 0) pmtQ->Fill(peForTube);
					else mpmtQ->Fill(peForTube);
					Qdata[tubeNumber - 1] += peForTube;
					totalQ += peForTube;

					if (verbose) printf("q, t, tubeid: %f %f %d \n", peForTube, wcsimrootcherenkovdigihit->GetT(), tubeNumber);

					WCSimRootPMT pmt;
					if (pmtType == 0) pmt = geo->GetPMT(tubeNumber - 1, false);
					else pmt = geo->GetPMT(tubeNumber - 1, true);

					double xpos = pmt.GetPosition(0);
					double ypos = pmt.GetPosition(1);
					double zpos = pmt.GetPosition(2);
					std::vector<double> pos{ xpos, ypos, zpos };
					Pos position(pos);
					positions.push_back(position);
					if (verbose) printf("{X: %f, Y: %f, Z: %f}, ", xpos, ypos, zpos);
					if (pmtType != 0) {
						mpmtPosX->Fill(xpos);
						mpmtPosY->Fill(ypos);
						mpmtPosZ->Fill(zpos);
					}
					else {
						pmtPosX->Fill(xpos);
						pmtPosY->Fill(ypos);
						pmtPosZ->Fill(zpos);
					}
				} // End of loop over Cherenkov digihits
				txtfile << "Total Hits: " << ncherenkovdigihits << ", " << "Q: " << totalQ << std::endl;
			}

			// reinitialize super event between loops.
			wcsimrootsuperevent->ReInitialize();
			if (hybrid) wcsimrootsuperevent2->ReInitialize();
		} // End of loop over events

		float win_scale = 0.75;
		int n_wide(3);
		int n_high(1);

		std::vector<double> filterQ1, filterQ2;
		for (int i = 0; i < positions.size(); i++) {
			for (int j = i; j < positions.size(); j++) {
				try {
					if (IDs.at(i) != IDs.at(j)) {
						Pos pos1, pos2; //pos1 = mPMT if hybrid
						if (hybrid) {
							if (is_mPMT.at(i) && !is_mPMT.at(j)) {
								pos1 = positions.at(i);
								pos2 = positions.at(j);
							}
							else if (!is_mPMT.at(i) && is_mPMT.at(j)) {
								pos1 = positions.at(j);
								pos2 = positions.at(i);
							}
							else continue;

						}
						else {
							pos1 = positions.at(i);
							pos2 = positions.at(j);
						}
						double x1 = pos1.Get(0);
						double y1 = pos1.Get(1);
						double z1 = pos1.Get(2);

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
								double Q1, Q2;
								if (!hybrid) {
									Q1 = Qdata.at(IDs.at(i));
									Q2 = Qdata.at(IDs.at(j));
								}
								else {
									if (is_mPMT.at(i) && !is_mPMT.at(j)) {
										Q1 = Qdata.at(IDs.at(j));
										Q2 = Qdata.at(IDs.at(i));
									}
									else if (!is_mPMT.at(i) && is_mPMT.at(j)) {
										Q1 = Qdata.at(IDs.at(i));
										Q2 = Qdata.at(IDs.at(j));
									}
									else {
										Q1 = Qdata.at(IDs.at(i));
										Q2 = Qdata.at(IDs.at(j));
									}
								} // Force Q2 to be mPMT if hybrid
								if (Q1 > maxQ1) maxQ1 = Q1;
								if (Q1 < minQ1) minQ1 = Q1;
								if (Q2 > maxQ2) maxQ2 = Q2;
								if (Q2 < minQ2) minQ2 = Q2;
								filterQ1.push_back(Q1);
								filterQ2.push_back(Q2);
							}
						}

						if ((i % 5 == 0) && (j % (positions.size() - 1) == 0)) {
							std::system("clear");
							auto now = std::chrono::system_clock::now();
							if (idx > 0) for (auto st : consoletxt) std::cout << st << std::endl;
							std::cout << filename.c_str() << "| " << i << ", " << j << "| Time Running: " << chrono::duration<double>(now - start).count() << " s" << std::endl;
						}
					}
				}
				catch (std::out_of_range& err) {
					std::cout << i << ", " << j << " (" << IDs.at(i) << ", " << IDs.at(j) << ")" << std::endl;
					return 1;
				}
			}
			if (i == (positions.size() - 1)) {
				auto now = std::chrono::system_clock::now();
				std::string s = filename + "| Finished at: " + std::to_string(chrono::duration<double>(now - start).count()) + " s";
				consoletxt.push_back(s);
			}
		}
		std::system("clear");
		if (idx > 0) for (auto st : consoletxt) std::cout << st << std::endl;


		TCanvas* c1 = new TCanvas("c1", "First canvas", 1920 * win_scale, 1080 * win_scale);
		c1->Draw();
		c1->Divide(n_high, n_wide);
		c1->cd(1); pmtPosX->GetXaxis()->SetTitle("PMT Position, x [cm]"); pmtPosX->GetYaxis()->SetTitle("# of PMTs"); pmtPosX->Draw();
		c1->cd(2); pmtPosY->GetXaxis()->SetTitle("PMT Position, y [cm]"); pmtPosY->GetYaxis()->SetTitle("# of PMTs"); pmtPosY->Draw();
		c1->cd(3); pmtPosZ->GetXaxis()->SetTitle("PMT Position, z [cm]"); pmtPosZ->GetYaxis()->SetTitle("# of PMTs"); pmtPosZ->Draw();


		TCanvas* c2 = new TCanvas("c2", "Second canvas", 1920 * win_scale, 1080 * win_scale);
		c2->Draw();
		c2->Divide(n_high, n_wide);
		c2->cd(1); vertX->GetXaxis()->SetTitle("Event Vertex, x [cm]"); vertX->GetYaxis()->SetTitle("Count"); vertX->Draw();
		c2->cd(2); vertY->GetXaxis()->SetTitle("Event Vertex, y [cm]"); vertY->GetYaxis()->SetTitle("Count"); vertY->Draw();
		c2->cd(3); vertZ->GetXaxis()->SetTitle("Event Vertex, z [cm]"); vertZ->GetYaxis()->SetTitle("Count"); vertZ->Draw();


		std::cout << "num_trig " << num_trig << "\n";


		TCanvas* c3 = new TCanvas("c3", "Third canvas", 1920 * win_scale, 1080 * win_scale);
		c3->Draw();
		c3->Divide(n_high, n_wide);
		c3->cd(1); h1->GetXaxis()->SetTitle("Hits on a PMT"); h1->GetYaxis()->SetTitle("Count"); h1->Draw();
		c3->cd(2); pmtQ->GetXaxis()->SetTitle("Charge on hit PMT [eV]"); pmtQ->GetYaxis()->SetTitle("Count"); pmtQ->Draw();
		c3->cd(3); dists->GetXaxis()->SetTitle("Distance between hit PMTs [cm]"); dists->GetYaxis()->SetTitle("Count"); dists->Draw();


		std::string fn(filename);
		fn = fn.substr(0, fn.find(".root"));

		int nbins = 200;
		TH2D* pmtQvQ = new TH2D("QvQ", ((std::string)"Q2 vs Q1 for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1 * 1.1, nbins, minQ2, maxQ2 * 1.1);
		nbins = nbins / 2;
		TProfile* QvQProfile = new TProfile("QvQ Profile", ((std::string)"Q2 vs Q1 for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1 * 1.1, minQ2, maxQ2 * 1.1);

		for (int i = 0; i < filterQ1.size(); i++) {
			pmtQvQ->Fill(filterQ1.at(i), filterQ2.at(i));
			QvQProfile->Fill(filterQ1.at(i), filterQ2.at(i));
		}

		TCanvas* c4 = new TCanvas("c4", "Fourth canvas", 1920 * 2 * win_scale, 1080 * win_scale);
		c4->Draw();
		c4->Divide(2, 1);
		c4->cd(1); pmtQvQ->SetContour(100); pmtQvQ->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); pmtQvQ->GetYaxis()->SetTitle("Charge on other hit PMT, Q2"); pmtQvQ->Draw("COLZ");
		c4->cd(2); QvQProfile->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); QvQProfile->GetYaxis()->SetTitle("Charge on other hit PMT, Q2"); QvQProfile->Draw();



		for (int i = 0; i < filterQ1.size(); i++) {
			filterQ1.at(i) = log(filterQ1.at(i));
		}
		for (int i = 0; i < filterQ2.size(); i++) {
			filterQ2.at(i) = log(filterQ2.at(i));
		}

		nbins = 200;
		TH2D* pmtLogQvQ = new TH2D("LogQvQ", ((std::string)"ln(Q2) vs ln(Q1) for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1 * 1.1, nbins, minQ2, maxQ2 * 1.1);
		nbins = nbins / 2;
		TProfile* LogQvQProfile = new TProfile("Log QvQ Profile", ((std::string)"ln(Q2) vs ln(Q1) for hit PMTs (" + fn.substr(25, fn.find(")") - 1)).c_str(), nbins, minQ1, maxQ1 * 1.1, minQ2, maxQ2 * 1.1);

		for (int i = 0; i < filterQ1.size(); i++) {
			pmtLogQvQ->Fill(filterQ1.at(i), filterQ2.at(i));
			LogQvQProfile->Fill(filterQ1.at(i), filterQ2.at(i));
		}

		TCanvas* c5 = new TCanvas("c5", "Fifth canvas", 1920 * 2 * win_scale, 1080 * win_scale);
		c5->Divide(2, 1);
		c5->cd(1); pmtLogQvQ->SetContour(100); pmtLogQvQ->GetXaxis()->SetTitle("Log Q for hit PMT, ln(Q1)"); pmtLogQvQ->GetYaxis()->SetTitle("Log Q for other hit PMT, ln(Q2)"); pmtLogQvQ->Draw("COLZ");
		c5->cd(2); LogQvQProfile->GetXaxis()->SetTitle("Log Q for hit PMT, ln(Q1)"); LogQvQProfile->GetYaxis()->SetTitle("Log Q for other hit PMT, ln(Q2)"); LogQvQProfile->Draw();
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



		txtfile << "------------------------------------------\n" << std::endl;
	}
	txtfile.close();


	return 0;
}