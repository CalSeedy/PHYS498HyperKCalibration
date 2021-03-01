#include <ctime>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <TApplication.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
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

#include "./getopt.h"


double poly2(double val, double A, double B, double C) {
	return A * val * val + B * val + C;
}


int main(int argc, char** argv) {
	char* func1 = "A";
	char* func2 = "m*x + c";
	char* filename = "./Output/test/wcsim_hybrid\(500MeVe-\)-out.root";
	char* outfilename = "fit.root";
	bool verbose = false;
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

	TFile* file;
	// Open the file
	if (filename == NULL) {
		filename = "./Output/test/wcsim_hybrid\(500MeVe-\)-out.root";
	}

	file = new TFile(filename, "read");

	if (!file->IsOpen()) {
		std::cout << "Error, could not open input file: " << filename << std::endl;
		return -1;
	}

	file->ls();

	TProfile* prof = (TProfile*)file->Get("Log QvQ Profile");
	/*
	std::vector<double> x, y, dy;
	int numBins = prof->GetNbinsX();
	for (int i = 0; i < numBins; i++) {
		x.push_back(prof->GetBinCenter(i));
		y.push_back(prof->GetBinContent(i));
		dy.push_back(prof->GetBinError(i));
	}
	*/

	double start = prof->GetBinCenter(0);
	double end = prof->GetBinCenter(-1);
	float w = prof->GetBarWidth();
	double xstart = start - (w / 2);
	double xend = end + (w / 2);

	auto fit = prof->Fit("pol2");
	TF1* g = (TF1*)prof->GetListOfFunctions()->FindObject("pol2");
	
	//TMatrixDSym cov = fit->GetCovarianceMatrix();

	double ex2 = g->GetParameter(0);
	double ex = g->GetParameter(1);
	double con = g->GetParameter(2);

	std::cout << "A: " << ex2 << ", B: " << ex << ", C: " << con << std::endl;


	int N = 200;
	double xs[N];
	double ys[N];
	double diff = (xend - xstart) / N;
	for (int i = 0; i < N; i++) {
		xs[i] = xstart + i * diff;
		ys[i] = poly2(xs[i], ex2, ex, con);
	}

	TGraph* graph = new TGraph(N, xs, ys);

	TCanvas* can = new TCanvas("can", "canvas", 3*1920, 2*1080);
	prof->Draw();
	graph->Draw("SAME");


	can->Print("fit.png");


	return 0;
}