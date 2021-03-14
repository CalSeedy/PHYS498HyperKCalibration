#include <ctime>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
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
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <Math/MinimizerOptions.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

#include "./getopt.h"


double func(double val, double m, double c) {
	return m * val + c;
}

int main(int argc, char** argv) {
	char* filename = NULL;
	char* outfilename = NULL;
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
		filename = "./Output/test/wcsim_hybrid\(10GeVmu\)-out.root";
	}

	file = new TFile(filename, "read");

	if (!file->IsOpen()) {
		std::cout << "Error, could not open input file: " << filename << std::endl;
		return -1;
	}

	TProfile* prof = (TProfile*)file->Get("QvQ Profile");
	TProfile* logProf = (TProfile*)file->Get("Log QvQ Profile");

	// Normal Charge Profile
	double ymin = prof->GetYmin();
	double ymax = prof->GetYmax();
	int under = prof->GetBinEntries(0);
	int over = prof->GetBinEntries(prof->GetNbinsX() + 1);
	std::cout << "Profile| Underflows: " << under << ", Overflows: " << over << std::endl;
	double start = prof->GetBinCenter(1);
	double end = prof->GetBinCenter(prof->GetNbinsX());
	double w = prof->GetBinCenter(2) - prof->GetBinCenter(1);
	std::cout << "Start: " << start << ", End: " << end << std::endl;
	std::cout << "Min: " << ymin << ", Max: " << ymax << std::endl;
	auto linearFit = new TF1("linear", "[0]*x + [1]", start, end + (w / 2));
	linearFit->SetLineColor(kBlack);

	auto linearConf95 = new TGraphErrors(prof->GetNbinsX());
	for (int i = 1; i <= logProf->GetNbinsX(); i++) {
		linearConf95->SetPoint(i, prof->GetBinCenter(i), 0.0);
	}
	auto linearConf99 = new TGraphErrors(prof->GetNbinsX());
	for (int i = 1; i <= logProf->GetNbinsX(); i++) {
		linearConf99->SetPoint(i, prof->GetBinCenter(i), 0.0);
	}

	prof->Fit(linearFit, "R", "SAME");
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(linearConf95);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(linearConf99, 0.99);

	// Log Charge Profile
	ymin = logProf->GetYmin();
	ymax = logProf->GetYmax();
	under = logProf->GetBinEntries(0);
	over = logProf->GetBinEntries(logProf->GetNbinsX() + 1);
	std::cout << "Profile| Underflows: " << under << ", Overflows: " << over << std::endl;
	start = logProf->GetBinCenter(1);
	end = logProf->GetBinCenter(logProf->GetNbinsX());
	w = logProf->GetBinCenter(2) - logProf->GetBinCenter(1);
	std::cout << "Start: " << start << ", End: " << end << std::endl;
	std::cout << "Min: " << ymin << ", Max: " << ymax << std::endl;
	auto piecewiseFit = new TF1("piecewise", "(x<[0]) ? [1] : [2]*x+[1]-[2]*[0]", start, end + (w / 2));
	piecewiseFit->SetLineColor(kBlack);
	piecewiseFit->SetParameter(0, end / 4.0);
	piecewiseFit->SetParLimits(0, 1.0, 4.0);
	piecewiseFit->SetParameter(1, 1.0);
	piecewiseFit->SetParLimits(1, 0.0, 2.0);
	/*
	piecewiseFit->SetParameter(2, 1.1);
	piecewiseFit->SetParLimits(2, 0.9, 1.1);
	*/
	auto piecewiseConf95 = new TGraphErrors(logProf->GetNbinsX());
	for (int i = 1; i <= logProf->GetNbinsX(); i++) {
		piecewiseConf95->SetPoint(i, logProf->GetBinCenter(i), 0.0);
	}
	auto piecewiseConf99 = new TGraphErrors(logProf->GetNbinsX());
	for (int i = 1; i <= logProf->GetNbinsX(); i++) {
		piecewiseConf99->SetPoint(i, logProf->GetBinCenter(i), 0.0);
	}

	logProf->Fit(piecewiseFit, "EWR", "SAME");
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(piecewiseConf95);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(piecewiseConf99, 0.99);

	// Plot everything on the canvas
	auto c1 = new TCanvas("can", "Fit", 1920, 2*1080);
	c1->Divide(1, 2);

	c1->cd(1);
	prof->Draw();
	linearConf99->SetFillStyle(3001);
	linearConf99->SetFillColor(kGreen);
	linearConf99->Draw("SAME E4");

	linearConf95->SetFillStyle(3001);
	linearConf95->SetFillColor(kRed);
	linearConf95->Draw("SAME E4");
	
	TLegend* leg1 = new TLegend(0.12, 0.7, 0.3, 0.85);
	leg1->AddEntry(prof, "Simulation Data", "lep");
	leg1->AddEntry(linearFit, "Fit", "l");
	leg1->AddEntry(linearConf95, "95% CI", "f");
	leg1->AddEntry(linearConf99, "99% CI", "f");
	leg1->Draw();
	
	c1->cd(2);

	logProf->Draw();
	piecewiseConf99->SetFillStyle(3001);
	piecewiseConf99->SetFillColor(kGreen);
	piecewiseConf99->Draw("SAME E4");

	piecewiseConf95->SetFillStyle(3001);
	piecewiseConf95->SetFillColor(kRed);
	piecewiseConf95->Draw("SAME E4");

	TLegend* leg2 = new TLegend(0.12, 0.7, 0.3, 0.85);
	leg2->AddEntry(logProf, "Simulation Data", "lep");
	leg2->AddEntry(piecewiseFit, "Fit", "l");
	leg2->AddEntry(piecewiseConf95, "95% CI", "f");
	leg2->AddEntry(piecewiseConf99, "99% CI", "f");
	leg2->Draw();


	c1->Modified();
	c1->Update();


	std::string fn(filename);
	size_t a = fn.rfind("/"); // path/to/folder/file.ext    -- find last slash (/)
	std::string fn1 = fn.substr(a + 1, fn.length() - a - 6); // attempt to grab "file"
	std::cout << fn1 << std::endl;

	c1->Print((fn1 + (std::string)"-fit.png").c_str());
	
	TFile* out;
	if (outfilename == NULL) out = new TFile((fn1 + (std::string)"-rootfile.root").c_str(), "RECREATE");
	else out = new TFile(outfilename, "RECREATE");

	prof->Write();
	logProf->Write();
	
	linearFit->Write();
	piecewiseFit->Write();
	piecewiseConf95->Write();
	piecewiseConf99->Write();
	linearConf95->Write();
	linearConf99->Write();

	out->Close();
	
	return 0;
}