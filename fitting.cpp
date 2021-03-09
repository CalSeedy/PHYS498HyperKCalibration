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

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

#include "./getopt.h"


double func(double val, double m, double c) {
	return m * val + c;
}

int main(int argc, char** argv) {
	char* filename = NULL;
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
		filename = "./Output/test/wcsim_hybrid\(10GeVmu\)-out.root";
	}

	file = new TFile(filename, "read");

	if (!file->IsOpen()) {
		std::cout << "Error, could not open input file: " << filename << std::endl;
		return -1;
	}


	TProfile* logProf = (TProfile*)file->Get("Log QvQ Profile");
	TProfile* prof = (TProfile*)file->Get("QvQ Profile");


	double start = logProf->GetBinCenter(0);
	double end = logProf->GetBinCenter(logProf->GetNbinsX());
	std::cout << "Start: " << start << ", End: " << end << std::endl;
	float w = logProf->GetBarWidth();
	double xstart = start - (w / 2);
	double xend = end + (w / 2);


	TF1* f1 = new TF1("f1", "(x<[0])*[1] + (x>[0])*([2]*x + [1] - [2]*[0])");
	auto tmpFit = logProf->Fit(f1);
	TGraphErrors* ge95 = new TGraphErrors(logProf->GetNbinsX());
	for (int i = 0; i < logProf->GetNbinsX(); i++) ge95->SetPoint(i, logProf->GetBinCenter(i), 0);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(ge95, 0.95);
	TGraphErrors* ge99 = new TGraphErrors(logProf->GetNbinsX());
	for (int i = 0; i < logProf->GetNbinsX(); i++) ge99->SetPoint(i, logProf->GetBinCenter(i), 0);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(ge99, 0.99);


	TF1* f2 = new TF1("f2", "(x>[0])*([1]*x) + (x<[0])*[3]");
	f2->SetParameter(0, exp(f1->GetParameter(0)));
	f2->SetParameter(3, exp(f1->GetParameter(1)));
	f2->SetParameter(1, exp(f1->GetParameter(1) - f1->GetParameter(0) * f1->GetParameter(2)));
	f2->SetParameter(2, f1->GetParameter(2));
	auto fit = prof->Fit(f2);
	TGraphErrors* ge295 = new TGraphErrors(prof->GetNbinsX());
	for (int i = 0; i < prof->GetNbinsX(); i++) ge295->SetPoint(i, prof->GetBinCenter(i), 0);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(ge295, 0.95);
	TGraphErrors* ge299 = new TGraphErrors(prof->GetNbinsX());
	for (int i = 0; i < prof->GetNbinsX(); i++) ge299->SetPoint(i, prof->GetBinCenter(i), 0);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(ge299, 0.99);

	//gStyle->SetHatchesSpacing(0.1);
	TCanvas* can = new TCanvas("can", "canvas", 3*1920, 2*1080);
	can->Divide(1, 2);

	prof->SetAxisRange(0, 500, "X");
	f2->SetRange(0, 500);
	f2->SetLineColor(kBlack);
	ge295->GetXaxis()->SetRangeUser(0, 500);
	ge299->GetXaxis()->SetRangeUser(0, 500);

	can->cd(1);

	ge299->SetFillColor(kRed);
	ge299->SetFillStyle(3002);
	ge299->Draw("a4+");

	ge295->SetFillColor(kYellow);
	ge295->SetFillStyle(3002);
	ge295->Draw("a4+ SAME");

	prof->Draw("SAME");
	f2->Draw("SAME");

	TLegend* legend = new TLegend(0.1, 0.7, 0.48, 0.9);
	legend->SetHeader("Legend");
	legend->AddEntry(prof, "Simulation Data", "pel");
	legend->AddEntry(f2, "Fit", "l");
	legend->AddEntry(ge295, "95% CI", "f");
	legend->AddEntry(ge299, "99% CI", "f");
	legend->Draw("SAME");

	can->cd(2);
	logProf->SetAxisRange(xstart, xend, "X");
	f1->SetRange(xstart, xend);
	f1->SetLineColor(kBlack);
	ge95->GetXaxis()->SetRangeUser(xstart, xend);
	ge99->GetXaxis()->SetRangeUser(xstart, xend);

	ge99->SetFillColor(kRed);
	ge99->SetFillStyle(3002);
	ge99->Draw("a4+");

	ge95->SetFillColor(kYellow);
	ge95->SetFillStyle(3002);
	ge95->Draw("a4+ SAME");

	logProf->Draw("SAME");
	f1->Draw("SAME");
	
	TLegend* legend2 = new TLegend(0.1, 0.7, 0.48, 0.9);
	legend2->SetHeader("Legend2"); 				
	legend2->AddEntry(logProf, "Simulation Data", "pel");
	legend2->AddEntry(f1, "Fit", "l");
	legend2->AddEntry(ge95, "95% CI", "f");
	legend2->AddEntry(ge99, "99% CI", "f");
	legend2->Draw("SAME");

	can->Print("fit.png");

	return 0;
}