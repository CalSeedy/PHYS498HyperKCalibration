#include <algorithm>
#include <ctime>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <stdexcept>
#include <tuple>
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
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TMultiGraph.h>
#include <TPaletteAxis.h>
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

std::vector<std::string> getFiles() {
	std::vector<std::string> out;
	
	//out.push_back("Output/HK/WCSim_hk_20bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	//out.push_back("Output/HK/WCSim_hk_20bl_4200hz_r14374_0hz_435nmCol_middle-analysis.root");
	//out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_3kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	//out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_3kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle-analysis.root");
	//out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_5kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	//out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_5kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle-analysis.root");
	//out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_10kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	//out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_10kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle-analysis.root");

	//out.push_back("Output/test/wcsim_hybrid\(10GeVe-\)_bias-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(10GeVe-\)-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(10GeVmu\)_bias-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(10GeVmu\)-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(500MeVe-\)_bias-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(500MeVe-\)-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(500MeVmu\)_bias-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(500MeVmu\)-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(10MeVe-\)_bias-analysis.root");
	//out.push_back("Output/test/wcsim_hybrid\(10MeVe-\)-analysis.root");

	//out.push_back("Output/atm_emu/wcsim_atm_emu_20kBL_43_bias-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_20kBL_43-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_40kBL_1011_bias-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_40kBL_1011-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_3kmPMT_1000_bias-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_3kmPMT_1000-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_5kmPMT_666_bias-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_5kmPMT_666-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_10kmPMT_888_bias-analysis.root");
	//out.push_back("Output/atm_emu/wcsim_atm_emu_10kmPMT_888-analysis.root");
	
	//out.push_back("Output/atm_tau/wcsim_atm_tau_20kBL_1000_bias-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_20kBL_1000-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_40kBL_1000_bias-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_40kBL_1000-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_3kmPMT_1000_bias-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_3kmPMT_1000-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_5kmPMT_1000_bias-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_5kmPMT_1000-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_10kmPMT_1000_bias-analysis.root");
	//out.push_back("Output/atm_tau/wcsim_atm_tau_10kmPMT_1000-analysis.root");
	
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_20kBL_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_20kBL_0-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_40kBL_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_40kBL_0-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_5kmPMT_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_5kmPMT_0-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_10kmPMT_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_10kmPMT_0-analysis.root");
	
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_20kBL_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_20kBL_0-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_40kBL_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_40kBL_0-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_5kmPMT_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_5kmPMT_0-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_10kmPMT_0_bias-analysis.root");
	//out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_10kmPMT_0-analysis.root");
	///////////////////////////////////////////////////////////////////////////////////////

	out.push_back("Charge_Data/wcsim_hybrid\(10GeVe-\)_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(10GeVe-\)-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(10GeVmu\)_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(10GeVmu\)-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(500MeVe-\)_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(500MeVe-\)-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(500MeVmu\)_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(500MeVmu\)-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(10MeVe-\)_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_hybrid\(10MeVe-\)-chargeData.root");
	
	out.push_back("Charge_Data/WCSim_hk_20bl_4200hz_r14374_0hz_435nmCol_middle_bias-chargeData.root");
	out.push_back("Charge_Data/WCSim_hk_20bl_4200hz_r14374_0hz_435nmCol_middle-chargeData.root");
	out.push_back("Charge_Data/WCSim_hkhybridmpmt10pc_20bl_3kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-chargeData.root");
	out.push_back("Charge_Data/WCSim_hkhybridmpmt10pc_20bl_3kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle-chargeData.root");
	out.push_back("Charge_Data/WCSim_hkhybridmpmt10pc_20bl_5kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-chargeData.root");
	out.push_back("Charge_Data/WCSim_hkhybridmpmt10pc_20bl_5kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle-chargeData.root");
	out.push_back("Charge_Data/WCSim_hkhybridmpmt10pc_20bl_10kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-chargeData.root");
	out.push_back("Charge_Data/WCSim_hkhybridmpmt10pc_20bl_10kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle-chargeData.root");
	
	out.push_back("Charge_Data/wcsim_atm_emu_20kBL_43_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_20kBL_43-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_40kBL_1011_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_40kBL_1011-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_3kmPMT_1000_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_3kmPMT_1000-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_5kmPMT_666_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_5kmPMT_666-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_10kmPMT_888_bias-chargeData.root");
	out.push_back("Charge_Data/wcsim_atm_emu_10kmPMT_888-chargeData.root");
	return out;
}

TGraphErrors* getInterval(TProfile* prof, TF1* func, double CI = 0.95) {

	auto confInt = new TGraphErrors(prof->GetNbinsX());
	for (int i = 1; i < prof->GetNbinsX(); i++) {
		confInt->SetPoint(i, prof->GetBinCenter(i), 0.0);
	}

	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(confInt, CI);
	return confInt;

}

TH1D divideProfiles(int n, TProfile* p1, TProfile* p2, bool verbose = false) {
	if (p1->GetNbinsX() != p2->GetNbinsX()) { 
		std::cerr << "Bins do not match.\n";
		return TH1D("out", "Failed", n, 0.0, 1.0);
	}

	double w1 = p1->GetBinCenter(n) - p1->GetBinCenter(n - 1);
	double w2 = p2->GetBinCenter(n) - p2->GetBinCenter(n - 1);
	double maxX = std::max(p1->GetBinCenter(n) + w1 / 2.0, p2->GetBinCenter(n) + w2 / 2.0);
	double minX = std::min(p1->GetBinCenter(0) - w1 / 2.0, p2->GetBinCenter(0) - w2 / 2.0);

	TH1D out("ratios", "Extracting Bias Function", n, minX, maxX);
	for (int i = 0; i <= n; i++) {
		double count1 = p1->GetBinContent(i);
		double count2 = p2->GetBinContent(i);
		double dx = sqrt(pow(w1 / 2.0, 2) + pow(w2 / 2.0, 2));
		std::cout << "X: " << count1 << ", Y: " << count2 << std::endl;
		if (count2 == 0.0000) {
			std::cout << "Division by zero... skipping bin." << std::endl;
			continue;
		}

		double y = count1 / count2;
		//double dy = y * sqrt(pow(p1->GetBinError(i) / count1, 2) + pow(p2->GetBinError(i) / count2, 2));
		double bErr = p1->GetBinError(i);
		double ubErr = p2->GetBinError(i);
		out.SetBinContent(i, y);
		out.SetBinError(i, sqrt(pow(bErr / count2, 2) + pow(ubErr * count1 / (count2 * count2), 2)));
		if (verbose) std::cout << "X: " << count1 << ", Y: " << count2 << " | X/Y: " << y << " +- (" << sqrt(pow(bErr / count2, 2) + pow(ubErr * count1 / (count2 * count2), 2)) << ")\n";
	}
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

std::tuple<int, double> chi2(TH1D* data, double xlow, double xhigh) {
	double chiSq = 0.0;
	int c = 0;
	for (int i = 1; i < data->GetNbinsX(); i++) {
		double x = data->GetBinCenter(i);
		if ((x <= xhigh) && (x >= xlow)) {
			double val = data->GetBinContent(i);
			double err = data->GetBinError(i);
			if (err < 1e-5) continue;
			chiSq += pow(val - bias_function(x), 2) / pow(err, 2);
			c++;
		}
	}
	std::tuple<int, double> out{ c, chiSq };
	return out;
}

int main(int argc, char* argv[]) {
	char* infilename = NULL;
	char* outfilename = NULL;
	bool verbose = false;
	int c = -1;
	while ((c = getopt(argc, argv, ":f:o:v")) != -1) {//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
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
		}
	}
	//"[&](double *x, double *p){ return p[0]*f1(x) + p[1]*f2(x); }"
	
	std::vector<std::string> files;
	// Open the file
	if (infilename == NULL) {
		files = getFiles();
	}
	else {
		files.push_back((std::string)infilename);
	}
	TFile *fileBiased, *fileUnbiased;
	TProfile *biased, *unbiased, *biasedLog, *unbiasedLog;
	TH1D *p1, *p2, *dataPtr;
	TAxis *axBias, *axUnbias;
	TLegend *leg1, *leg2;
	TCanvas *c1, *c2;

	for (int i = 0; i < files.size() / 2; i++) {
		fileBiased = new TFile(files.at(2 * i).c_str(), "read");
		TTree* treeB = (TTree*)fileBiased->Get("charge data");
		if (!fileBiased->IsOpen()) {
			std::cout << "Error, could not open input file: " << files.at(2 * i) << std::endl;
			return -1;
		}
		char* args;
		int nbins = 50;
		TCanvas* testCan = new TCanvas("test", "test", 1920,1080);

		double minY = treeB->GetMinimum("q1");
		double maxY = treeB->GetMaximum("q1");
		double minR = 0.25; // treeB->GetMinimum("ratio");
		double maxR = 4.0;  // treeB->GetMaximum("ratio");
		asprintf(&args, "ratio:q1>>altRatioQvQProfile(%i, %g, %g)", nbins, 0.0, 200.0);
		treeB->Draw(args, "q1<200", "PROF");
		biased = ((TProfile*)gPad->GetPrimitive("altRatioQvQProfile"))->Clone();
		testCan->Clear();
		delete args;
		asprintf(&args, "ratio:q1>>altRatioQvQProfile(%i, %g, %g)", nbins, 0.0, maxY);
		treeB->Draw(args, "", "PROF");
		biasedLog = ((TProfile*)gPad->GetPrimitive("altRatioQvQProfile"))->Clone();

		fileUnbiased = new TFile(files.at(2 * i + 1).c_str(), "read");
		if (!fileUnbiased->IsOpen()) {
			std::cout << "Error, could not open input file: " << files.at(2 * i + 1) << std::endl;
			return -1;
		}
		TTree* treeU = (TTree*)fileUnbiased->Get("charge data");
		testCan->Clear();
		delete args;
		asprintf(&args, "ratio:q1>>altRatioQvQProfile(%i, %g, %g)", nbins, 0.0, 200.0);
		treeU->Draw(args, "q1<200", "PROF");
		unbiased = ((TProfile*)gPad->GetPrimitive("altRatioQvQProfile"))->Clone();
		testCan->Clear();
		delete args;
		asprintf(&args, "ratio:q1>>altRatioQvQProfile(%i, %g, %g)", nbins, 0.0, maxY);
		treeU->Draw(args, "", "PROF");
		unbiasedLog = ((TProfile*)gPad->GetPrimitive("altRatioQvQProfile"))->Clone();
		testCan->Close();

		auto c1 = new TCanvas("can", "ratios", 1920*2, 1080);
		c1->Divide(2, 1);
		c1->cd(1);
		biased->Draw("E1 P"); biased->SetMarkerStyle(12); biased->SetLineWidth(5);
		unbiased->Draw("SAME E1 P"); unbiased->SetMarkerStyle(12); unbiased->SetMarkerColor(kRed); unbiased->SetLineWidth(5); unbiased->SetLineColor(kRed);
		biased->GetXaxis()->SetTitle("Charge on hit PMT, Q1"); biased->GetYaxis()->SetTitle("Charge Ratio, Q1/Q2");
		leg1 = new TLegend(0.12, 0.7, 0.3, 0.85);
		leg1->AddEntry(biased, "Biased Data", "EP");
		leg1->AddEntry(unbiased, "Unbiased Data", "EP");
		leg1->Draw();

		c1->cd(2);
		gPad->SetLogx();
		biasedLog->Draw("E1 P"); biasedLog->SetMarkerStyle(12); biasedLog->SetLineWidth(5);
		unbiasedLog->Draw("SAME E1 P"); unbiasedLog->SetMarkerStyle(12); unbiasedLog->SetMarkerColor(kRed); unbiasedLog->SetLineWidth(5);
		unbiasedLog->SetLineColor(kRed);
		biasedLog->GetXaxis()->SetTitle("Charge on hit PMT, ln(Q1)"); biasedLog->GetYaxis()->SetTitle("Charge Ratio, Q1/Q2");
		leg1->Draw();
		c1->Draw();

		TH1D* biasFunc = new TH1D("biasFunc", "bias Func", (int)(1.1*maxY) * 10, 0.0, maxY);
		for (int n = 1; n <= biasFunc->GetNbinsX(); ++n) {
			biasFunc->SetBinContent(n, bias_function(biasFunc->GetBinCenter(n)));
		}

		nbins = biased->GetNbinsX();
		TH1D data = divideProfiles(nbins, biased, unbiased, verbose);
		data.SetDrawOption("E1 P"); data.SetStats(0); data.SetMarkerStyle(20); data.SetMarkerSize(2); data.SetLineWidth(5);
		data.GetXaxis()->SetTitle("Charge on hit PMT, Q1"); data.GetYaxis()->SetTitle("biased / unbiased");
		auto linearFit = new TF1("linear", "[0]*x + [1]", 0, 200); linearFit->SetLineColor(kRed); linearFit->SetLineWidth(8);
		dataPtr = &data;
		dataPtr->Fit(linearFit, "EX0 R C ROB N");
		TH1D dataLog = divideProfiles(nbins, biasedLog, unbiasedLog, verbose);
		dataLog.SetDrawOption("E1 P"); dataLog.SetStats(0); dataLog.SetMarkerStyle(20); dataLog.SetMarkerSize(2); dataLog.SetLineWidth(5);
		dataLog.GetXaxis()->SetTitle("Charge on hit PMT, ln(Q1)"); dataLog.GetYaxis()->SetTitle("biased / unbiased");

		c2 = new TCanvas("can2", "biases", 1920*2, 1080);
		c2->Divide(2, 1);
		c2->cd(1);
		data.GetXaxis()->SetLimits(0.0, 200.0);
		data.GetYaxis()->SetRangeUser(0.25, 2.0);
		dataPtr->Draw();
		biasFunc->Draw("SAME L"); biasFunc->SetLineWidth(8); biasFunc->SetLineColor(kBlack);
		linearFit->Draw("SAME");
		leg2 = new TLegend(0.12, 0.75, 0.25, 0.85);
		leg2->AddEntry(dataPtr, "Data", "LEP");
		leg2->AddEntry(biasFunc, "Bias", "L");
		leg2->AddEntry(linearFit, "Linear Fit", "L");
		leg2->Draw("SAME");

		c2->cd(2);
		gPad->SetLogx();
		(&dataLog)->Draw();
		biasFunc->Draw("SAME L");
		TLegend* leg2L = new TLegend(0.12, 0.75, 0.25, 0.85);
		leg2L->AddEntry((&dataLog), "Data", "LEP");
		leg2L->AddEntry(biasFunc, "Bias", "L");
		leg2L->Draw("SAME");
		c2->Draw();

		std::string filename = files.at(2 * i + 1);

		std::string fn = filename.substr(filename.rfind("/") + 1, filename.rfind("-") - filename.rfind("/") - 1);
		c1->Print(((std::string)"./" + fn + (std::string)"-profiles.png").c_str());
		c2->Print(((std::string)"./" + fn + (std::string)"-biasComparison.png").c_str());

		TFile* outfile;
		std::ofstream fittingParams;
		if (outfilename == NULL) {
			outfile = new TFile((fn + (std::string)"-ratioPlots.root").c_str(), "RECREATE");
			fittingParams.open((fn + (std::string)"-params.txt").c_str());
		}
		else {
			outfile = new TFile(outfilename, "RECREATE");
			fittingParams.open(outfilename);
		}

		char* line;
		fittingParams << "~~~~~~~~~~~Fitting Parameters~~~~~~~~~~~\n--Linear Fit--\n";
		asprintf(&line, "\tGradient: %f +/- %f\n\tOffset: %f +/- %f\n\tChi2: %f (NDF=%i)\n", linearFit->GetParameter(0), linearFit->GetParError(0), linearFit->GetParameter(1), linearFit->GetParError(1), linearFit->GetChisquare(), linearFit->GetNDF());
		fittingParams << line << std::endl;
		delete line;
		std::tuple<int, double> chisq = chi2(dataPtr, 0.0, maxY);
		fittingParams << "--Data to Bias Stats--\n";
		asprintf(&line, "\tChi2: %f (NDF=%i)\n", std::get<1>(chisq), std::get<0>(chisq));
		fittingParams << line << std::endl;
		biased->Write();
		unbiased->Write();
		data.Write();
		biasFunc->Write();
		linearFit->Write();
		c1->Write();
		c2->Write();
		fittingParams.close();
		outfile->Close();

		delete line, biasFunc, fileBiased, fileUnbiased, biased, unbiased, biasedLog, unbiasedLog, leg1, leg2, c1, c2, axBias, axUnbias, p1, p2, dataPtr, testCan;
	}
	return 0;
}