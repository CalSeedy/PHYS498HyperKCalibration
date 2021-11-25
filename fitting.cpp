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

std::vector<std::string> getFiles() {
	std::vector<std::string> out;
	out.push_back("Output/HK/WCSim_hk_20bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_3kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_5kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");
	out.push_back("Output/HK/WCSim_hkhybridmpmt10pc_20bl_10kmpmt_bl_4200hz_r14374_0hz_435nmCol_middle_bias-analysis.root");

	out.push_back("Output/test/wcsim_hybrid\(10GeVe-\)_bias-analysis.root");
	out.push_back("Output/test/wcsim_hybrid\(10GeVmu\)_bias-analysis.root");
	out.push_back("Output/test/wcsim_hybrid\(500MeVe-\)_bias-analysis.root");
	out.push_back("Output/test/wcsim_hybrid\(500MeVmu\)_bias-analysis.root");
	out.push_back("Output/test/wcsim_hybrid\(10MeVe-\)_bias-analysis.root");
	
	out.push_back("Output/atm_emu/wcsim_atm_emu_20kBL_43_bias-analysis.root");
	out.push_back("Output/atm_emu/wcsim_atm_emu_40kBL_1011_bias-analysis.root");
	out.push_back("Output/atm_emu/wcsim_atm_emu_3kmPMT_1000_bias-analysis.root");
	out.push_back("Output/atm_emu/wcsim_atm_emu_5kmPMT_666_bias-analysis.root");
	out.push_back("Output/atm_emu/wcsim_atm_emu_10kmPMT_888_bias-analysis.root");
	
	out.push_back("Output/atm_tau/wcsim_atm_tau_20kBL_1000_bias-analysis.root");
	out.push_back("Output/atm_tau/wcsim_atm_tau_40kBL_1000_bias-analysis.root");
	out.push_back("Output/atm_tau/wcsim_atm_tau_3kmPMT_1000_bias-analysis.root");
	out.push_back("Output/atm_tau/wcsim_atm_tau_5kmPMT_1000_bias-analysis.root");
	out.push_back("Output/atm_tau/wcsim_atm_tau_10kmPMT_1000_bias-analysis.root");
	
	out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_20kBL_0_bias-analysis.root");
	out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_40kBL_0_bias-analysis.root");
	out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_5kmPMT_0_bias-analysis.root");
	out.push_back("Output/Rn/scalingA/wcsim_Rn_scalingA_10kmPMT_0_bias-analysis.root");
	
	out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_20kBL_0_bias-analysis.root");
	out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_40kBL_0_bias-analysis.root");
	out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_5kmPMT_0_bias-analysis.root");
	out.push_back("Output/Rn/scalingB/wcsim_Rn_scalingA_10kmPMT_0_bias-analysis.root");

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


int main(int argc, char** argv) {
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

	std::vector<std::string> files;
	// Open the file
	if (infilename == NULL) {
		files = getFiles();
	}
	else {
		files.push_back((std::string)infilename);
	}
	for (std::string filename : files) {
		TFile* file = new TFile(filename.c_str(), "read");

		if (!file->IsOpen()) {
			std::cout << "Error, could not open input file: " << filename << std::endl;
			return -1;
		}

		TProfile* prof = (TProfile*)file->Get("QvQ Profile");
		TProfile* logProf = (TProfile*)file->Get("Log QvQ Profile");
		

		// Log Charge Profile
		double ymin = logProf->GetYmin();
		double ymax = logProf->GetYmax();
		int under = logProf->GetBinEntries(0);
		int over = logProf->GetBinEntries(logProf->GetNbinsX() + 1);
		std::cout << "Profile| Underflows: " << under << ", Overflows: " << over << std::endl;
		double start = logProf->GetBinCenter(1);
		double end = logProf->GetBinCenter(logProf->GetNbinsX());
		double w = logProf->GetBinCenter(2) - logProf->GetBinCenter(1);
		std::cout << "Start: " << start << ", End: " << end << std::endl;
		std::cout << "Min: " << ymin << ", Max: " << ymax << std::endl;
		auto piecewiseFit = new TF1("piecewise", "(x<[0]) ? [1] : [2]*x+[1]-[2]*[0]", start, end - (5*w / 2));
		piecewiseFit->SetLineColor(kBlack);
		piecewiseFit->SetLineWidth(5.0);
		piecewiseFit->SetParameter(0, end / 6.0);
		piecewiseFit->SetParLimits(0, 1.0, 4.0);
		piecewiseFit->SetParameter(1, 1.0);
		piecewiseFit->SetParLimits(1, 0.0, 2.0);
		piecewiseFit->SetParameter(2, 1.0);
		piecewiseFit->SetParLimits(2, 0.9, 1.1);

		logProf->Fit(piecewiseFit, "IR", "SAME");
		//auto piecewiseConf99 = getInterval(logProf, piecewiseFit, 0.99);
		//auto piecewiseConf95 = getInterval(logProf, piecewiseFit);

		// Normal Charge Profile
		ymin = prof->GetYmin();
		ymax = prof->GetYmax();
		under = prof->GetBinEntries(0);
		over = prof->GetBinEntries(prof->GetNbinsX() + 1);
		std::cout << "Profile| Underflows: " << under << ", Overflows: " << over << std::endl;
		start = prof->GetBinCenter(1);
		end = prof->GetBinCenter(prof->GetNbinsX());
		w = prof->GetBinCenter(2) - prof->GetBinCenter(1);
		std::cout << "Start: " << start << ", End: " << end << std::endl;
		std::cout << "Min: " << ymin << ", Max: " << ymax << std::endl;
		auto linearFit = new TF1("linear", "[0]*x + [1]", start, end + (5* w / 2));
		double pwf0 = piecewiseFit->GetParameter(0);
		double pwf1 = piecewiseFit->GetParameter(1);
		double pwf2 = piecewiseFit->GetParameter(2);
		double pwfe0 = piecewiseFit->GetParError(0);
		double pwfe1 = piecewiseFit->GetParError(1);
		double pwfe2 = piecewiseFit->GetParError(2);
		linearFit->SetLineColor(kBlack);
		linearFit->SetLineWidth(5.0);
		linearFit->SetParameter(0, pwf1/(pwf2*pwf0));
		double pwfError = (pwf1 / (pwf2 * pwf0 * pwf0)) * pwfe0 + (1 / (pwf2 * pwf0 )) * pwfe1 + (pwf1 / (pwf2 * pwf2 * pwf0)) * pwfe2;
		linearFit->SetParLimits(0, pwf1 / (pwf2 * pwf0) - pwfError , pwf1 / (pwf2 * pwf0) + pwfError);
		linearFit->SetParameter(1, 1.0);
		linearFit->SetParLimits(1, 0.0, 2.0);

		prof->Fit(linearFit, "CER", "SAME");
		//auto linearConf99 = getInterval(prof,linearFit, 0.99);
		//auto linearConf95 = getInterval(prof,linearFit);

		

		// Plot everything on the canvas
		auto c1 = new TCanvas("can", "Fit", 3840, 2 * 2160);
		c1->Divide(1, 2);

		c1->cd(1);
		prof->Draw();
		//linearConf99->SetFillStyle(3001);
		//linearConf99->SetFillColor(kGreen);
		//linearConf99->Draw("SAME E4");
		//
		//linearConf95->SetFillStyle(3001);
		//linearConf95->SetFillColor(kRed);
		//linearConf95->Draw("SAME E4");

		TLegend* leg1 = new TLegend(0.12, 0.7, 0.3, 0.85);
		leg1->AddEntry(prof, "Simulation Data", "lep");
		leg1->AddEntry(linearFit, "Fit", "l");
		//leg1->AddEntry(linearConf95, "95% CI", "f");
		//leg1->AddEntry(linearConf99, "99% CI", "f");
		leg1->Draw();

		c1->cd(2);

		logProf->Draw();
		//piecewiseConf99->SetFillStyle(3001);
		//piecewiseConf99->SetFillColor(kGreen);
		//piecewiseConf99->Draw("SAME E4");
		//
		//piecewiseConf95->SetFillStyle(3001);
		//piecewiseConf95->SetFillColor(kRed);
		//piecewiseConf95->Draw("SAME E4");

		TLegend* leg2 = new TLegend(0.12, 0.7, 0.3, 0.85);
		leg2->AddEntry(logProf, "Simulation Data", "lep");
		leg2->AddEntry(piecewiseFit, "Fit", "l");
		//leg2->AddEntry(piecewiseConf95, "95% CI", "f");
		//leg2->AddEntry(piecewiseConf99, "99% CI", "f");
		leg2->Draw();

		c1->Modified();
		c1->Update();

		std::string fn(filename);
		size_t a = fn.rfind("/"); // path/to/folder/file.ext    -- find last slash (/)
		size_t b = fn.substr(0, a - 1).rfind("/"); // find second to last slash (/)
		std::string fn1 = fn.substr(b + 1, fn.rfind(".") - b - 1); // attempt to grab "file"
		std::cout << fn1 << std::endl;

		c1->Print((fn1.substr(fn1.rfind("/") + 1, fn1.size() - fn1.rfind("/")) + (std::string)"-fit.png").c_str());

		TFile* out;
		if (outfilename == NULL) out = new TFile((fn1 + (std::string)"-fit-rootfile.root").c_str(), "RECREATE");
		else out = new TFile(outfilename, "RECREATE");

		prof->Write();
		logProf->Write();

		linearFit->Write();
		piecewiseFit->Write();
		//piecewiseConf95->Write();
		//piecewiseConf99->Write();
		//linearConf95->Write();
		//linearConf99->Write();

		out->Close();
		file->Close();
		outfilename = NULL;

		ofstream txtfile((fn1 + (std::string)"-fit.txt").c_str());
		if (txtfile.is_open()) {
			txtfile << "------------FITTING PARAMETERS-------------\n";
			txtfile << "Piecewise fit -> (x<[0]) ? [1] : [2]*x+[1]-[2]*[0]\n";
			txtfile << "(knee/turning point) p0: " << pwf0 << " +/- " << pwfe0 << std::endl;
			txtfile << "(vertical offset)    p1: " << pwf1 << " +/- " << pwfe1 << std::endl;
			txtfile << "(gradient)           p2: " << pwf0 << " +/- " << pwfe0 << std::endl;

			txtfile << "\nLinear fit -> [0]*x+[1]\n";
			txtfile << "(gradient)           p0: " << linearFit->GetParameter(0) << " +/- " << linearFit->GetParError(0) << std::endl;
			txtfile << "(vertical offset)    p1: " << linearFit->GetParameter(1) << " +/- " << linearFit->GetParError(1) << std::endl;
		}
		txtfile.close();
	}
	return 0;
}