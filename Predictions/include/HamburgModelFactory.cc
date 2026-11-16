#include "HamburgModelFactory.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// Read fluence/temperature scenario
//-------------------------------

void HamburgModelFactory::readFluenceTempScenario(
    std::string filename = "scenario.txt") {

    Int_t nper = 0;

    Int_t step;
    Float_t fluence, temperature;
    Bool_t LVOn;
    std::string def, var;

    std::ifstream lumifile;
    std::string line;
    lumifile.open(filename.c_str());
    if (lumifile.is_open()) {
        std::cout << "Reading file : " << filename << std::endl;
        while (getline(lumifile, line)) {
            std::istringstream iss(line);
            // std::cout << line << std::endl;
            if (!(iss >> step >> fluence >> temperature >> LVOn)) {

                // Check if a variable is set
                std::istringstream iss2(line);
                if ((iss2 >> def >> var >> step) && def == "#define") {
                    std::cout << " Defining: " << var << " : " << step
                              << std::endl;
                    std::cout << " ----------------------- " << std::endl;
                    if (var == "StepInHours")
                        periodInHours = step;
                } else
                    continue;
            } else {
                Fluence.push_back(fluence);
                Temp.push_back(temperature);
                if (LVOn)
                    Active.push_back(1);
                else
                    Active.push_back(0);
                nper++;
            }
        }
        lumifile.close();
    } else
        std::cerr << "File '" << filename << "' not available." << std::endl;

    /*for (unsigned int i = 0; i < nper; i++) {
        std::cout << " " << i << " " << Fluence[i] << " " << Temp[i] <<
    std::endl;
    }*/

    Nperiods = nper;
    std::cout << " " << Nperiods << " periods found" << std::endl;

    // initialize containers
    for (unsigned int i = 0; i < Nperiods; i++) {
        Feq.push_back(0);
        Temp_evol.push_back(0);
        I_leak.push_back(0);
        Na.push_back(0);
        NY.push_back(0);
        Nc.push_back(0);
        N.push_back(0);
        U.push_back(0);
    }

    return;
}

// Loop over time steps
//---------------------

void HamburgModelFactory::simulateSensorEvolution(SiSensor sensor) {

    sensor.Print();

    // clear and initialize containers
    for (Int_t i = 0; i < Nperiods; i++) {
        Feq[i] = Fluence[i];
        Temp_evol[i] = ZeroC + Temp[i];
        I_leak[i] = sensor.GetIniLeak() / 1000 *
                    LeakCorrection(293.15, Temp_evol[i]); // A
        // charge carriets concentration
        Na[i] = 0;
        NY[i] = 0;
        Nc[i] = 0;
        intFeq.clear();
        std::cout << i << " " << Fluence[i] << " " << Temp[i] << std::endl;
    }

    // Ileak
    Double_t alpha1, alpha2, alpha3, delIleak;

    // Vdep
    TF1 *fff =
        new TF1("fff", "-0.00111/(x*x)+0.0586/x+0.240-0.651*x+0.355*x*x");
    Float_t pitch, w2p, strip_length, thickness;
    Int_t nstrips;
    sensor.GetSensorGeometry(pitch, w2p, strip_length, nstrips, thickness);
    Double_t volume = strip_length * pitch/10000 * nstrips * thickness/10000; // in cm^3
    std::cout<< " Si volume : " <<volume<<" cm^3"<<std::endl;
    Int_t sign = -1;
    if (sensor.GetBulkDopingType()==0) sign = +1;
    Double_t Ufd_diode =
        sensor.GetIniVdep() / (1. + 2. * pitch / thickness * fff->Eval(w2p));
    // constants eps, q_el defined in HamburgModel.h
    Double_t Neff0 = sign * Ufd_diode * 2. * eps * eps0 / q_el /
                     (thickness * thickness * 1e-8); // initial Neff
    std::cout << " Neff corr: "
              << (1. + 2. * pitch / thickness * fff->Eval(w2p)) << std::endl;

    Double_t Na_temp, NY_temp;
    Double_t temperature;

    Float_t total_Feq = 0;

    // Container for Ileak contributions from each previous time periods
    std::vector<std::vector<double>> A1;
    std::vector<std::vector<double>> A2;
    std::vector<std::vector<double>> A3;

    //----------------------
    int TempAlgo = 1; // 1 from Christian
    //----------------------

    // Loop over periods for ILEAK
    for (Int_t i = 0; i < Nperiods; i++) {

        // integrated fluence
        total_Feq += Feq[i];
        intFeq.push_back(total_Feq);

        // std::cout<<i<<" Feq: "<<Feq[i]<<" temp: "<<Temp_evol[i]<<std::endl;

        alpha1 = a10;
        alpha2 = I_a2(Temp_evol[i]); // Term with no annealing, use temperature
                                     // in the time period
        alpha3 = I_a3(1, periodInHours);

        // Annealing: consider a flash irradiation at beginning of the period
        // Update leakage current in current and future periods
        // Update consequent temperature increase

        if (TempAlgo == 1) // Christian's algorithm
        {
            // Loop on current and next periods to compute contributions with
            // annealing
            for (Int_t j = i; j < Nperiods;
                 j++) { /// Start from i : compute leakage current for the end
                        /// of the period i

                    
                alpha1 *=
                    I_a1(Temp_evol[j], 1,
                         periodInHours); /// *= is equivalent to add the
                                         /// additionnal time period in the
                                         /// exponential. Takes into account the
                                         /// temperature in the current time
                                         /// period (can be different than the
                                         /// one at the time of the irradiation)

                // use lowest value from current and previous periods
                if (I_a2(Temp_evol[j]) < alpha2)
                    alpha2 =
                        I_a2(Temp_evol[j]); /// intercept (t=0) is temperature
                                            /// dependent ; can only decrease

                alpha3 = I_a3(j - i + 1, periodInHours);
                delIleak = (alpha1 + alpha2 + alpha3) * Feq[i] * volume * 1000;

                //if(j==i) std::cout<<" Ileak "<<j<<" "<<alpha1<<" "<<alpha2<<" "<<alpha3<<std::endl;
                //if(j==i) std::cout<<" "<<j<<" "<<Feq[i]<<" "<<Temp_evol[j]<<" "<<delIleak<<std::endl;

                // if detector active, could simulated selfHeating
                // not done here
                
                I_leak[j]+=delIleak; //in mA per module

            }

        } else // alternative algo. Take into account summed Ileak for computing
               // increased temperature
        {
            // Container preparation
            A1.push_back(std::vector<double>(i + 1, 0));
            A2.push_back(std::vector<double>(i + 1, 0));
            A3.push_back(std::vector<double>(i + 1, 0));

            //  Current contribution:
            Double_t I_a2_i = I_a2(Temp_evol[i]);
            A1[i][0] = a10 * I_a1(Temp_evol[i], 1, periodInHours);
            A2[i][0] = I_a2_i;
            A3[i][0] = I_a3(1, periodInHours);

            // Loop on time previous periods to update their contribution with
            // annealing (using current temperature) in A[i][j]  j=0 is current
            // one, j=1 is previous one, ...
            if (i > 0)
                for (Int_t j = 1; j <= i; j++) {
                    A1[i][j] =
                        A1[i - 1][j - 1] * I_a1(Temp_evol[i], 1, periodInHours);
                    A2[i][j] = A2[i - 1][j - 1];
                    if (I_a2_i < A2[i][j])
                        A2[i][j] = I_a2_i;
                    A3[i][j] = I_a3(j + 1, periodInHours);
                }

            /*std::cout<<"A1 "<<std::endl;
            for(Int_t j=0; j<=i; j++){
                for(Int_t k=0; k<=j; k++) std::cout<<" "<<A1[j][k];
                std::cout<<std::endl;
            }*/

            // Sum contributions from previous periods
            alpha1 = 0;
            alpha2 = 0;
            alpha3 = 0;
            for (Int_t j = 0; j <= i; j++) {
                alpha1 += A1[i][j] *
                          Feq[i - j]; // take into account fluence at that time
                alpha2 += A2[i][j] * Feq[i - j];
                alpha3 += A3[i][j] * Feq[i - j];
            }
            delIleak = (alpha1 + alpha2 + alpha3) * volume * 1000;
            // std::cout<<" "<<i<<" "<<alpha1<<" "<<alpha2<<" "<<alpha3<<"
            // "<<delIleak<<std::endl;

            I_leak[i] += delIleak; // in mA per module

        } // end of alternative algo


        // std::cout<<i<<" "<<I_leak[i]<<std::endl;

    } // end loop over periods ILEAK

    // Loop over periods for VDEP
    for (Int_t i = 0; i < Nperiods; i++) {

        Nc[i] = -N_stable(Neff0, intFeq[i]);
        Na_temp = ga * Feq[i];
        NY_temp = 1;
        // annealing
        for (Int_t j = i; j < Nperiods; j++) {
            Na_temp = N_short(1., Temp_evol[j], Na_temp, periodInHours);
            Na[j] += -Na_temp;
            NY_temp *= N_reverse2(1., Temp_evol[j], periodInHours);
            NY[j] += -gY * Feq[i] * (1. - NY_temp);
        }
        // std::cout<<i<<" "<<Neff0<<" "<<Nc[i]<<" "<<Na[i]<<"
        // "<<NY[i]<<std::endl; std::cout<<i<<" "<<Nc[i]<<" "<<Na[i]<<"
        // "<<NY[i]<<std::endl;

        // filling Neff and full depletion voltage
        N[i] = Neff0 + Nc[i] + Na[i] + NY[i];
        U[i] = (1. + 2. * pitch / thickness * fff->Eval(w2p)) * fabs(N[i]) *
               (thickness * thickness * 1e-8) * q_el / (2. * eps * eps0);

    } // end loop over periods VDEP

    std::cout << " total Feq " << total_Feq << std::endl;
}

// Fill, draw and save plots
//---------------------------

void HamburgModelFactory::drawSaveSensorSimu(std::string detname = "",
                                             bool print_plots = false) {

    // Could decide to average over periods

    Int_t maxTime = Nperiods;
    fout->cd();

    // Module temperature & Ileak
    TH1D *h_T =
        new TH1D("h_T", "Module Temperature", maxTime, 0.5, maxTime + 0.5);
    TH1D *h_I_leak = new TH1D("I_leak", "Leakage Current (not corrected)",
                              maxTime, 0.5, maxTime + 0.5);
    TH1D *h_I_leak_corr = new TH1D("I_leak_corr", "Leakage Current", maxTime,
                                   0., periodInHours * 60. * 60. * maxTime);

    // Vdep related plots
    TH1D *h_Nc = new TH1D("Nc", "Stable Damage", maxTime, 0.5, maxTime + 0.5);
    TH1D *h_Na =
        new TH1D("Na", "Short Term Annealing", maxTime, 0.5, maxTime + 0.5);
    TH1D *h_NY =
        new TH1D("NY", "Reverse Annealing", maxTime, 0.5, maxTime + 0.5);
    TH1D *h_N = new TH1D("N", "Neff", maxTime, 0.5, maxTime + 0.5);
    TH1D *h_U = new TH1D("U", "V_dep", maxTime, 0.,
                         periodInHours * 60. * 60. * maxTime);

    TGraph *lumigr = new TGraph(maxTime);
    TGraph *feqgr = new TGraph(maxTime);

    for (Int_t i = 0; i < Nperiods; i++) {

        // if(i>200) std::cout << i+1 <<" "<<Temp_evol[i]<<"
        // "<<I_leak[i]<<std::endl;
        h_T->SetBinContent(i + 1, Temp_evol[i]);
        h_I_leak->SetBinContent(i + 1, I_leak[i]);
        h_I_leak_corr->SetBinContent(
            i + 1, I_leak[i] * LeakCorrection(Temp_evol[i], ZeroC + 20));

        h_Nc->SetBinContent(i + 1, Nc[i]);
        h_Na->SetBinContent(i + 1, Na[i]);
        h_NY->SetBinContent(i + 1, NY[i]);
        h_N->SetBinContent(i + 1, N[i]);
        if (debug && (i % 100 == 0))
            std::cout << "     U " << U[i] << std::endl;
        h_U->SetBinContent(i + 1, U[i]);

        feqgr->SetPoint(i, intFeq[i], U[i]);
        if (i == Nperiods - 1)
            std::cout << " total Feq " << intFeq[i] << std::endl;
    }

    // gStyle->SetTitleX(0.1f);
    // gStyle->SetTitleW(0.8f);
    gStyle->SetLineColor(kWhite);
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas();
    TString title, name;
    if (print_plots) {

        TCanvas *c2 = new TCanvas();
        c2->cd();
        if (detname!="") title.Form("Module temperature - %s", detname.c_str());
        else title="Module temperature";
        h_T->SetTitle(title);
        h_T->GetXaxis()->SetTitle("Time (steps)");
        h_T->GetYaxis()->SetTitle("Temperature (#circK)");
        h_T->DrawCopy();
        name.Form("VdepEvol_T_%s.png", detname.c_str());
        c2->SaveAs(name);
        h_T->SetName(Form("T_%s", detname.c_str()));
        h_T->Write();

        TCanvas *c3 = new TCanvas();
        c3->cd();
        if (detname!="") title.Form("Corrected leakage current - %s", detname.c_str());
        else title="Corrected leakage current";
        h_I_leak_corr->SetTitle(title);
        h_I_leak_corr->GetXaxis()->SetTimeDisplay(1);
        // h_I_leak_corr->GetXaxis()->SetTimeFormat("%d\/%m\/%y%F2010-01-01
        // 00:00:00");
        h_I_leak_corr->GetXaxis()->SetTitle("Time");
        h_I_leak_corr->GetYaxis()->SetTitle("Leakage current (mA)");
        h_I_leak_corr->DrawCopy();
        name.Form("VdepEvol_Ileak_%s_.png", detname.c_str());
        c3->SaveAs(name);
        h_I_leak_corr->SetName(Form("Ileak_%s", detname.c_str()));
        h_I_leak_corr->Write();

        TCanvas *c4 = new TCanvas();
        c4->cd();
        if (detname!="") title.Form("Full depletion voltage - %s", detname.c_str());
        else title="Full depletion voltage";
        h_U->SetTitle(title);
        h_U->GetXaxis()->SetTimeDisplay(1);
        // h_U->GetXaxis()->SetTimeFormat("%d\/%m\/%y%F2010-01-01 00:00:00");
        h_U->GetXaxis()->SetTitle("Time");
        h_U->GetYaxis()->SetTitle("Full depletion voltage (V)");
        h_U->DrawCopy();
        name.Form("VdepEvol_Vdep_%s.png", detname.c_str());
        c4->SaveAs(name);
        h_U->SetName(Form("Vdep_%s", detname.c_str()));
        h_U->Write();

        TCanvas *c5 = new TCanvas();
        c5->cd();
        if (detname!="") title.Form("Effective space charge contributions - %s", detname.c_str());
        else title="Effective space charge contributions";
        h_Nc->SetTitle(title);
        h_Nc->SetLineColor(kBlue);
        h_Nc->DrawCopy();
        h_Na->SetLineColor(kGreen);
        h_Na->DrawCopy("SAME");
        h_NY->SetLineColor(kRed);
        h_NY->DrawCopy("SAME");
        name.Form("VdepEvol_Ncontrib_%s.png", detname.c_str());
        c5->SaveAs(name);
        h_N->SetName(Form("N_%s", detname.c_str()));
        h_N->Write();
        h_Nc->SetName(Form("Nc_%s", detname.c_str()));
        h_Nc->Write();
        h_Na->SetName(Form("Na_%s", detname.c_str()));
        h_Na->Write();
        h_NY->SetName(Form("NY_%s", detname.c_str()));
        h_NY->Write();

        c5->cd();
        if (detname!="") title.Form("Effective space charge (cm-3) - %s",
                   detname.c_str());
        else title = "Effective space charge (cm-3)";
        h_N->SetTitle(title);
        h_N->GetXaxis()->SetTitle("Time (steps)");
        h_N->GetYaxis()->SetTitle("Effective space charge (cm-3)");
        h_N->DrawCopy();
        name.Form("VdepEvol_Neff_%s.png", detname.c_str());
        c5->SaveAs(name);
    }

    TCanvas *cgr = new TCanvas();
    cgr->cd();

    if (detname!="") title.Form("Vdep vs Fluence - %s", detname.c_str());
    else title = "Vdep vs Fluence";
    name.Form("VdepEvol_VdepvsFeq_%s.png", detname.c_str());
    feqgr->GetHistogram()->GetYaxis()->SetTitle("Full depletion voltage (V)");
    feqgr->GetHistogram()->GetXaxis()->SetTitle("Fluence (1MeV neq/cm^{-2})");
    feqgr->SetTitle(title);
    feqgr->SetName(name);
    feqgr->SetLineColor(1);
    feqgr->Draw("alp");
    if (print_plots)
        cgr->SaveAs(name);
    feqgr->SetName(Form("feqgr_%s", detname.c_str()));
    feqgr->Write();

    h_T->Delete();
    h_I_leak->Delete();
    h_I_leak_corr->Delete();

    h_Nc->Delete();
    h_Na->Delete();
    h_NY->Delete();
    h_N->Delete();
    h_U->Delete();
}

void HamburgModelFactory::drawFluenceTempScenario() {

    Int_t maxTime = Nperiods;
    fout->cd();

    // fluence/temperature history
    TH1D *h_Feq = new TH1D("h_Feq", "Fluence", maxTime, 0.5,
                           periodInHours * 60. * 60. * maxTime + 0.5);
    TH1D *h_Feq_total = new TH1D("h_Feq_total", "Total Fluence", maxTime, 0.5,
                                 periodInHours * 60. * 60. * maxTime + 0.5);
    TH1D *h_T_base = new TH1D("h_T_base", "Temperature set points", maxTime,
                              0.5, periodInHours * 60 * 60 * maxTime + 0.5);

    for (Int_t i = 0; i < Nperiods; i++) {
        h_Feq->SetBinContent(i + 1, Feq[i]);
        h_Feq_total->SetBinContent(i + 1, intFeq[i]);
        h_T_base->SetBinContent(i + 1, Temp[i]);
    }

    TString title;
    TCanvas *cF = new TCanvas();
    cF->cd();
    //	  title.Form("Fluence vs time");
    //	  h_Feq_total->SetTitle(title);
    h_Feq_total->SetName("intFeq");
    h_Feq_total->GetXaxis()->SetTimeDisplay(1);
    // h_Feq_total->GetXaxis()->SetTimeFormat("%H");
    h_Feq_total->GetXaxis()->SetTitle("Time");
    h_Feq_total->GetYaxis()->SetTitle("Fluence (1MeV neq/cm^{-2})");
    h_Feq_total->DrawCopy();
    cF->SaveAs("VdepEvol_Feq_total.png");
    h_Feq_total->Write();

    TCanvas *cT = new TCanvas();
    cT->cd();
    title.Form("Temperature without corrections");
    h_T_base->SetTitle(title);
    h_T_base->SetName("Temp");
    h_T_base->GetXaxis()->SetTimeDisplay(1);
    // h_T_base->GetXaxis()->SetTimeFormat("%H");
    h_T_base->GetXaxis()->SetTitle("Time");
    h_T_base->GetYaxis()->SetTitle("Temperature (#circC)");
    h_T_base->DrawCopy();
    cT->SaveAs("VdepEvol_T_base.png");
    h_T_base->Write();
}

// Save results in a tree
//------------------------

void HamburgModelFactory::saveSensorSimuInTree(
    std::string filename = "SimuTree.root") {

    /* if (!foutTree) {
         foutTree = new TFile(filename.c_str(), "recreate");
         simutree = new TTree("simutree", "");

         simutree->Branch("detid", &detid, "detid/I");
         simutree->Branch("period", &iper, "period/I");
         simutree->Branch("fluence", &flu_out, "fluence/F");
         simutree->Branch("lumi", &lumi_out, "lumi/F");
         simutree->Branch("T", &T_out, "T/F");
         simutree->Branch("Ileak", &Ileak_out, "Ileak/F");
         simutree->Branch("Vdep", &Vdep_out, "Vdep/F");
     }

     foutTree->cd();

     for (Int_t i = 0; i < Nperiods; i++) {
         flu_out = intFeq[i];
         T_out = Temp_evol[i];
         Ileak_out = I_leak[i] * LeakCorrection(Temp_evol[i], ZeroC + 20);
         Vdep_out = U[i];
         iper = i + 1;
         simutree->Fill();
     }*/
}
