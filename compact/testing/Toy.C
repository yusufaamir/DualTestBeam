/*
 * This is an example to generate readout of FrontEnd
 *
 * Copy this file "Toy.C" and the root file "pulse_shape.root" in the same directory
 *
 * To run first example
 * $> root -l
 * [] .L Toy.C+
 * [] PlotOneEvent()
 *
 *
 * To run second example:
 * $> root -l
 * [] .L Toy.C+
 * [] GenerateFewEvents()
 * and examine TGraphs in "output.root"
 * 
 * To run third analysis example:
 * $> root -l
 * [] .L Toy.C+
 * [] TFile *f = new TFile("output.root");
 * [] PlotSummedSignal();        // for summed signal
 * [] PlotSummedSignal(true);    // for normalized signal
 * 
 * See more explanations in PlotOneEvent() below
 */

#include <vector>
#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

TRandom3 rnd;

// Global number of events to simulate and analyze
const int nEvents = 10;

//-------------------------------------------------------------------

class SPR {
public:
    SPR(int);
    ~SPR();
    double Evaluate(double);
//    double tMin()   { return tMin_; };
private:
    std::vector<double> vspr_;
    double dt_;
    double tMin_;
    double tMax_;
};



SPR::SPR(int iopt)
{

    vspr_.clear();
    if(iopt==0){
        /*
        * This example of SPR corresponds to Calvision TB at FNAL in 2023
        */
        dt_    = 0.1;
        tMin_  = 0.0;
        tMax_  = 1000.0;

        double tRise       = 0.853;
        double tDecay      = 6.538;
        double tUnderShoot = 101.7;
        double norm        = 0.111051;

        double tNow = tMin_;
        while(tNow <= tMax_){
            double a = 1./ tRise;
            double b = 1./ tDecay;
            double A = -a * b / (a - b);
            double B = -A;
            double result = A * exp(-a*tNow) + B * exp(-b*tNow);

            double g = 1./ tUnderShoot;
            double Atmp = -A * g / ( a - g);
            double Btmp = -B * g / ( b - g);
            double G = - Atmp - Btmp ;
            A = Atmp;
            B = Btmp;
            result -= A * exp(-a*tNow) + B * exp(-b*tNow) + G * exp(-g*tNow);

            vspr_.push_back(result / norm);
            tNow += dt_;
        }
    }else{
        /*
        * This example of SPR corresponds to Calvision TB at DESY in 2024
        */
        dt_    = 0.1;
        tMin_  = 0.0;
        tMax_  = 900.0;
        TFile *fin = new TFile("pulse_shape.root");
        TGraph *gr = (TGraph*)fin->Get("gr");
        double tNow = tMin_;
        while(tNow <= tMax_){
            vspr_.push_back(gr->Eval(tNow));
            tNow += dt_;
        }
        fin->Close();
        delete gr;
    }

}



SPR::~SPR()
{
}


double SPR::Evaluate(double x)
{
    double result = 0.;
    if(x > tMin_){
        unsigned int indx = (x - tMin_) / dt_;
        if(indx < vspr_.size() - 1){
            result = vspr_[indx] + (vspr_[indx+1] - vspr_[indx]) * (x - (tMin_ + dt_*indx)) / dt_;
        }
    }
    return result;
}


//-----------------------------------------------------------------

class PhotoElectrons {
public:
    PhotoElectrons(double,double);
    ~PhotoElectrons();
    void Generate(double,double);
    int Nc()   { return Nc_; };
    int Ns()   { return Ns_; };
    double t(int i) { return vpe_[i]; };
private:
    std::vector<double> vpe_;
    int Nc_;
    int Ns_;
    TF1 *fDispersion_;
    TF1 *fScintillation_;
    double tauDispersion_;
    double tauScintillation_;
};



PhotoElectrons::PhotoElectrons(double tau1, double sigma)
{
    tauScintillation_ = tau1;
    tauDispersion_ = sigma / 1.40;
    fDispersion_ = new TF1("fDispersion_","x*exp(-x/[0])",0,10*tauDispersion_);
    fDispersion_->SetParameter(0,tauDispersion_);
    fScintillation_ = new TF1("fScintillation_","exp(-x/[0])",0,10*tauScintillation_);
    fScintillation_->SetParameter(0,tauScintillation_);
}



PhotoElectrons::~PhotoElectrons()
{
    delete fDispersion_;
    delete fScintillation_;
}




void PhotoElectrons::Generate(double nS, double nC)
{
    vpe_.clear();
    Nc_ = 0;
    Ns_ = 0;
    if(nC>0){
        Nc_ = rnd.Poisson(nC);
    }
    if(nS>0){
        Ns_ = rnd.Poisson(nS);
    }
    for(int i=0; i<Nc_; i++){
        double t = 0;
        t += fDispersion_->GetRandom();
        vpe_.push_back(t);
    }
    for(int i=0; i<Ns_; i++){
        double t = 0;
        t += fDispersion_->GetRandom();
        t += fScintillation_->GetRandom();
        vpe_.push_back(t);
    }
}


//-----------------------------------------------------

class FrameFE {
public:
    FrameFE(int,double,double);
    ~FrameFE();
    void Clear();
    void Add(double, SPR*);
    int    n()      { return n_; };
    double x(int i) { return x_[i]; };
    double y(int i) { return y_[i]; };
private:
    int n_;
    std::vector<double> x_;
    std::vector<double> y_;
    double a1pe_;
};



FrameFE::FrameFE(int n, double dt, double a1pe)
{
    n_ = n;
    x_.clear();
    y_.clear();
    for(int i=0; i<n_; i++){
        x_.push_back(i*dt);
        y_.push_back(0.);
    }
    a1pe_ = a1pe;
}


FrameFE::~FrameFE()
{
}


void FrameFE::Clear()
{
    for(int i=0; i<n_; i++){
        y_[i] = 0.;
    }
}


void FrameFE::Add(double t, SPR *spr)
{
    for(int i=0; i<n_; i++){
        double x = x_[i] - t;
        y_[i] += a1pe_ * spr->Evaluate(x);
    }
}

//-----------------------------------------------------


void PlotOneEvent()
{
    /*
     * This example generates one event and plots it.
     * It uses SiPM + amplifier parameters from CALVISION TB 2023 at FNAL
     */


    /*
     * Define Single Photo-Electron Response function (SPR)
     * It is a waveform from a single PE at the FrontEnd with amplitude = 1
     *
     * There are two options here
     *  0 - SPR is available as analytic function. This one corresponds to SPR in CALVISION TB 2023 at FNAL
     *  1 - SPR is available as TGraph. This one corresponds to SPR in CALVISION TB 2024 at DESY
     */

    int iOption = 0;
    SPR *spr = new SPR(iOption);

    /*
     * Define FrontEnd Frame for readout from one channel:
     * - the number of samples
     * - the time interval between samples
     * - the amplitude of 1-pe pulse. It depends on gain of the SiPM and FrontEnd amplification
     */

    int nSamples        = 1000;
    double timeInterval = 0.2;          // nanoseconds
    double amp1pe       = 0.60;         // millivolts;
    FrameFE *fe = new FrameFE(nSamples, timeInterval, amp1pe);


    /*
     * Now we need a vector of values for time of arrival of photo-electrons in the event. We use a dummy generator that produces Nc Cherenkov photo-electrons with "prompt" timing and Ns scintillation photo-electrons with additional random delays with scintillation decay time. "Prompt" time is a random value with "prompt" resolution that corresponds to fluctuations due to light collection, SiPM avalanche developement etc
     */

    double tauScintillation = 300.;   // nanoseconds
    double sigmaPrompt      = 0.30;   // nanoseconds
    PhotoElectrons *pe = new PhotoElectrons( tauScintillation, sigmaPrompt);


    /*
     * Generate FrontEnd Frame with expected rate of Cherenkov and scintillation photo-electrons nC and nS, respectively
     * Start of the pulses is at T0
     */

    double nC = 20.;     // photo-electrons
    double nS = 200.;    // photo-electrons
    double T0 = 20.0;    // nanoseconds

    fe->Clear();
    pe->Generate( nS,nC);
    double ntot = pe->Nc()+pe->Ns();
    for(int i=0; i<ntot; i++){
        fe->Add(pe->t(i) + T0, spr);
    }
    cout << "Ns = " << pe->Ns() << "     Nc = " << pe->Nc() << endl;



    /*
     * Make TGraph with FE readout
     */

    TGraph *gr = new TGraph();
    for(int i=0; i<fe->n(); i++){
        gr->SetPoint(i, fe->x(i), fe->y(i));
    }
    gr->GetXaxis()->SetTitle("t (ns)");
    gr->GetYaxis()->SetTitle("a(t) (mV)");
    gr->Draw("APL");

}




void GenerateFewEvents()
{
    /*
     * This examples generates 10 events. Readout for each event is saved in root file "output.root" as TGraph
     * See detailed comments in PlotOneEvent()
     */

    int iOption = 1;
    SPR *spr = new SPR(iOption);

    int nSamples        = 1000;
    double timeInterval = 1.0;          // nanoseconds
    double amp1pe       = 0.2;          // millivolts;
    FrameFE *fe = new FrameFE(nSamples, timeInterval, amp1pe);

    double tauScintillation = 300.;   // nanoseconds
    double sigmaPrompt      = 0.30;   // nanoseconds
    PhotoElectrons *pe = new PhotoElectrons( tauScintillation, sigmaPrompt);

    double nC = 20.;     // photo-electrons
    double nS = 200.;    // photo-electrons
    double T0 = 100.;    // nanoseconds


    TGraph *gr[nEvents];
    for(int ievt=0; ievt<nEvents; ievt++){
        fe->Clear();
        pe->Generate( nS,nC);
        double ntot = pe->Nc()+pe->Ns();
        for(int i=0; i<ntot; i++){
            fe->Add(pe->t(i) + T0, spr);
        }
        gr[ievt] = new TGraph();
        for(int i=0; i<fe->n(); i++){
            gr[ievt]->SetPoint(i, fe->x(i), fe->y(i));
        }
        gr[ievt]->GetXaxis()->SetTitle("t (ns)");
        gr[ievt]->GetYaxis()->SetTitle("a(t) (mV)");
        gr[ievt]->SetName(Form("gr_%d",ievt));
    }

    TFile *fout = new TFile("output.root","recreate");
    for(int i=0; i<nEvents; i++){
        gr[i]->Write();
    }
    
    TCanvas *c1 = new TCanvas("c1", "All events", 800, 600);
	int colors[10] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet, kBlack, kGray+1, kPink+6};

	// Plotting the first 10 event signals

	for (int i = 0; i < 10; i++) {
    	TString name = Form("gr_%d", i);
    	TGraph *g = (TGraph*)gDirectory->Get(name);
    	g->SetLineColor(colors[i]);
    	if (i == 0)
        	g->Draw("AL");
    	else
        	g->Draw("L SAME");
	}

	// Save the canvas as a PNG image
	c1->SaveAs(Form("all_events_overlay_%d_events.png", nEvents));
    
    fout->Close();
}


void PlotSummedSignal(bool normalize = false)
{
    TGraph *gr[nEvents];

    // Load graphs from the current ROOT directory
    for (int i = 0; i < nEvents; i++) {
        TString name = Form("gr_%d", i);
        gr[i] = (TGraph*)gDirectory->Get(name);
        if (!gr[i]) {
            std::cerr << "Error: TGraph " << name << " not found in memory.\n";
            return;
        }
    }

    // This assumes all graphs have same time points
    int nPoints = gr[0]->GetN();
    double *x = gr[0]->GetX();
    std::vector<double> sumY(nPoints, 0.0);

    // Sum (or average) across all events
    for (int i = 0; i < nEvents; i++) {
        double *y = gr[i]->GetY();
        for (int j = 0; j < nPoints; j++) {
            if (normalize)
                sumY[j] += y[j] / nEvents;
            else
                sumY[j] += y[j];
        }
    }

    // Create summed TGraph
    TGraph *grSum = new TGraph(nPoints);
    for (int i = 0; i < nPoints; i++) {
        grSum->SetPoint(i, x[i], sumY[i]);
    }

    grSum->SetTitle(normalize ? "Average Signal from 10 Events" : "Summed Signal from 10 Events");
    grSum->GetXaxis()->SetTitle("t (ns)");
    grSum->GetYaxis()->SetTitle(normalize ? "Average a(t) (mV)" : "Total a(t) (mV)");
    grSum->SetLineColor(kRed+1);
    grSum->SetLineWidth(3);

    TCanvas *c4 = new TCanvas("c4", "Summed Signal", 800, 600);
    grSum->Draw("AL");
    c4->SaveAs(normalize 
    	? Form("average_signal_%d_events.png", nEvents) 
    	: Form("summed_signal_%d_events.png", nEvents));
}
