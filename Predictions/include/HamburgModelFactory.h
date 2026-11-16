#ifndef Hamburg_Model_Factory
#define Hamburg_Model_Factory

#include "HamburgModel.h"
//#include "SiSensor.h"
#include "TFile.h"
#include "TTree.h"

class SiSensor {
  public:
    SiSensor() {

        pitch_ = 90; // um
        w2p_ = 0.25;
        strip_length_ = 5*2; // cm
        nstrips_ = 1016;
        thickness_ = 150; // um
        
        bulk_is_n_type_ = true;

        ini_leak_ = 0;
        ini_vdep_ = 200; // V
        
        detname_="";

        debug = false;
    }
    ~SiSensor() {};

    // Setters
    void SetSensorGeometry(Float_t pitch, Float_t w2p, Float_t strip_length,
                           Int_t nstrips, Float_t thickness) {
        pitch_ = pitch;
        w2p_ = w2p;
        strip_length_ = strip_length;
        nstrips_ = nstrips;
        thickness_ = thickness;
    }
    
    void SetBulkDopingType(Bool_t bulk_is_n_type) {
        bulk_is_n_type_ = bulk_is_n_type;
    }

    void SetInitialValues(Float_t ini_leak, Float_t ini_vdep) {
        ini_leak_ = ini_leak;
        ini_vdep_ = ini_vdep;
    }
    
    void SetDetName(std::string detname){
        detname_ = detname;
    }

    // Getters
    void GetSensorGeometry(Float_t &pitch, Float_t &w2p, Float_t &strip_length, Int_t                        &nstrips, Float_t &thickness) {
        pitch = pitch_;
        w2p = w2p_;
        strip_length = strip_length_;
        nstrips = nstrips_;
        thickness = thickness_;
        return;
    }

    Bool_t GetBulkDopingType() { return bulk_is_n_type_; }
    
    Float_t GetIniLeak() { return ini_leak_; }
    Float_t GetIniVdep() { return ini_vdep_; }
    
    std::string GetDetName() { return detname_; }
    
    // Printer
    void Print(){
        std::cout << "--------------------" << std::endl;
        std::cout << "Sensor data :" << std::endl;
        std::cout << "   pitch : "<< pitch_ << " um" << std::endl;
        std::cout << "   w2p : " << w2p_ << std::endl;
        std::cout << "   strip_length : " << strip_length_ << " cm" << std::endl;
        std::cout << "   nstrips : " << nstrips_ << std::endl;
        std::cout << "   thickness : " << thickness_ << " um " << std::endl;
        std::cout << "   ini_leak : "<< ini_leak_ << std::endl;
        std::cout << "   ini_vdep : "<< ini_vdep_ << " V " << std::endl;
        std::cout << "--------------------" << std::endl;
    }
    

  private:
    // data members
    // -------------

    // Sensor geometry
    Float_t pitch_;        // strip pitch in um
    Float_t w2p_;          // width of pitch ratio
    Float_t strip_length_; // strip length in cm
    Int_t nstrips_;
    Float_t thickness_; // sensor effective thickness in um

    // Doping type
    Bool_t bulk_is_n_type_;
    
    // Sensor name
    std::string detname_;

    // initial values;
    Float_t ini_leak_; // leakage current in A
    Float_t ini_vdep_; // full depletion voltage in V

    Bool_t debug;
};


class HamburgModelFactory {
    
  public:
    HamburgModelFactory(){

        periodInHours = 1;
        debug=true;
        
        // file
        fout = new TFile("Graph_new.root", "recreate");
        foutTree = 0;
    
        total_Feq=0;
    
    };
    ~HamburgModelFactory(){
        if(fout) fout->Close();
        if(foutTree) {foutTree->cd(); simutree->Write(); foutTree->Close();}
    };

    void readFluenceTempScenario(std::string filename);
    void drawFluenceTempScenario();
    void simulateSensorEvolution(SiSensor sensor);
    void drawSaveSensorSimu(std::string detname, bool print_plots);
    void saveSensorSimuInTree(std::string filename);
    
  private:
    
    
    
    // data members
    // -------------
    
    Int_t periodInHours;
    Bool_t debug;
    
    
    // containers for fluence/temp scenario
    
    std::vector< int > Step;
    std::vector< float > Fluence;
    std::vector< float > Temp;
    std::vector< int > Active;
    unsigned int Nperiods;
    
    
    // containers for Hamburg model (simu over time)

    std::vector< double > Feq;
    std::vector< double > Temp_evol;
    std::vector< double > I_leak;
    std::vector< double > Na;
    std::vector< double > NY;
    std::vector< double > Nc;
    std::vector< double > N;
    std::vector< double > U;
    
    std::vector< double > intFeq;
    Float_t total_Feq;
    
    // File to save graphs
    TFile* fout;


    // for output tree writing
    
    TFile *foutTree;
    TTree *simutree;
    
    Int_t iper;
    Float_t flu_out, T_out, Ileak_out, Vdep_out;
    
};

#endif
