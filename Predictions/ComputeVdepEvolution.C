#include "include/HamburgModelFactory.h"
#include "include/HamburgModelFactory.cc"

int main(){

    SiSensor sen; // Default parameters set for 2S module sensor
    // Phase1 TIB sensor
    //sen.SetSensorGeometry(80, 0.25, 11.7, 768, 290);
    //sen.SetBulkDopingType(0);
    //sen.SetInitialValues(0, 270); //ileak, vdep
    
    HamburgModelFactory factory;
    factory.readFluenceTempScenario("Inputs/scenario_BeamTest.txt");
    factory.simulateSensorEvolution(sen);
    factory.drawFluenceTempScenario();
    factory.drawSaveSensorSimu("", true);
}


