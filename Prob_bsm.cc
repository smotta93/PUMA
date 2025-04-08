#include <iostream>
#include <cmath>
#include <string.h>
#include<float.h>
#include<complex.h>
#include <vector>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include <globes/globes.h>
#include<fstream>

#include <algorithm>

extern "C"
{
	#include "bsm.h"
}

using namespace std;

char AEDLFILE[] = "/home/smotta93/globes-3.2.18/PUMA/DUNE_GLoBES.glb";


int main(int argc, char * argv[])
{

	glbInit(argv[0]);
	glbInitExperiment(AEDLFILE, &glb_experiment_list[0], &glb_num_of_exps);

    ofstream outProb_nue_app, outProb_nuebar_app, outProb_numu_dis, outProb_numubar_dis;

    outProb_nue_app.open("prob/prob_nue_STD_bsm.dat");
    outProb_nuebar_app.open("prob/prob_nuebar_STD_bsm.dat");

    outProb_numu_dis.open("prob/prob_numu_STD_bsm.dat");
    outProb_numubar_dis.open("prob/prob_numubar_STD_bsm.dat");

/* 	outProb_nue_app.open("prob/prob_nue_STD_bsm.dat");
    outProb_nuebar_app.open("prob/prob_nuebar_STD_bsm.dat");

    outProb_numu_dis.open("prob/prob_numu_STD_bsm.dat");
    outProb_numubar_dis.open("prob/prob_numubar_STD_bsm.dat"); */

	double dm21 = 7.49e-5; //7.55e-5;
	double dm31 = 2.513e-3; //2.50e-3;
	double theta12 = 33.68 * (M_PI/180); //asin(sqrt(0.308)); 
	double theta23 = 43.30 * (M_PI/180); //asin(sqrt(0.574));
	double theta13 = 8.56 * (M_PI/180); // asin(sqrt(0.022));
	double deltacp = 212.0 * (M_PI/180); //1.08 * M_PI; 

	bsm_init_probability_engine_3();

	glbRegisterProbabilityEngine(8 * 9 - 3,
                               &bsm_probability_matrix,
							   &bsm_set_oscillation_parameters,
  							   &bsm_get_oscillation_parameters,
  							   NULL);
	/* Define "true" oscillation parameter vector */
	glb_params true_values = glbAllocParams();
  
	for(unsigned int i=0; i < 69; i++)
	{
	glbSetOscParams(true_values, 0.0, i);
	}

	glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);

	double abs_a_ee = 0;
	double abs_a_mue = 0;
	double arg_a_mue = 0;
	double abs_a_etau = 0;
	double arg_a_etau = 0; // 270 * (M_PI/180);
	double abs_a_mumu = 0;
	double abs_a_mutau = 0;
	double arg_a_mutau = 0;
	double abs_a_tautau = 0;

	double abs_c_ee = 0;
	double abs_c_mue = 0;
	double arg_c_mue = 0;
	double abs_c_etau = 0;
	double arg_c_etau = 0;
	double abs_c_mumu = 0;
	double abs_c_mutau = 0; //1.0e-32 / 1.0e-8;
	double arg_c_mutau = 0 ; // 0 * (M_PI/180);
	double abs_c_tautau = 0;

	//############ LIV Parameter #################################//
	glbSetOscParams(true_values, abs_a_ee, 51);  // a_ee 
	glbSetOscParams(true_values, abs_a_mue, 52);  // a_mue magnitude
    glbSetOscParams(true_values, arg_a_mue, 53);  // a_mue phase
    glbSetOscParams(true_values, abs_a_etau, 54);  // a_etau 
    glbSetOscParams(true_values, arg_a_etau, 55);  // a_etau phase
    glbSetOscParams(true_values, abs_a_mumu, 56);  // a_mumu
    glbSetOscParams(true_values, abs_a_mutau, 57);  // a_mutau
    glbSetOscParams(true_values, arg_a_mutau, 58);  // a_mutau phase
    glbSetOscParams(true_values, abs_a_tautau, 59);  // a_tautau

	glbSetOscParams(true_values, abs_c_ee, 60);  // c_ee 
	glbSetOscParams(true_values, abs_c_mue, 61);  // c_mue magnitude
    glbSetOscParams(true_values, arg_c_mue, 62);  // c_mue phase
    glbSetOscParams(true_values, abs_c_etau, 63);  // c_etau 
    glbSetOscParams(true_values, arg_c_etau, 64);  // c_etau phase
    glbSetOscParams(true_values, abs_c_mumu, 65);  // c_mumu
    glbSetOscParams(true_values, abs_c_mutau, 66);  // c_mutau
    glbSetOscParams(true_values, arg_c_mutau, 67);  // c_mutau phase
    glbSetOscParams(true_values, abs_c_tautau, 68);  // c_tautau

	glbSetDensityParams(true_values,1.0,GLB_ALL);

	glb_params input_errors = glbAllocParams();
  	
	glbSetDensityParams(input_errors, 0.05, GLB_ALL);
	glbSetInputErrors(input_errors);

	glbSetOscillationParameters(true_values);
	glbSetRates();
	
    double energy, Prob_nue_app, Prob_nuebar_app, Prob_numu_dis, Prob_numubar_dis;
	double emin= 0.25 ; //GeV
	double emax= 10.1 ; //GeV
	double step= 1000;
	//double L = 1300; // km
  
    for (energy = emin; energy <= emax; energy += (emax - emin) / step){

        // Probabilidade para o canal nue appearance
        Prob_nue_app = glbProfileProbability(0, 2, 1, +1, energy);
        outProb_nue_app << energy << "  " << Prob_nue_app << endl;

        // Probabilidade para o canal nuebar appearance
        Prob_nuebar_app = glbProfileProbability(0, 2, 1, -1, energy);
        outProb_nuebar_app << energy << "  " << Prob_nuebar_app << endl;
		
		// Probabilidade para o canal numu disappearance
        Prob_numu_dis = glbProfileProbability(0, 2, 2, +1, energy);
        outProb_numu_dis << energy << "  " << Prob_numu_dis << endl;
		
		// Probabilidade para o canal numubar disappearance
        Prob_numubar_dis = glbProfileProbability(0, 2, 2, -1, energy);
        outProb_numubar_dis << energy << "  " << Prob_numubar_dis << endl;
    }
	
	outProb_nue_app.close();
    outProb_nuebar_app.close();
    outProb_numu_dis.close();
    outProb_numubar_dis.close();

    glbFreeParams(true_values);
	
	cout << "done" << endl;

    return 0;
	
}