#include <iostream>
#include <cmath>
#include <string.h>
#include <float.h>
#include <complex.h>
#include <vector>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <globes/globes.h>
#include <fstream>

using namespace std;

char AEDLFILE[] = "/home/smotta93/globes-3.2.18/PUMA/DUNE_GLoBES.glb";

int main(int argc, char * argv[])
{

    glbInit(argv[0]);
    glbInitExperiment(AEDLFILE, &glb_experiment_list[0], &glb_num_of_exps);

    ofstream outProb_nue_app, outProb_nuebar_app, outProb_numu_dis, outProb_numubar_dis;

    outProb_nue_app.open("prob/prob_nue_STD.dat");
    outProb_nuebar_app.open("prob/prob_nuebar_STD.dat");

    outProb_numu_dis.open("prob/prob_numu_STD.dat");
    outProb_numubar_dis.open("prob/prob_numubar_STD.dat");
	

	double dm21 = 7.49e-5; //7.55e-5;
	double dm31 = 2.513e-3; //2.50e-3;
	double theta12 = 33.68 * (M_PI/180); //asin(sqrt(0.308)); 
	double theta23 = 43.30 * (M_PI/180); //asin(sqrt(0.574));
	double theta13 = 8.56 * (M_PI/180); // asin(sqrt(0.022));
	double deltacp = 212.0 * (M_PI/180); //1.08 * M_PI; 
	
	
/* 	//DUNE TDR Vol. 2 - Table 5.1 data
	double dm21 = 7.39e-5;
    double dm31 = 2.451e-3;
    double theta12 = 0.5903;//33.45 * (M_PI/180) (se necessário converter para radianos)
    double theta23 = 0.866;
    double theta13 = 0.150; */

    cout << "dm21 = " << dm21 << endl;
    cout << "dm31 = " << dm31 << endl;
    cout << "theta12 = " << theta12 << " (" << pow(sin(theta12),2) << ")" << endl;
    cout << "theta23 = " << theta23 << " (" << pow(sin(theta23),2) << ")" << endl;
    cout << "theta13 = " << theta13 << " (" << pow(sin(theta13),2) << ")" << endl;
    cout << "Pi = " << M_PI << endl;

	glb_params true_values = glbAllocParams();
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, dm21, dm31);
	/*glbPrintParams(stdout, true_values);*/

    glbSetDensityParams(true_values, 1.0, GLB_ALL);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    double energy, Prob_nue_app, Prob_nuebar_app, Prob_numu_dis, Prob_numubar_dis;
	double emin= 0.25 ; //GeV
	double emax= 100.1 ; //GeV
	double step= 1000;
	
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
    return 0;
}
