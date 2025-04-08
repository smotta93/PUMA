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
							   
	glb_params true_values = glbAllocParams();

	// Inicializa todos os parâmetros a zero
	for (unsigned int i = 0; i < 69; i++)
	{
		glbSetOscParams(true_values, 0.0, i);
	}

	// Define os parâmetros padrão do modelo padrão
	glbDefineParams(true_values, theta12, theta13, theta23, deltacp, dm21, dm31);

/* 	// Define os parâmetros LIV 
	double abs_a_ee = 0;
	double abs_a_mue = 0;
	double arg_a_mue = 0;
	double abs_a_etau = 0; //2.0e-23 / 1.0e-9;
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
	double abs_c_tautau = 0; */

	//############ LIV Parameter #################################//
	glbSetOscParams(true_values, 0, 51);  // a_ee 
	glbSetOscParams(true_values, 0, 52);  // a_mue magnitude 
    glbSetOscParams(true_values, 0, 53);  // a_mue phase
	glbSetOscParams(true_values, 0, 54);  // a_etau
    glbSetOscParams(true_values, 0, 55);  // a_etau phase
    glbSetOscParams(true_values, 0, 56);  // a_mumu
    glbSetOscParams(true_values, 0, 57);  // a_mutau
    glbSetOscParams(true_values, 0, 58);  // a_mutau phase
    glbSetOscParams(true_values, 0, 59);  // a_tautau
	/* glbSetOscParams(true_values, 0, 60);  // c_ee  */
	glbSetOscParams(true_values, 0, 61);  // c_mue magnitude
    glbSetOscParams(true_values, 0, 62);  // c_mue phase
    glbSetOscParams(true_values, 0, 63);  // c_etau 
    glbSetOscParams(true_values, 0, 64);  // c_etau phase
    glbSetOscParams(true_values, 0, 65);  // c_mumu
    glbSetOscParams(true_values, 0, 66);  // c_mutau
    glbSetOscParams(true_values, 0, 67);  // c_mutau phase
    glbSetOscParams(true_values, 0, 68);  // c_tautau
	
	// Parâmetros que variam
	double num_steps = 4;
	
	double a_test_value_start = -1e-27 * 10e-9;   //se estiver usando a divida pelo parâmetro 1.0e-9
	double a_test_value_end = -4e-26 * 10e-9;  //se estiver usando a divida pelo parâmetro 1.0e-9
	
/* 	double c_test_value_start = 1e-25 * 10e-36;   //se estiver usando a divida pelo parâmetro 1.0e-9
	double c_test_value_end = 1e-27 * 10e-36;  //se estiver usando a divida pelo parâmetro 1.0e-9 */
	
	// Calcula o fator multiplicativo para cobrir os passos desejados.
 	double a_test_value_step = pow(a_test_value_end / a_test_value_start, 1.0 / (num_steps - 1));
 	/* double a_test_value_step = pow(c_test_value_end / c_test_value_start, 1.0 / (num_steps - 1)); */

	double emin = 0.25; // GeV
	double emax = 100.1; // GeV
	double step = 10000;

	int file = 0;

	for (double a_test_value = a_test_value_start; a_test_value >= a_test_value_end; a_test_value *= a_test_value_step)
	/* for (double c_test_value = c_test_value_start; c_test_value >= c_test_value_end; c_test_value *= a_test_value_step)  */
	{
		double a_abs_par = a_test_value;
		/* double c_abs_par = c_test_value; */
/* 		cout << "valor do par: " << abs_c_mue << endl;*/

		// Atualiza os parâmetros para a fase atual
		glbSetOscParams(true_values, a_abs_par, 51); // a_ee
		glbSetOscParams(true_values, a_abs_par, 52);  // a_mue magnitude
		glbSetOscParams(true_values, a_abs_par, 54);  // a_etau 
/* 		glbSetOscParams(true_values, c_abs_par, 60);  // c_ee 
 */
		glbSetOscillationParameters(true_values);
		glbSetRates();
		
		// Extrai a fase em graus para o nome do arquivo
		string file_number = "file" + to_string(static_cast<int>(file));
		file++;

		// Criação dos nomes de arquivos
		string filename_nue_app = "prob/data/a/prob_nue_" + file_number + ".dat";
		string filename_nuebar_app = "prob/data/a/prob_nuebar_" + file_number + ".dat";
		string filename_numu_dis = "prob/data/a/prob_numu_" + file_number + ".dat";
		string filename_numubar_dis = "prob/data/a/prob_numubar_" + file_number + ".dat";

		ofstream outProb_nue_app(filename_nue_app);
		ofstream outProb_nuebar_app(filename_nuebar_app);
		ofstream outProb_numu_dis(filename_numu_dis);
		ofstream outProb_numubar_dis(filename_numubar_dis);
		
		// Verifica se os arquivos foram abertos corretamente
		if (!outProb_nue_app.is_open() || !outProb_nuebar_app.is_open() || !outProb_numu_dis.is_open() || !outProb_numubar_dis.is_open()) {
			cout << "Erro ao abrir os arquivos para a fase " << file << endl;
		} else {
			cout << "Arquivo para a fase " << file << " aberto com sucesso!" << endl;
		}


		// Calcula probabilidades e grava nos arquivos
		for (double energy = emin; energy <= emax; energy += (emax - emin) / step)
		{
			double Prob_nue_app = glbProfileProbability(0, 2, 1, +1, energy);
			double Prob_nuebar_app = glbProfileProbability(0, 2, 1, -1, energy);
			double Prob_numu_dis = glbProfileProbability(0, 2, 2, +1, energy);
			double Prob_numubar_dis = glbProfileProbability(0, 2, 2, -1, energy);

			outProb_nue_app << energy << " " << Prob_nue_app << " " << a_test_value / 10e-9 << endl; //se estiver usando "a" multiplique pelo parâmetro 1.0e-9
			outProb_nuebar_app << energy << " " << Prob_nuebar_app << " " << a_test_value / 10e-9<< endl; //se estiver usando "a" multiplique pelo parâmetro 1.0e-9
			outProb_numu_dis << energy << " " << Prob_numu_dis << " " << a_test_value / 10e-9<< endl; //se estiver usando "a" multiplique pelo parâmetro 1.0e-9
			outProb_numubar_dis << energy << " " << Prob_numubar_dis << " " << a_test_value / 10e-9 << endl; //se estiver usando "a" multiplique pelo parâmetro 1.0e-9
		}

		outProb_nue_app.close();
		outProb_nuebar_app.close();
		outProb_numu_dis.close();
		outProb_numubar_dis.close();
	}

	glbFreeParams(true_values);

	cout << "done" << endl;
	return 0;
}