/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Imperial College London                                              ///
///  m.white08@imperial.ac.uk                                             ///
///                                                                       ///
///  Please feel free to use and modify if you wish. However,             ///
///  please provide appropriate acknowledgement and get in touch          ///
///  if you have any questions. This is not necessarily the               ///
///  final, canonical version of the code - contact me to get that.       ///
///                                                                       ///
///  There is a prize of a pint for reporting any serious bugs or         ///
///  finding something new that results in >20% speed up.                 ///
///                                                                       ///
///                                                                       ///
///  The code below is structured as follows:                             ///
///                                                                       ///
///  0. SETTING UP STRUCTURES AND CLASSES                                 ///
///     All parameter values are stored in a structure called params.     ///
///     A class is created which stores all the information of a          ///
///     single individual.                                                ///
///     A structure called population stores all individuals.             ///
///     The time-dependent output of the model is stored in a             ///
///     structure called simulation.                                      ///
///                                                                       ///
///  1. MAIN - SIMULATION                                                 ///
///     Here we read in parameters from files (model parameters,          ///
///     mosquito parameters, intervention parameters).                    ///
///     A population of individuals is created at equilibrium and         ///
///     then simulated.                                                   ///
///                                                                       ///
///  2. FUNCTIONS                                                         ///
///     Useful functions needed for model simulations.                    ///
///     Mosquitoes are simulated using a deterministic ODE solver.        ///
///     The functions for simulating a population of individuals are      ///
///     provided here - note that the model itself appears in Section 4.  ///
///                                                                       ///
///  3. EQUILIBRIUM SETUP                                                 ///
///     This set of functions calculates the equilibrium set up of the    ///
///     population. It is only called once while the population is        ///
///     being initialised.                                                ///
///                                                                       ///
///  4. INDIVIDUAL-BASED MODEL                                            ///
///     Details of the stochastic individual-based model for each         ///
///     person. Transitions occur with a fixed time step according to     ///
///     compting hazards                                                  ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include "randlib.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

#include "Params.h"
#include "Intervention.h"
#include "Population.h"


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//          //                                      //
//   ####   //   ###  ##### ######   ##  ## #####   //
//  ##  ##  //  ##    ##      ##     ##  ## ##  ##  //
//  ##  ##  //   ###  ####    ##     ##  ## #####   //
//  ##  ##  //     ## ##      ##     ##  ## ##      //
//   ####   //   ###  #####   ##      ####  ##      //
//          //                                      //
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.5. Define a structure for storing the output of a simulation                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

struct simulation
{
    //////////////////////////////////////////
    // 0.5.1. Vector of simulation times

    int N_time;

    vector<double> t_vec;


    //////////////////////////////////////////
    // 0.5.2. Tracking output

    vector<vector<int>> yH_t;
    vector<vector<vector<double>>> yM_t;


    vector<double> EIR_t;

    vector<vector<int>> prev_all;   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T} 
    vector<vector<int>> prev_U5;    // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T} 
    vector<vector<int>> prev_U10;   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T} 


    /////////////////////////////////////////
    // 0.5.3. Tracking coverage over time

    vector<int> LLIN_cov_t;
    vector<int> IRS_cov_t;
    vector<int> ACT_treat_t;
    vector<int> PQ_treat_t;
    vector<int> pregnant_t;


    //////////////////////////////////////////
    // 0.5.4. Tracking immunity over time

    vector<double> A_par_mean_t;
    vector<double> A_clin_mean_t;
};


////////////////////////////////////////////////////////////
//                                                        //
// 0.6. Function declarations                             //
//                                                        //
////////////////////////////////////////////////////////////

void mosq_derivs(const double t, double(&yM)[N_spec][N_M_comp], double(&dyMdt)[N_spec][N_M_comp], Params& theta, Population& POP);
void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_spec][N_M_comp], Params& theta, Population& POP);
void mosquito_step(double t, Params& theta, Population& POP);
void human_step(Params& theta, Population& POP);
void intervention_dist(double t, Params& theta, Population& POP, Intervention& INTVEN);
void POP_summary(Population& POP);
void model_simulator(Params& theta, Population& POP, Intervention& INTVEN, simulation& SIM);
double phi_inv(double pp, double mu, double sigma);


////////////////////////////////////////////
////////////////////////////////////////////
//        //                              //
//   ##   //  #     #  ####  #### #   ##  //
//  ###   //  ##   ## ##  ##  ##  ##  ##  //
//   ##   //  ####### ######  ##  ### ##  //
//   ##   //  ## # ## ##  ##  ##  ## ###  //
//  ####  //  ##   ## ##  ## #### ##  ##  //
//        //                              //
////////////////////////////////////////////
////////////////////////////////////////////

int main(int argc, char** argv)
{
    setall(time(NULL), 7);

    clock_t clock_time;
    clock_time = clock();


    ////////////////////////////////////////////
    //                                        //  
    //  1.1 Read in file names                //
    //                                        //
    ////////////////////////////////////////////

    // do we have the correct command line?
    if (argc != 4 + N_spec_max)
    {
        std::cout << "Incorrect command line.\n";
        return 0;
    }

    const char* parameter_File = argv[1];

    const char* mosquito_File[N_spec_max];
    for (int g = 0; g < N_spec_max; g++)
    {
        mosquito_File[g] = argv[2 + g];
    }

    const char* coverage_File = argv[5];
    const char* output_File = argv[6];


    /*
    char* parameter_File = "C:/U/Pv_mod/IB_mod/Mod_4/model_parameters.txt";

    char* mosquito_File[N_spec_max];
    mosquito_File[0] = "C:/U/Pv_mod/IB_mod/Mod_4/farauti_parameters.txt";
    mosquito_File[1] = "C:/U/Pv_mod/IB_mod/Mod_4/punctulatus_parameters.txt";
    mosquito_File[2] = "C:/U/Pv_mod/IB_mod/Mod_4/koliensis_parameters.txt";

    char* coverage_File = "C:/U/Pv_mod/IB_mod/Mod_4/intervention_coverage.txt";
    char* output_File = "C:/U\Pv_mod/IB_mod\Mod_4/Output/model_output.txt";
    */


    ////////////////////////////////////////////
    //                                        //
    // 1.2. Initialise objects                //
    //                                        //
    ////////////////////////////////////////////

    Population PNG_pop;
    Params Pv_mod_par;

    Pv_mod_par.set(parameter_File, mosquito_File);
    PNG_pop.N_pop = Pv_mod_par.N_pop;


    /////////////////////////////////////////////////////////////////////////
    //                                                                     //
    // 1.7. Read in intervention coverage                                  //
    //                                                                     //
    /////////////////////////////////////////////////////////////////////////

    Intervention PNG_intven;
    PNG_intven.set(coverage_File);


    ///////////////////////////////////////////////////////////////////////////
    //                                                                       //
    // 1.8. Initialise Population of individuals                             //
    //      Note that they begin with exponential age distribution           //
    //      and susceptible without immunity                                 //
    //                                                                       // 
    ///////////////////////////////////////////////////////////////////////////

    cout << "Initialise population of individuals for simulation at equilbirium EIR of " << 365.0*Pv_mod_par.EIR_equil << endl;
    cout << endl;

    PNG_pop.setup_equilibrium(Pv_mod_par);

    cout << "Population of size " << PNG_pop.N_pop << " initialised!" << endl;
    cout << endl;


    /////////////////////////////////////////////////////////////////////////
    //                                                                     //
    // 1.9. Create simulation object                                       //
    //                                                                     //
    /////////////////////////////////////////////////////////////////////////

    simulation PNG_sim;


    /////////////////////////////////////////////////////////////////////////
    // 1.9.1. Vector of simulation times

    int N_time = Pv_mod_par.N_time;
    PNG_sim.N_time = N_time;

    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.t_vec.push_back((double)(
              Pv_mod_par.time_start * 365
            - Pv_mod_par.burnin_time * 365 + i*t_step));
    }


    /////////////////////////////////////////////////////////////////////////
    // 1.9.2. Create storage for output

    PNG_sim.yH_t.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.yH_t[i].resize(N_H_comp);
    }


    PNG_sim.yM_t.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.yM_t[i].resize(N_spec);
        for (int g = 0; g < N_spec; g++)
        {
            PNG_sim.yM_t[i][g].resize(N_M_comp);
        }
    }


    PNG_sim.prev_all.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.prev_all[i].resize(11);
    }

    PNG_sim.prev_U5.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.prev_U5[i].resize(11);
    }

    PNG_sim.prev_U10.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.prev_U10[i].resize(11);
    }

    PNG_sim.EIR_t.resize(N_time);

    PNG_sim.LLIN_cov_t.resize(N_time);
    PNG_sim.IRS_cov_t.resize(N_time);
    PNG_sim.ACT_treat_t.resize(N_time);
    PNG_sim.PQ_treat_t.resize(N_time);
    PNG_sim.pregnant_t.resize(N_time);

    PNG_sim.A_par_mean_t.resize(N_time);
    PNG_sim.A_clin_mean_t.resize(N_time);


    //////////////////////////////////////////////////////
    //                                                  //
    // 1.10. Begin stochastic simulations               //
    //                                                  //
    ////////////////////////////////////////////////////// 

    cout << "Starting model simulations......." << endl;

    model_simulator(Pv_mod_par, PNG_pop, PNG_intven, PNG_sim);

    cout << "Model simulations completed....." << endl;
    cout << endl;


    //////////////////////////////////////////////////////
    //                                                  //
    // 1.11. Output to file                             //
    //                                                  //
    ////////////////////////////////////////////////////// 

    cout << "Start writing output to file......" << endl;
    cout << endl;

    ofstream output_Stream(output_File);

    for (int i = 0; i<N_time; i++)
    {
        output_Stream << PNG_sim.t_vec[i] << "\t";

        for (int k = 0; k<N_H_comp; k++)
        {
            output_Stream << PNG_sim.yH_t[i][k] << "\t";
        }

        for (int g = 0; g < N_spec; g++)
        {
            for (int k = 0; k < N_M_comp; k++)
            {
                output_Stream << PNG_sim.yM_t[i][g][k] << "\t";
            }
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << PNG_sim.prev_all[i][k] << "\t";
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << PNG_sim.prev_U5[i][k] << "\t";
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << PNG_sim.prev_U10[i][k] << "\t";
        }

        output_Stream << PNG_sim.EIR_t[i] << "\t";
        output_Stream << PNG_sim.LLIN_cov_t[i] << "\t";
        output_Stream << PNG_sim.IRS_cov_t[i] << "\t";
        output_Stream << PNG_sim.ACT_treat_t[i] << "\t";
        output_Stream << PNG_sim.PQ_treat_t[i] << "\t";
        output_Stream << PNG_sim.pregnant_t[i] << "\t";

        output_Stream << PNG_sim.A_par_mean_t[i] << "\t";
        output_Stream << PNG_sim.A_clin_mean_t[i] << "\t";

        output_Stream << endl;
    }

    output_Stream.close();


    cout << "Output successfully written to file......" << endl;
    cout << endl;


    cout << "Time taken: " << ( (double) clock() - clock_time)/( (double) CLOCKS_PER_SEC ) << " seconds" << endl;


    return 0;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//          //                                                              //
//   ####   //  ##### ##  ## #   ##  ####  ###### ####  ####  #   ##  ###   //
//  ##  ##  //  ##    ##  ## ##  ## ##  ##   ##    ##  ##  ## ##  ## ##     // 
//     ##   //  ####  ##  ## ### ## ##       ##    ##  ##  ## ### ##  ###   //
//    ##    //  ##    ##  ## ## ### ##  ##   ##    ##  ##  ## ## ###    ##  //
//   #####  //  ##     ####  ##  ##  ####    ##   ####  ####  ##  ##  ###   //
//          //                                                              //
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.1. Derivatives of mosquito ODE model                                  //
//                                                                          //
//  0.01721421 = 2*pi/365                                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_derivs(const double t, double (&yM)[N_spec][N_M_comp], double (&dyMdt)[N_spec][N_M_comp], Params& theta, Population& POP)
{
    double Karry_seas_inv[N_spec];

    for (int g = 0; g < N_spec; g++)
    {
        Karry_seas_inv[g] = 1.0 / (theta.Karry[g] * ( theta.dry_seas[g] + (1 - theta.dry_seas[g])*pow(0.5 + 0.5*cos(0.01721421*(t - theta.t_peak_seas[g])), theta.kappa_seas[g])/ theta.denom_seas[g] ) );

        //Karry_seas_inv[g] = 1.0/theta.Karry[g];

        dyMdt[g][0] = POP.beta_VC[g] * (yM[g][3] + yM[g][4] + yM[g][5]) - yM[g][0] / theta.d_E_larvae - yM[g][0] * theta.mu_E0*(1.0 + (yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
        dyMdt[g][1] = yM[g][0] / theta.d_E_larvae - yM[g][1] / theta.d_L_larvae - yM[g][1] * theta.mu_L0*(1.0 + theta.gamma_larvae*(yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
        dyMdt[g][2] = yM[g][1] / theta.d_L_larvae - yM[g][2] / theta.d_pupae - yM[g][2] * theta.mu_P;
        dyMdt[g][3] = 0.5*yM[g][2] / theta.d_pupae - theta.lam_M[g] * yM[g][3] - POP.mu_M_VC[g] * yM[g][3];
        dyMdt[g][4] = + theta.lam_M[g] * yM[g][3] - theta.lam_S_M_track[g][0] * POP.exp_muM_tauM_VC[g] - POP.mu_M_VC[g] * yM[g][4];
        dyMdt[g][5] =                              + theta.lam_S_M_track[g][0] * POP.exp_muM_tauM_VC[g] - POP.mu_M_VC[g] * yM[g][5];
    }
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.2. Runge-Kutta 4 step updater for mosquito model                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_rk4(const double t, const double t_step_mosq, double (&yM)[N_spec][N_M_comp], Params& theta, Population& POP)
{
    double k1_yM[N_spec][N_M_comp], k2_yM[N_spec][N_M_comp], k3_yM[N_spec][N_M_comp], k4_yM[N_spec][N_M_comp], yM_temp[N_spec][N_M_comp];


    //////////////////////////
    // step 1

    mosq_derivs(t, yM, k1_yM, theta, POP);


    //////////////////////////
    // step 2

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k < N_M_comp; k++)
        {
            yM_temp[g][k] = yM[g][k] + 0.5*t_step_mosq*k1_yM[g][k];
        }
    }

    mosq_derivs(t + 0.5*t_step_mosq, yM_temp, k2_yM, theta, POP);


    //////////////////////////
    // step 3

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k < N_M_comp; k++)
        {
            yM_temp[g][k] = yM[g][k] + 0.5*t_step_mosq*k2_yM[g][k];
        }
    }

    mosq_derivs(t + 0.5*t_step_mosq, yM_temp, k3_yM, theta, POP);


    //////////////////////////
    // step 4

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k < N_M_comp; k++)
        {
            yM_temp[g][k] = yM[g][k] + t_step_mosq*k3_yM[g][k];
        }
    }

    mosq_derivs(t + t_step_mosq, yM_temp, k4_yM, theta, POP);


    //////////////////////////
    // output

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k < N_M_comp; k++)
        {
            yM[g][k] = yM[g][k] + t_step_mosq*k1_yM[g][k] / 6.0 + t_step_mosq*k2_yM[g][k] / 3.0 + t_step_mosq*k3_yM[g][k] / 3.0 + t_step_mosq*k4_yM[g][k] / 6.0;
        }
    }

}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  2.3. Update step for mosquitoes                                          //
//                                                                           // 
//       For every human step we take mosq_steps (=10) steps for mosquitoes  //
//       The smaller step size ensures that the ODE solver works smoothly.   //
//       Especially an issue for the larval stages                           //   
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void mosquito_step(double t, Params& theta, Population& POP)
{
    //////////////////////////////////
    // Set up mosquito state vector

    double yM[N_spec][N_M_comp];

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k<N_M_comp; k++)
        {
            yM[g][k] = POP.yM[g][k];
        }
    }


    double t_step_mosq = (double(t_step)) / (double(mosq_steps));


    //////////////////////////////////
    // Force of infection on mosquitoes

    for (int g = 0; g < N_spec; g++)
    {
        theta.lam_M[g] = 0.0;
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            theta.lam_M[g] = theta.lam_M[g] + POP.lam_n[n][g] * (theta.c_PCR*POP.people[n].I_PCR + theta.c_LM*POP.people[n].I_LM +
                                                                    theta.c_D*POP.people[n].I_D + theta.c_T*POP.people[n].T);
        }
    }


    //////////////////////////////////////
    // Carry out the mosq_steps

    for (int j = 0; j<mosq_steps; j++)
    {
        mosq_rk4(t, t_step_mosq, yM, theta, POP);

        for (int g = 0; g < N_spec; g++)
        {
            theta.lam_S_M_track[g].push_back(theta.lam_M[g]* yM[g][3]);
            theta.lam_S_M_track[g].erase(theta.lam_S_M_track[g].begin());
        }
    }

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k < N_M_comp; k++)
        {
            POP.yM[g][k] = yM[g][k];
        }
    }

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.4. Update the vector of human classes                                 //
//                                                                          // 
//       THINK CAREFULLY ABOUT THE ORDERING OF EVENTS                       //
//////////////////////////////////////////////////////////////////////////////

void human_step(Params& theta, Population& POP)
{

    //////////////////////////////////////////////////////////////////////////
    // 2.4.1 Temporary objects for setting up individuals' intervention   
    //       access characteristics 

    float GMN_parm[(N_int)*(N_int + 3) / 2 + 1];
    float GMN_work[N_int];
    float GMN_zero[N_int];
    float zz_GMN[N_int];

    for (int k = 0; k<N_int; k++)
    {
        GMN_zero[k] = 0.0;
    }


    ///////////////////////////////////////////////
    // 2.4.2. Apply ageing

    for (int n = 0; n<POP.N_pop; n++)
    {
        POP.people[n].ager(theta);
    }


    ///////////////////////////////////////////////
    // 2.4.3. Deaths
    //
    // Look again at how things are erased from vectors.

    int N_dead = 0;

    for (size_t n = 0; n < POP.people.size(); n++)
    {
        /////////////////////////////////////////////
        // Everyone has an equal probability of dying

        if (theta.P_dead > genunf(0, 1))
        {
            POP.people.erase(POP.people.begin() + n);
            
            POP.pi_n.erase(POP.pi_n.begin() + n);
            POP.lam_n.erase(POP.lam_n.begin() + n);

            N_dead = N_dead + 1;
            n = n - 1;      // If we erase something, the next one moves into it's place so we don't want to step forward.
        }
        else {

            ///////////////////////////////////////////
            // People die once they reach the maximum age

            if (POP.people[n].age > theta.age_max)
            {
                POP.people.erase(POP.people.begin() + n);

                POP.pi_n.erase(POP.pi_n.begin() + n);
                POP.lam_n.erase(POP.lam_n.begin() + n);

                N_dead = N_dead + 1;
                n = n - 1;       // If we erase something, the next one moves into it's place so we don't want to step forward.
            }
        }
    }

    /////////////////////////////////////////////////////////////
    // 2.4.4. Births - set up to ensure balanced population.
    //        Can be adjusted to account for changing demography.

    double zeta_start, het_dif_track, q_rand;

    vector<double> zero_push(N_spec);
    for (int g = 0; g < N_spec; g++)
    {
        zero_push[g] = 0.0;
    }


    for (int n = 0; n<N_dead; n++)
    {
        zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));

        while (zeta_start > theta.het_max)
        {
            zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));
        }

        Human HH(0.0, zeta_start);

        HH.S = 1;
        HH.I_PCR = 0;
        HH.I_LM = 0;
        HH.I_D = 0;
        HH.T = 0;
        HH.P = 0;


        HH.A_par = 0.0;
        HH.A_clin = 0.0;

        HH.A_par_boost = 0;
        HH.A_clin_boost = 0;

        HH.A_par_timer = -1.0;
        HH.A_clin_timer = -1.0;

        HH.PQ_proph = 0;
        HH.PQ_proph_timer = -1.0;

        HH.Hyp = 0;

        if (genunf(0.0, 1.0) < 0.5)
        {
            HH.gender = 0;
        }
        else {
            HH.gender = 1;
        }

        if (HH.gender == 0)
        {
            if (genunf(0.0, 1.0) < theta.G6PD_prev)
            {
                HH.G6PD_def = 1;
            } else {
                HH.G6PD_def = 0;
            }
        } else {

            q_rand = genunf(0.0, 1.0);

            if(q_rand <= theta.G6PD_prev*theta.G6PD_prev)
            {
                HH.G6PD_def = 2;
            }

            if ((q_rand > theta.G6PD_prev*theta.G6PD_prev) && (q_rand <= theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
            {
                HH.G6PD_def = 1;
            }

            if (q_rand > (theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
            {
                HH.G6PD_def = 0;
            }
        }

        if (genunf(0.0, 1.0) < theta.CYP2D6_prev)
        {
            HH.CYP2D6 = 1;
        }
        else {
            HH.CYP2D6 = 0;
        }

        HH.preg_age = 0;
        HH.pregnant = 0;
        HH.preg_timer = 0.0;


        /////////////////////////////////
        // Assign levels of maternally-acquired immunity
        // by finding women of child-bearing age with the
        // closest level of heterogeneity

        HH.A_par_mat = 0.0;
        HH.A_clin_mat = 0.0;

        het_dif_track = 1e10;

        for (size_t j = 0; j<POP.people.size(); j++)
        {
            if (POP.people[j].preg_age == 1)
            {
                if (abs(HH.zeta_het - POP.people[j].zeta_het) < het_dif_track)
                {
                    HH.A_par_mat = theta.P_mat*POP.people[j].A_par_mat;
                    HH.A_clin_mat = theta.P_mat*POP.people[j].A_clin_mat;

                    het_dif_track = (HH.zeta_het - POP.people[j].zeta_het)*(HH.zeta_het - POP.people[j].zeta_het);
                }
            }
        }


        ///////////////////////////////////////////////////
        // Lagged exposure equals zero - they're not born yet!

        for (int k = 0; k<theta.H_track; k++)
        {
            HH.lam_bite_track.push_back(0.0);
            HH.lam_rel_track.push_back(0.0);
        }


        ///////////////////////////////////////////////////
        // Assign intervention access scores

        for (int p = 0; p<N_int; p++)
        {
            for (int q = 0; q<N_int; q++)
            {
                theta.V_int_dummy[p][q] = theta.V_int[p][q];
            }
        }

        setgmn(GMN_zero, *theta.V_int_dummy, N_int, GMN_parm);

        genmn(GMN_parm, zz_GMN, GMN_work);

        for (int k = 0; k<N_int; k++)
        {
            HH.zz_int[k] = zz_GMN[k];
        }


        ///////////////////////////////////////////////////
        // Born with no interventions

        HH.LLIN = 0;
        HH.IRS = 0;

        for (int g = 0; g < N_spec; g++)
        {
            HH.w_VC[g] = 1.0;
            HH.y_VC[g] = 1.0;
            HH.z_VC[g] = 0.0;
        }


        /////////////////////////////////////////////////////////////////
        // 2.4.5. Push the created individual onto the vector of people

        POP.people.push_back(move(HH));

        POP.pi_n.push_back(zero_push);
        POP.lam_n.push_back(zero_push);
    }



    ///////////////////////////////////////////////////
    // 2.4.6. Update individual-level vector control

    for (int n = 0; n<POP.N_pop; n++)
    {
        POP.people[n].intervention_updater(theta);
    }


    ///////////////////////////////////////////////////
    // 2.4.7. Update proportion of bites
    //
    //        Note the ordering of n and g loops. Need to 
    //        check if this makes a difference for speed.
    //
    //        Should be able to make this quicker


    for (int n = 0; n<POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.pi_n[n][g] = POP.people[n].zeta_het*(1.0 - theta.rho_age*exp(-POP.people[n].age*theta.age_0_inv));

            //POP.pi_n[n][g] = POP.people[n].zeta_het - (POP.people[n].zeta_het - POP.people[n].zeta_het)*POP.P_age_bite;   // Slightly quicker - no calling of exponentials
        }
    }

    double SIGMA_PI[N_spec];
    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 0.0;
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            SIGMA_PI[g] = SIGMA_PI[g] + POP.pi_n[n][g];
        }
    }

    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 1.0 / SIGMA_PI[g];
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.pi_n[n][g] = POP.pi_n[n][g] * SIGMA_PI[g];
        }
    }


    ///////////////////////////////////////////////////
    // 2.4.8 Update population-level vector control quantities

    for (int g = 0; g < N_spec; g++)
    {
        POP.SUM_pi_w[g] = 0;
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.SUM_pi_w[g] = POP.SUM_pi_w[g] + POP.pi_n[n][g] * POP.people[n].w_VC[g];
        }
    }


    for (int g = 0; g < N_spec; g++)
    {
        POP.W_VC[g] = 1.0 - theta.Q_0[g] + theta.Q_0[g] * POP.SUM_pi_w[g];
        POP.Z_VC[g] = theta.Q_0[g] * POP.SUM_pi_z[g];

        POP.delta_1_VC[g] = theta.delta_1 / (1.0 - POP.Z_VC[g]);
        POP.delta_VC[g] = POP.delta_1_VC[g] + theta.delta_2;

        POP.p_1_VC[g] = theta.p_1[g] * POP.W_VC[g] / (1.0 - POP.Z_VC[g] * theta.p_1[g]);

        POP.mu_M_VC[g] = -log(POP.p_1_VC[g] * theta.p_2[g]) / POP.delta_VC[g];

        POP.Q_VC[g] = 1.0 - (1.0 - theta.Q_0[g]) / POP.W_VC[g];

        POP.aa_VC[g] = POP.Q_VC[g] / POP.delta_VC[g];

        POP.exp_muM_tauM_VC[g] = exp(-POP.mu_M_VC[g] * theta.tau_M[g]);
        POP.beta_VC[g] = theta.eps_max[g] * POP.mu_M_VC[g] / (exp(POP.delta_VC[g] * POP.mu_M_VC[g]) - 1.0);
    }


    ///////////////////////////////////////////////////
    // 2.4.9. Update individual-level force of infection on humans

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.lam_n[n][g] = POP.aa_VC[g] * POP.pi_n[n][g] * POP.people[n].w_VC[g] / POP.SUM_pi_w[g];
        }
    }


    ///////////////////////////////////////////////////
    // 2.4.10. Implement moves between compartments
    //
    // TO DO: Can take some multiplications out of the loop.

    double lam_bite_base[N_spec];
    double lam_bite_n;     // better notation (this is lam_bite)

    for (int g = 0; g < N_spec; g++)
    {
        lam_bite_base[g] = (double(POP.N_pop))*theta.bb*POP.yM[g][5];
    }

    for (int n = 0; n<POP.N_pop; n++)
    {
        lam_bite_n = 0.0;

        for (int g = 0; g < N_spec; g++)
        {
            lam_bite_n = lam_bite_n + POP.lam_n[n][g] * lam_bite_base[g];
        }

        POP.people[n].state_mover(theta, lam_bite_n);
    }

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.5. Summarise the output from the population                           //
//                                                                          // 
//////////////////////////////////////////////////////////////////////////////

void POP_summary(Population& POP)
{
    for (int k = 0; k<N_H_comp; k++)
    {
        POP.yH[k] = 0.0;
    }

    for (int k = 0; k<10; k++)
    {
        POP.prev_all[k] = 0.0;
        POP.prev_U5[k] = 0.0;
        POP.prev_U10[k] = 0.0;
    }


    for (int n = 0; n<POP.N_pop; n++)
    {
        ////////////////////////////////////////
        // Numbers in each compartment

        POP.yH[0] = POP.yH[0] + POP.people[n].S;
        POP.yH[1] = POP.yH[1] + POP.people[n].I_PCR;
        POP.yH[2] = POP.yH[2] + POP.people[n].I_LM;
        POP.yH[3] = POP.yH[3] + POP.people[n].I_D;
        POP.yH[4] = POP.yH[4] + POP.people[n].T;
        POP.yH[5] = POP.yH[5] + POP.people[n].P;


        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - full population

        ////////////////////////////////////////
        // Prevalence

        POP.prev_all[0] = POP.prev_all[0] + 1;                                                                        // Numbers - denominator
        POP.prev_all[1] = POP.prev_all[1] + POP.people[n].I_PCR + POP.people[n].I_LM +
                                            + POP.people[n].I_D + POP.people[n].T;                                      // PCR detectable infections
        POP.prev_all[2] = POP.prev_all[2] + POP.people[n].I_LM + POP.people[n].I_D + POP.people[n].T;                // LM detectable infections
        POP.prev_all[3] = POP.prev_all[3] + POP.people[n].I_D + POP.people[n].T;                                      // Clinical episodes

        if (POP.people[n].Hyp > 0)
        {
            POP.prev_all[4] = POP.prev_all[4] + 1;                     // Hypnozoite positive

            POP.prev_all[5] = POP.prev_all[5] + POP.people[n].Hyp;    // Number of batches of hypnozoites
        }


        ////////////////////////////////////////
        // Incidence

        POP.prev_all[6]  = POP.prev_all[6]  + POP.people[n].I_PCR_new;
        POP.prev_all[7]  = POP.prev_all[7]  + POP.people[n].I_LM_new;
        POP.prev_all[8]  = POP.prev_all[8]  + POP.people[n].I_D_new;
        POP.prev_all[9]  = POP.prev_all[9]  + POP.people[n].ACT_new;
        POP.prev_all[10] = POP.prev_all[10] + POP.people[n].PQ_new;


        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - under 5's

        if (POP.people[n].age < 1825.0)
        {
            ////////////////////////////////////////
            // Prevalence

            POP.prev_U5[0] = POP.prev_U5[0] + 1;                                                                // Numbers - denominator
            POP.prev_U5[1] = POP.prev_U5[1] + POP.people[n].I_PCR + POP.people[n].I_LM
                                              + POP.people[n].I_D + POP.people[n].T;                              // PCR detectable infections
            POP.prev_U5[2] = POP.prev_U5[2] + POP.people[n].I_LM + POP.people[n].I_D + POP.people[n].T;        // LM detectable infections
            POP.prev_U5[3] = POP.prev_U5[3] + POP.people[n].I_D + POP.people[n].T;                              // Clinical episodes

            if (POP.people[n].Hyp > 0)
            {
                POP.prev_U5[4] = POP.prev_U5[4] + 1;                     // Hypnozoite positive

                POP.prev_U5[5] = POP.prev_U5[5] + POP.people[n].Hyp;    // Number of batches of hypnozoites
            }


            ////////////////////////////////////////
            // Incidence

            POP.prev_U5[6]  = POP.prev_U5[6]  + POP.people[n].I_PCR_new;
            POP.prev_U5[7]  = POP.prev_U5[7]  + POP.people[n].I_LM_new;
            POP.prev_U5[8]  = POP.prev_U5[8]  + POP.people[n].I_D_new;
            POP.prev_U5[9]  = POP.prev_U5[9]  + POP.people[n].ACT_new;
            POP.prev_U5[10] = POP.prev_U5[10] + POP.people[n].PQ_new;
        }

        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - under 10's

        if (POP.people[n].age < 3650.0)
        {
            ////////////////////////////////////////
            // Prevalence

            POP.prev_U10[0] = POP.prev_U10[0] + 1;                                                            // Numbers - denominator
            POP.prev_U10[1] = POP.prev_U10[1] + POP.people[n].I_PCR + POP.people[n].I_LM
                                                + POP.people[n].I_D + POP.people[n].T;                          // PCR detectable infections
            POP.prev_U10[2] = POP.prev_U10[2] + POP.people[n].I_LM + POP.people[n].I_D + POP.people[n].T;    // LM detectable infections
            POP.prev_U10[3] = POP.prev_U10[3] + POP.people[n].I_D + POP.people[n].T;                          // Clinical episodes

            if (POP.people[n].Hyp > 0)
            {
                POP.prev_U10[4] = POP.prev_U10[4] + 1;                     // Hypnozoite positive

                POP.prev_U10[5] = POP.prev_U10[5] + POP.people[n].Hyp;    // Number of batches of hypnozoites
            }


            ////////////////////////////////////////
            // Incidence

            POP.prev_U10[6]  = POP.prev_U10[6]  + POP.people[n].I_PCR_new;
            POP.prev_U10[7]  = POP.prev_U10[7]  + POP.people[n].I_LM_new;
            POP.prev_U10[8]  = POP.prev_U10[8]  + POP.people[n].I_D_new;
            POP.prev_U10[9]  = POP.prev_U10[9]  + POP.people[n].ACT_new;
            POP.prev_U10[10] = POP.prev_U10[10] + POP.people[n].PQ_new;
        }
    }


    //////////////////////////////
    // Intervention coverage

    POP.LLIN_cov_t = 0;
    POP.IRS_cov_t = 0;
    POP.ACT_treat_t = 0;
    POP.PQ_treat_t = 0;
    POP.pregnant_t = 0;

    for (int n = 0; n<POP.N_pop; n++)
    {
        POP.LLIN_cov_t  = POP.LLIN_cov_t  + POP.people[n].LLIN;
        POP.IRS_cov_t   = POP.IRS_cov_t   + POP.people[n].IRS;
        POP.ACT_treat_t = POP.ACT_treat_t + POP.people[n].ACT_new;
        POP.PQ_treat_t  = POP.PQ_treat_t + POP.people[n].PQ_new;
        POP.pregnant_t  = POP.pregnant_t  + POP.people[n].pregnant;
    }


    //////////////////////////////
    // Immunity

    double A_par_mean = 0.0, A_clin_mean = 0.0;

    for (int n = 0; n<POP.N_pop; n++)
    {
        A_par_mean = A_par_mean + POP.people[n].A_par;
        A_clin_mean = A_clin_mean + POP.people[n].A_clin;
    }

    POP.A_par_mean_t = A_par_mean / ((double)POP.N_pop);
    POP.A_clin_mean_t = A_clin_mean / ((double)POP.N_pop);
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.6. Simulate the model and store the output in SIM                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void model_simulator(Params& theta, Population& POP, Intervention& INTVEN, simulation& SIM)
{

    for (int i = 0; i<SIM.N_time; i++)
    {
        if (SIM.t_vec[i] / 365.0 - floor(SIM.t_vec[i] / 365.0) < 0.5*t_step / 365.0)
        {
            cout << "time = " << SIM.t_vec[i] / 365.0 << "\t" << 100.0*(SIM.t_vec[i] - SIM.t_vec[0]) / (double(t_step*SIM.N_time)) << "% complete" << endl;
        }

        human_step(theta, POP);

        mosquito_step(SIM.t_vec[i], theta, POP);

        intervention_dist(SIM.t_vec[i], theta, POP, INTVEN);

        POP_summary(POP);

        //////////////////////////////////////
        // Fill out simulation object

        for (int k = 0; k<N_H_comp; k++)
        {
            SIM.yH_t[i][k] = POP.yH[k];
        }

        for (int k = 0; k<N_M_comp; k++)
        {
            for (int g = 0; g < N_spec; g++)
            {
                SIM.yM_t[i][g][k] = POP.yM[g][k];
            }
        }

        for (int k = 0; k<11; k++)
        {
            SIM.prev_all[i][k] = POP.prev_all[k];
            SIM.prev_U5[i][k]  = POP.prev_U5[k];
            SIM.prev_U10[i][k] = POP.prev_U10[k];
        }


        SIM.LLIN_cov_t[i]  = POP.LLIN_cov_t;
        SIM.IRS_cov_t[i]   = POP.IRS_cov_t;
        SIM.ACT_treat_t[i] = POP.ACT_treat_t;
        SIM.PQ_treat_t[i]  = POP.PQ_treat_t;
        SIM.pregnant_t[i]  = POP.pregnant_t;

        SIM.EIR_t[i] = 0.0;
        for (int g = 0; g < N_spec; g++)
        {
            SIM.EIR_t[i] = SIM.EIR_t[i] + POP.aa_VC[g] * POP.yM[g][5];
        }

        SIM.A_par_mean_t[i] = POP.A_par_mean_t;
        SIM.A_clin_mean_t[i] = POP.A_clin_mean_t;

    }

}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.7. Vector control distribution                                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


void intervention_dist(double t, Params& theta, Population& POP, Intervention& INTVEN)
{
    double QQ;

    //////////////////////////////////////////////////////////
    // Intervention 1: LLINS

    for (size_t m = 0; m<INTVEN.LLIN_year.size(); m++)
    {
        if ((t > INTVEN.LLIN_year[m] - 0.5*t_step) &&
            (t < INTVEN.LLIN_year[m] + 0.51*t_step))
        {
            cout << "LLIN distribution" << endl;

            QQ = phi_inv(INTVEN.LLIN_cover[m], 0.0, sqrt(1.0 + theta.sig_round_LLIN*theta.sig_round_LLIN));

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[0], theta.sig_round_LLIN) < QQ)
                {
                    POP.people[n].LLIN = 1;
                    POP.people[n].LLIN_age = 0.0;

                    for (int g = 0; g < N_spec; g++)
                    {
                        POP.people[n].d_LLIN[g] = theta.d_LLIN_0[g];
                        POP.people[n].r_LLIN[g] = theta.r_LLIN_0[g];
                        POP.people[n].s_LLIN[g] = 1.0 - POP.people[n].d_LLIN[g] - POP.people[n].r_LLIN[g];
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 2: IRS

    for (size_t m = 0; m<INTVEN.IRS_year.size(); m++)
    {
        if ((t > INTVEN.IRS_year[m] - 0.5*t_step) &&
            (t < INTVEN.IRS_year[m] + 0.51*t_step))
        {
            cout << "IRS distribution" << endl;

            QQ = phi_inv(INTVEN.IRS_cover[m], 0.0, sqrt(1.0 + theta.sig_round_IRS*theta.sig_round_IRS));

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[1], theta.sig_round_IRS) < QQ)
                {
                    POP.people[n].IRS = 1;
                    POP.people[n].IRS_age = 0.0;

                    for (int g = 0; g < N_spec; g++)
                    {
                        POP.people[n].d_IRS[g] = theta.d_IRS_0[g];
                        POP.people[n].r_IRS[g] = theta.r_IRS_0[g];
                        POP.people[n].s_IRS[g] = 1.0 - POP.people[n].d_IRS[g] - POP.people[n].r_IRS[g];
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 3: MDA (blood-stage)

    for (size_t m = 0; m<INTVEN.MDA_BS_year.size(); m++)
    {
        if ((t > INTVEN.MDA_BS_year[m] - 0.5*t_step) &&
            (t < INTVEN.MDA_BS_year[m] + 0.51*t_step))
        {
            cout << "MDA (BS) distribution" << endl;

            theta.MDA_BS_cover   = INTVEN.MDA_BS_cover[m];
            theta.MDA_BS_BSeff   = INTVEN.MDA_BS_BSeff[m];
            theta.MDA_BS_BSproph = INTVEN.MDA_BS_BSproph[m];

            QQ = phi_inv(theta.MDA_BS_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[2], theta.sig_round_MDA) < QQ)
                {
                    POP.people[n].ACT_new = 1;

                    if (genunf(0.0, 1.0) < theta.MDA_BS_BSeff)
                    {
                        if (gennor(POP.people[n].zz_int[2], theta.sig_round_MDA) < QQ)
                        {
                            if (POP.people[n].S == 1    ) { POP.people[n].S = 0;     POP.people[n].P = 1; }
                            if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1;  }
                            if (POP.people[n].I_LM == 1 ) { POP.people[n].I_LM = 0;  POP.people[n].P = 1;  }
                            if (POP.people[n].I_D == 1  ) { POP.people[n].I_D = 0;   POP.people[n].T = 1;  }
                        }
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 4: MDA (blood-stage and liver-stage)

    for (size_t m = 0; m<INTVEN.MDA_PQ_year.size(); m++)
    {
        if ((t > INTVEN.MDA_PQ_year[m] - 0.5*t_step) &&
            (t < INTVEN.MDA_PQ_year[m] + 0.51*t_step))
        {
            cout << "MDA (BS+PQ) distribution" << endl;

            theta.MDA_PQ_cover   = INTVEN.MDA_PQ_cover[m];
            theta.MDA_PQ_BSeff   = INTVEN.MDA_PQ_BSeff[m];
            theta.MDA_PQ_PQeff   = INTVEN.MDA_PQ_PQeff[m];
            theta.MDA_PQ_BSproph = INTVEN.MDA_PQ_BSproph[m];
            theta.MDA_PQ_PQproph = INTVEN.MDA_PQ_PQproph[m];
            theta.MDA_PQ_CYP2D6  = INTVEN.MDA_PQ_CYP2D6[m];

            QQ = phi_inv(theta.MDA_PQ_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[3], theta.sig_round_MDA) < QQ)
                {
                    POP.people[n].ACT_new = 1;

                    if (genunf(0.0, 1.0) < theta.MDA_PQ_BSeff)
                    {
                        if (POP.people[n].S == 1    ) { POP.people[n].S = 0;     POP.people[n].P = 1; }
                        if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
                        if (POP.people[n].I_LM == 1 ) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
                        if (POP.people[n].I_D == 1  ) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
                    }

                    if( (POP.people[n].G6PD_def == 0) && (POP.people[n].pregnant == 0) && (POP.people[n].age > 180.0) )
                    {
                        POP.people[n].PQ_new = 1;

                        if (theta.MDA_PQ_CYP2D6 == 0)    // Is CYP2D6 low metabolization a problem? No = 0, e.g. TQ; Otherwise Yes = 1, e.g. PQ
                        {
                            if (genunf(0.0, 1.0) < theta.MDA_PQ_PQeff)
                            {
                                POP.people[n].Hyp = 0;

                                POP.people[n].PQ_proph = 1;
                                POP.people[n].PQ_proph_timer = theta.PQ_treat_PQproph;
                            }
                        }else{
                            if (POP.people[n].CYP2D6 == 0)          // If CYP2D6 low metabolization is a problem - it only effects the low metabolizers
                            {
                                if (genunf(0.0, 1.0) < theta.MDA_PQ_PQeff)
                                {
                                    POP.people[n].Hyp = 0;

                                    POP.people[n].PQ_proph = 1;
                                    POP.people[n].PQ_proph_timer = theta.PQ_treat_PQproph;
                                }
                            }
                        }

                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 7: first-line treatment (blood-stage)

    for (size_t m = 0; m<INTVEN.BS_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.BS_treat_year_on[m] - 0.5*t_step) &&
            (t < INTVEN.BS_treat_year_on[m] + 0.51*t_step))
        {
            cout << "New front-line BS treatment" << endl;

            theta.BS_treat_cover   = INTVEN.BS_treat_cover[m];
            theta.BS_treat_BSeff   = INTVEN.BS_treat_BSeff[m];
            theta.BS_treat_BSproph = INTVEN.BS_treat_BSproph[m];

            theta.treat_cov = theta.BS_treat_cover;
            theta.treat_eff = theta.BS_treat_BSeff;
            theta.r_P = 1.0 / theta.BS_treat_BSproph;
        }
    }


    //////////////////////////////////////////////////////////
    // Switching back to baseline.

    for (size_t m = 0; m<INTVEN.BS_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.BS_treat_year_off[m] - 0.5*t_step) &&
            (t < INTVEN.BS_treat_year_off[m] + 0.51*t_step))
        {
            cout << "End of changing front-line BS treatment" << endl;

            theta.treat_cov = theta.BS_treat_cov_base;
            theta.treat_eff = theta.BS_treat_eff_base;
            theta.r_P = 1.0 / theta.BS_treat_BSproph_base;
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 8: first-line treatment (primaquine)

    for (size_t m = 0; m<INTVEN.PQ_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.PQ_treat_year_on[m] - 0.5*t_step) &&
            (t < INTVEN.PQ_treat_year_on[m] + 0.51*t_step))
        {
            cout << "New front-line PQ treatment" << endl;

            theta.PQ_treat_cover   = INTVEN.PQ_treat_cover[m];
            theta.PQ_treat_PQcover = INTVEN.PQ_treat_PQcover[m];
            theta.PQ_treat_BSeff   = INTVEN.PQ_treat_BSeff[m];
            theta.PQ_treat_PQeff   = INTVEN.PQ_treat_PQeff[m];
            theta.PQ_treat_BSproph = INTVEN.PQ_treat_BSproph[m];
            theta.PQ_treat_PQproph = INTVEN.PQ_treat_PQproph[m];
            theta.PQ_treat_CYP2D6  = INTVEN.PQ_treat_CYP2D6[m];

            theta.treat_cov = theta.PQ_treat_cover;
            theta.treat_eff = theta.PQ_treat_BSeff;
            theta.r_P = 1.0 / theta.PQ_treat_BSproph;
        }
    }


    //////////////////////////////////////////////////////////
    // Switching back to baseline.

    for (size_t m = 0; m<INTVEN.PQ_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.PQ_treat_year_off[m] - 0.5*t_step) &&
            (t < INTVEN.PQ_treat_year_off[m] + 0.51*t_step))
        {
            cout << "End of changing front-line PQ treatment" << endl;

            theta.PQ_treat_cover = 0.0;
            theta.PQ_treat_PQeff = 0.0;
            theta.PQ_treat_BSproph = 10.0;
            theta.PQ_treat_PQproph = 10.0;
            theta.PQ_treat_CYP2D6 = 1;

            theta.treat_cov = theta.BS_treat_cov_base;
            theta.treat_eff = theta.BS_treat_eff_base;
            theta.r_P = 1.0 / theta.BS_treat_BSproph_base;
        }
    }

}


//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 
//        //                                                                            // 
//  2.9.  //  Inverse of the cumulative normal distribution function required           //
//        //  for implementing correlated intervention coverage.                        //
////////////                                                                            //
////////////  The code is based on the following website                                //
////////////  http://www.johndcook.com/blog/cpp_phi_inverse/                            //
////////////  which is based the algorithm in Abramowitz and Stegun formula 26.2.23.    // 
////////////                                                                            //
////////////  The absolute value of the error should be less than 4.5 e-4 and it tests  //
////////////  out nicely in R.                                                          //
////////////                                                                            //
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double phi_inv(double pp, double mu, double sigma)
{
    if (pp <= 0.0 || pp >= 1.0)
    {
        throw("bad vlaue of pp (coverage) in erfinv");
    }

    double cc[] = { 2.515517, 0.802853, 0.010328 };
    double dd[] = { 1.432788, 0.189269, 0.001308 };

    double tt, temp;

    if (pp < 0.5)
    {
        tt = sqrt(-2.0*log(pp));

        temp = -(tt - ((cc[2] * tt + cc[1])*tt + cc[0]) / (((dd[2] * tt + dd[1])*tt + dd[0])*tt + 1.0));

    }
    else {

        tt = sqrt(-2.0*log(1.0 - pp));

        temp = (tt - ((cc[2] * tt + cc[1])*tt + cc[0]) / (((dd[2] * tt + dd[1])*tt + dd[0])*tt + 1.0));
    }

    return mu + sigma*temp;
}
