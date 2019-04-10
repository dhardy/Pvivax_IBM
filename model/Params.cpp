/////////////////////////////////////////////////////////////////////////////
///  Individual-based Plasmodium vivax transmission model.                ///
///  Model paramters                                                      ///
///  Dr Michael White, Imperial College London, m.white08@imperial.ac.uk  ///
/////////////////////////////////////////////////////////////////////////////


#include "Params.h"

void Params::set(const char *parameter_File)
{
    ////////////////////////////////////////////
    //                                        //
    // 1.3. Read in model parameters          //
    //                                        //
    ////////////////////////////////////////////
    
    cout << "Reading in parameter file............." << endl;
    cout << endl;

    string discard;

    std::ifstream parameter_Stream(parameter_File);

    if (parameter_Stream.fail())
    {
        std::cout << "Failure reading in data." << endl;
    }


    //////////////////////////////////////////////////
    // Population size and simulation time

    parameter_Stream >> discard >> PNG_pop.N_pop >> discard;               // Number of participants

    parameter_Stream >> discard >> Pv_mod_par.EIR_equil >> discard;        // EIR at equilibrium

    for (int g = 0; g < N_spec_max; g++)
    {
        parameter_Stream >> discard >> Pv_mod_par.Prop_mosq[g] >> discard; // Proportions of each mosquito species
    }

    parameter_Stream >> discard >> time_start >> discard;                  // Start time for simulation
    parameter_Stream >> discard >> time_end >> discard;                    // End time for simulation
    parameter_Stream >> discard >> burnin_time >> discard;                 // End time for simulation

    N_time = (1 / t_step)*(burnin_time + time_end - time_start) * 365;     // Number of time steps for simulation


    /////////////////////////////////
    // Human demography
    // Note that the mean age isn't technically the mean 
    // due to the effects of truncation at maximum age

    parameter_Stream >> discard >> Pv_mod_par.age_mean >> discard;    // mean age of human population
    parameter_Stream >> discard >> Pv_mod_par.age_max >> discard;     // maximum age in human population

    Pv_mod_par.mu_H = 1.0 / Pv_mod_par.age_mean;                      // human death rate

    Pv_mod_par.het_max = 100.0;                                       // Maximum relative heterogeneity value


    ////////////////////////////////////////////////////////
    // Heterogeneity in exposure and age-dependent biting

    parameter_Stream >> discard >> Pv_mod_par.age_0 >> discard;       // age-dependent biting parameter
    parameter_Stream >> discard >> Pv_mod_par.rho_age >> discard;     // age-dependent biting parameter
    parameter_Stream >> discard >> Pv_mod_par.sig_het >> discard;     // heterogeneity in exposure - standard deviation on a log scale


    /////////////////////////////////
    // Transmission probabilities

    parameter_Stream >> discard >> Pv_mod_par.bb >> discard;          // mosquito -> human transmission probability
    parameter_Stream >> discard >> Pv_mod_par.c_PCR >> discard;       // human -> mosquito transmission probability (PCR detectable)
    parameter_Stream >> discard >> Pv_mod_par.c_LM >> discard;        // human -> mosquito transmission probability (LM detectable)
    parameter_Stream >> discard >> Pv_mod_par.c_D >> discard;         // human -> mosquito transmission probability (clinical disease)
    parameter_Stream >> discard >> Pv_mod_par.c_T >> discard;         // human -> mosquito transmission probability (during treatment)


    ////////////////////////////////
    // Human recovery paramters

    parameter_Stream >> discard >> Pv_mod_par.d_latent >> discard;             // latent period in liver    
    parameter_Stream >> discard >> Pv_mod_par.r_LM >> discard;                 // rate of recovery from LM detectable infection
    parameter_Stream >> discard >> Pv_mod_par.r_D >> discard;                  // rate of recovery from symptomatic disease
    parameter_Stream >> discard >> Pv_mod_par.r_T >> discard;                  // rate of progression through treatment 
    parameter_Stream >> discard >> Pv_mod_par.d_PCR_min >> discard;            // minimum duration of PCR-detectable infection - full immunity
    parameter_Stream >> discard >> Pv_mod_par.d_PCR_max >> discard;            // maximum duration of PCR-detectable infection - no immunity

    parameter_Stream >> discard >> Pv_mod_par.BS_treat_cov_base >> discard;      // proportion of episodes of symptomatic disease treated (baseline)
    parameter_Stream >> discard >> Pv_mod_par.BS_treat_eff_base >> discard;      // efficacy of front-line treatment (baseline)
    parameter_Stream >> discard >> Pv_mod_par.BS_treat_BSproph_base >> discard;  // duration of prophylaxis of front-line treatment (baseline)

    parameter_Stream >> discard >> Pv_mod_par.A_PCR_50pc >> discard;           // PCR_detectable infection scale parameter
    parameter_Stream >> discard >> Pv_mod_par.K_PCR >> discard;                // PCR_detectable infection shape parameter


    Pv_mod_par.H_track = int(Pv_mod_par.d_latent / t_step);                                 // Number of time steps for duration of latency

    Pv_mod_par.treat_cov = Pv_mod_par.BS_treat_cov_base;                       // Treatment coverage
    Pv_mod_par.treat_eff = Pv_mod_par.BS_treat_eff_base;                       // Efficacy of treatment
    Pv_mod_par.r_P = 1.0 / Pv_mod_par.BS_treat_BSproph_base;                   // rate of recovery from prophylaxis

    Pv_mod_par.PQ_treat_cover   = 0.0;
    Pv_mod_par.PQ_treat_PQcover = 0.0;
    Pv_mod_par.PQ_treat_BSeff   = 0.0;
    Pv_mod_par.PQ_treat_PQeff   = 0.0;
    Pv_mod_par.PQ_treat_BSproph = 10.0;

    //////////////////
    // temporary setting treatment cov to zero
    //
    //Pv_mod_par.treat_cov_eff = 0.0;


    /////////////////////////////////
    // Blood-stage immunity paramters

    parameter_Stream >> discard >> Pv_mod_par.u_par >> discard;         // scale paramter for acquisition of blood-stage immunity
    parameter_Stream >> discard >> Pv_mod_par.r_par >> discard;         // rate of decay of blood-stage immunity

    parameter_Stream >> discard >> Pv_mod_par.phi_LM_max >> discard;    // probability of blood-stage infection with no immunity 
    parameter_Stream >> discard >> Pv_mod_par.phi_LM_min >> discard;    // probability of blood-stage infection with maximum immunity
    parameter_Stream >> discard >> Pv_mod_par.A_LM_50pc >> discard;     // blood-stage immunity scale parameter
    parameter_Stream >> discard >> Pv_mod_par.K_LM >> discard;          // blood-stage immunity shape parameter


    /////////////////////////////////
    // Clinical immunity paramters

    parameter_Stream >> discard >> Pv_mod_par.u_clin >> discard;       // scale paramter for acquisition of blood-stage immunity
    parameter_Stream >> discard >> Pv_mod_par.r_clin >> discard;       // rate of decay of clinical immunity

    parameter_Stream >> discard >> Pv_mod_par.phi_D_max >> discard;    // probability of clinical episode with no immunity
    parameter_Stream >> discard >> Pv_mod_par.phi_D_min >> discard;    // probability of clinical episode with maximum immunity
    parameter_Stream >> discard >> Pv_mod_par.A_D_50pc >> discard;     // clinical immunity scale parameter
    parameter_Stream >> discard >> Pv_mod_par.K_D >> discard;          // clinical immunity shape parameter


    /////////////////////////////////////////////
    // maternal immunity

    parameter_Stream >> discard >> Pv_mod_par.P_mat >> discard;        // New-born immunity relative to mother's
    parameter_Stream >> discard >> Pv_mod_par.d_mat >> discard;        // Inverse of decay rate of maternal immunity


    /////////////////////////////////
    // Relapse paramters

    parameter_Stream >> discard >> Pv_mod_par.ff >> discard;           // relapse rate
    parameter_Stream >> discard >> Pv_mod_par.gamma_L >> discard;      // liver clearance rate


    ////////////////////////////////
    // Human genotype prevalences

    parameter_Stream >> discard >> Pv_mod_par.G6PD_prev >> discard;    // prevalence of G6PD deficiency
    parameter_Stream >> discard >> Pv_mod_par.CYP2D6_prev >> discard;  // prevalence of CYP2D6 phenotype 


    Pv_mod_par.PQ_treat_CYP2D6 = 1;


    ////////////////////////////////////////////////////////
    // Intervention distribution parameters

    parameter_Stream >> discard >> Pv_mod_par.rho_round_LLIN >> discard;
    parameter_Stream >> discard >> Pv_mod_par.rho_round_IRS >> discard;
    parameter_Stream >> discard >> Pv_mod_par.rho_round_MDA >> discard;

    parameter_Stream >> discard >> Pv_mod_par.rho_LLIN_IRS >> discard;
    parameter_Stream >> discard >> Pv_mod_par.rho_MDA_IRS >> discard;
    parameter_Stream >> discard >> Pv_mod_par.rho_MDA_LLIN >> discard;


    Pv_mod_par.sig_round_LLIN = sqrt((1.0 - Pv_mod_par.rho_round_LLIN) / Pv_mod_par.rho_round_LLIN);
    Pv_mod_par.sig_round_IRS = sqrt((1.0 - Pv_mod_par.rho_round_IRS) / Pv_mod_par.rho_round_IRS);
    Pv_mod_par.sig_round_MDA = sqrt((1.0 - Pv_mod_par.rho_round_MDA) / Pv_mod_par.rho_round_MDA);


    Pv_mod_par.V_int[0][0] = 1.0;
    Pv_mod_par.V_int[0][1] = Pv_mod_par.rho_LLIN_IRS;
    Pv_mod_par.V_int[0][2] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[0][3] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[0][4] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[0][5] = Pv_mod_par.rho_MDA_LLIN;

    Pv_mod_par.V_int[1][0] = Pv_mod_par.rho_LLIN_IRS;
    Pv_mod_par.V_int[1][1] = 1.0;
    Pv_mod_par.V_int[1][2] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[1][3] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[1][4] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[1][5] = Pv_mod_par.rho_MDA_IRS;

    Pv_mod_par.V_int[2][0] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[2][1] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[2][2] = 1.0;
    Pv_mod_par.V_int[2][3] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[2][4] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[2][5] = Pv_mod_par.rho_round_MDA;

    Pv_mod_par.V_int[3][0] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[3][1] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[3][2] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[3][3] = 1.0;
    Pv_mod_par.V_int[3][4] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[3][5] = Pv_mod_par.rho_round_MDA;

    Pv_mod_par.V_int[4][0] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[4][1] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[4][2] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[4][3] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[4][4] = 1.0;
    Pv_mod_par.V_int[4][5] = Pv_mod_par.rho_round_MDA;

    Pv_mod_par.V_int[5][0] = Pv_mod_par.rho_MDA_LLIN;
    Pv_mod_par.V_int[5][1] = Pv_mod_par.rho_MDA_IRS;
    Pv_mod_par.V_int[5][2] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[5][3] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[5][4] = Pv_mod_par.rho_round_MDA;
    Pv_mod_par.V_int[5][5] = 1.0;



    /////////////////////////////////////////////////////////
    // We need to make a dummy covariance matrix for intervention
    // distribution because of the way genmn works

    for (int p = 0; p<N_int; p++)
    {
        for (int q = 0; q<N_int; q++)
        {
            Pv_mod_par.V_int_dummy[p][q] = Pv_mod_par.V_int[p][q];
        }
    }

    parameter_Stream.close();

    cout << "Parameter values read in from file!" << endl;
    cout << endl;


    ////////////////////////////////////////////
    //                                        //
    // 1.4. Read in mosquito parameters       //
    //                                        //
    ////////////////////////////////////////////

    cout << "Reading in mosquito files............." << endl;
    cout << endl;

    for (int g = 0; g < N_spec; g++)
    {

        std::ifstream mosquito_Stream(mosquito_File[g]);

        if (mosquito_Stream.fail())
        {
            std::cout << "Failure reading in mosquito parameters." << endl;
        }


        /////////////////////////////////
        // Death rate and duration of sporogony

        mosquito_Stream >> discard >> Pv_mod_par.mu_M[g] >> discard;      // mosquito death rate
        mosquito_Stream >> discard >> Pv_mod_par.tau_M[g] >> discard;     // duration of sporogony

        Pv_mod_par.M_track = int(Pv_mod_par.tau_M[g] * mosq_steps / t_step);


        /////////////////////////////////
        // Larval paramters    
        // From White et al (2011) P&V

        mosquito_Stream >> discard >> Pv_mod_par.d_E_larvae >> discard;      // Development time of early larval instars
        mosquito_Stream >> discard >> Pv_mod_par.d_L_larvae >> discard;      // Development time of late larval instars
        mosquito_Stream >> discard >> Pv_mod_par.d_pupae >> discard;         // Development time of pupae
        mosquito_Stream >> discard >> Pv_mod_par.mu_E0 >> discard;           // Mortality rate of early larval instars (low density)
        mosquito_Stream >> discard >> Pv_mod_par.mu_L0 >> discard;           // Mortality rate of late larval instars (low density)
        mosquito_Stream >> discard >> Pv_mod_par.mu_P >> discard;            // Mortality rate of pupae
        mosquito_Stream >> discard >> Pv_mod_par.beta_larvae >> discard;     // Number of eggs laid per day per mosquito
        mosquito_Stream >> discard >> Pv_mod_par.gamma_larvae >> discard;    // Effect of density dependence on late instars relative to early instars

        Pv_mod_par.omega_larvae[g] = Pv_mod_par.gamma_larvae*Pv_mod_par.mu_L0 / Pv_mod_par.mu_E0 - Pv_mod_par.d_E_larvae / Pv_mod_par.d_L_larvae + (Pv_mod_par.gamma_larvae - 1.0)*Pv_mod_par.mu_L0*Pv_mod_par.d_E_larvae;
        Pv_mod_par.omega_larvae[g] = -0.5*Pv_mod_par.omega_larvae[g] + sqrt(0.25*Pv_mod_par.omega_larvae[g] * Pv_mod_par.omega_larvae[g] + 0.5*Pv_mod_par.gamma_larvae*Pv_mod_par.beta_larvae*Pv_mod_par.mu_L0*Pv_mod_par.d_E_larvae /
                                       (Pv_mod_par.mu_E0*Pv_mod_par.mu_M[g] * Pv_mod_par.d_L_larvae*(1.0 + Pv_mod_par.d_pupae*Pv_mod_par.mu_P)));


        /////////////////////////////////
        // Seasonality paramters
        // Denominator for seasonality - see Griffin (2015) PLoS Comp Biol

        mosquito_Stream >> discard >> Pv_mod_par.dry_seas[g] >> discard;      // Proportion of dry season transmission compared to mean 
        mosquito_Stream >> discard >> Pv_mod_par.kappa_seas[g] >> discard;    // Shape parameter for seasonality
        mosquito_Stream >> discard >> Pv_mod_par.t_peak_seas[g] >> discard;   // Timing of peak for seasonal transmission

        Pv_mod_par.denom_seas[g] = exp(gammln(0.5) + gammln(Pv_mod_par.kappa_seas[g] + 0.5) - gammln(Pv_mod_par.kappa_seas[g] + 1.0)) / 3.14159265359;


        //////////////////////////////////////////////
        // Entomology paramters

        mosquito_Stream >> discard >> Pv_mod_par.Q_0[g] >> discard;           // Human Blood Index (proportion of blood meals taken on humans)
        mosquito_Stream >> discard >> Pv_mod_par.CHI_endo[g] >> discard;      // Endophily - proportion of mosquitoes resting indoors after feeding (no intervention)
        mosquito_Stream >> discard >> Pv_mod_par.PSI_indoors[g] >> discard;   // Proportion of bites taken on humans indoors
        mosquito_Stream >> discard >> Pv_mod_par.PSI_bed[g] >> discard;       // Proportion of bites taken on humans in bed

        mosquito_Stream >> discard >> Pv_mod_par.delta_1 >> discard;          // Time spent foraging for a blood meal
        mosquito_Stream >> discard >> Pv_mod_par.delta >> discard;            // Duration of gonotrophic cycle

        Pv_mod_par.delta_2 = Pv_mod_par.delta - Pv_mod_par.delta_1;

        Pv_mod_par.p_1[g] = exp(-Pv_mod_par.mu_M[g] * Pv_mod_par.delta_1);
        Pv_mod_par.p_2[g] = exp(-Pv_mod_par.mu_M[g] * Pv_mod_par.delta_2);

        Pv_mod_par.aa[g] = Pv_mod_par.Q_0[g] / (Pv_mod_par.delta_1 + Pv_mod_par.delta_2);

        Pv_mod_par.eps_max[g] = Pv_mod_par.beta_larvae*(exp(Pv_mod_par.delta*Pv_mod_par.mu_M[g]) - 1.0) / Pv_mod_par.mu_M[g];


        //////////////////////////////////////////////
        // LLIN paramters

        mosquito_Stream >> discard >> Pv_mod_par.LLIN_half_life >> discard;

        mosquito_Stream >> discard >> Pv_mod_par.PYR_half_life >> discard;

        mosquito_Stream >> discard >> Pv_mod_par.r_LLIN_0[g] >> discard;       // Probability mosquito repelled (with full insecticide activity)
        mosquito_Stream >> discard >> Pv_mod_par.r_LLIN_net[g] >> discard;     // Probability mosquito repelled due to barrier effect of net (no insecticide)
        mosquito_Stream >> discard >> Pv_mod_par.d_LLIN_0[g] >> discard;       // Probability mosquito dies during feeding attempt


        Pv_mod_par.P_LLIN_loss = 1.0 - exp(-t_step*log(2.0) / Pv_mod_par.LLIN_half_life);   // Probability of losing LLIN in a time step
        Pv_mod_par.PYR_decay = log(2.0) / Pv_mod_par.PYR_half_life;                            // Rate of pyrethroid decay

        Pv_mod_par.s_LLIN_0[g] = 1.0 - Pv_mod_par.r_LLIN_0[g] - Pv_mod_par.d_LLIN_0[g];     // Probability mosquito feeds successfully


        ////////////////////////////////////////////////////////
        // IRS parameters

        mosquito_Stream >> discard >> Pv_mod_par.IRS_half_life >> discard;    // IRS insecticide half-life

        mosquito_Stream >> discard >> Pv_mod_par.r_IRS_0[g] >> discard;          // IRS repellency
        mosquito_Stream >> discard >> Pv_mod_par.d_IRS_0[g] >> discard;          // IRS death

        Pv_mod_par.IRS_decay = log(2.0) / Pv_mod_par.IRS_half_life;         // IRS decay rate

        Pv_mod_par.s_IRS_0[g] = 1.0 - Pv_mod_par.d_IRS_0[g] - Pv_mod_par.s_IRS_0[g];   // Feeding success of mosquito on IRS protected person
    }

    cout << "Mosquito parameter values read in from file!" << endl;
    cout << endl;


    //////////////////////////////////////////////
    // Normalise relative proprotions of 
    // different mosquito species

    double Prop_mosq_denom = 0.0;

    for (int g = 0; g < N_spec; g++)
    {
        Prop_mosq_denom = Prop_mosq_denom + Pv_mod_par.Prop_mosq[g];
    }

    for (int g = 0; g < N_spec; g++)
    {
        Pv_mod_par.Prop_mosq[g] = Pv_mod_par.Prop_mosq[g] / Prop_mosq_denom;
    }


    ///////////////////////////////////////////////////////////
    //                                                       //
    // 1.5. Pre-multiplication of quantities for efficiency  //
    //                                                       //
    ///////////////////////////////////////////////////////////

    Pv_mod_par.A_par_decay = exp(-Pv_mod_par.r_par*t_step);
    Pv_mod_par.A_clin_decay = exp(-Pv_mod_par.r_clin*t_step);
    Pv_mod_par.mat_decay = exp(-Pv_mod_par.d_mat*t_step);

    Pv_mod_par.age_0_inv = 1.0 / Pv_mod_par.age_0;                 // Inverse of age-dependent biting parameter

    Pv_mod_par.A_PCR_50pc_inv = log2 / Pv_mod_par.A_PCR_50pc;      // Immune scalar for clearance of infection
    Pv_mod_par.A_LM_50pc_inv = 1.0 / Pv_mod_par.A_LM_50pc;        // Immune scalar for BS infection
    Pv_mod_par.A_D_50pc_inv = 1.0 / Pv_mod_par.A_D_50pc;         // Immune scalar for clinical disease

    Pv_mod_par.P_dead = 1.0 - exp(-t_step*Pv_mod_par.mu_H);
    Pv_mod_par.P_preg = 0.0014189;


    Pv_mod_par.P_PYR_decay = exp(-Pv_mod_par.PYR_decay*t_step);
    Pv_mod_par.P_IRS_decay = exp(-Pv_mod_par.IRS_decay*t_step);


    ///////////////////////////////////////////////////////////
    //                                                       //
    // 1.6. Fill out hypnozoite transition matrices          //
    //                                                       //
    ///////////////////////////////////////////////////////////

    for (int k1 = 0; k1<(K_max + 1); k1++)
    {
        for (int k2 = 0; k2 < (K_max + 1); k2++)
        {
            Pv_mod_par.D_MAT[k1][k2] = 0.0;
            Pv_mod_par.OD_MAT[k1][k2] = 0.0;
            Pv_mod_par.K_MAT[k1][k2] = 0.0;
            Pv_mod_par.L_MAT[k1][k2] = 0.0;
            Pv_mod_par.H_MAT[k1][k2] = 0.0;
        }
    }

    for (int k = 0; k < (K_max + 1); k++)
    {
        Pv_mod_par.D_MAT[k][k] = 1.0;
    }

    for (int k = 0; k < K_max; k++)
    {
        Pv_mod_par.OD_MAT[k + 1][k] = 1.0;
    }
    Pv_mod_par.OD_MAT[K_max][K_max] = 1.0;

    for (int k = 0; k < (K_max + 1); k++)
    {
        Pv_mod_par.K_MAT[k][k] = (double)(k);
    }

    for (int k = 0; k < K_max; k++)
    {
        Pv_mod_par.L_MAT[k][k + 1] = +(double)(k + 1);
        Pv_mod_par.L_MAT[k + 1][k + 1] = -(double)(k + 1);
    }

    for (int k = 0; k < K_max; k++)
    {
        Pv_mod_par.H_MAT[k][k] = -1.0;
        Pv_mod_par.H_MAT[k + 1][k] = +1.0;
    }
}
