/////////////////////////////////////////////////////////////////////////////
///  Individual-based Plasmodium vivax transmission model.                ///
///  Human individual                                                     ///
///  Dr Michael White, Imperial College London, m.white08@imperial.ac.uk  ///
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_HUMAN
#define PVIVAX_MODEL_HUMAN

#include "Params.h"


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Human class
//
/////////////////////////////////////////////////////////////////////////////////////////

struct Human
{
    //////////////////////////////////////////////////////////////////////////
    //  Functions
    //////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////
    // 0.2.1. Class constructor

    Human(double a, double zeta)
    {
        age = a;
        zeta_het = zeta;
    }
    
    // Delete unwanted copy constructors
    Human(Human&) =delete;
    void operator= (Human&) =delete;
    // Allow default move constructors
    Human(Human&&) = default;
    Human& operator= (Human&&) = default;


    ////////////////////////////////////////////////////
    // 0.2.2. Function declarations within the human class

    void state_mover(Params theta, double lam_bite);
    void ager(Params theta);
    void intervention_updater(Params theta);


    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // 0.2.3. Person-specific age and exposure classifiers

    double age;                      // Person's age
    double zeta_het;                 // Heterogeneity in exposure to mosquitoes

    bool gender;                     // 0 = male; 1 = female
    bool G6PD_def;                   // Is the person G6PD deficient? 0 = normal; 1 = deficient (homozygous); 2 = deficient (heterozygous - female only)  
    bool CYP2D6;                     // Does the person have CYP2D6 phenotype? 0 = normal; 1 = low metabolizer


    ////////////////////////////////////////////////////
    //  0.2.4. Child-bearing age.
    //         Indicator for people of 18-22 years of age. New-born
    //         children acquire a proportion of the immunity of a
    //         women aged 20 years.

    bool preg_age;
    bool pregnant;
    double preg_timer;


    double lam_bite_lag;             // Lagged force of infection due to moquito bites
    vector<double> lam_bite_track;   // Tracking lagged force of infection due to moquito bites

    double lam_rel_lag;              // Lagged force of infection due to relapses 
    vector<double> lam_rel_track;    // Tracking lagged force of infection due to relapses

    double lam_H_lag;                // Lagged total force of infection


    ////////////////////////////////////////////////////
    // 0.2.5. Indicators for compartments

    bool S;
    bool I_PCR;
    bool I_LM;
    bool I_D;
    bool T;
    bool P;

    /////////////////////////////////////////////////////////////////
    //  0.2.7. Number of batches of hypnozoites. Must be an integer. 

    int Hyp;


    ////////////////////////////////////////////////////
    // Indicator for competing hazards move 

    int CH_move;


    ////////////////////////////////////////////////////
    //  0.2.6. Indicators for new events

    bool I_PCR_new;    // New PCR-detectable infection (I_LM, I_D & T included here)
    bool I_LM_new;     // New LM-detectable infection (I_D & T included here)
    bool I_D_new;      // New clinical episode (treated or untreated)  
    bool ACT_new;      // 
    bool PQ_new;       // 


    ////////////////////////////////////////////////////
    //  0.2.8. Person-specific levels of immunity and indicators 
    //         for whether immunity is suppressed

    double A_par;
    double A_clin;

    double A_par_mat;
    double A_clin_mat;

    bool A_par_boost;
    bool A_clin_boost;

    double A_par_timer;
    double A_clin_timer;

    bool PQ_proph;
    double PQ_proph_timer;


    //////////////////////////////////////////////////////////
    // 0.2.9. Person-specific intervention access parameter

    double zz_int[N_int];


    ///////////////////////////////
    // LLINs

    bool LLIN;                 // Does the person have an LLIN?
    double LLIN_age;           // How old is the LLIN

    double r_LLIN[N_spec];     // repellency
    double d_LLIN[N_spec];     // killing effect
    double s_LLIN[N_spec];     // survival


    ///////////////////////////////
    // IRS

    bool IRS;                 // Is the person protected by IRS
    double IRS_age;           // How long ago was their house sprayed?

    double r_IRS[N_spec];     // repellency
    double d_IRS[N_spec];     // killing effect
    double s_IRS[N_spec];     // survival


    ///////////////////////////////////////////////////
    // Individual-level effect of vector control

    double z_VC[N_spec];      // probability of mosquito being repelled from this individual during a single feeding attempt
    double y_VC[N_spec];      // probability of mosquito feeding on this individual during a single attempt
    double w_VC[N_spec];      // probability of mosquito feeding and surviving on this individual during a single feeding attempt
};

#endif
