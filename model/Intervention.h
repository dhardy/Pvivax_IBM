/////////////////////////////////////////////////////////////////////////////
///  Individual-based Plasmodium vivax transmission model.                ///
///  Intervention data                                            ///
///  Dr Michael White, Imperial College London, m.white08@imperial.ac.uk  ///
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_INTERV
#define PVIVAX_MODEL_INTERV

#include <vector>

using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////
//
// Intervention details
//
/////////////////////////////////////////////////////////////////////////////////////////

struct Intervention
{
    //////////////////////////////////////////////////////////////////////////
    //  Functions
    //////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////
    // Set parameters from coverage file
    void set(const char *coverage_File);


    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////
    // LLINs

    vector<double> LLIN_year;
    vector<double> LLIN_cover;


    ////////////////////////////////
    // IRS

    vector<double> IRS_year;
    vector<double> IRS_cover;


    ////////////////////////////////
    // MDA - blood-stage drugs

    vector<double> MDA_BS_year;
    vector<double> MDA_BS_cover;
    vector<double> MDA_BS_BSeff;
    vector<double> MDA_BS_BSproph;


    ////////////////////////////////
    // MDA - blood-stage drugs + primaquine

    vector<double> MDA_PQ_year;
    vector<double> MDA_PQ_cover;
    vector<double> MDA_PQ_BSeff;
    vector<double> MDA_PQ_PQeff;
    vector<double> MDA_PQ_BSproph;
    vector<double> MDA_PQ_PQproph;
    vector<int>    MDA_PQ_CYP2D6;


    ////////////////////////////////
    // First-line treatment - blood-stage drugs

    vector<double> BS_treat_year_on;
    vector<double> BS_treat_year_off;
    vector<double> BS_treat_cover;
    vector<double> BS_treat_BSeff;
    vector<double> BS_treat_BSproph;


    ////////////////////////////////
    // First-line treatment - blood-stage drugs
    // plus primaquine

    vector<double> PQ_treat_year_on;
    vector<double> PQ_treat_year_off;
    vector<double> PQ_treat_cover;
    vector<double> PQ_treat_PQcover;
    vector<double> PQ_treat_BSeff;
    vector<double> PQ_treat_PQeff;
    vector<double> PQ_treat_BSproph;
    vector<double> PQ_treat_PQproph;
    vector<int>    PQ_treat_CYP2D6;
};

#endif
