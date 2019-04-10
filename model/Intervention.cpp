/////////////////////////////////////////////////////////////////////////////
///  Individual-based Plasmodium vivax transmission model.                ///
///  Intervention data                                            ///
///  Dr Michael White, Imperial College London, m.white08@imperial.ac.uk  ///
/////////////////////////////////////////////////////////////////////////////

#include "Intervention.h"
#include <iostream>
#include <fstream>


void Intervention::set(const char *coverage_File)
{
    /////////////////////////////////////////////////////////////////////////
    // 1.7.1. Read in matrix of years and coverages
    //        Note that the matrix we read in may have variable size

    int N_cov_rounds = 0;

    cout << "Read in intervention coverage file............" << endl;

    std::ifstream coverage_Stream(coverage_File);

    if (coverage_Stream.fail())
    {
        std::cout << "Failure reading in data." << endl;
    }


    vector<vector<double>> coverage;
    coverage.resize(0);

    vector<double> cov_read(24);

    do {
        for (int j = 0; j<24; j++)
        {
            coverage_Stream >> cov_read[j];
        }
        if (cov_read[0] > -0.5)
        {
            N_cov_rounds = N_cov_rounds + 1;
            coverage.push_back(cov_read);
        }
    } while (cov_read[0] > -0.5);

    cout << "Intervention coverage file successfully read in!" << endl;


    /////////////////////////////////////////////////////////////////////////
    // 1.7.2. Fill out intervention structure

    for (int i = 0; i<N_cov_rounds; i++)
    {
        //////////////////////////////////////////////////////////////
        // LLINs

        if ((coverage[i][0] > -0.5) && (coverage[i][1] > -0.5))
        {
            LLIN_year.push_back(coverage[i][0] * 365.0);
            LLIN_cover.push_back(coverage[i][1]);
        }


        //////////////////////////////////////////////////////////////
        // IRS

        if ((coverage[i][0] > -0.5) && (coverage[i][2] > -0.5))
        {
            IRS_year.push_back(coverage[i][0] * 365.0);
            IRS_cover.push_back(coverage[i][2]);
        }


        //////////////////////////////////////////////////////////////
        // MDA - blood-stage drugs

        if ((coverage[i][0] > -0.5) && (coverage[i][3] > -0.5))
        {
            MDA_BS_year.push_back(coverage[i][0] * 365.0);
            MDA_BS_cover.push_back(coverage[i][3]);
            MDA_BS_BSeff.push_back(coverage[i][4]);
            MDA_BS_BSproph.push_back(coverage[i][5]);
        }


        //////////////////////////////////////////////////////////////
        // MDA - blood-stage drugs plus primaquine

        if ((coverage[i][0] > -0.5) && (coverage[i][6] > -0.5))
        {
            MDA_PQ_year.push_back(coverage[i][0] * 365.0);
            MDA_PQ_cover.push_back(coverage[i][6]);
            MDA_PQ_BSeff.push_back(coverage[i][7]);
            MDA_PQ_PQeff.push_back(coverage[i][8]);
            MDA_PQ_BSproph.push_back(coverage[i][9]);
            MDA_PQ_PQproph.push_back(coverage[i][10]);
            MDA_PQ_CYP2D6.push_back((int) (coverage[i][11]));
        }


        //////////////////////////////////////////////////////////////
        // Front-line treatment - blood-stage drugs

        if ((coverage[i][0] > -0.5) && (coverage[i][12] > -0.5))
        {
            BS_treat_year_on.push_back(coverage[i][0] * 365.0);
            BS_treat_cover.push_back(coverage[i][12]);
            BS_treat_BSeff.push_back(coverage[i][13]);
            BS_treat_BSproph.push_back(coverage[i][14]);
            BS_treat_year_off.push_back(coverage[i][15] * 365.0);
        }


        //////////////////////////////////////////////////////////////
        // Front-line treatment - primaquine

        if ((coverage[i][0] > -0.5) && (coverage[i][16] > -0.5))
        {
            PQ_treat_year_on.push_back(coverage[i][0] * 365.0);
            PQ_treat_cover.push_back(coverage[i][16]);
            PQ_treat_PQcover.push_back(coverage[i][17]);
            PQ_treat_BSeff.push_back(coverage[i][18]);
            PQ_treat_PQeff.push_back(coverage[i][19]);
            PQ_treat_BSproph.push_back(coverage[i][20]);
            PQ_treat_PQproph.push_back(coverage[i][21]);
            PQ_treat_CYP2D6.push_back((int) (coverage[i][22]));
            PQ_treat_year_off.push_back(coverage[i][23] * 365.0);
        }
    }
}
