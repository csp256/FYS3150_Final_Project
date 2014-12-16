#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "mc_sampling.h"
#include "E_pot.h"

using namespace std;
using namespace arma;

ofstream ofile;

void  output(vec&, mat&, double, double, double, double, int, int, int, vec, vec, mat, mat, vec, vec, mat, mat, int, double);

int main(int argc, char* argv[]) {
    char *outfilename;
    int number_cycles=0, max_variations=0, charge=0, dimension=0, number_particles=0, J=0;
    double timestep=0, omega=0, stepa=0, stepb=0, starta=0, startb=0;
    vector<Random*> randoms;
    randoms.push_back(new Random(-1));
    randoms.push_back(new Random(-1));
    randoms.push_back(new Random(-1));
    randoms.push_back(new Random(-1));

    if(argc != 2) {
        cout << "Call this program with:   ./monte_carlo output_file seed" << endl;
        exit(1);
    } else {
        outfilename = argv[1];
    }
    ofile.open(outfilename);

    initialise(stepa, stepb, omega, dimension, starta, startb, number_particles, charge, max_variations, number_cycles, timestep);

    vec cumulative_e    (max_variations+1);
    vec cumulative_e2   (max_variations+1);
    mat m_cumulative    (max_variations+1,  max_variations +1);
    mat m_cumulative2   (max_variations+1,  max_variations +1);
    vec e_kv            (max_variations+1);
    vec e_potv          (max_variations+1);
    mat e_km            (max_variations+1,  max_variations+1);
    mat e_potm          (max_variations+1,  max_variations+1);
    mat r_mean2         (max_variations+1,  max_variations+1);
    vec r_mean          (max_variations+1);
    cumulative_e.zeros();
    cumulative_e2.zeros();
    m_cumulative.zeros();
    m_cumulative2.zeros();
    e_kv.zeros();
    e_potv.zeros();
    e_km.zeros();
    e_potm.zeros();
    r_mean.zeros();
    r_mean2.zeros();


  //  Do the mc sampling
    mc_sampling(stepa, stepb, omega, dimension, number_particles, charge,
                max_variations, starta, startb,
                number_cycles, timestep,
                r_mean, r_mean2,
                e_kv, e_potv, e_km, e_potm, cumulative_e, cumulative_e2, m_cumulative, m_cumulative2, J, randoms);

    // Print out results
    output(r_mean, r_mean2, starta, startb, stepa, stepb, max_variations, number_cycles, charge, e_kv, e_potv, e_km, e_potm, cumulative_e, cumulative_e2, m_cumulative, m_cumulative2, J, number_particles);
    cumulative_e.reset();
    cumulative_e2.reset();
    m_cumulative.reset();
    m_cumulative2.reset();
    e_kv.reset();
    e_potv.reset();
    e_km.reset();
    e_potm.reset();
    r_mean.reset();
    ofile.close();  // close output file
    return 0;
}



void output(vec& r_mean, mat& r_mean2, double starta, double startb, double stepa, double stepb, int max_variations, int number_cycles, int charge,
            vec e_kv, vec e_potv, mat e_km, mat e_potm,
            vec cumulative_e, vec cumulative_e2, mat m_cumulative, mat m_cumulative2, int J, double number_particles) {
    int i, j;
    double alpha=0, beta=0, variance=0, error=0;
    alpha = starta;
    for( i=1; i <= max_variations; i++){
        alpha += stepa;
        if(J == 1) {
            beta = startb;
            for( j=1; j<=max_variations; j++) {
                beta  += stepb;
                variance = m_cumulative2(i,j)-m_cumulative(i,j)*m_cumulative(i,j);
                error=sqrt(variance/number_cycles);
                ofile << setiosflags(ios::showpoint | ios::uppercase);
                ofile << setw(15) << setprecision(8) << alpha;
                ofile << setw(15) << setprecision(8) << beta;
                ofile << setw(15) << setprecision(8) << m_cumulative(i,j);
                ofile << setw(15) << setprecision(8) << e_km(i,j);
                ofile << setw(15) << setprecision(8) << e_potm(i,j);
                ofile << setw(15) << setprecision(8) << variance;
                ofile << setw(15) << setprecision(8) << error << endl;
                //if (number_particles == 2) {
                //    ofile << setw(15) << setprecision(8) << r_mean2(i,j) << endl;
                //}
            }
        } else {
            variance = cumulative_e2[i]-cumulative_e[i]*cumulative_e[i];
            error=sqrt(variance/number_cycles);
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setw(15) << setprecision(8) << alpha;
            ofile << setw(15) << setprecision(8) << 1.0; 
            ofile << setw(15) << setprecision(8) << cumulative_e[i];
            ofile << setw(15) << setprecision(8) << e_kv(i);
            ofile << setw(15) << setprecision(8) << e_potv(i);
            ofile << setw(15) << setprecision(8) << variance;
            ofile << setw(15) << setprecision(8) << error << endl;
            //if( number_particles == 2) {
            //    ofile << setw(15) << setprecision(8) << r_mean(i) << endl;
            //}
        }
    }
}