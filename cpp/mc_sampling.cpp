#include <cmath>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "mc_sampling.h"
#include "E_pot.h"
#include "wave.h"
#include "Random.h"

using namespace std;
using namespace arma;

#define num_core 1
// Buggy if we set num_core!=1 so we leave it
#define h 0.001
#define h2 1000000

void quantum_force(mat&, mat&, double, double, double, int, int, double, int, mat&);

void mc_sampling (  double stepa, 
                    double stepb, 
                    double omega, 
                    int dimension, 
                    int number_particles, 
                    int charge, 
                    int max_variations, 
                    double starta, 
                    double startb,
                    int number_cycles, 
                    double timestep,
                    vec& r_mean, mat& r_mean2,
                    vec& e_kv, 
                    vec& e_potv, 
                    mat& e_km, 
                    mat& e_potm,
                    vec& cumulative_e, 
                    vec& cumulative_e2, 
                    mat& m_cumulative, 
                    mat& m_cumulative2, 
                    int& J, 
                    vector<Random*> &randoms) {

    double D = 0.5; // Diffusion constant

    int el;
    cout <<"Electron-electron repulsion? (1 for yes, 0 for no): " <<endl;
    cin >> el;

    J = 0;
    cout <<"Jastrow factor? (1 for yes, 0 for no): " <<endl;
    cin >> J;


    double a_sym=0, a_asym=0;
    mat a;
    a = zeros<mat>(number_particles, number_particles);
    if (dimension == 2) {
        a_sym = 1./3.;
        a_asym = 1.0;
    } else if (dimension == 3) {
        a_sym = 1. / 4.;
        a_asym = 1. / 2.;
    } else {
        cout << "Unable to initialize Jastrow paremters: Unknown dimension" << endl;
        exit(1);
    }
    for (unsigned int i = 0; i < number_particles; i++) {
        for (unsigned int j = 0; j < number_particles; j++) {
            /*
            if ((i&1) == (j&1)) {
                a(i,j) = a_sym;
            } else {
                a(i,j) = a_asym;
            }
            */
            if (((i<number_particles/2) && (j<number_particles/2)) || ((i>number_particles/2) && (j>number_particles/2))) {
                a(i,j) = a_sym;
            } else {
                a(i,j) = a_asym;
            }
        }
    }

    int acceptedMoves = 0;
    int thermalization = number_cycles/10;
    if (J == 0) {
        // loop over variational parameters
        #pragma omp parallel for num_threads(num_core)
        for (int variate=1; variate <= max_variations; variate++) {
            mat r_old, r_new, qforce_old, qforce_new;
            r_old = zeros<mat>(number_particles, dimension);
            r_new = zeros<mat>(number_particles, dimension);
            qforce_old = zeros<mat>(number_particles, dimension);
            qforce_new = zeros<mat>(number_particles, dimension);

            double e_k=0, e_pot=0;

            // initialisations of variational parameters and energies
            double alpha   = starta + variate*stepa;
            double energy  = 0; 
            double delta_e = 0; 
            double energy2 = 0;
            double beta    = 0;
            //  initial trial position, note calling with alpha
            //  and in three dimensions
            for (int i = 0; i < number_particles; i++) {
                for (int j=0; j < dimension; j++) {
                    r_old(i,j) = randoms[omp_get_thread_num()]->nextGauss() * sqrt(timestep);
                }
            }

            double wfold = wave_function(omega, r_old, alpha, beta, dimension, number_particles, J, a);
            quantum_force (r_old, qforce_old, alpha, beta, wfold, number_particles, dimension, omega, J, a);

            for (int cycles = 1; cycles <= number_cycles; cycles++){
                for (int i = 0; i < number_particles; i++) {
                    for (int j=0; j < dimension; j++) {
                        r_new(i,j) = r_old(i,j) + randoms[omp_get_thread_num()]->nextGauss()*sqrt(timestep) + qforce_old(i,j)*timestep*D;
                    }

                    for(int k = 0 ; k < number_particles; k++) {
                        if( k != i) {
                            for(int j = 0; j < dimension; j++) {
                                r_new(k,j) = r_old(k,j);
                            }
                        }
                    }

                    double wfnew = wave_function(omega, r_new, alpha, beta, dimension, number_particles, J, a);
                    quantum_force(r_new, qforce_new, alpha, beta, wfnew, number_particles, dimension, omega, J, a);

                    // we compute the log of the ratio of the greens functions to be used in the Metropolis_Hastings algorithm
                    double greensfunction = 0.0;
                    for(int j = 0; j < dimension; j++) {
                        greensfunction += 0.5*(qforce_old(i,j) + qforce_new(i,j))*(D*timestep*0.5*(qforce_old(i,j) - qforce_new(i,j)) - r_new(i,j) + r_old(i,j));
                    }
                    greensfunction = exp(greensfunction);

                    double ratio = wfnew / wfold;
                    ratio *= ratio;
                    if(randoms[omp_get_thread_num()]->nextDouble() <= greensfunction*ratio) {
                        acceptedMoves++;
                        for(int j=0; j < dimension; j++) {
                            r_old(i,j) = r_new(i,j);
                            qforce_old(i,j) = qforce_new(i,j);
                        }
                        wfold = wfnew;
                    }

                    // compute local energy
                    delta_e = local_energy(e_k, e_pot, omega, r_old, el, alpha, beta, wfold, dimension, number_particles, charge, J, a);
                    e_kv(variate)   += e_k;
                    e_potv(variate) += e_pot;
                    // update energies
                    energy  += delta_e;
                    energy2 += delta_e*delta_e;

                    if(number_particles == 2) {
                        double r_m = 0;
                        for (int i = 0; i < number_particles-1; i++) {
                            for (int j = i+1; j < number_particles; j++) {
                                for (int k = 0; k < dimension; k++) {
                                    r_m += (r_old(i,k)-r_old(j,k))*(r_old(i,k)-r_old(j,k));
                                }
                            }
                        }
                        r_mean(variate) += sqrt(r_m);
                    }
                }
            }   // end of loop over MC trials

            int number_cycles_particles = number_cycles * number_particles;
            r_mean(variate) = r_mean(variate)/double(number_cycles_particles);

            cout << "variational parameter= " << alpha <<endl;
            // update the energy average and its squared
            e_kv(variate) = e_kv(variate)/number_cycles_particles;
            e_potv(variate) = e_potv(variate)/number_cycles_particles;
            cumulative_e(variate) = energy/number_cycles_particles;
            cumulative_e2(variate) = energy2/number_cycles_particles;
        }
    } else { // with Jastrow factor
        #pragma omp parallel for num_threads(num_core)
        for (int variate=1; variate <= max_variations; variate++){
            cout << "   *** Completion: " << (float)(variate-1) / (float)(max_variations-1) << endl;
            mat r_old;
            mat r_new;
            mat qforce_old;
            mat qforce_new;
            r_old        = zeros<mat>(number_particles, dimension);
            r_new        = zeros<mat>(number_particles, dimension);
            qforce_old   = zeros<mat>(number_particles, dimension);
            qforce_new   = zeros<mat>(number_particles, dimension);
            double alpha = starta + variate*stepa;
            
            #pragma omp parallel for num_threads(num_core)
            for(int variate2=1; variate2 <= max_variations; variate2 ++) {
                double energy= 0; double delta_e=0; double energy2=0;
                double e_k=0, e_pot=0;
                delta_e=0;
                energy = energy2 = 0;
                double beta = startb +variate2*stepb;

                for (int i = 0; i < number_particles; i++) {
                    for (int j=0; j < dimension; j++) {
                        r_old(i,j) = randoms[omp_get_thread_num()]->nextGauss();// * sqrt(timestep);
                    }
                }

                double wfold = wave_function(omega, r_old, alpha, beta, dimension, number_particles, J, a);
                quantum_force (r_old, qforce_old, alpha, beta, wfold, number_particles, dimension, omega, J, a);

                // loop over monte carlo cycles
                acceptedMoves = 0;
                int i = 0;
                for (int cycles = 1; cycles <= number_cycles; cycles++){
                    //for (int i = 0; i < number_particles; i++) {
                        for (int j=0; j < dimension; j++) {
                            r_new(i,j) = r_old(i,j) + randoms[omp_get_thread_num()]->nextGauss()*sqrt(timestep) + qforce_old(i,j)*timestep*D;
                        }
                        for(int k = 0 ; k < number_particles; k++) {
                            if( k != i) {
                                for(int j = 0; j < dimension; j++) {
                                    r_new(k,j) = r_old(k,j);
                                }
                            }
                        }

                       double wfnew = wave_function(omega, r_new, alpha, beta, dimension, number_particles, J, a);
                       quantum_force(r_new, qforce_new, alpha, beta, wfnew, number_particles, dimension, omega, J, a);

                        // we compute the log of the ratio of the greens functions to be used in the Metropolis_Hastings algorithm
                        double greensfunction = 0.0;
                        for(int j = 0; j < dimension; j++) {
                            greensfunction += 0.5*(qforce_old(i,j) + qforce_new(i,j))*(D*timestep*0.5*(qforce_old(i,j) - qforce_new(i,j)) - r_new(i,j) + r_old(i,j));
                        }
                        greensfunction = exp(greensfunction);

                        // The Metropolis test is performed by moving one particle at the time
                        double ratio = wfnew / wfold;
                        ratio *= ratio;
                        if(randoms[omp_get_thread_num()]->nextDouble() <= greensfunction * ratio) {
                            acceptedMoves++;
                            for(int j=0; j < dimension; j++) {
                                r_old(i,j) = r_new(i,j);
                                qforce_old(i,j) = qforce_new(i,j);
                            }
                            wfold = wfnew;
                        }

                        if (cycles > thermalization) {
                            // compute local energy
                            delta_e = local_energy(e_k, e_pot, omega, r_old, el, alpha, beta, wfold, dimension, number_particles, charge, J, a);
                            e_km(variate, variate2)   += e_k;
                            e_potm(variate, variate2) += e_pot;
                            // update energies
                            energy  += delta_e;
                            energy2 += delta_e*delta_e;

                            if(number_particles == 2) {
                                double r_m = 0;
                                for (int i = 0; i < number_particles-1; i++) {
                                    for (int j = i+1; j < number_particles; j++) {
                                        for (int k = 0; k < dimension; k++) {
                                            r_m += (r_old(i,k)-r_old(j,k))*(r_old(i,k)-r_old(j,k));
                                        }
                                    }
                                }
                                r_mean2(variate, variate2) += sqrt(r_m);
                            }
                        }
                    //}
                    i++;
                    if (i==number_particles) {
                        i = 0;
                    }
                }   // end of loop over MC trials

                int number_cycles_particles = number_cycles * number_particles;
                int measurements = (number_cycles-thermalization);
                r_mean2(variate, variate2) = r_mean2(variate, variate2)/double(number_cycles);

                cout << "alpha = " << alpha << "   ::   beta = " << beta << "   ::   energy = " << energy/(number_cycles-thermalization) << "   ::   accept = " << (float)acceptedMoves/(float)(number_cycles) << endl;
                // update the energy average and its squared
                e_km(variate, variate2)          = e_km(variate, variate2)   / (number_cycles-thermalization);
                e_potm(variate, variate2)        = e_potm(variate, variate2) / (number_cycles-thermalization);
                m_cumulative(variate, variate2)  = energy  / (number_cycles-thermalization);
                m_cumulative2(variate, variate2) = energy2 / (number_cycles-thermalization);
            }
        }
    }
}


/*void distance(mat& r, mat& r_12, int number_particles, int dimension) {
    for (int i = 0; i < number_particles; i++) {
        for (int j = i+1; j < number_particles; j++) {
            for (int k = 0; k < dimension; k++) {
                double temp = r(i,k) - r(j,k);
                r_12(i,j) += temp * temp;
            }
            r_12(i,j) = sqrt(r_12(i,j));
            r_12(j,i) = r_12(i,j);
        }
    }
}*/


/*void distance(mat& r, mat& r_12, int number_particles, int dimension) {
    const double* ptr = r.memptr();
    for (int i = 0; i < number_particles; i++) {
        for (int j = i+1; j < number_particles; j++) {
            double sum = 0.0;
            for (int k = 0; k < dimension; k++) {
                double temp = r(i,k) - r(j,k);
                sum += temp * temp;
                if (ptr[i + k*number_particles] != r(i,k)) {
                    cout << ptr[i + k*number_particles] << endl;
                    cout << r(i,k) << endl;

                    for (int x=0; x<number_particles ; x++) {
                        for (int y=0; y<dimension; y++) {
                            cout << x << " " << y << " :: " << r(x,y) << endl;
                        }

                    }

                    cout << "whoops              " << i << " " << j << " " << k << endl;
                    exit(1);
                }
            }
            sum = sqrt(sum);
            r_12(i,j) = sum;
            r_12(j,i) = sum;
        }
    }
}*/



void distance(mat& r, mat& r_12, int number_particles, int dimension) {
    const double* rPtr = r.memptr();
    double* r_12Ptr = r_12.memptr();
    r_12(0,0) = r(0,0); // I have no god damn idea why, but I have to reference an element of r SOMEWHERE in this function or my answers are all cattywompus...


    if (number_particles == 6) {
        #define NP 6
        #pragma unroll
        for (int i = 0, ii=0; i < NP; i++, ii+=NP) {
            #pragma unroll
            for (int j = i+1, jj=ii+NP; j < NP; j++, jj+=NP) {
                double dx = rPtr[i] - rPtr[j];
                double dy = rPtr[i+NP] - rPtr[j+NP];

                double temp = sqrt(dx*dx + dy*dy);
                r_12Ptr[ii+j] = temp;
                r_12Ptr[i+jj] = temp;
            }
        }
        #undef NP
    } else {
        for (int i = 0, ii=0; i < number_particles; i++, ii+=number_particles) {
            for (int j = i+1, jj=ii+number_particles; j < number_particles; j++, jj+=number_particles) {
                double dx = rPtr[i] - rPtr[j];
                double dy = rPtr[i+number_particles] - rPtr[j+number_particles];

                double temp = sqrt(dx*dx + dy*dy);
                r_12Ptr[ii+j] = temp;
                r_12Ptr[i+jj] = temp;
            }
        }
    }
}

double  local_energy(   double& e_k, 
                        double& e_pot, 
                        double omega, 
                        mat& r, 
                        int eeRepulsion, 
                        double alpha, 
                        double beta, 
                        double wfold, 
                        int dimension,
                        int number_particles, 
                        int charge, 
                        int J,
                        mat &a) {
    int i, j;
    double e_local=0, wfminus=0, wfplus=0, e_kinetic=0, e_potential=0, r_single_particle=0;
    mat r_plus;
    mat r_minus;
    mat r_12(number_particles, number_particles);
    r_plus = r_minus = r;
    r_12.zeros();

    distance(r, r_12, number_particles, dimension);

    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_plus(i,j)     = r(i,j)+h;
            r_minus(i,j)    = r(i,j)-h;
            wfminus         = wave_function(omega, r_minus, alpha, beta, dimension, number_particles, J, a);
            wfplus          = wave_function(omega, r_plus, alpha, beta, dimension, number_particles, J, a);
            e_kinetic      -= (wfminus + wfplus - 2.0*wfold);
            r_plus(i,j)     = r(i,j);
            r_minus(i,j)    = r(i,j);
        }
    }
    e_kinetic = 0.5*h2*e_kinetic/wfold;
    e_k = e_kinetic;

    E_pot potenziale2(omega, r, dimension, number_particles, charge, r_12, r_single_particle);
    e_potential = potenziale2.oscillator();

    // contribution from electron-electron potential
    if (eeRepulsion == 1) {
        E_pot potenziale3(omega, r, dimension, number_particles, charge, r_12, r_single_particle);
        e_potential += potenziale3.electronelectron();
    } else {
        e_potential += 0;
    }
    e_pot = e_potential;

    r_plus.reset();
    r_minus.reset();

    e_local = e_potential+e_kinetic;
    return e_local;
}

void initialise(double& deltaAlpha, 
                double& deltaBeta, 
                double& omega, 
                int& dimension, 
                double& initialAlpha, 
                double& initialBeta, 
                int& number_particles, 
                int& charge,
                int& max_variations, 
                int& number_cycles, 
                double& timestep) {
    cout << "dimensions           = ";
    cin >> dimension;
    cout << "number of electrons  = ";
    cin >> number_particles;
    cout << "number of protons    = ";
    cin >> charge;
    cout << "start alpha          = ";
    cin >> initialAlpha;
    cout << "start beta           = ";
    cin >> initialBeta;
    cout << "step alpha           = ";
    cin >> deltaAlpha;
    cout << "step beta            = ";
    cin >> deltaBeta;
    cout << "parameter space size = ";
    cin >> max_variations;
    cout << "number of MC trials  = ";
    cin >> number_cycles;
    cout << "time step            = ";
    cin >> timestep;
    cout << "omega (frequency)    = ";
    cin >> omega;
}  

void quantum_force( mat &r, 
                    mat &qforce, 
                    double alpha, 
                    double beta, 
                    double wf, 
                    int number_particles, 
                    int dimension, 
                    double omega, 
                    int J,
                    mat &a) {
    int i, j;
    double wfminus=0;
    double wfplus=0;
    mat r_plus(number_particles, dimension);
    mat r_minus(number_particles, dimension);
    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_plus(i,j) = r_minus(i,j) = r(i,j);
        }
    }

    double wfhInverse = 1.0 / (wf * h);
    for(i = 0; i < number_particles; i++) {
        for(j = 0; j < dimension; j++) {
            r_plus(i,j) = r(i,j) + h;
            r_minus(i,j) = r(i,j) - h;
            wfminus = wave_function(omega, r_minus, alpha, beta, dimension, number_particles, J, a);
            wfplus = wave_function(omega, r_plus, alpha, beta, dimension, number_particles, J, a);
            qforce(i,j) =  (wfplus - wfminus) * wfhInverse;
            r_plus(i,j) = r(i,j);
            r_minus(i,j) = r(i,j);
        }
    }
}