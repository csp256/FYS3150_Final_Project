#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "mc_sampling.h"
#include "E_pot.h"
#include "wave.h"

using namespace std;
using namespace arma;

double Jas(mat &r_12, double beta, int dimension, int number_particles, int J, mat &a) {
    if(J == 0) {
        return 1.0;
    }

    double argument;

    const double* aPtr = a.memptr();
    const double* rPtr = r_12.memptr();
    for (int i = 0, ii=0; i < number_particles; i++, ii+=number_particles) {
        for (int j = i + 1; j < number_particles; j++) {
            double temp = rPtr[ii+j];
            argument += aPtr[ii+j] * temp / (1.0 + beta * temp);
        }
    }

    return exp(argument);
}
