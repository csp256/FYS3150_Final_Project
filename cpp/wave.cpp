#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "mc_sampling.h"
#include "Jas.h"
#include "Slater1.h"

using namespace std;
using namespace arma;

// Function to compute the squared wave function, simplest form

double  wave_function(double omega, mat &r, double alpha, double beta, int dimension, int number_particles, int J, mat &a) {
    mat r_12;
    r_12 = zeros <mat> (number_particles, number_particles);
    double jas;
    double slater;

    distance(r, r_12, number_particles, dimension);
    jas 	= Jas(r_12, beta, dimension, number_particles, J, a);
    slater 	= Slater(r, dimension, number_particles, omega, alpha);
    return jas * slater;
}
