#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
//#include "lib.h"
#include "Hermite1.h"

using namespace std;
using namespace arma;

double Hermite(double coor, int q_num) {
    if (q_num == 0) {
        return 1;
    } else {
        return coor;// + coor; // we should be able to get away with just "coor" because only ratios of determinants matter
    }
}