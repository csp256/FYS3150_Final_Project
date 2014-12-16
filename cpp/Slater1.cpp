#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "Slater1.h"
#include "Hermite1.h"

using namespace std;
using namespace arma;

/*
double Slater(  mat& r, 
                int dimension, 
                int number_particles, 
                double omega, 
                double alpha) {
    int half_number_particles = number_particles / 2;
    int i, ii, j, jj, k;
    vec q_num(dimension);
    q_num.zeros();

    mat psi1(half_number_particles, half_number_particles);
    mat psi2(half_number_particles, half_number_particles);
    psi1.zeros();
    psi2.zeros();

    double argument1 = 0.0;
    double argument2 = 0.0;
    const double halfNegativeAlphaOmega = -0.5*alpha*omega; // Insert bad joke here.

    for(i=0, ii=half_number_particles; i<number_particles/2; i++, ii++) {
        switch (i) {
            //case 0: 
            //    q_num(0) = 0;
            //    q_num(1) = 0;
            //    break;
            case 1:
                q_num(0) = 1;
            //    q_num(1) = 0;
                break;
            case 2:
                q_num(0) = 0;
                q_num(1) = 1;
                break;
            //default:
            //    cout << "Oops! Why is your quantum number " << i << "?" << endl;
            //    exit(1);
        }

        for(j=0, jj=half_number_particles; j<half_number_particles; j++, jj++) {
            double product1 = 1.0;
            double product2 = 1.0;
            for(k=0; k<dimension; k++) {
                double r1 = r(j,k);
                double r2 = r(jj,k);
                argument1 += r1*r1;
                argument2 += r2*r2;
                product1 *= Hermite(r1, q_num(k));
                product2 *= Hermite(r2, q_num(k));
            }
            psi1(i,j) = product1 * exp(halfNegativeAlphaOmega*argument1);
            psi2(i,j) = product2 * exp(halfNegativeAlphaOmega*argument2);
        }
        argument1 = 0.0;
        argument2 = 0.0;
    }

    double det1, det2;
    if (number_particles == 2) {
        Det result(r, omega, alpha);
        det1 = result.det2();
        return det1;
    } else if(number_particles == 6) {
        Det result1(psi1, omega, alpha);
        Det result2(psi2, omega, alpha);
        det1 = result1.det6();
        det2 = result2.det6();
        return det1*det2;
    }
}*/



// double Slater(  mat& r, 
//                 int dimension, 
//                 int number_particles, 
//                 double omega, 
//                 double alpha) {
//     double det1, det2;
//     if (number_particles == 2) {
//         Det result(r, omega, alpha);
//         det1 = result.det2();
//         return det1;
//     } else if (number_particles == 6 && dimension == 2) {
//         //#define NP 6
//         //#define DIM 2

// //        int half_number_particles = NP / 2;
//        // int i, ii, j, jj, k;
//     //    int* quantumNumbers[dimension];

//   //      mat psi1(half_number_particles, half_number_particles);
//       //  mat psi2(half_number_particles, half_number_particles);
// //        psi1.zeros();
// //        psi2.zeros();
// //        double* psiPtr1 = psi1.memptr();
//   //      double* psiPtr2 = psi2.memptr();

//         const double* rPtr = r.memptr();

//   //      double argument1;
// //        double argument2;
//         double product;// = 1.0;
// //        double product2 = 1.0;

//         const double halfNegativeAlphaOmega = -0.5*alpha*omega; // Insert bad joke here.



//         // wait if we assume that same spin particles are adjacent in memory here, is that convention followed in the a(i,j) matrix for the jastrow factor?

//         // double arg = 0.0;
//         // #pragma unroll
//         // for (i=0; i<DIM*NP; i++) {
//         //     arg += rPtr[i]*rPtr[i];
//         // }

//         const double sumOfMagnitudesSquared = rPtr[0]*rPtr[0] + rPtr[1]*rPtr[1] + rPtr[2]*rPtr[2] + rPtr[3]*rPtr[3] + rPtr[4]*rPtr[4] + rPtr[5]*rPtr[5] + rPtr[6]*rPtr[6] + rPtr[7]*rPtr[7] + rPtr[8]*rPtr[8] + rPtr[9]*rPtr[9] + rPtr[10]*rPtr[10] + rPtr[11]*rPtr[11];

// //         #pragma unroll
// //         for(i=0, ii=half_number_particles; i<half_number_particles; i++, ii++) {
// // //            argument1 = 0.0;
// // //            argument2 = 0.0;
// //   //          for(k=0; k<DIM; k++) {
// // //                double r1 = r(i,k);
// //                 double r1 = rPtr[i];
// //                 double r2 = rPtr[ii];//r(ii,k);
// //                 argument1 += r1*r1;
// //                 argument2 += r2*r2;


// //                 r1 = rPtr[i  + NP];
// //                 r2 = rPtr[ii + NP];//r(ii,k);
// //                 argument1 += r1*r1;
// //                 argument2 += r2*r2;

// //     //        }
// //   //          argument1 = exp(halfNegativeAlphaOmega * argument1);
// // //            argument2 = exp(halfNegativeAlphaOmega * argument2);
// // //            product1 *= exp(halfNegativeAlphaOmega * (argument1+argument2));
// //             //product2 *= exp(halfNegativeAlphaOmega * argument2);

// // //            product1 *= argument1;
// // //            product2 *= argument2;

// //             //for(j=0, jj=half_number_particles; j<half_number_particles; j++, jj++) { 
// //               //  switch(j) {
// //                 //    case 0:
// //                       //  #define J 0
// //                         //psi1(i,j) = argument1; // should be able to remove use of these function calls alltogether.
// //                       //  psi2(i,j) = argument2;
// //                   //      psiPtr1[i] = 1.0;//argument1;
// //                      //   psiPtr2[i] = 1.0;//argument2;
// //                         //#undef J 
// //                     //break;
// //                   //  case 1:
// //                 //        psiPtr1[i+half_number_particles] = /*argument1 * */ rPtr[i+number_particles]; //r(i,1);
// //                    //     psiPtr2[i+half_number_particles] = /*argument2 * */ rPtr[ii+number_particles]; // r(ii,1);
// // //                        psi1(i,j) = argument1 * r(i,1);
// // //                        psi2(i,j) = argument2 * r(ii,1);
// //                     //break;
// //                     //case 2:
// //               //          psiPtr1[i+2*half_number_particles] = /*argument1 * */ rPtr[i];//r(i,0);
// //                  //       psiPtr2[i+2*half_number_particles] = /*argument2 * */rPtr[ii];//r(ii,0);
// // //                        psi1(i,j) = argument1 * r(i,0);
// // //                        psi2(i,j) = argument2 * r(ii,0);
// //                     //break;
// //                 //}
// //             //}
// //         }
//         // #undef NP
//         // #undef DIM

//        // product1 *= exp(halfNegativeAlphaOmega * (argument1+argument2));
//         product = exp(halfNegativeAlphaOmega * sumOfMagnitudesSquared);

//     //    Det result1(psi1, omega, alpha);
// //        Det result2(psi2, omega, alpha);
//         product *= rPtr[7]*(rPtr[2]-rPtr[0]) + rPtr[6]*(rPtr[1]-rPtr[2]) + rPtr[8]*(rPtr[0]-rPtr[1]);
//         product *= rPtr[7+3]*(rPtr[2+3]-rPtr[0+3]) + rPtr[6+3]*(rPtr[1+3]-rPtr[2+3]) + rPtr[8+3]*(rPtr[0+3]-rPtr[1+3]);
// //        det1 = result1.det6() * product1;
// //        det2 = result2.det6() * product2;
//         return product;///*det1*/product1*product2;//det2;
//     } else {
//         cout << "We do not have support for that many particles/dimensions at this time. (but it should be pretty easy to fix)" << endl;
//         exit(1);
//     }
// }








double Slater(  mat& r, 
                int dimension, 
                int number_particles, 
                double omega, 
                double alpha) {
    double det1, det2;
    if (number_particles == 2) {
        Det result(r, omega, alpha);
        det1 = result.det2();
        return det1;
    } else if (number_particles == 6 && dimension == 2) {
        const double* rPtr = r.memptr();
        const double sumOfMagnitudesSquared = rPtr[0]*rPtr[0] + rPtr[1]*rPtr[1] + rPtr[2]*rPtr[2] + rPtr[3]*rPtr[3] + rPtr[4]*rPtr[4] + rPtr[5]*rPtr[5] + rPtr[6]*rPtr[6] + rPtr[7]*rPtr[7] + rPtr[8]*rPtr[8] + rPtr[9]*rPtr[9] + rPtr[10]*rPtr[10] + rPtr[11]*rPtr[11];
        double product = exp(-0.5*alpha*omega * sumOfMagnitudesSquared);
        product *= rPtr[7]*(rPtr[2]-rPtr[0]) + rPtr[6]*(rPtr[1]-rPtr[2]) + rPtr[8]*(rPtr[0]-rPtr[1]);
        product *= rPtr[7+3]*(rPtr[2+3]-rPtr[0+3]) + rPtr[6+3]*(rPtr[1+3]-rPtr[2+3]) + rPtr[8+3]*(rPtr[0+3]-rPtr[1+3]);
        return product;
    } else {
        cout << "We do not have support for that many particles/dimensions at this time. (but it should be pretty easy to fix)" << endl;
        exit(1);
    }
}




double Det::det6() {
//    return det(m_r);
    const double* a = m_r.memptr();
    return a[0]*(a[4]*a[8]-a[5]*a[7])-a[3]*(a[1]*a[8]-a[2]*a[7])+(a[1]*a[5]-a[2]*a[4])*a[6];

    // Below is the 6x6 case, O(n!), incase you need it.
    // return a[0]*(-a[8]*(-a[15]*(-a[22]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[34]*a[25]-a[28]*a[31])+(a[28]*a[35]-a[29]*a[34])*a[19])+a[16]*(-a[21]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[35]-a[29]*a[33])*a[19])-a[17]*(-a[21]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[34]-a[28]*a[33])*a[19])+a[13]*((a[27]*a[34]-a[28]*a[33])*a[23]-(a[27]*a[35]-a[29]*a[33])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[21]))+a[9]*(-a[14]*(-a[22]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[34]*a[25]-a[28]*a[31])+(a[28]*a[35]-a[29]*a[34])*a[19])+a[16]*(-a[20]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[35]-a[29]*a[32])*a[19])-a[17]*(-a[20]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[34]-a[28]*a[32])*a[19])+a[13]*((a[26]*a[34]-a[28]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[20]))-a[10]*(-a[14]*(-a[21]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[35]-a[29]*a[33])*a[19])+a[15]*(-a[20]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[35]-a[29]*a[32])*a[19])-a[17]*(-a[20]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[33]-a[27]*a[32])*a[19])+a[13]*((a[26]*a[33]-a[27]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[21]+(a[27]*a[35]-a[29]*a[33])*a[20]))+a[11]*(-a[14]*(-a[21]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[34]-a[28]*a[33])*a[19])+a[15]*(-a[20]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[34]-a[28]*a[32])*a[19])-a[16]*(-a[20]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[33]-a[27]*a[32])*a[19])+a[13]*((a[26]*a[33]-a[27]*a[32])*a[22]-(a[26]*a[34]-a[28]*a[32])*a[21]+(a[27]*a[34]-a[28]*a[33])*a[20]))+a[7]*(a[14]*((a[27]*a[34]-a[28]*a[33])*a[23]-(a[27]*a[35]-a[29]*a[33])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[21])-a[15]*((a[26]*a[34]-a[28]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[20])+a[16]*((a[26]*a[33]-a[27]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[21]+(a[27]*a[35]-a[29]*a[33])*a[20])-a[17]*((a[26]*a[33]-a[27]*a[32])*a[22]-(a[26]*a[34]-a[28]*a[32])*a[21]+(a[27]*a[34]-a[28]*a[33])*a[20])))+a[2]*(a[6]*(-a[15]*(-a[22]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[34]*a[25]-a[28]*a[31])+(a[28]*a[35]-a[29]*a[34])*a[19])+a[16]*(-a[21]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[35]-a[29]*a[33])*a[19])-a[17]*(-a[21]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[34]-a[28]*a[33])*a[19])+a[13]*((a[27]*a[34]-a[28]*a[33])*a[23]-(a[27]*a[35]-a[29]*a[33])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[21]))+a[9]*(a[12]*(-a[22]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[34]*a[25]-a[28]*a[31])+(a[28]*a[35]-a[29]*a[34])*a[19])+a[16]*(a[18]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[35]*a[24]-a[29]*a[30]))-a[17]*(a[18]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[34]*a[24]-a[28]*a[30]))-a[13]*(-a[22]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[34]*a[24]-a[28]*a[30])+(a[28]*a[35]-a[29]*a[34])*a[18]))-a[10]*(a[12]*(-a[21]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[35]-a[29]*a[33])*a[19])+a[15]*(a[18]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[35]*a[24]-a[29]*a[30]))-a[17]*(a[18]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[33]*a[24]-a[27]*a[30]))-a[13]*(-a[21]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[35]-a[29]*a[33])*a[18]))+a[11]*(a[12]*(-a[21]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[34]-a[28]*a[33])*a[19])+a[15]*(a[18]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[34]*a[24]-a[28]*a[30]))-a[16]*(a[18]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[33]*a[24]-a[27]*a[30]))-a[13]*(-a[21]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[34]-a[28]*a[33])*a[18]))-a[7]*(-a[15]*(-a[22]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[34]*a[24]-a[28]*a[30])+(a[28]*a[35]-a[29]*a[34])*a[18])+a[16]*(-a[21]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[35]-a[29]*a[33])*a[18])-a[17]*(-a[21]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[34]-a[28]*a[33])*a[18])+a[12]*((a[27]*a[34]-a[28]*a[33])*a[23]-(a[27]*a[35]-a[29]*a[33])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[21])))-a[3]*(a[6]*(-a[14]*(-a[22]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[34]*a[25]-a[28]*a[31])+(a[28]*a[35]-a[29]*a[34])*a[19])+a[16]*(-a[20]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[35]-a[29]*a[32])*a[19])-a[17]*(-a[20]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[34]-a[28]*a[32])*a[19])+a[13]*((a[26]*a[34]-a[28]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[20]))+a[8]*(a[12]*(-a[22]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[34]*a[25]-a[28]*a[31])+(a[28]*a[35]-a[29]*a[34])*a[19])+a[16]*(a[18]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[35]*a[24]-a[29]*a[30]))-a[17]*(a[18]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[34]*a[24]-a[28]*a[30]))-a[13]*(-a[22]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[34]*a[24]-a[28]*a[30])+(a[28]*a[35]-a[29]*a[34])*a[18]))-a[10]*(a[12]*(-a[20]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[35]-a[29]*a[32])*a[19])+a[14]*(a[18]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[35]*a[24]-a[29]*a[30]))-a[17]*(a[18]*(a[32]*a[25]-a[26]*a[31])+a[20]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[32]*a[24]-a[26]*a[30]))-a[13]*(-a[20]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[35]-a[29]*a[32])*a[18]))+a[11]*(a[12]*(-a[20]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[34]-a[28]*a[32])*a[19])+a[14]*(a[18]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[34]*a[24]-a[28]*a[30]))-a[16]*(a[18]*(a[32]*a[25]-a[26]*a[31])+a[20]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[32]*a[24]-a[26]*a[30]))-a[13]*(-a[20]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[34]-a[28]*a[32])*a[18]))-a[7]*(-a[14]*(-a[22]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[34]*a[24]-a[28]*a[30])+(a[28]*a[35]-a[29]*a[34])*a[18])+a[16]*(-a[20]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[35]-a[29]*a[32])*a[18])-a[17]*(-a[20]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[34]-a[28]*a[32])*a[18])+a[12]*((a[26]*a[34]-a[28]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[20])))+a[4]*(a[6]*(-a[14]*(-a[21]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[35]-a[29]*a[33])*a[19])+a[15]*(-a[20]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[35]-a[29]*a[32])*a[19])-a[17]*(-a[20]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[33]-a[27]*a[32])*a[19])+a[13]*((a[26]*a[33]-a[27]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[21]+(a[27]*a[35]-a[29]*a[33])*a[20]))+a[8]*(a[12]*(-a[21]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[35]-a[29]*a[33])*a[19])+a[15]*(a[18]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[35]*a[24]-a[29]*a[30]))-a[17]*(a[18]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[33]*a[24]-a[27]*a[30]))-a[13]*(-a[21]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[35]-a[29]*a[33])*a[18]))-a[9]*(a[12]*(-a[20]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[35]-a[29]*a[32])*a[19])+a[14]*(a[18]*(a[35]*a[25]-a[29]*a[31])+a[23]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[35]*a[24]-a[29]*a[30]))-a[17]*(a[18]*(a[32]*a[25]-a[26]*a[31])+a[20]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[32]*a[24]-a[26]*a[30]))-a[13]*(-a[20]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[35]-a[29]*a[32])*a[18]))+a[11]*(a[12]*(-a[20]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[33]-a[27]*a[32])*a[19])+a[14]*(a[18]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[33]*a[24]-a[27]*a[30]))-a[15]*(a[18]*(a[32]*a[25]-a[26]*a[31])+a[20]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[32]*a[24]-a[26]*a[30]))-a[13]*(-a[20]*(a[33]*a[24]-a[27]*a[30])+a[21]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[33]-a[27]*a[32])*a[18]))-a[7]*(-a[14]*(-a[21]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[35]-a[29]*a[33])*a[18])+a[15]*(-a[20]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[35]-a[29]*a[32])*a[18])-a[17]*(-a[20]*(a[33]*a[24]-a[27]*a[30])+a[21]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[33]-a[27]*a[32])*a[18])+a[12]*((a[26]*a[33]-a[27]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[21]+(a[27]*a[35]-a[29]*a[33])*a[20])))-a[5]*(a[6]*(-a[14]*(-a[21]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[34]-a[28]*a[33])*a[19])+a[15]*(-a[20]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[34]-a[28]*a[32])*a[19])-a[16]*(-a[20]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[33]-a[27]*a[32])*a[19])+a[13]*((a[26]*a[33]-a[27]*a[32])*a[22]-(a[26]*a[34]-a[28]*a[32])*a[21]+(a[27]*a[34]-a[28]*a[33])*a[20]))+a[8]*(a[12]*(-a[21]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[33]*a[25]-a[27]*a[31])+(a[27]*a[34]-a[28]*a[33])*a[19])+a[15]*(a[18]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[34]*a[24]-a[28]*a[30]))-a[16]*(a[18]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[33]*a[24]-a[27]*a[30]))-a[13]*(-a[21]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[34]-a[28]*a[33])*a[18]))-a[9]*(a[12]*(-a[20]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[34]-a[28]*a[32])*a[19])+a[14]*(a[18]*(a[34]*a[25]-a[28]*a[31])+a[22]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[34]*a[24]-a[28]*a[30]))-a[16]*(a[18]*(a[32]*a[25]-a[26]*a[31])+a[20]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[32]*a[24]-a[26]*a[30]))-a[13]*(-a[20]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[34]-a[28]*a[32])*a[18]))+a[10]*(a[12]*(-a[20]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[32]*a[25]-a[26]*a[31])+(a[26]*a[33]-a[27]*a[32])*a[19])+a[14]*(a[18]*(a[33]*a[25]-a[27]*a[31])+a[21]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[33]*a[24]-a[27]*a[30]))-a[15]*(a[18]*(a[32]*a[25]-a[26]*a[31])+a[20]*(a[31]*a[24]-a[30]*a[25])-a[19]*(a[32]*a[24]-a[26]*a[30]))-a[13]*(-a[20]*(a[33]*a[24]-a[27]*a[30])+a[21]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[33]-a[27]*a[32])*a[18]))-a[7]*(-a[14]*(-a[21]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[34]-a[28]*a[33])*a[18])+a[15]*(-a[20]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[34]-a[28]*a[32])*a[18])-a[16]*(-a[20]*(a[33]*a[24]-a[27]*a[30])+a[21]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[33]-a[27]*a[32])*a[18])+a[12]*((a[26]*a[33]-a[27]*a[32])*a[22]-(a[26]*a[34]-a[28]*a[32])*a[21]+(a[27]*a[34]-a[28]*a[33])*a[20])))-a[1]*(-a[8]*(-a[15]*(-a[22]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[34]*a[24]-a[28]*a[30])+(a[28]*a[35]-a[29]*a[34])*a[18])+a[16]*(-a[21]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[35]-a[29]*a[33])*a[18])-a[17]*(-a[21]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[34]-a[28]*a[33])*a[18])+a[12]*((a[27]*a[34]-a[28]*a[33])*a[23]-(a[27]*a[35]-a[29]*a[33])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[21]))+a[9]*(-a[14]*(-a[22]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[34]*a[24]-a[28]*a[30])+(a[28]*a[35]-a[29]*a[34])*a[18])+a[16]*(-a[20]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[35]-a[29]*a[32])*a[18])-a[17]*(-a[20]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[34]-a[28]*a[32])*a[18])+a[12]*((a[26]*a[34]-a[28]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[20]))-a[10]*(-a[14]*(-a[21]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[35]-a[29]*a[33])*a[18])+a[15]*(-a[20]*(a[35]*a[24]-a[29]*a[30])+a[23]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[35]-a[29]*a[32])*a[18])-a[17]*(-a[20]*(a[33]*a[24]-a[27]*a[30])+a[21]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[33]-a[27]*a[32])*a[18])+a[12]*((a[26]*a[33]-a[27]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[21]+(a[27]*a[35]-a[29]*a[33])*a[20]))+a[11]*(-a[14]*(-a[21]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[33]*a[24]-a[27]*a[30])+(a[27]*a[34]-a[28]*a[33])*a[18])+a[15]*(-a[20]*(a[34]*a[24]-a[28]*a[30])+a[22]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[34]-a[28]*a[32])*a[18])-a[16]*(-a[20]*(a[33]*a[24]-a[27]*a[30])+a[21]*(a[32]*a[24]-a[26]*a[30])+(a[26]*a[33]-a[27]*a[32])*a[18])+a[12]*((a[26]*a[33]-a[27]*a[32])*a[22]-(a[26]*a[34]-a[28]*a[32])*a[21]+(a[27]*a[34]-a[28]*a[33])*a[20]))+a[6]*(a[14]*((a[27]*a[34]-a[28]*a[33])*a[23]-(a[27]*a[35]-a[29]*a[33])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[21])-a[15]*((a[26]*a[34]-a[28]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[22]+(a[28]*a[35]-a[29]*a[34])*a[20])+a[16]*((a[26]*a[33]-a[27]*a[32])*a[23]-(a[26]*a[35]-a[29]*a[32])*a[21]+(a[27]*a[35]-a[29]*a[33])*a[20])-a[17]*((a[26]*a[33]-a[27]*a[32])*a[22]-(a[26]*a[34]-a[28]*a[32])*a[21]+(a[27]*a[34]-a[28]*a[33])*a[20])));
}

double Det::det2() {
    double r1, r2;
    r1 = m_r(0,0)*m_r(0,0) + m_r(0,1)*m_r(0,1);
    r2 = m_r(1,0)*m_r(1,0) + m_r(1,1)*m_r(1,1);
    return exp(-m_alpha*m_omega*(r1 + r2)*0.5);
}

Det::Det(   mat& r, 
            double omega, 
            double alpha) {
    m_r = r;
    m_alpha = alpha;
    m_omega = omega;
}