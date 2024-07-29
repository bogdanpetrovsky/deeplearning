//
//  elemGeneration.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 03.06.2024.
//

#include "elemgen.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define rnode drand48()*(cmax-cmin)+cmin

using namespace std;

int main() {
    int i,j,k,nel;
    double node0[8][3], node[8][3], cmin, cmax, rseed;
    /*---Part 1--*/
    scanf ("%le %le %d %lf", &cmin, &cmax, &nel, &rseed) ;
    srand48(rseed) ;
    node[0][0] = 0.0; node [0][1] = 0.0; node [0][2] = 0.0;
    node[1][0] = 1.0; node[1][1] = 0.0; node[1][2] = 0.0;
    node0[2][0] = 1.0; node0[2][1] = 1.0; node0[2][2] = 0.0;
    node0[3][0] = 0.0; node0[3][1] = 1.0; node0[3][2] = 0.0;
    node0[4][0]= 0.0; node0[4][1] = 0; node0[4][2] = 1;
    node0[5][0] = 1.0; node0[5][1] = 0.0; node0[5][2] = 1.0;
    node0[6][0] = 1.0; node0[6][1] = 1.0; node0[6][2] = 1.0;
    node0[7][0] = 0.0; node0[7][1] = 1.0; node0[7][2] = 1.0;
    /*---Part 2ーーーーーーーー*/
    printf ("%d\n", nel);
    for(i=0;i<nel;i++) {
        for (j=2;j<8;j++) {
            for (k=0;k<3;k++) node[j][k] = node0[j][k] + rnode;
        }
    }
    node [3] [2] = 0.0;
    for (j=0;j<8;j++) {
        printf("%d %d", i, j);
        for (k=0;k<3;k++) printf (" %e" ,node[j][k]);
            printf ("\n");
    }
    return 0;
}
