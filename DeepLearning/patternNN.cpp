//
//  patternNN.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 15.09.2024.
//

#include "patternNN.hpp"
/* patternNN.c */
#include <stdio.h>
#include <stdlib.h>
int main() {
    int i,j,k, ia,ib,nel,n_in,n_out;
    float f1, f2, f3;
    FILE *fp1,*fp2;
    scanf ("%d %d", &n_in,&n_out) ;
    /*---Part 1---*/
    fp1 = fopen ("elem_shapeNN.dat", "r") ;
    fscanf (fp1, "%d" ,&nel) ;
    fp2 = fopen ("elem_npgNN.dat","r") ;
    fscanf (fp2, "%d" ,&nel);
    printf(" %d", n_in);
    
    /*---Part 2-----*/
    for (i=0;i<nel;i++) {
        printf ("%d", i);
        fscanf (fp1, "%d", &ia) ;
        for(j=0;j<n_in;j++) {
            fscanf (fp1, "%e", &f1);
            printf(" %e", f1);
        }
        fscanf (fp2, "%d",&ia);
        for(j=0;j<n_out;j++) {
            fscanf (fp2, "%d", &ib);
            printf(" %d", ib);
        }
        printf ("\n");
    }
    /*---Part 3----- */
    fclose(fp1);
    fclose(fp2);
    
    return 0;
}
