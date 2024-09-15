//
//  shapeNN.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 31.07.2024.
//

#include <stdio.h>
#include <stdlib.h>

int main (void) {
    int i, j,k,ia, ib,nel;
    double **shape, smin [7], smax [7], swidth [7];
    FILE *fp;
    /*---Part 1----*/
    fp = fopen ("elem_shape.dat", "r") ;
    fscanf (fp, "%d",&nel) ;
    shape = (double **)malloc (nel*sizeof (double *) ) ;
    for (i=0;i<nel;i++) shape[i] = (double *)malloc(7*sizeof(double)) ;
    for (i=0; i<nel; i++) {
        fscanf (fp, "%d",&ia);
        for (j=0;j<7;j++) fscanf(fp, "%le", shape[i]+j);
    }
    fclose (fp) ;
    /*---Part 2------*/
    for (i=0;i<7;i++) smax[i] = -1.0e30;
    for (i=0;i<7;i++) smin[i] = 1.0e30;
    for (i=0;i<nel; i++){
        for (j=0;j<7;j++) {
            if (shape [i] [j] > smax[j]) smax[j] = shape[i][j];
            if(shape[i][j] < smin[j]) smin[j] = shape[i][j];
        }
    }
    
    for (j=0;j<7;j++) swidth[j]= smax[j] - smin[j];
    /*---Part 3ーーーー--- -*/
    for (i=0;i<nel;i++) {
        for (j=0;j<7;j++) shape[i][j] =(shape[i][j] - smin[j])/swidth[j];
    }
    /*---Part 4------*/
    printf ("%d\n",nel);
    for (i=0;i<nel;i++){
        printf ("%d",i);
        for (j=0;j<7;j++) printf(" %e", shape[i][j]);
        printf ("\n");
    }
    for (j=0;j<7;j++) printf("%d %e %e\n",j, smin[j], smax[j]);
    
    return 0;
}
