//
//  npgNN.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 15.09.2024.
//

#include "npgNN.hpp"
/* ngpNN.c */
#include <stdio.h>
#include <stdlib.h>
int main (void) {
    int i,j,k,ia, ib,nel, *ngp, **t, gmin, gmax, gcat;
    FILE *fp;
    /*---Part 1-----*/
    fp = fopen ("elem_ngp.dat","r");
    fscanf (fp, "%d", &nel) ;
    ngp = (int *)malloc (nel*sizeof (int));
    for (i=0;i<nel;i++) fscanf (fp, "%d %d",&ia, ngp+i) ;
    fclose (fp);
    /*---Part 2ーー*/
    gmin = 1000; gmax = -1000 ;
    for (i=0;i<nel;i++) {
        if(ngp[i] > gmax) gmax = ngp[i] ;
        if(ngp[i] < gmin) gmin = ngp[i] ;
    }
    gcat = gmax - gmin + 1 ;
    /*---Part 3---*/
    t = (int **)malloc (nel*sizeof(int *));
    for (i=0;i<nel;i++) t[i] = (int *)malloc (gcat*sizeof (int));
    for (i=0;i<nel;i++) {
        for (j=0;j<gcat;j++) t[i] [3] = 0;
    }
    for (i=0;i<nel; i++) t[i][ngp[i]-gmin] = 1 ;
    /*---Part 4----*/
    printf ("%d\n",nel);
    for (i=0;i<nel; i++){
        printf ("%d", i);
        for(j=0;j<gcat;j++) printf(" %d", t[i][j]);
        printf ("\n") ;
    }
    printf ("%d %d\n", gmin, gmax) ;
    return 0;
}
