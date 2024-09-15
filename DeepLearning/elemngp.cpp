//
//  elemngp.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 31.07.2024.
//

/* elemngp.c */
#include <stdio.h>
#include <stdlib.h>

int main (void) {
    int i,j,ia, ib,nel;
    double er[30], th=1.0e-7;
    scanf ("%d", &nel) ;
    printf ("%d\n" ,nel) ;
    for (i=0;i<nel;i++) {
        for (j=2; j<30;j++) scanf ("%d %d %le", &ia, &ib, er+j) ;
        for (j=2;j<30;j++)if(er[j]<th) break;
        printf("%d %d\n",i,j);
    }
    
    return 0;
}
