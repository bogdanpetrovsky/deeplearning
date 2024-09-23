//
//  DAneuro.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 15.09.2024.
//

#include "DAneuro.hpp"
#include <math.h>
#define FNAMELENGTH 100
#define NHU_V 1
#define NHU_C 0
#define Momi 0.1
#define Mom2 0.1
#include "DLneuro.cpp"
/*ーーー*/
int main()
{
    int i, j,k, il, j1, rseed, MaxPattern, nIU, nOU, *nHU, nHL, nHUO,nhflag;
    double *zOU, **zIU, **zHU, ***w, **bias, **zdHU, *zdOU;
    char fname1 [FNAMELENGTH] , fname2 [FNAMELENGTH] ;
    FILE *fp;
    
    scanf ("%d %d %d %d %d %d %s %s", &MaxPattern, &nIU, &nHUO, &nOU, &nHL, &nhflag, fname1, fname2) ;
    
    nHU = ivector (0, nHL+1) ;
    if (nhflag == NHU_V) {
        for (i=1; i<=nHL;i++) scanf ("%d" ,nHU+i);
    } else {
        for(i = 1; i<=nHL;i++) nHU[i] = nHUO;
    }
    nHU [0] = nIU ;
    nHU [nHL+1] = nOU ;
    printf(fname1);
    printf(fname2);
    
    zIU = matrix (0,MaxPattern-1,0,nIU-1) ;
    zHU = (double **)malloc ((nHL+2)*sizeof(double *)) ;
    for (i=0;i<nHL+2;i++) zHU[i] = vector (0,nHU[i]-1) ;
    zdHU = (double **)malloc ((nHL+2) *sizeof (double *));
    for (i=0;i<nHL+2;i++) zdHU[i] = vector (0,nHU[i]-1) ;
    zOU = vector (0, nOU-1) ;
    zdOU = vector (0, nOU-1) ;
    w = (double ***)malloc ((nHL+1) *sizeof (double **)) ;
    for(i=0;i<=nHL;i++) w[i] = matrix(0, nHU[i+1]-1,0,nHU[i]-1);
    bias = (double **)malloc ((nHL+2) *sizeof (double *));
    for (i=0;i<=nHL+1;i++) bias [i] = vector (0,nHU[i]-1) ;
    read_fileA (fname1, zIU,nIU,MaxPattern);
    printf("%d %d %d \n", nIU, nOU, nHL);
    load_weight (fname2,w, bias, nIU, nHU, nOU, nHL) ;
    for (i = 0; i< nHL+1;i++) printf ("%f ", w[i][1]);
    printf("\n");
    for (i = 0; i< nHL+2;i++) printf ("%f ", bias[i][1]);
    printf("\n");
    
    for (i=0; i<MaxPattern; i++){
        propagation (i, zIU, zHU, zdHU, zOU, zdOU, w, bias, nIU, nHU,nOU, nHL);
        printf ("%d", i);
        
        for(j=0;j<nOU;j++) printf(" %e", zOU[j]);
        printf ("\n" ) ;
    }
        
    return 0;
}
