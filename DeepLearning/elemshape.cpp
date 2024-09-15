//
//  elemshape.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 31.07.2024.
//


#include "elemgen.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define rnode drand48()*(cmax-cmin)+cmin

/* ElementShapeMetric.c */
void inv_mat3 (double invA[][3], double A[][3]) {
    double det;
    det = A[0] [0]*A[1] [1]*A[2] [2] + A[0] [1]*A[1][2]*A[2] [0]
    + A[0] [2]*A[1][0]*A[2][1] - A[0] [0]*A[1][2]*A[2] [1]
    - A[0] [1]*A[1] [0]*A[2] [2] - A[0] [2]*A[1] [1]*A[2] [0] ;
    
    invA[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1]) /det;
    invA[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2]) /det;
    invA[0][2] = (A[0][1]*A[1][2] - A[1][1]*A[0][2]) /det;
    invA[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2]) /det;
    invA[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0]) /det;
    invA[1][2] = (A[1][0]*A[0][2] - A[0][0]*A[1][2]) /det;
    invA[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0]) /det;
    invA[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1]) /det;
    invA[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0]) /det;
}

void matMulti3 (double C[][3], double A[][3], double B[][3]) {
    int i,j,k;
    for(i=0;i<3;i++) {
        for (j = 0;j < 3;j++) {
            C[i][j] = 0.0;
            for (k=0;k<3; k++) C[i][j] += A[i][k]*B[k][j];
        }
    }
}
double f_norm3 (double A[][3] ) {
    int i,j;
    double dsum;
    dsum = 0.0;
    for (i=0;i<3;i++)
        for(j = 0;j < 3;j++)
            dsum += A[i][j]*A[i][j];
    
    return dsum ;
}

double shape_metric( int *elem, double **node, int nnpe, int nfpn) {
    int i, j,k, ia, ib, ic, id, inode;
    int idx[8][3] = { {1,3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2,7}, {7,5,0}, {4,6, 1}, {5,7, 2}, {6, 4,3}};
    double kp[8],W[3][3],A[3][3],T[3][3],invW[3][3],
    invT[3][3],cond;
    double cube[8][3] = {
        {0.0, 0.0, 0.0},{1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}
    };
    
    for (inode=0; inode<nnpe; inode++) {
        ia = elem[inode] ;
        ic = inode ;
        for (j=0;j<nfpn;j++) {
            ib = elem[idx[inode] [j]] ;
            id = idx[inode] [j] ;
            for (i=0;i<nfpn;i++)
                A[i][j] = node[ib][i] - node[ia][i] ;
            for (i=0;i<nfpn;i++)
                W[i][j] = cube[id][i] - cube[ic][i] ;
        }
        
        inv_mat3 (invW,W) ;
        matMulti3 (T,A, invW) ;
        inv_mat3 (invT,T) ;
        kp [inode] = f_norm3 (T) *f_norm3 (invT) ;
    }
    for (inode=0,cond=0.0;inode<nnpe;inode++)
        cond += kp[inode] /9;
    
    return 8.0/cond ;
}

/* ElementShape.c */
void check_shape (double *s_data,int *elem,double **node,int nfpn) {
    
    int i,j,k,i0,i1,i2,n,ne,nnpe=8;
    double dl[12], da[48], de[24] , e[12][3],d1,d2,dx0[3], dx1[3], dx2[3], nv[8][3][3], v1[3],v2[3],v3[3], v4[3], ndata[8][3] , dl_min, dl_max, da_min, da_max, de_max, de_min;
    int idx[8][3] ={
        {0,3,8}, {1,0,9}, {2,1, 10}, {3,2,11},
        {7,4,8}, {4,5,9}, {5,6,10}, {6,7,11}
    };
    int sg [8] [3] = {
        { 1, -1, 1}, { 1, -1, 1}, { 1,-1, 1}, { 1,-1, 1},
        {-1, 1, -1}, {-1, 1, -1}, {-1, 1, -1}, {-1, 1, -1}
    } ;
    int id[12][6] = {
        {0,2, 1, 1,0, 1}, {1,2, 1, 2, 0, 1},{2,2, 1, 3, 0, 1}, {3, 2, 1, 0, 0, 1},
        {4,1,0, 5,1,2}, {5,1, 0, 6,1,2}, {6,1,0, 7,1,2}, {7,1,0, 4,1,2},
        {0,0,2, 4,2, 0}, {1,0,2, 5,2,0}, {2,0, 2, 6,2, 0}, {3, 0, 2, 7,2, 0}
    };
    for (i=0;i<nnpe;i++){
        k = elem[i] ;
        for(j=0;j<nfpn;j++) ndata[i][j] = node[k][j];
    }
    
    for (i = 0; i<=3; i++) {
        for (j=0;j<nfpn;j++)
            e[i][j] = ndata[(i+1)%4][j] - ndata[i][j];
        
        for (j=0,d1=0.0;j<nfpn;j++) d1 += e[i][j] *e[i][j];
        dl[i] = sqrt(d1) ;
        for(j=0;j<nfpn;j++) e[i][j] /= dl[i];
    }
    for(i=4; i<=7;i++) {
        for(j=0;j<nfpn;j++)
            e[i][j] = ndata[(1+1) % 4 + 4][j] - ndata[i][j];
        for(j=0,d1=0.0;j<nfpn;j++)
            d1 += e[i][j]*e[i][j];
        dl[i] = sqrt(d1);
        for (j=0;j<nfpn;j++)
            e[i][j] /= dl[i];
    }
    for(i = 0;i<=3;i++) {
        for(j=0;j<nfpn;j++)
            e[i+8][j] = ndata[i+4][j] - ndata[i][j];
        for(j=0,d1=0.0;j<nfpn;j++)
            d1 += e[i+8][j] *e[i+8][j];
        dl[i+8] = sqrt(d1);
        for (j=0;j<nfpn;j++)
            e[i+8][j] /=dl[i+8];
    }
    /*- anglel --*/
    for (n=0;n<nnpe;n++) {
        for (i=0;i<nfpn;i++) dx0[i] = sg[n][0]*e[idx[n][0]][i];
        for (i=0;i<nfpn;i++) dx1[i] = sg[n][1]*e[idx[n][1]][i];
        for (i=0;i<nfpn;i++) dx2[i] = sg[n][2]*e[idx[n][2]][i];
        
        de [n*3+0] = dx0[0] * dx1[0] + dx0[1] * dx1[1] + dx0[2] * dx1[2];
        de [n*3+1] = dx1[0] * dx2[0] + dx1[1] * dx2[1] + dx1[2] * dx2[2];
        de [n*3+2] = dx2[0] * dx0[0] + dx2[1] * dx0[1] + dx2[2] * dx0[2];
        nv[n][0][0] = dx2[1] * dx1[2] - dx2[2] * dx1[1];
        nv[n][0][1] = dx2[2] * dx1[0] - dx2[0] * dx1[2];
        nv[n][0][2] = dx2[0] * dx1[1] - dx2[1] * dx1[0];
        
        nv[n][1][0] = dx1[1] * dx0[2] - dx1[2] * dx0[1];
        nv[n][1][1] = dx1[2] * dx0[0] - dx1[0] * dx0[2];
        nv[n][1][2] = dx1[0] * dx0[1] - dx1[1] * dx0[0];
        
        nv[n][2][0] = dx0[1] * dx2[2] - dx0[2] * dx2[1];
        nv[n][2][1] = dx0[2] * dx2[0] - dx0[0] * dx2[2];
        nv[n][2][2] = dx0[0] * dx2[1] - dx0[1] * dx2[0];
        for (j = 0;j <3; j++) {
            for (i=0,d1=0.0;i<nfpn;i++) d1 += nv[n][j][i] * nv[n][j][i];
            d2 = 1.0/sqrt(d1);
            for (i=0; i<nfpn;i++) nv[n] [j][i] *= d2 ;
        }
    }
    /* angle2 */
    for (ne=0;ne<12;ne++) {
        for (i=0;i<nfpn;i++) v1[i] = nv[id[ne][0]][id[ne][1]][i];
        for (i=0;i<nfpn;i++) v2[i] = nv[id[ne][0]][id[ne][2]][i];
        for (i=0;i<nfpn;i++) v3[i] = nv[id[ne][3]][id[ne][4]][i];
        for (i=0;i<nfpn;i++) v4[i] = nv[id[ne][3]][id[ne][5]][i];
        
        da[ne*4+0] = v1[0] *v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        da[ne*4+1] = v3[0] *v4[0] + v3[1] * v4[1] + v3[2] * v4[2];
        da[ne*4+2] = v1[0] *v4[0] + v1[1] * v4[1] + v1[2] * v4[2];
        da[ne*4+3] = v2[0] *v3[0] + v2[1] * v3[1] + v2[2] * v3[2];
    }
    /* -minmax- */
    dl_max = -1.0;
    dl_min = 1.0e30;
    for(i=0;i<12;i++) {
        if(dl[i] > dl_max) dl_max = dl[i];
        if(dl[i] < dl_min) dl_min = dl[i];
    }
    de_max = -1.0e30 ;
    de_min = 1.0e30 ;
    for (i=0;i<24;i++) {
        if (de[i] > de_max) de_max = de[i];
        if(de[i] < de_min) de_min = de[i];
    }
    da_max = -1.0e30;
    da_min = 1.0e30;
    for (i=0; i<48;i++) {
        if (da[i] > da_max) da_max = da[i];
        if(da[i] < da_min) da_min = da[i];
    }
    s_data[0] = dl_min;
    s_data[1] = dl_max;
    s_data[2] = acos(da_min) * 180 / 3.1415926;
    s_data[3] = acos(da_max) * 180 / 3.1415926;
    s_data[4] = acos(de_min) * 180 / 3.1415926;
    s_data[5] = acos(de_max) * 180 / 3.1415926;
}


int main(void) {
    int i,j,k,ia, ib, nel, nfpn=3, nnpe=8, elem[8]={0,1,2,3,4,5, 6,7} ;
    double eshape [7], **node;
    scanf ("%d",&nel) ;
    printf ("%d\n", nel) ;
    node = (double **)malloc (nnpe*sizeof (double *) );
    for (i=0;i<nnpe;i++) node[i] = (double *)malloc (nfpn*sizeof(double)) ;
    for (i=0;i<nel;i++) {
        for (j=0;j<nnpe;j++) {
            scanf ("%d %d",&ia, &ib);
            for (k=0;k<nfpn;k++) scanf ("%le", node[j]+k);
        }
        check_shape(eshape, elem, node, nfpn) ;
        eshape [6] = shape_metric (elem, node, nnpe, nfpn);
        printf ("%d", i);
        for(j=0;j<7;j++) printf(" %e", eshape[j]);
        printf ("\n");
    }

    return 0;
}
