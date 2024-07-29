//
//  main.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 02.06.2024.
//

#include "elemgen.hpp"
#include<cmath>
#include<math.h>

using namespace std;

void esm3D08 (int *elem, double **node, double *mate, double **esm, int ngauss, double *gc, double *gw, int nfpn) {
    int i,j,k,ii,jj,kk,counter,necm=6,nnpe=8,kdim=24;
    double e, v, ee, det, coord[8][3],J[3][3], invJ[3][3], D[6][6],s, t,u, ra, rs, sa,ss, ua, us,N[8][7],B[6][24], DB[6][24] ;
    double dtmp, ra2, rs2, ssus, ssua, saus, saua;
    
    for (ii=0;ii<kdim; ii++) {
        for(jj=0;jj < kdim;jj++) esm[ii][jj] = 0.0;
    }
    for(i=0;i<nnpe; i++) {
        ii = elem[i];
        for (j=0;j<nfpn; j++) coord [i][j] = node[ii][j];
    }
    for (i=0;i<necm; i++)
        for (j=0;j<kdim;j++) B[i][j] = 0.0;
    
    e = mate[0];
    v = mate[1];
    for (i=0;i<necm;i++)
        for (j=0;j<necm;j++) D[i][j] = 0.0;
    ee = e* (1.0 - v) / (1.0 + v) / (1.0 - 2.0*v);
    
    D[0][0] = ee;
    D[1][1] = ee;
    D[2][2] = ee;
    D[3][3] = ee*(1.0 - 2.0*v)/2.0/(1.0 - v);
    D[4][4] = D[3][3];
    D[5][5] = D[3][3];
    D[0][1] = ee*v/(1.0 - v);
    D[0][2] = D[0][1];
    D[1][2] = D[0][1];
    D[1][0] = D[0] [1];
    D[2][0] = D[0][2];
    D[2][1] = D[1][2];
    
    for (i=0; i<ngauss; i++){
        s = gc[i];
        for (j=0;j<ngauss; j++) {
            t= gc[j];
            for (k=0;k<ngauss; k++) {
                u = gc [k];
                ra = (1.0 + s)*0.5;
                rs = (1.0 - s)*0.5;
                sa = (1.0 + t)*0.5;
                ss = (1.0 - t)*0.5;
                
                ua = (1.0 + u) *0.5;
                us = (1.0 - u) *0.5;
                ssus = ss*us;
                saus = sa*us;
                ssua = ss*ua;
                saua = sa*ua;
                N[0][6] = rs*ssus;
                N[1][6] = ra*ssus;
                N[2][6] = ra*saus;
                N[3][6] = rs*saus;
                N[4][6] = rs*ssua;
                N[5][6] = ra*ssua;
                N[6][6] = ra*saua;
                N[7][6] = rs*saua;
                N[0][0] = -0.5*ssus;
                N[1][0] = 0.5*ssus;
                N[2][0] = 0.5*saus;
                N[3][0] = -0.5*saus;
                N[4][0] = -0.5*ssua;
                N[5][0] = 0.5*ssua;
                N[6][0] = 0.5*saua;
                N[7][0] = -0.5*saua;
                rs2 = 0.5*rs;
                ra2 = 0.5*ra;
                N[0][1] = -rs2*us;
                N[1][1] = -ra2*us;
                N[2][1] = ra2*us;
                N[3][1] = rs2*us;
                N[4][1] = -rs2*ua;
                N[5][1] = -ra2*ua;
                N[6][1] = ra2*ua;
                N[7][1] = rs2*ua;
                N[0][2] = -rs2*ss;
                N[1][2] = -ra2*ss;
                N[2][2] = -ra2*sa;
                N[3][2] = -rs2*sa;
                N[4][2] = rs2*ss;
                N[5][2] = ra2*ss;
                N[6][2] = ra2*sa;
                N[7][2] = rs2*sa;
                
                for(ii=0;ii<nfpn;ii++) {
                    for (jj=0;jj<nfpn;jj++) {
                        J[ii][jj] = 0.0;
                        for (kk=0; kk<nnpe; kk++) {
                            J[ii][jj] += N[kk][ii]*coord[kk][jj];
                        }
                    }
                }
            }
            
            det = J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1] - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0];
            
            invJ[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1])/det;
            invJ[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2])/det;
            invJ[0][2] = (J[0][1]*J[1][2] - J[1][1]*J[0][2])/det;
            invJ[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2])/det;
            
            invJ[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0])/det;
            invJ[1][2] = (J[1][0]*J[0][2] - J[0][0]*J[1][2])/det;
            invJ[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0])/det;
            invJ[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1])/det;
            invJ[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0])/det;
            
            for (ii=0;ii<nnpe; ii++){
                N[ii][3] = 0.0;
                N[ii][4] = 0.0;
                N[ii][5] = 0.0;
                for(jj=0;jj<3;jj++) {
                    N[ii][3] += invJ[0][jj]*N[ii][jj];
                    N[ii][4] += invJ[1][jj]*N[ii][jj];
                    N[ii][5] += invJ[2][jj]*N[ii][jj];
                }
            }
            
            for (ii=0;ii<nnpe;ii++){
                jj = ii*nfpn;
                B[0][jj] = N[ii][3];
                B[1][1+jj] = N[ii][4];
                B[2][2+jj] = N[ii][5];
                B[3][jj] = N[ii][4];
                B[3][1+jj] = N[ii][3];
                B[4][1+jj] = N[ii][5];
                B[4][2+jj] = N[ii][4];
                B[5][jj] = N[ii][5];
                B[5][2+jj] = N[ii][3];
            }
            
            for (ii=0;ii<necm;ii++) {
                for(jj=0;jj <kdim;jj++) {
                    DB[ii][jj] = 0.0;
                    for (kk=0; kk<necm; kk++) {
                        DB[ii][jj] += D[ii][kk]*B[kk][jj];
                    }
                }
            }
            dtmp = gw[i]*gw[j]*gw[k]*det;
            for (ii=0;ii<kdim;ii++) {
                for(jj=0;jj<kdim;jj++) {
                    for (kk=0; kk<necm; kk++) {
                        esm[ii][jj] += B[kk][ii]*DB[kk][jj]*dtmp;
                    }
                }
            }
        }
    }
}

/* ElementShapeMetric.c */
void inv_mat3 (double invA[][3], double A[][3]) {
    double det;
    det = A[0] [0]*A[1] [1]*A[2] [2] + A[0] [1]*A[1][2]*A[2] [0]
    + A[0] [2]*A[1][0]*A[2][1] - A[0] [0]*A[1][2]*A[2] [1]
    - A[0] [1]*A[1] [0]*A[2] [2] - A[0] [2]*A[1] [1]*A[2] [0] ;
    invA[0] [0] = (A[1] [1]*A[2] [2] - A[1] [2] *A[2] [1]) /det;
    invA[0] [1] = (A[0] [2]*A[2] [1] - A[0] [1] *A[2] [2]) /det;
    invA[0] [2] = (A[0] [1]*A[1] [2] - A[1] [1] *A[0] [2]) /det;
    invA[1][0] = (A[1][2]*A[2] [0] - A[1] [0]*A[2] [2])/det;
    invA[1][1] = (A[0] [0]*A[2] [2] - A[0] [2]*A[2] [0])/det;
    invA[1] [2] = (A[1] [0]*A[0] [2] - A[0] [0]*A[1] [2]) /det;
    invA[2] [0] = (A[1] [0]*A[2] [1] - A[1] [1] *A[2] [0]) /det;
    invA[2] [1] = (A[0] [1]*A[2] [0] - A[0] [0]*A[2] [1]) /det;
    invA[2] [2] = (A[0] [0]*A[1] [1] - A[0] [1] *A[1] [0]) /det;
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
        for(j = 0;i < 3;j++)
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
            e[i][j] = ndata[(i+1)%4] [j] - ndata[i][j];
        
        for (j=0,d1=0.0;j<nfpn;j++) d1 += e[i][j] *e[i][j];
        dl[i] = sqrt(d1) ;
        for(j=0;j<nfpn;j++) e[i][j] /= dl[i];
    }
    for(i=4; 1<=7;i++) {
        for(j=0;j<nfpn;j++)
            e[i][j] = ndata[(1+1) % 4 + 4][j] - ndata[i][j];
        for(j=0,d1=0.0;j<nfpn;j++) d1 += e[i][j]*e[i][j];
        dl[i] = sqrt(d1);
        for (j=0;j<nfpn;j++) e[i][j] /= dl[i];
    }
    for(i=0;i<=3;i++) {
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
        for (i=0;i<nfpn;i++) dx0 [i] = sg[n][0]*e[idx[n][0]][i];
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
    dl_min = 1.030;
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

/* bspline.c */
#define NMAX_KV 100
#define KV_EPS 1.0e-5
#define NURBS_EPS 1.0e-20

/*-----------*/
void Bspline00 (double *N, double *dN, double xi, int ni,int p,int nkv, double *KV) {
    int i,j,k;
    double d, e,d_d, d_e, temp[NMAX_KV], dtemp [NMAX_KV], dw;
    for (i=ni-p;i<=ni+p;i++){
        temp[i] = 0.0;
        dtemp[i] = 0.0;
    }
    for(i=ni-p;i<=ni;i++) {
        if ((xi>=KV[i]) && (xi<KV[i+1])) {
            temp[i] = 1.0;
        } else {
            temp[i] =0.0;
        }
        dtemp[i] = 0.0 ;
    }
    for (k=1; k<=p;k++) {
        for (i=ni-k;i<=ni;i++) {
            if (fabs (temp [i]) > NURBS_EPS) {
                dw = KV[i+k] - KV[i];
                if (fabs (dw) > KV_EPS) {
                    d = ((xi - KV[i]) *temp[i]) /dw ;
                    d_d = k*temp[i]/dw ;
                } else {
                    d = 0.0 ;
                    d_d = 0.0 ;
                }
            } else {
                d = 0.0;
                d_d = 0.0;
            }
            if (fabs (temp [i+1]) > NURBS_EPS) {
                dw = KV[i+k+1] - KV[i+1] ;
                if (fabs (dw) > KV_EPS) {
                    e = ( (KV[i+k+1] - xi) * temp[i+1]) /dw ;
                    d_e = k*temp[i+1] / dw ;
                } else {
                    e = 0.0;
                    d_e = 0.0 ;
                }
            } else {
                e = 0.0;
                d_e =0.0;
            }
            temp [i] = d + e;
            dtemp [i] = d_d - d_e;
        }
    }
    for (i=0;i<=p;i++) {
        N[i] = temp[ni-p+i] ;
        dN[i] = dtemp [ni-p+i] ;
    }
}
/*ーー*/
void BsplinelD(double *N, double *dN, int n, double xi, int p, int nkv, double *KV) {
    int i,j, k, nb;
    double d1, a2, d3;
    for (i=0;i<n;i++){
        N[2] = 0.0;
        dN[i] = 0.0;
    } for (i=0;i<nkv;i++)
        if ((xi >= KV[i]) && (xi < KV[i+1])) break ;
    Bspline00 (N+ (i-p) ,dN+ (i-p) ,xi, i,p,nkv,KV);
}
/*------------------------*/
void Bspline3D(double ***B, double ***dB_xi, double ***dB_et, double ***dB_ze, double *N, double *dN, int n, double *M, double *aM, int m, double *L, double *dL, int l) {
    
    int i,j,k;
    double d1, d2, d3, dd1, dd2, dd3 ;
    for(i=0;i<n;i++) {
        d1 = N[i] ;
        dd1 = dN [i];
        for(j=0;j<m;j++) {
            d2 = M[j] ;
            dd2 = aM [j] ;
            for (k=0;k<1;k++) {
                d3 = L[k] ;
                dd3 = dL[k] ;
                B[i][j][k] = d1*d2*d3;
                dB_xi[i][j][k] = dd1*d2*d3;
                dB_et[i][j][k] = d1*dd2*d3;
                dB_ze[i][j][k] = d1*d2*dd3;
            }
        }
    }
}
/*------------------------*/
/* nurbs.c */
/*---*/
void Nurbs1D(double *R, double *dR, int n, double *N, double *dN, double xi,int p, int nkv,
             double *KV, double *w) {
    int i;
    double dsum, dsum2, ddsum;
    BsplinelD(N, dN, n, xi, p, nkv, KV) ;
    for (i=0,dsum=0.0;i<n;i++) dsum += N[i]*w[i];
    dsum2 = 1.0/dsum/dsum ;
    for (i=0;i<n; i++) R[i] = N[i] *w[i]/dsum ;
    for (i=0, ddsum=0.0;i<n; i++) ddsum += dN[i]*w[i] ;
    for (i=0;i<n;i++) dR [i] = (dN[i]*w[i]*dsum - N[i]*w[i]*ddsum) *dsum2 ;
}
/*--*/
void Nurbs3D(double ***R3,double ***dR3_xi,double ***dR3_et, double ***dR3_ze,double *N,double *dN,int n,double *M,double *dM,int m,double *L,double *dL,int l,double *wl,double *w2, double *w3) {
    int i,j,k;
    double d1, d2, d3,dd1,dd2, dd3, dsum, dsum2, ddsum1, ddsum2 , ddsum3 ;
    for (i=0,dsum=0.0,ddsum1=0.0,ddsum2=0.0,ddsum3=0.0;i<n;j++) {
        d1 = N[i]*wl[i] ;
        dd1 = dN[i] * wl[i] ;
        for(j=0;j<m;j++) {
            d2 = M[j]*w2[j] ;
            dd2 = dM[j]*w2[j];
            for (k=0; k<l;k++) {
                d3 = L[k]*w3[k];
                dd3 = dL[k] *w3 [k];
                dsum += d1*d2*d3;
                ddsum1 += dd1*d2*d3;
                ddsum2 += d1*dd2*d3;
                ddsum3 += d1*d2*dd3;
            }
        }
    }
    dsum2 = 1.0/dsum/dsum;
    for (i=0;i<n;i++){
        d1 = N[i]*wl[i] ;
        dd1 = dN[i] *wl[i];
        for(j=0;j<m;j++) {
            d2 = M[j]*w2[j] ;
            dd2 = dM[j]*w2[j];
            for (k=0;k<1; k++) {
                d3 = L[k] *w3 [k] ;
                dd3 = dL[k] *w3 [k] ;
                R3 [i][j] [k] = d1*d2*d3/dsum ;
                dR3_xi [i][j][k] = (dd1*dsum - d1*ddsum1) *dsum2 ;
                dR3_et [i][j] [k] = (dd2*dsum - d2*ddsum2) *dsum2 ;
                dR3_ze [i] [j][k] = (dd3*dsum - d3*ddsum3) *dsum2;
            }
        }
    }
}
/*---------*/
/* DLcommon.c */
/*-ー*/
#define rnd() (drand48() * (Wmax - Wmin) + Wmin)
#define noise () ((drand48 ()-0.5f) *2.0f*NoiseLevel)
#define FNAMELENGTH 100
#define NHU_V 1
#define NHU_C 0
float Mom1=0.05f; float Mom2=0.05f;
float Wmin = -0.10f ;
float Wmax = 0.10f ;

void s_shuffle(int *a,int n) {
    int i, ic,tsize, itemp;
    for (i=0,tsize=n;i<n-1;i++,tsize--) {
        ic = floor (drand48 ()*(tsize-1) + 0.1) ;
        itemp = a[tsize-1] ;
        a[tsize-1] = a[ic] ;
        a [ic] = itemp ;
    }
}
/*--*/
void a0f (float *fv, float *fvd, float x) {
    float dd;
    dd = (1.0f+ (float)tanh(x/2.0f))/2.0f;
    *fv = dd;
    *fvd = dd* (1.0 - dd) ;
}
/*--*/
void alf ( float *fv, float *fvd, float x) {
    float dd;
    dd = (1.0f+ (float) tanh (x/2.0f))/2.0f;
    *fv = dd;
    *fvd = dd* (1.0 - dd);
}
/*--*/
void read_file( char *name, float **o,float **t,int nIU, int nOU, int npattern)
{
    int i,j,k;
    FILE *fp;
    fp = fopen( name, "r" ) ;
    for (i=0; i<npattern; i++) { fscanf (fp, "%d", &k) ;
        for (j=0;j<nIU;j++) fscanf (fp, "%e", o[i]+j);
        for (j=0;j<nOU;j++) fscanf (fp, "%e", t[1]+j);
    }
    fclose(fp);
}
/*-*/
void initialize(float ***w, float **bias, int nIU, int *nHU, int nOU, int nHL) {
    int i,j,k;
    for (i=0;i<=nHL;i++)
        for (j=0;j<nHU[i+1];j++)
            for (k=0;k<nHU[i];k++) w[i][j][k] = rnd();
    for(j=1;j<=nHL+1;j++)
        for(i=0;i<nHU[j];i++) bias[j][i] = rnd();
}
/**/
void store_weight (float ***w, float **bias, float ***w_min, float **bias_min, int nIU, int *nHU, int noU, int nHL) {
int i,j,k;
for (i=0;i<=nHL;i++)
    for(j=0;j<nHU[i+1];j++)
        for (k=0;k<nHU[i];k++) w_min [i][j][k] =w[i][j][k];
    for(j = 1;j<=nHL+1;j++)
        for (i=0;i<nHU[j];i++) bias_min[j][i] =bias[j][i] ;
}
/*--*/
void show_results ( float ***w, float **bias, float ***w_m, float **bias_m, int nIU, int *nHU, int nou, int nHL) {
    int i,j,k,iL;
    for (iL=0;iL<=nHL;iL++) {
        for (i=0; i<nHU[iL]; i++) {
            printf ("%5d", i);
            for(j=0;j<nHU[iL+1];j++)
                printf("%e", w_m[iL][j][i]);
            printf ("\n");
        }
    }
    for (iL=1;iL<=nHL+1;iL++){
        for (j=0; j<nHU[iL];j++) printf("%e ",bias_m[iL][j]);
        printf ("\n") ;
    }
    for (iL=0;iL<=nHL;iL++) {
        for (i=0;i<nHU[iL];i++) {
            printf ("%5d", i);
            for (j=0;j<nHU[iL+1];j++)
                printf(" %e",w[iL][j][i]);
            printf ("\n");
        }
    }
    for (iL=1; iL<=nHL+1;iL++) {
        for(j=0;j<nHU[iL];j++) printf("%e ", bias[iL][j]);
        printf ("\n");
    }
}
/*.*/
void clear_dweight ( float ***dw, float **dbias, int nIU, int *nHU, int nou, int nHL) {
    int i,j,k;
    for (i=0;i<=nHL;i++)
        for(j=0;j<nHU[i+1];j++)
            for (k=0;k<nHU[i];k++) dw[i][j][k] = 0.0 ;
    for (j=1;j<=nHL+1;j++)
        for (i=0;i<nHU[j];i++) dbias[j][i] = 0.0 ;
}
/**/
/* DLebp.c */
void propagation ( int p, float **zIU, float **zHU, float **zdHU, float *zOU, float *zdOU, float ***w, float **bias, int nIU,int *nHU, int nOU, int nHL)
{
    int i,j, iL;
    float net;
    for (i=0;i<nHU[1];i++) {
        for(net=0,j=0;j<nIU;j++)
            net += w[0][i][j] * zIU[p][j];
        net += bias [1][i];
        a0f (zHU[1]+i, zdHU[1]+i, net);
    }
    for (iL=2;iL<=nHL;iL++) {
        for (i=0;i<nHU[iL];i++) {
            for (net=0, j=0;j<nHU[iL-1];j++)
                net += w[iL-1] [i][j] * zHU[iL-1] [j];
            net += bias [iL][i];
            a0f(zHU[iL]+i, zdHU[iL]+i, net);
        }
    }
    for (i=0;i<nOU;i++) {
        for(net=0, j=0;j<nHU[nHL];j++)
            net += w[nHL][i][j]* zHU[nHL][j];
            net += bias[nHL+1] [i];
        alf(zOU+i, zdOU+i, net);
    }
}
    /*-*/
void back_propagation ( int p,float **t,float **zIU, float **zHU, float **zdHU, float *zOU, float *zdOU, float ***w, float **bias, float ***dw, float **dbias, float **d, int nIU, int *nHU, int nOU, int nHL, float Alpha, float Beta) {
    int i,j, iL; float sum;
    for (i=0;i<nOU;i++)
        d[nHL+1][i] = (t[p] [i] - zOU[i]) * zdOU[i];
    for (i=0;i<nHU[nHL];i++) {
        for (sum=0.0f, j=0;j<nOU; j++)
            sum += d[nHL+1] [j] * w[nHL] [j] [i];
        d [nHL] [i] = zHU [nHL] [i] * sum;
    }
    for (iL=nHL-1;iL>=1;iL--){
        for (j=0;j<nHU[iL];j++) {
            for (sum=0.0f,i=0;i<nHU[iL+1];i++)
                sum += d[iL+1] [i] * w[iL][i][j];
            d[iL] [j] = zdHU[iL] [j] * sum;
        }
    }
    for (iL=nHL+1; iL>=1; iL--) {
        for (i=0;i<nHU[iL];i++){
            dbias [iL][i] = Beta*d[iL] [i] + Mom1*dbias[iL][i];
            bias [iL][i] += dbias [iL] [i];
        }
    }
    for (iL=nHL; iL>=1; iL--){
        for (j=0;j<nHU[iL+1];j++) {
            for (i=0;i<nHU[iL];i++) {
                dw[iL] [j] [i] = Mom2*dw[iL] [j][i];
                dw[iL][j][i] += Alpha*d[iL+1][j]*zHU[iL][i];
                w[iL][j][i] += dw[iL][j][i];
            }
        }
    }
    for (j=0; j<nHU[1];j++) {
        for (i=0;i<nIU;i++) {
            dw[0][j][i] = Alpha*d[1][j]*zIU[p][i] + Mom2*dw[0][j][i];
            + Mom2*dw[0][j][i];
            w[0][j][i] += dw[0][j][i];
        }
    }
}
/*--*/

int main(int argc, const char * argv[]) {
//    elemgen();
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
        for(j=0;j<7;j++) printf("%e", eshape[j]);
        printf ("\n");
    }

    return 0;
}
