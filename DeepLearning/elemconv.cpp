//
//  elemconv.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 31.07.2024.
//

/* elemconv.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
}

/*-*/
void gausslg (double *weight, double *pos, int ngauss) {
    int i, j,k, ia, ib, ic;
    double w, a, ww,v,p,q,r,d0, d1, d2,p1,p2, u, df;
    ia = ngauss ;
    ib = ngauss/2 ;
    w = 1.0* ngauss;
    a = 3.1415926535897932/ (w+w) ;
    ww = w* (w + 1.0)*0.5 ;
    for (k=1;k<=ib;k++) {
        v = cos (a* (2*k-1)) ;
    Loop1: ;
        p = 1.0 ;
        q = v;
        for (j=2;j<=ngauss;j++) {
            r =((2*j-1) *v*q -(j-1) *p)/j;
            p = q;
            q = r;
        }
        u = (1.0 - v) * (1.0 + v);
        d0 = (p - v*q) *w/u ;
        d1 = (v*d0 - ww*q) /u ;
        d2 = (q*d1/ (d0*d0) +1.0) *q/d0 ;
        v -= d2 ;
        if (fabs (d2) >= 1.0e-16) goto Loop1 ;
        df = d2*v/u ;
        weight [k-1] = 2.0/ (w*d0* (1.0 - df*w) *p* (1.0 - df*2.0)) ;
        pos [k-1] = v ;
    }
    if (ib*2 < ngauss) {
        d0=1.0;
        for(j=1;j<=ib;j++) d0 = (1.0 + 0.5/j)*d0;
        weight [ib] = 2.0/ (d0*d0);
        pos [ib] = 0.0 ;
    }
    for(i=0;i<ib;i++){
        weight[ngauss-1-i] = weight[i];
        pos [ngauss-1-i] = pos [i] ;
        pos[i]*= -1.0;
    }
}
    /*.*/
double set_refdata ( double **esm, int edim) {
    int i,j;
    double d1 ;
    d1 = esm[0][0];
    for (i=0; i<edim; i++) {
        for (j=0; j<edim;j++) {
            if (esm[i][j] > d1) d1 = esm[i][j];
        }
    }
    return d1;
}
/**/
double check_esm(double **esm, double **esm0, double ref_value, int edim) {
    int i,j,k;
    double sum ;
    
    for (i=0, sum=0.0;i<edim; i++) {
        for (j=0;j<edim;j++) sum += fabs(esm[i][j] - esm0[i][j]);
    }
    return sum/ref_value ;
}

int main (){
    int i, j,k, ia, ib, ig, nel, nfpn=3, elem [8] ={0,1,2,3,4,5,6,7};
    double **esm, **esm0, *gc, *gw, **node, mate [2] ={2.0e11, 0.3} ;
    int max_ngp=30,nnpe=8, edim=24;
    double ref_value, chk_data;
    node = (double **)malloc (nnpe*sizeof (double *)) ;
    for (i=0;i<nnpe; i++) node[i] = (double *)malloc (nfpn*sizeof(double));
    esm = (double **)malloc (edim*sizeof (double *)) ;
    for (i=0; i<edim; i++)esm[i] = (double *)malloc (edim*sizeof(double));
    esm0 = (double **)malloc (edim*sizeof (double *)) ;
    for (i=0;i<edim; i++)esm0 [i] = (double *)malloc (edim*sizeof(double));
    gc = (double *)malloc (max_ngp*sizeof (double)) ;
    gw = (double *)malloc (max_ngp*sizeof (double)) ;
    scanf ("%d", &nel) ;
    printf ("%d\n" ,nel);
    for(i=0;i<nel;i++) {
        for(j=0;j<8;j++) {
            scanf ("%d %d" ,&ia, &ib);
            for (k=0;k<nfpn;k++) scanf ("%le", node[j]+k) ;
        }
        ig = max_ngp ;
        gausslg (gw,gc, ig) ;
        esm3D08 (elem, node, mate, esm0, ig, gc, gw, nfpn) ;
        ref_value = set_refdata (esm0, edim) ;
        for (ig=2; ig<max_ngp;ig++) {
            gausslg (gw,gc, ig) ;
            esm3D08 (elem, node, mate, esm, ig, gc, gw,nfpn);
            chk_data = check_esm (esm, esm0, ref_value, edim) ;
            printf ("%d %d %e\n", i, ig, chk_data);
        }
    }
}
