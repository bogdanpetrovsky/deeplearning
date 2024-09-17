//
//  DLneuro.cpp
//  DeepLearning
//
//  Created by Bogdan Petrovsky on 15.09.2024.
//

#include "DLneuro.hpp"
#include "nrutil.cpp"
#include<cmath>
#define rnode drand48()*(cmax-cmin)+cmin
#define rnd() (drand48() * (Wmax - Wmin) + Wmin)
#define noise() ((drand48 ()-0.5f) *2.0f*NoiseLevel)
#define FNAMELENGTH 100
#define NHU_V 1
#define NHU_C 0

float Mom3=0.05f;
float Mom1=0.05f;
float Wmin = -0.10f;
float Wmax = 0.10f;

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
void a0f (double *fv, double *fvd, double x) {
    float dd;
    dd = (1.0f+ (float)tanh(x/2.0f))/2.0f;
    *fv = dd;
    *fvd = dd* (1.0 - dd) ;
}
/*--*/
void alf ( double *fv, double *fvd, double x) {
    float dd;
    dd = (1.0f+ (float) tanh (x/2.0f))/2.0f;
    *fv = dd;
    *fvd = dd* (1.0 - dd);
}
/*--*/
void read_file( char *name, double **o,double **t,int nIU, int nOU, int npattern)
{
    int i,j,k;
    FILE *fp;
    fp = fopen( name, "r" ) ;
    for (i=0; i<npattern; i++) { fscanf (fp, "%d", &k) ;
        for (j=0;j<nIU;j++) fscanf (fp, "%le", o[i]+j);
        for (j=0;j<nOU;j++) fscanf (fp, "%le", t[1]+j);
    }
    fclose(fp);
}

void read_fileA( char *name, double **o, int nIU, int npattern)
{
    int i,j,k;
    FILE *fp;
    fp = fopen( name, "r" ) ;
    for (i=0; i<npattern; i++) { fscanf (fp, "%d", &k) ;
        for (j=0;j<nIU;j++) fscanf (fp, "%le", o[i]+j);
    }
    fclose(fp);
}

void load_weight (char *fname,double ***w, double **bias, int nIU, int *nHU, int noU, int nHL)
{
    int i,j,k,iL;
    FILE * fp;
    fp = fopen (fname, "r") ;
    for (iL=0;iL<=nHL; iL++) {
        for (i=0;i<nHU[iL]; i++){
            fscanf (fp, "%d" ,&k) ;
            for(j=0;j<nHU[iL+1];j++) fscanf (fp, "%le", w[iL][j]+i);
        }
    }
    for (iL=1; iL<=nHL+1; iL++) {
        for (j=0;j<nHU[iL];j++) fscanf (fp, "%le" ,bias[iL]+j);
    }
    printf("%f \n", bias[iL-1][j-1]);
    fclose (fp) ;
}

void initialize(double ***w, double **bias, int nIU, int *nHU, int nOU, int nHL) {
    int i,j,k;
    for (i=0;i<=nHL;i++)
        for (j=0;j<nHU[i+1];j++)
            for (k=0;k<nHU[i];k++) w[i][j][k] = rnd();
    for(j=1;j<=nHL+1;j++)
        for(i=0;i<nHU[j];i++) bias[j][i] = rnd();
}
void store_weight (double ***w, double **bias, double ***w_min, double **bias_min, int nIU, int *nHU, int noU, int nHL) {
    printf("Stored weight:\n");
int i,j,k;
for (i=0;i<=nHL;i++)
    for(j=0;j<nHU[i+1];j++)
        for (k=0;k<nHU[i];k++) w_min [i][j][k] =w[i][j][k];
    for(j = 1;j<=nHL+1;j++)
        for (i=0;i<nHU[j];i++) bias_min[j][i] =bias[j][i] ;
    
    printf("\nStored weight end\n");
}
/*--*/
void show_results ( double ***w, double **bias, double ***w_m, double **bias_m, int nIU, int *nHU, int nou, int nHL) {
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
    printf ("\nWeights: \n");
    for (iL=0;iL<=nHL;iL++) {
        for (i=0;i<nHU[iL];i++) {
            printf ("%5d", i);
            for (j=0;j<nHU[iL+1];j++)
                printf(" %e",w[iL][j][i]);
            printf ("\n");
        }
    }
    printf ("\nBiases: \n");
    for (iL=1; iL<=nHL+1;iL++) {
        for(j=0;j<nHU[iL];j++) printf("%e ", bias[iL][j]);
        printf ("\n");
    }
}
/*.*/
void clear_dweight ( double ***dw, double **dbias, int nIU, int *nHU, int nou, int nHL) {
    int i,j,k;
    for (i=0;i<=nHL;i++)
        for(j=0;j<nHU[i+1];j++)
            for (k=0;k<nHU[i];k++) dw[i][j][k] = 0.0 ;
    for (j=1;j<=nHL+1;j++)
        for (i=0;i<nHU[j];i++) dbias[j][i] = 0.0 ;
}
/**/
/* DLebp.c */
void propagation ( int p, double **zIU, double **zHU, double **zdHU, double *zOU, double *zdOU, double ***w, double **bias, int nIU,int *nHU, int nOU, int nHL)
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
void back_propagation ( int p,double **t,double **zIU, double **zHU, double **zdHU, double *zOU, double *zdOU, double ***w, double **bias, double ***dw, double **dbias, double **d, int nIU, int *nHU, int nOU, int nHL, double Alpha, double Beta) {
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
                dw[iL] [j] [i] = Mom3*dw[iL] [j][i];
                dw[iL][j][i] += Alpha*d[iL+1][j]*zHU[iL][i];
                w[iL][j][i] += dw[iL][j][i];
            }
        }
    }
    for (j=0; j<nHU[1];j++) {
        for (i=0;i<nIU;i++) {
            dw[0][j][i] = Alpha*d[1][j]*zIU[p][i] + Mom3*dw[0][j][i];
//            + Mom3*dw[0][j][i];
            w[0][j][i] += dw[0][j][i];
        }
    }
}

int main1 () {
    int i,j,k, iteration_min, i1, j1,rseed, o_freq, MaxPattern, MaxEpochs, lp_no, tp_no,nIU, nOU, *nHU, nHL, nHU0, nhflag,
    *idx1;
    double *zOU, **zOU_min, **zHU, **zIU, **zIUor, ***w, **bias,***dw, **dbias, ***w_min, **bias_min, **dtemp, **zdHU, *zdOU, ef1, ef2, ef2_min=16,NoiseLevel, Alpha, Beta, **t ; char fname1[FNAMELENGTH];
    FILE *fp;
    /*ー ー*/
    scanf ("%d %d %d %d %d %d %d %d %d %s %d %le %le %le", &MaxPattern, &lp_no, &nIU, &nHU0, &nOU, &nHL, &nhflag, &MaxEpochs, &o_freq, fname1, &rseed, &Alpha, &Beta,&NoiseLevel) ;
    tp_no = MaxPattern - lp_no ;
    /*.-*/
    nHU = ivector (0, nHL+1) ;
    if (nhflag == NHU_V) {
        for (i=1; i<=nHL;i++) scanf ("%d" ,nHU+i) ;
    } else {
        for (i=1;i<=nHL;i++) nHU[i]= nHU0;
    }
    nHU [0] = nIU;
    nHU [nHL+1] = nOU;

    /*-----*/
    t = matrix(0,MaxPattern-1,0, nOU-1);
    zIU = matrix (0,MaxPattern-1,0,nIU-1) ;
    zIUor = matrix(0,MaxPattern-1,0,nIU-1) ;
    zHU = (double **)malloc ((nHL+2)*sizeof(double *)) ;
    for (i=0; i<nHL+2;i++) zHU[i] = vector (0, nHU[i]-1) ;
    zdHU = (double **)malloc ((nHL+2) *sizeof (double *)) ;
    for (i=0;i<nHL+2;i++) zdHU[i] = vector (0,nHU[i]-1) ;
    zOU = vector (0, nOU-1) ;
    zdOU = vector (0, nOU-1) ;
    zOU_min = matrix (0, MaxPattern-1, 0, nOU-1) ;
    w = (double ***)malloc ((nHL+1)*sizeof(double **));
    for(i=0;i<=nHL;i++)
        w[i] = matrix(0, nHU[i+1]-1,0, nHU[i]-1);
    w_min = (double ***)malloc ( (nHL+1) *sizeof(double **));
    for(i=0;i<=nHL;i++)
        w_min[i] = matrix(0,nHU[i+1]-1,0,nHU[i]-1);
    dw = (double ***)malloc ((nHL+1) *sizeof (double **));
    for(i=0;i<=nHL;i++)
        dw[i] = matrix(0,nHU[i+1]-1,0,nHU[i]-1) ;
    bias = (double **)malloc ((nHL+2) *sizeof (double *));
    for(i=0;i<=nHL+1;i++) bias[i] = vector(0, nHU[i]-1);
    bias_min = (double **)malloc ((nHL+2) *sizeof(double *));
    for (i=0; i<=nHL+1;i++) bias_min[i] = dvector (0, nHU [i]-1);
    dbias = (double **)malloc ((nHL+2) *sizeof (double *));
    for(i=0;i<=nHL+1;i++) dbias[i] = dvector(0, nHU[i]-1);
    dtemp = (double **)malloc ((nHL+2)*sizeof(double *));
    for (i=0; i<nHL+2;i++) dtemp[i] = dvector (0, nHU [i]-1);
    
    idx1 = (int *)malloc(lp_no*sizeof (int));
    for (i=0;i<lp_no;i++) idx1[i]=i ;
    /*ーーーーー*/
    read_file(fname1, zIUor, t, nIU, nOU, MaxPattern);
    srand48 (rseed) ;
    initialize (w,bias,nIU, nHU, nOU, nHL) ;
    clear_dweight (dw, dbias, nIU, nHU, nOU, nHL) ;
    /*-----*/
    for (i=0; i<=MaxEpochs; i++) {
        for (i1=0;i1<lp_no;i1++)
            for (j=0;j<nIU;j++)
                zIU[i1][j] = (1.0 + noise()) * zIUor[i1][j];
        s_shuffle(idx1, lp_no) ;
        for (j=0;j<lp_no;j++) {
            propagation(idx1[j], zIU, zHU, zdHU, zOU, zdOU,w,bias, nIU, nHU, nOU, nHL) ;
            back_propagation (idx1[j], t, zIU, zHU, zdHU, zOU, zdOU,w, bias, dw, dbias, dtemp, nIU, nHU, nOU, nHL, Alpha, Beta);
        }
        if(i%o_freq==0) {
            for (ef1=0.0f,j=0;j<lp_no;j++) {
                propagation(j, zIUor, zHU, zdHU, zOU, zdOU, w, bias, nIU, nHU, nOU, nHL) ;
                for (k=0;k<nOU;k++) ef1 += fabs (t [j] [k] - zOU[k]) ;
            }
            for (ef2=0.0f,j=lp_no;j<lp_no+tp_no;j++) {
                propagation (j, zIUor, zHU, zdHU, zOU, zdOU, w, bias,nIU, nHU, nOU, nHL) ;
                for (k=0;k<nOU; k++) ef2 += fabs (t[j][k] - zOU[k]) ;
            }
            printf ("%d th Error : %.5f %.5f\n",i,ef1/lp_no,ef2/tp_no) ;
            if(ef2<=ef2_min) {
                ef2_min = ef2;
                iteration_min = i;
                store_weight (w, bias,w_min, bias_min,nIU, nHU, nOU, nHL) ;
            }
        }
    }

    show_results (w, bias, w_min, bias_min, nIU, nHU, nOU, nHL);
    return 0;
}

