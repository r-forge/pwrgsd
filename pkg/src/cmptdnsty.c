#include<R.h>
#include<Rmath.h>

#define normut 8.20953615160139
#define logrt2pi 0.9189385332046727
#define dnormGI(x) exp(-x*x-logrt2pi)

double ratiodnorm(double x,double y);
void printmat(double *,int,int,char *);

void cmptdnsty(int *dofu, int *pnlook, double *pfrac_k, double *pfrac_kp1, 
	       double *mu_o, double *mu_n, double *Psiab_o, double *Psiab_n, 
	       double *Psiminfa_o, double *Psiminfa_n, double *Psibinf_o, 
	       double *Psibinf_n, double *gqx, double *gqw, int *pngqnodes, 
	       double *bold, double *bnew, double *store, double *xgrid, 
	       double *dxgrid, double *dens)
{
  double eps, sqrf_o, sqrf_n, sqrdf, a_o, a_n, b_o, b_n, Phia_o, Phia_n;
  double Phib_o, Phib_n, sab_i, xab_i, sab_j, xab_j, dsab_j, sminfa_i;
  double xminfa_i, dsminfa_i, sbinf_i, xbinf_i;
  double dsbinf_i, psiab, psiminfa, psibinf, prob0, prob1, del;
  int nlook,i,j;
  int ngqnodes;
  int zero=0, one=1;

  eps = pnorm5(-normut, 0.0, 1.0, one, zero);

  ngqnodes=*pngqnodes;
  nlook = *pnlook;

  sqrf_o = pow(*pfrac_k,0.5);
  sqrf_n = pow(*pfrac_kp1,0.5);
  sqrdf = pow(*pfrac_kp1 - *pfrac_k,0.5);

  if(*dofu==1){
    a_o = *(bold+1);
    Phia_o = pnorm5(sqrf_o * a_o - (*mu_o),0.0,1.0,one,zero);
    a_n = *(bnew+1);
    Phia_n = pnorm5(sqrf_n * a_n - (*mu_n),0.0,1.0,one,zero);
  }
  else {
    Phia_o = eps;
    Phia_n = eps;
  }

  b_o = *bold;
  Phib_o = pnorm5(sqrf_o * b_o - (*mu_o),0.0,1.0,one,zero);
  b_n = *bnew;
  Phib_n = pnorm5(sqrf_n * b_n - (*mu_n),0.0,1.0,one,zero);

  for(i=0;i<ngqnodes;i++){
    if(nlook==1) {
      sab_i = (1.0 - *(gqx+i))*Phia_o/2.0 + (1.0 + *(gqx+i))*Phib_o/2.0;
      xab_i = qnorm5(sab_i, 0.0, 1.0, one, zero);
      if(*dofu==1){
	sminfa_i = (1.0+ *(gqx+i))*Phia_o/2.0;
	dsminfa_i = *(gqw+i)*Phia_o/2.0;
	xminfa_i = qnorm5(sminfa_i, 0.0, 1.0, one, zero);
      }
      sbinf_i = (1.0 - *(gqx+i))*Phib_o/2.0 + (1.0 + *(gqx+i))/2.0;
      dsbinf_i = *(gqw + i) * (1.0 - Phib_o)/2.0;
      xbinf_i = qnorm5(sbinf_i, 0.0, 1.0, one, zero);

      *(Psiab_o + i) = exp(-xab_i*xab_i/(2.0*(*pfrac_k)) - logrt2pi - log(sqrf_o));
      *(Psiminfa_o + i) = exp(-xminfa_i*xminfa_i/(2.0*(*pfrac_k)) - logrt2pi - log(sqrf_o));
      *(Psibinf_o + i) = exp(-xbinf_i*xbinf_i/(2.0*(*pfrac_k)) - logrt2pi - log(sqrf_o));
      *store = 1.0;
    }

    if(nlook>1){
      *(Psiab_o + i) = *(Psiab_n + i);
      *(Psibinf_o + i) = *(Psibinf_n + i);
      if(*dofu==1){
	*(Psiminfa_o + i) = *(Psiminfa_n + i);
	sminfa_i = (1.0+ *(gqx+i))*Phia_o/2.0;
	dsminfa_i = *(gqw+i)*Phia_o/2.0;
	xminfa_i = qnorm5(sminfa_i, 0.0, 1.0, one, zero);
      }
      sbinf_i = (1.0 - *(gqx+i))*Phib_o/2.0 + (1.0 + *(gqx+i))/2.0;
      dsbinf_i = *(gqw + i) * (1.0 - Phib_o)/2.0;
      xbinf_i = qnorm5(sbinf_i, 0.0, 1.0, one, zero);
    }

    *(xgrid + i) = xminfa_i;
    *(xgrid + ngqnodes + i) = xbinf_i;
    *(dxgrid + i) = dsminfa_i/2.0;
    *(dxgrid + ngqnodes + i) = dsbinf_i/2.0;
    *(dens + i) = exp(log(*(Psiminfa_o + i)) + xminfa_i*xminfa_i/2.0 + logrt2pi-log(*store));
    *(dens + ngqnodes + i) = exp(log(*(Psibinf_o + i)) + xbinf_i*xbinf_i/2.0 + logrt2pi-log(*store));
  }

  for(i=0;i<ngqnodes;i++) {
    sab_i = (1.0 - *(gqx+i))*Phia_n/2.0 + (1.0 + *(gqx+i))*Phib_n/2.0;
    xab_i = qnorm5(sab_i, 0.0, 1.0, one, zero);
    if(*dofu==1){
      sminfa_i = (1.0+ *(gqx+i))*Phia_n/2.0;
      xminfa_i = qnorm5(sminfa_i, 0.0, 1.0, one, zero);
    }
    sbinf_i = (1.0 - *(gqx+i))*Phib_n/2.0 + (1.0 + *(gqx+i))/2.0;
    xbinf_i = qnorm5(sbinf_i, 0.0, 1.0, one, zero);
    psiab = 0.0;
    psiminfa = 0.0;
    psibinf = 0.0;
    for(j=0;j<ngqnodes;j++){
      sab_j = (1.0 - *(gqx+j))*Phia_o/2.0 + (1.0 + *(gqx+j))*Phib_o/2.0;
      xab_j = qnorm5(sab_j, 0.0, 1.0, one, zero);
      dsab_j = *(gqw + j) * (Phib_o - Phia_o)/2.0;
      del = (xab_i-xab_j)/sqrdf;
      psiab +=exp(log(*(Psiab_o+j)) - del*del/2.0 + xab_j*xab_j/2.0 - log(sqrdf))*dsab_j;
      if(*dofu==1){
	del = (xminfa_i-xab_j)/sqrdf;
	psiminfa+=exp(log(*(Psiab_o+j))-del*del/2.0+xab_j*xab_j/2.0 - log(sqrdf))*dsab_j;
      }
      del = (xbinf_i-xab_j)/sqrdf;
      psibinf+=exp(log(*(Psiab_o+j))-del*del/2.0+xab_j*xab_j/2.0 - log(sqrdf))*dsab_j;
    }
    *(Psiminfa_n + i) = psiminfa;
    *(Psiab_n + i) = psiab;
    *(Psibinf_n + i) = psibinf;
    *store = 1.0;
      //sqrdf;
  }
}
