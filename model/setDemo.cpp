#include "population.h"
#include "inputData.h"

double setDemoParameters (double param1, double param2) {
        /*
        
        For input values for alpha1 and mu1, the fit to field data (residualSS) is calculated and returned
        function called iteratively by estimateRemovalRates
        
        */

    double midpt;
    double M_inf;
    double M_nn;
    double ppred;
    double residual;
    double sumpred;
    double perc_inf;
    double L_inf;
    double L1;
    double IMR;
    int i;
    double valsetDemoParameters;
    rho=get_growthrate()/(100.0*(intervalsPerYear));
    IMR=0.1;
    M_inf=-log(1-IMR);
    mu1=exp(param1)/100;
    alpha1=exp(param2)/100;
    alpha0=4.0;
    mu0=(M_inf-mu1*(exp(alpha1*0.5)-1)*alpha0)/(alpha1*(1-exp(-alpha0*0.5)));
    valsetDemoParameters=0.0;
    sumpred=0.0;
    ppred=1.0;
    for ( i=1;i<=ngroups-1; i++) {
        midpt=(a1[i - 1]+a0[i - 1])*0.5;
        M1[i - 1]=mu0*(1-exp(-alpha0*midpt))/alpha0;
        M2[i - 1]=mu1*(exp(alpha1*midpt)-1)/alpha1;
        M[i - 1]=M1[i - 1]+M2[i - 1];
        pred[i - 1]=(a1[i - 1]-a0[i - 1])*exp(-rho*midpt-M[i - 1]);
        sumpred=sumpred+pred[i - 1];
    }
    for ( i=1;i<=ngroups-1; i++) {
        pred[i - 1]=pred[i - 1]/sumpred*100.0;
    }
    L_inf=exp(-rho*0.5-M[2 - 1]);
    M_nn=-log(1-0.4*(1-exp(-M[2 - 1])));
    L1=1.0/12.0*exp(-rho/24.0-M_nn);
    perc_inf=perc[1 - 1]+perc[2 - 1];
    perc[1 - 1]=perc_inf*L1/L_inf;
    perc[2 - 1]=perc_inf-perc[1 - 1];
    for ( i=1;i<=ngroups-1; i++) {
        residual=log(pred[i - 1])-log(perc[i - 1]);
        valsetDemoParameters=valsetDemoParameters+residual*residual;
    }
    return valsetDemoParameters;
}

