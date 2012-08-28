
/*
 * //#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>
double I_d(double V_c){
    return V_c*V_c*V_c - 1.5*V_c*V_c + 0.6*V_c ;
}

int TunnelDiodeOscillatorDynamic(double t, const double y[], double f[], void *params){
    double C = 1e-9 ;
    double G = 5  ;
    double L = 1e-6 ;
    //double V0 = 0.8 ;
    //double I0 = 0.04e-3 ;

    double V = *(double *)params;
    f[0] = (1/C)*(-I_d(y[0])+y[1]) ;
    f[1] = (1/L)*(-y[0] - (1/G)*y[1] + V);

    return GSL_SUCCESS;
}


int VanderpolOscillatorDynamic(double t, const double y[], double f[], void *params){
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

int VanderpolOscillatorJacobian (double t, const double y[], double *dfdy, double dfdt[], void *params){
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}
*/
