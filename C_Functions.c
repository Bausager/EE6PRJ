/* Includes ------------------------------------------------------------------*/
#include "C_Functions.h"

double MAF(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[25000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}

double MAF_Error_1(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[25000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}

double MAF_Error_2(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[25000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}

double MAF_1(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[45000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}
double MAF_2(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[45000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}
double MAF_3(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[45000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}
double MAF_4(double sample, unsigned long len){

    float a = (float)1/(float)len;

    static float x[45000] = {0};
    static float y = 0;

    static unsigned long N = 0;


    y = y + a*sample - a*x[N];

    x[N] = sample;

    N = (N + 1) % len;

    return y;
}


double Clarke_alpha(double x_a, double x_b, double x_c){
    static float a = 2.0f/3.0f;
    static float b = -(1.0f/3.0f);
    static float c = sqrt(3.0f)/3.0f;

    return (a*x_a + b*x_b + b*x_c);
}

double Clarke_beta(double x_a, double x_b, double x_c){

    static float a = sqrt(3.0f)/3.0f;
    static float b = -1.0f*(sqrt(3.0f)/3.0f);

    return (a*x_b + b*x_c);
}

double Clarke_inv_a(double x_alpha, double x_beta){

    return (x_alpha);
}
double Clarke_inv_b(double x_alpha, double x_beta){
    static double a = -0.5f;
    static double b = (sqrt(3.0f))/2.0f;

    double y = ((a*x_alpha) + (b*x_beta));

    return y;
}
double Clarke_inv_c(double x_alpha, double x_beta){
    static double a = -0.5f;
    static double b = -1.0f*(sqrt(3.0f)/2.0f);

    double y = ((a*x_alpha) + (b*x_beta));

    return y;
}

double Park_d(double x_alpa, double x_beta, double theta){

    float a = cosf(theta);
    float b = sinf(theta);

    return(a*x_alpa + b*x_beta);
}

double Park_q(double x_alpa, double x_beta, double theta){

    float a = -sinf(theta);
    float b = cosf(theta);
    
    return(a*x_alpa + b*x_beta);
}

double Park_inv_alpha(double x_d, double x_q, double theta){

    float a = cosf(theta);
    float b = -sinf(theta);

    return(a*x_d + b*x_q);
}
double Park_inv_beta(double x_d, double x_q, double theta){

    float a = sinf(theta);
    float b = cosf(theta);

    return(a*x_d + b*x_q);
}

double Clarke_Park_d(double x_a, double x_b, double x_c, double theta){

    static float a = 2.0f/3.0f;
    static float b = 2.0f*M_PI/3.0f;

    return(a*(cosf(theta)*x_a + cosf(theta - b)*x_b + cosf(theta + b)*x_c));
}

double Clarke_Park_q(double x_a, double x_b, double x_c, double theta){

    static float a = 2.0f/3.0f;
    static float b = 2.0f*M_PI/3.0f;

    return(a*(-sinf(theta)*x_a + -sinf(theta - b)*x_b + -sinf(theta + b)*x_c));
}

double Clarke_Park_inv_a(double x_d, double x_q, double theta){
    double a = cosf(theta);
    double b = -sinf(theta);

    double y = (a*x_d) + (b*x_q);

    return y;
}
double Clarke_Park_inv_b(double x_d, double x_q, double theta){
    static double z = (2.0f*M_PI/3.0f);
    double a = cosf(theta - z);
    double b = -sinf(theta - z);

    double y = (a*x_d) + (b*x_q);

    return y;
}
double Clarke_Park_inv_c(double x_d, double x_q, double theta){
    static double z = (2.0f*M_PI/3.0f);
    double a = cosf(theta + z);
    double b = -sinf(theta + z);
    
    double y = (a*x_d) + (b*x_q);

    return y;
}
double PID_Type_1(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;

    if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    
    
    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}
double PID_Type_1_1(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;


    if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    
    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}
double PID_Type_1_2(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;

        if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    
    
    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}
double PID_Type_1_3(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;
    
    if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    

    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}
double PID_Type_1_4(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;

    if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    
    
    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}
double PID_Type_1_5(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;

    if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    
    
    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}
double PID_Type_1_6(double T_s, double T_settling, double error, double time, double I_sat, double out_sat){

    //static double K_p = 92.0f;
    double K_p = 9.2f/T_settling;
    //static double T_i = 0.000235f;
    double T_i = (1.0f/(sqrt(2.0f)) * (pow(T_settling, 2.0f))) / 21.16f;

    static double old_error = 0.0f;
    static double I_part = 0.0f;

    if (time <= 0.02)
    {
        old_error = 0.0f;
        I_part = 0.0f;
    }
    
    
    double P_part = K_p * error;
    I_part += (T_s/2.0f*T_i)*(error + old_error);

    if (I_part > I_sat){
        I_part = I_sat;
    }
    else if (I_part < -I_sat){
        I_part = -I_sat;
    }

    double PI = P_part + I_part;

    if (PI > out_sat){
        PI = out_sat;
    }
    else if (PI < -out_sat){
        PI = -out_sat;
    }

    old_error = error;

    return PI;}


double integrator(double error, double T_s){

    static float old_error = 0.0f;
    static double intetral = 0.0f;

    intetral += (T_s/2.0f) * (error + old_error);

    old_error = error;

    return intetral;}

double measure_P(double u_a, double u_b, double u_c, double i_a, double i_b, double i_c){
    double P = (u_a*i_a) + (u_b*i_b) + (u_c*i_c);
    return P;
}

double measure_Q(double u_a, double u_b, double u_c, double i_a, double i_b, double i_c){
    double Q = (1.0f/sqrt(3.0f)) * ((i_a*(u_c - u_b)) + (i_b*(u_a - u_c)) + (i_c*(u_b - u_a)));
    //double Q =  (3 * ((i_a*(u_b - u_c))) - ((u_a - u_b - u_c + u_a)*(i_b - i_c)))/(2.0f*sqrt(3.0f));
    
    return Q;
}