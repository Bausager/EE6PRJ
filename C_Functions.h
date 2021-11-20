
#ifndef INC_C_FUNCTIONS_H_
#define INC_C_FUNCTIONS_H_
#include <stdio.h>
#include <math.h>


double MAF(double sample, unsigned long len);

double MAF_Error_1(double sample, unsigned long len);
double MAF_Error_2(double sample, unsigned long len);

double MAF_1(double sample, unsigned long len);
double MAF_2(double sample, unsigned long len);
double MAF_3(double sample, unsigned long len);
double MAF_4(double sample, unsigned long len);

double Clarke_alpha(double x_a, double x_b, double x_c);
double Clarke_beta(double x_a, double x_b, double x_c);


double Clarke_inv_a(double x_alpha, double x_beta);
double Clarke_inv_b(double x_alpha, double x_beta);
double Clarke_inv_c(double x_alpha, double x_beta);

double Park_d(double x_alpa, double x_beta, double theta);
double Park_q(double x_alpa, double x_beta, double theta);

double Park_inv_alpha(double x_d, double x_q, double theta);
double Park_inv_beta(double x_d, double x_q, double theta);

double Clarke_Park_d(double x_a, double x_b, double x_c, double theta);
double Clarke_Park_q(double x_a, double x_b, double x_c, double theta);

double Clarke_Park_inv_a(double x_d, double x_q, double theta);
double Clarke_Park_inv_b(double x_d, double x_q, double theta);
double Clarke_Park_inv_c(double x_d, double x_q, double theta);

double PID_Type_1(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);

double PID_Type_1_1(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);
double PID_Type_1_2(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);
double PID_Type_1_3(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);
double PID_Type_1_4(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);
double PID_Type_1_5(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);
double PID_Type_1_6(double T_s, double T_settling, double error, double time, double I_sat, double out_sat);


double integrator(double error, double T_s);


double measure_P(double u_a, double u_b, double u_c, double i_a, double i_b, double i_c);
double measure_Q(double u_a, double u_b, double u_c, double i_a, double i_b, double i_c);

#endif 