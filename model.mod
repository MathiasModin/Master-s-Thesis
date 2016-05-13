%//////////////////////////////////////////////////////////////////////////
%// Estimation of DSGE model based on Smets and Wouters (2003) (modified)
%// Master's Thesis: Mathias Concha Modin, 2016-05-12
%//////////////////////////////////////////////////////////////////////////

// stable: addpath c:\dynare\4.4.3\matlab

%//////////////////////////////////////////////////////////////////////////
%// Endogenous Variables
%//////////////////////////////////////////////////////////////////////////

%// Endogenous variables
var Y C W L R K I Q PIE r Yf Kf Lf Wf Cf If Qf rf Rf gap

%// Shock processes
eps_B %// Preference shock
eps_L %// Labour supply shock
eps_I %// Investment/cost adjustment shock
eps_a %// Productivity shock
eps_R %// Inflation objective MP shock
eps_G; %// Government spending shock


%//////////////////////////////////////////////////////////////////////////
%// Exogenous Variables
%//////////////////////////////////////////////////////////////////////////

varexo
eta_Q %// Non-structural shock to external finance premium
eta_P %// Constant in time varying mark-up process 
eta_w %// Constant in shock to wage mark-up process
eta_R %// Monetary policy shock
eta_B %// Preference shock
eta_L %// Labour supply shock
eta_I %// Investment shock
eta_a %// Productivity shock
eta_G %// Government spendings shock
;
  

%//////////////////////////////////////////////////////////////////////////
%// Parameters and calibration
%//////////////////////////////////////////////////////////////////////////

parameters 
    cbeta sig_c sig_l h cdelta gamma_p phi_y phi_i theta_p calpha gamma_w 
    theta_w rho_L rho_G rho_B rho_I rho_a rho_R r_pi r_dy r_y crho lambda_w 
    cpsi ccs cinvs;

calpha  = 0.3;
cbeta   = 0.99;
cdelta  = 0.025;
lambda_w= 0.5;
phi_y   = 1.408;
sig_c   = 1.353; 
h       = 0.573; 
sig_l   = 2.400;
gamma_w = 0.763;
gamma_p = 0.469;
theta_w = 0.7337;
theta_p = 0.908;
r_y     = 0.099;
r_dy    = 0.159;
crho    = 0.961;
phi_i   = 1/6.771; %// S=6.711 
r_pi    = 1.684;
rho_a   = 0.823;
rho_B   = 0.855;
rho_L   = 0.889;
rho_I   = 0.927;
rho_R   = 0.924;
cpsi    = (1/0.169); 
rho_G   = 0.949;
ccs     = 0.6;
cinvs   = 0.22;



%//////////////////////////////////////////////////////////////////////////
%// Model Block
%//////////////////////////////////////////////////////////////////////////
model(linear);

%// Consumption equation
C = h/(1+h)*C(-1) + 1/(1+h)*C(+1) - ((1-h)/((1+h)*sig_c))*(R-PIE(+1)) + 
    ((1-h)/((1+h)*sig_c))*(eps_B-eps_B(+1));

%// Investment equation
I = 1/(1+cbeta)*I(-1) + cbeta/(1+cbeta)*I(+1) + phi_i/(1+cbeta)*Q + eps_I;

%// Q-equation
Q = -(R-PIE(+1)) + (1-cdelta)*cbeta*Q(+1) + (1-cbeta+cdelta*cbeta)*r(+1) + eta_Q;

%// Capital accumulation equation
K = (1-cdelta)*K(-1)+cdelta*I(-1);

%// Inflation equation / Phillips curve
PIE=cbeta*PIE(+1)/(1+cbeta*gamma_p) + gamma_p*PIE(-1)/(1+cbeta*gamma_p) + 
    (1/(1+cbeta*gamma_p))*((1-cbeta*theta_p)*(1-theta_p))/(theta_p)*
    (calpha*r+(1-calpha)*W-eps_a+eta_P);

%// Wage equation
W = cbeta*W(+1)/(1+cbeta) + W(-1)/(1+cbeta) + cbeta*PIE(+1)/(1+cbeta) - 
    (1+cbeta*gamma_w)*PIE(+1) + gamma_w*PIE(-1)/(1+cbeta) - 1/(1+cbeta)*
    ((1-cbeta*theta_w)*(1-theta_w))/(((1+(1+lambda_w)*sig_l)/lambda_w)*theta_w)
    *(W-sig_l*L-sig_c/(1-h)*(C-h*C(-1))-eps_L-eta_w);

%// Labour demand equation
L = -W + (1+cpsi)*r + K(-1);

%// Goods marked equilibrium condition
Y = phi_y*eps_a + phi_y*calpha*K(-1) + phi_y*calpha*cpsi*r + phi_y*(1-calpha)*L;

%// Gov equation
Y = (ccs*C+cinvs*I)+ eps_G;

%// Monetary policy reaction function
R = crho*R(-1) + (1-crho)*(r_pi*PIE+r_y*(Y-Yf))+r_dy*((Y-Yf)-(Y(-1)-Yf(-1)))+eps_R;



%// Frictionless model
Yf = phi_y*eps_a + phi_y*calpha*Kf(-1) + phi_y*calpha*cpsi*rf + phi_y*(1-calpha)*Lf;
Yf = (ccs*Cf+cinvs*If)+ eps_G;
Lf = -Wf + (1+cpsi)*rf + Kf(-1);
Wf = (sig_c/(1-h))*(Cf-h*Cf(-1)) + sig_l*Lf - eps_L;
Cf = h/(1+h)*Cf(-1) + 1/(1+h)*Cf(+1) - ((1-h)/((1+h)*sig_c))*Rf + 
    ((1-h)/((1+h)*sig_c))*(eps_B-eps_B(+1));
Kf = (1-cdelta)*Kf(-1)+cdelta*If(-1);
If = 1/(1+cbeta)*If(-1) + cbeta/(1+cbeta)*If(+1) + phi_i/(1+cbeta)*Qf + eps_I;
Qf = -Rf + (1-cdelta)*cbeta*Qf(+1) + (1-cbeta+cdelta*cbeta)*rf(+1) + eta_Q;
0 =  calpha*rf+(1-calpha)*Wf -eps_a ;
    


%// Shock Processes
  
%// Preference shock process
eps_B = rho_B*eps_B(-1) + eta_B;
%// Labour supply shock process
eps_L = rho_L*eps_L(-1) + eta_L;
%// Investment shock process
eps_I = rho_I*eps_I(-1) + eta_I;
%// Productivity shock
eps_a = rho_a*eps_a(-1) + eta_a;
%// Monetary Policy Shock
eps_R = rho_R*eps_R(-1) + eta_R; 
%// Government Spendings shock
eps_G = rho_G*eps_G(-1) + eta_G ;



    
%// Auxiliary output gap
gap = Y - Yf;  
end;    


%//////////////////////////////////////////////////////////////////////////
%// Shocks
%//////////////////////////////////////////////////////////////////////////
shocks;
var eta_a; stderr 0.598;
var eta_B; stderr 0.336;
var eta_I; stderr 0.085;
var eta_L; stderr 3.520;
var eta_P; stderr 0.160;
var eta_w; stderr 0.289;
var eta_R; stderr 0.081;
var eta_Q; stderr 0.604;
var eta_G; stderr 0.325;
end;

%//////////////////////////////////////////////////////////////////////////
steady;
check;


%//////////////////////////////////////////////////////////////////////////
%// Priors
%//////////////////////////////////////////////////////////////////////////
estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF

%// shocks
stderr eta_a,0.543,0.01,4,  INV_GAMMA_PDF,0.4,2;
stderr eta_B,0.2694,0.01,4, INV_GAMMA_PDF,0.2,2;
stderr eta_L,1.4575,0.1,6,  INV_GAMMA_PDF,1,2;
stderr eta_I,0.1318,0.01,4, INV_GAMMA_PDF,0.1,2;
stderr eta_R,0.1363,0.01,4, INV_GAMMA_PDF,0.1,2;
stderr eta_Q,0.4842,0.01,4, INV_GAMMA_PDF,0.4,2;
stderr eta_P,0.1731,0.01,4, INV_GAMMA_PDF,0.15,2;
stderr eta_w,0.2462,0.1,4,  INV_GAMMA_PDF,0.25,2;
stderr eta_G,0.3052,0.01,4, INV_GAMMA_PDF,0.3,2;


%// Persistence parameters
rho_a,.9722,.1,.9999,       BETA_PDF,0.85,0.1;
rho_R,.85,.1,.999,          BETA_PDF,0.85,0.1;
rho_B,.7647,.1,.99,         BETA_PDF,0.85,0.1;
rho_G,.9502,.1,.9999,       BETA_PDF,0.85,0.1;
rho_L,.9542,.1,.9999,       BETA_PDF,0.85,0.1;
rho_I,.6705,.1,.99,         BETA_PDF,0.85,0.1;

%// Deep parameters
sig_c,0.9817,0.25,3,        NORMAL_PDF,1.5,1;
sig_l,1.7526,0.5,5,         NORMAL_PDF,2,2.25;
h,0.5612,0.3,0.95,          BETA_PDF,0.65,0.1;
theta_w,0.7661,0.3,0.9,     BETA_PDF,0.75,0.05;
theta_p,0.8684,0.3,0.95,    BETA_PDF,0.75,0.05;
gamma_w,0.6202,0.1,0.99,    BETA_PDF,0.5,0.15;
gamma_p,0.6638,0.1,0.99,    BETA_PDF,0.5,0.15;
phi_i,0.1477,0,2,           NORMAL_PDF,0.1477,0.03;
cpsi,                       NORMAL_PDF,5,2;
phi_y,1.3011,1.001,2,       NORMAL_PDF,1.45,0.75;
r_pi,1.4616,1.2,2,          NORMAL_PDF,1.7,0.15;
crho,0.8865,0.5,0.99,       BETA_PDF,0.8,0.10;
r_y,0.0571,0.01,0.2,        NORMAL_PDF,0.125,0.05;
r_dy,0.2228,0.05,0.5,       NORMAL_PDF,0.0625,0.05;
end;


%// Obeserved variables
varobs Y C I L PIE W R;


%//////////////////////////////////////////////////////////////////////////
%// Estimation Options
%//////////////////////////////////////////////////////////////////////////

estimation(datafile=modeldata, lik_init=2, mh_nblocks=6, 
mh_jscale=0.1, mode_compute=6, mh_replic=1000000, mh_drop=0.5, plot_priors=0,
bayesian_irf, irf=40, moments_varendo, 
conditional_variance_decomposition=[1 2 4 8 12 16 20 40]) Y C W L R K I PIE r;
shock_decomposition;


