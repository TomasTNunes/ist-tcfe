close all
clear all

%%%%%%%%%%%%%%%%%%%SYMBOLIC COMPUTATIONS%%%%%%%%%%%%%%%%%%%

pkg load symbolic;

R1 = sym ('1.02943118797');
R2 = sym ('2.01395929215');
R3 = sym ('3.04865352258');
R4 = sym ('4.00813564187');
R5 = sym ('3.14877356154');
R6 = sym ('2.0324721756') ;
R7 = sym ('1.02032682556');
Va = sym ('5.05481864136');
Id = sym ('1.02872173547');
Kb = sym ('7.12052712169');
Kc = sym ('8.12923323408');
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G5 = 1/R5;
G6 = 1/R6;
G7 = 1/R7;


%%%%%Mesh Method%%%%%

syms I_A I_B I_C I_D;
Eqd_A = R1*I_A + R3*(I_A+I_B) + R4*(I_A+I_C) == Va;
Eqd_C = R4*(I_A+I_C)+R6*I_C+R7*I_C == Kc*I_C;
Eqd_B = I_B == Kb*R3*(I_A+I_B);
Eqd_D = I_D == Id;
sm = solve(Eqd_A,Eqd_B,Eqd_C,Eqd_D);


%%%%%Nodal method%%%%%

syms V1n V2n V3n V4n V5n V6n V7n 
Eqd_2 = (V2n-V4n)*G3 + (V2n-V3n)*G2 + (V2n-V1n)*G1 == 0;
Eqd_0_G6 = (V1n-V2n)*G1 - V6n*G6 - V4n*G4 == 0; 
Eqd_5 = (V5n-V4n)*G5 - Id + Kb*(V2n-V4n) == 0;
Eqd_0_G7 = -V4n*G4 + (V1n-V2n)*G1 + (V6n-V7n)*G7 == 0;
Eqd_Va = V1n == Va;
Eqd_Vc = V4n - V7n == -Kc*V6n*G6;
Eqd_3 = Kb*(V2n-V4n) == (V3n-V2n)*G2;
sn = solve(Eqd_2,Eqd_0_G6,Eqd_5,Eqd_0_G7,Eqd_Va,Eqd_Vc,Eqd_3);


%%%%%%%%%%%%%%%%%%%NUMERIC COMPUTATIONS%%%%%%%%%%%%%%%%%%%
R1 = 1.02943118797;
R2 = 2.01395929215;
R3 = 3.04865352258;
R4 = 4.00813564187;
R5 = 3.14877356154;
R6 = 2.0324721756;
R7 = 1.02032682556;
Va = 5.05481864136;
Id = 1.02872173547;
Kb = 7.12052712169;
Kc = 8.12923323408;
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G5 = 1/R5;
G6 = 1/R6;
G7 = 1/R7;


%%%%%Mesh Method%%%%%

I_A = double(sm.I_A)
I_B = double(sm.I_B)
I_C = double(sm.I_C)
I_D = double(sm.I_D) 
Ibm = I_B;
Idm = Id;
IR1m = I_A; 
IR2m = -I_B;
IR3m = I_A + I_B;
IR4m = I_A + I_C ;
IR5m = I_B - I_D ;
IR6m = I_C ;
IR7m = I_C ;
V1m = Va;
V2m = V1m-R1*I_A;
V3m = V2m+R2*I_B;
V4m = R4*(I_A+I_C);
V5m = V4m-R5*(I_B-I_D);
V6m = -R6*I_C;
V7m = V6m-R7*I_C;
V8m = V6m;

%%%%%Nodal method%%%%%

V1n = double(sn.V1n);
V2n = double(sn.V2n);
V3n = double(sn.V3n);
V4n = double(sn.V4n);
V5n = double(sn.V5n);
V6n = double(sn.V6n);
V7n = double(sn.V7n);
V8n = V6n;
Ibn = (V2n-V4n)*Kb;
Idn = Id;
IR1n = (V1n-V2n)*G1;
IR2n = (V2n-V3n)*G2;
IR3n = (V2n-V4n)*G3;
IR4n = (V4n)*G4;
IR5n = (V4n-V5n)*G5;
IR6n = (-V6n)*G6;
IR7n = (V6n-V7n)*G7;


%%%%%%%%%%%%%%%%%%%TEXT WITH TABLE DATA%%%%%%%%%%%%%%%%%%%

diary "Mesh_tab.tex"
diary on
printf("$I_b$ & %f\n", Ibm);
printf("$I_d$ & %f\n", Idm);
printf("$I_R$$_1$ & %f\n", IR1m);
printf("$I_R$$_2$ & %f\n", IR2m);
printf("$I_R$$_3$ & %f\n", IR3m);
printf("$I_R$$_4$ & %f\n", IR4m);
printf("$I_R$$_5$ & %f\n", IR5m);
printf("$I_R$$_6$ & %f\n", IR6m);
printf("$I_R$$_7$ & %f\n", IR7m);
printf("$V_1$ & %f\n", V1m);
printf("$V_2$ & %f\n", V2m);
printf("$V_3$ & %f\n", V3m);
printf("$V_4$ & %f\n", V4m);
printf("$V_5$ & %f\n", V5m);
printf("$V_6$ & %f\n", V6m);
printf("$V_7$ & %f\n", V7m);
printf("$V_8$ & %f\n", V8m);
diary off


diary "Nodal_tab.tex"
diary on
printf("$I_b$ & %f\n", Ibn);
printf("$I_d$ & %f\n", Idn);
printf("$I_R$$_1$ & %f\n", IR1n);
printf("$I_R$$_2$ & %f\n", IR2n);
printf("$I_R$$_3$ & %f\n", IR3n);
printf("$I_R$$_4$ & %f\n", IR4n);
printf("$I_R$$_5$ & %f\n", IR5n);
printf("$I_R$$_6$ & %f\n", IR6n);
printf("$I_R$$_7$ & %f\n", IR7n);
printf("$V_1$ & %f\n", V1n);
printf("$V_2$ & %f\n", V2n);
printf("$V_3$ & %f\n", V3n);
printf("$V_4$ & %f\n", V4n);
printf("$V_5$ & %f\n", V5n);
printf("$V_6$ & %f\n", V6n);
printf("$V_7$ & %f\n", V7n);
printf("$V_8$ & %f\n", V8n);
diary off

