close all
clear all

format long;
pkg load symbolic;

%%%%%%%%%%%%%%%%%%%READ AND WRITE DATA FILES%%%%%%%%%%%%%%%%%%%

dataf = fopen('data.txt','r');
DATA = fscanf(dataf,'%*s = %f');
#DATA
fclose(dataf);

data_ngf = fopen('data_ng.txt','w');
fprintf(data_ngf,'R1 1 2 %.11fk\n', DATA(1));
fprintf(data_ngf,'R2 2 3 %.11fk\n', DATA(2));
fprintf(data_ngf,'R3 2 4 %.11fk\n', DATA(3));
fprintf(data_ngf,'R4 4 0 %.11fk\n', DATA(4));
fprintf(data_ngf,'R5 4 5 %.11fk\n', DATA(5));
fprintf(data_ngf,'R6 0 8 %.11fk\n', DATA(6));
fprintf(data_ngf,'R7 6 7 %.11fk\n', DATA(7));
fprintf(data_ngf,'C 5 7 %.11fu\n', DATA(9));
fprintf(data_ngf,'Vs 1 0 %.11f\n', DATA(8));
fprintf(data_ngf,'Vaux 8 6 DC 0\n');
fprintf(data_ngf,'Hd 4 7 vaux %.11fk\n', DATA(11));
fprintf(data_ngf,'Gb 5 3 (2,4) %.11fm\n', DATA(10));
fclose(data_ngf);

%%%%%%%%%%%%%%%%%%%SYMBOLIC COMPUTATIONS%%%%%%%%%%%%%%%%%%%

R1 = sym (sprintf('%.11f',DATA(1)));
R2 = sym (sprintf('%.11f',DATA(2)));
R3 = sym (sprintf('%.11f',DATA(3)));
R4 = sym (sprintf('%.11f',DATA(4)));
R5 = sym (sprintf('%.11f',DATA(5)));
R6 = sym (sprintf('%.11f',DATA(6)));
R7 = sym (sprintf('%.11f',DATA(7)));
Vs = sym (sprintf('%.11f',DATA(8)));
C = sym (sprintf('%.11f',DATA(9)));
Kb = sym (sprintf('%.11f',DATA(10)));
Kd = sym (sprintf('%.11f',DATA(11)));

%%%%%Nodal method 1)%%%%%

syms V1n V2n V3n V4n V5n V6n V7n 
Eqd_2 = (V2n-V1n)/R1 + (V2n-V4n)/R3 + (V2n-V3n)/R2 == 0;
Eqd_0_R6 = (V1n-V2n)/R1 + (-V4n)/R4 + (-V6n)/R6 == 0;
Eqd_0_R7 = (V1n-V2n)/R1 + (-V4n)/R4 + (V6n-V7n)/R7 == 0;
Eqd_Vs_1 = V1n == Vs;
Eqd_Vd = V4n-V7n == Kd*(-V6n)/R6;
Eqd_5_1 = (V4n-V5n)/R5 == Kb*(V2n-V4n);
Eqd_3 = (V3n-V2n)/R2 == Kb*(V2n-V4n);
sn_1 = solve(Eqd_2,Eqd_0_R6,Eqd_0_R7,Eqd_Vs_1,Eqd_Vd,Eqd_5_1,Eqd_3);

%%%%%Nodal method 1) NUMERIC%%%%%

V1n = double(sn_1.V1n);
V2n = double(sn_1.V2n);
V3n = double(sn_1.V3n);
V4n = double(sn_1.V4n);
V5n = double(sn_1.V5n);
V6n = double(sn_1.V6n);
V7n = double(sn_1.V7n);
V8n = V6n;
IR1n = (V1n-V2n)/double(R1);
IR2n = (V2n-V3n)/double(R2);
IR3n = (V2n-V4n)/double(R3);
IR4n = (V4n)/double(R4);
IR5n = (V4n-V5n)/double(R5);
IR6n = (-V6n)/double(R6);
IR7n = (V6n-V7n)/double(R7);
Ibn = double(Kb)*(V2n-V4n);
Icn = 0

%%%%%Nodal method 1) TABLE FILE%%%%%

diary "Nodal1_tab.tex"
diary on
printf("$I_b$ & %e\n", Ibn/1000);
printf("$I_c$ & %d\n", Icn);
printf("$I_R$$_1$ & %e\n",IR1n/1000);
printf("$I_R$$_2$ & %e\n",IR2n/1000);
printf("$I_R$$_3$ & %e\n",IR3n/1000);
printf("$I_R$$_4$ & %e\n",IR4n/1000);
printf("$I_R$$_5$ & %e\n",IR5n/1000);
printf("$I_R$$_6$ & %e\n",IR6n/1000);
printf("$I_R$$_7$ & %e\n",IR7n/1000);
printf("$V_1$ & %f\n",V1n);
printf("$V_2$ & %f\n",V2n);
printf("$V_3$ & %f\n",V3n);
printf("$V_4$ & %f\n",V4n);
printf("$V_5$ & %f\n",V5n);
printf("$V_6$ & %f\n",V6n);
printf("$V_7$ & %f\n",V7n);
printf("$V_8$ & %f\n",V8n);
diary off

%%%%%Nodal method 2)%%%%%

Vs = sym ('0');
Vx = sym (sprintf('%.11f',V5n-V7n));
syms V1n V2n V3n V4n V5n V6n V7n 

Eqd_Vx = V5n-V7n == Vx;
Eqd_Vs_2 = V1n == Vs;
sn_2 = solve(Eqd_2,Eqd_0_R6,Eqd_0_R7,Eqd_Vs_2,Eqd_Vd,Eqd_Vx,Eqd_3);

%%%%%Nodal method 2) NUMERIC%%%%%

V1n = double(sn_2.V1n);
V2n = double(sn_2.V2n);
V3n = double(sn_2.V3n);
V4n = double(sn_2.V4n);
V5n = double(sn_2.V5n);
V6n = double(sn_2.V6n);
V7n = double(sn_2.V7n);
V8n = V6n;
IR1n = (V1n-V2n)/double(R1);
IR2n = (V2n-V3n)/double(R2);
IR3n = (V2n-V4n)/double(R3);
IR4n = (V4n)/double(R4);
IR5n = (V4n-V5n)/double(R5);
IR6n = (-V6n)/double(R6);
IR7n = (V6n-V7n)/double(R7);
Ibn = double(Kb)*(V2n-V4n);

Ix = IR5n - Ibn;
Req = double(Vx)/Ix;
tau = double(C) * 10^-6 * Req*1000;

%%%%%Nodal method 2) TABLE FILE%%%%%

diary "Req_tau_tab.tex"
diary on
printf("$R_e$$_q$ & %f\n",Req);
printf("$tau$ & %e\n",tau);
diary off

%%%%%NATURAL SOLUTION 3)%%%%%

syms vc_n(t)
syms A s

vc_n(t) = A*exp(s*t);

%%%%%NATURAL SOLUTION 3) NUMERIC%%%%%

A = double(Vx); %initial condition V6(0) = Vx
s = 1/tau; %pois Req e negativo logo tau e negativo

t=0:1e-6:20e-3;
vc_n = A*exp(s*t);
v6_n = vc_n;

%%%%%NATURAL SOLUTION 3) PLOT%%%%%

hf = figure (1);
plot (t*1000, v6_n, "g");
xlabel ("t [ms]");
ylabel ("v6_n [V]");
legend('v6_n(t)','Location','northeast');
print (hf, "v6_n.eps", "-depsc");

%%%%%FORCED SOLUTION 4)%%%%%

Vs_p = sym ('1');
f = 1000
w = 2*pi*f



















