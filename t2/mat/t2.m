close all
clear all

format long;
pkg load symbolic;

%%%%%%%%%%%%%%%%%%%READ AND WRITE DATA FILE DATA_1_NG.TXT%%%%%%%%%%%%%%%%%%%

dataf = fopen('data.txt','r');
DATA = fscanf(dataf,'%*s = %f');
fclose(dataf);

data_ngf1 = fopen('data_1_ng.txt','w');
fprintf(data_ngf1,'R1 1 2 %.11fk\n', DATA(1));
fprintf(data_ngf1,'R2 2 3 %.11fk\n', DATA(2));
fprintf(data_ngf1,'R3 2 4 %.11fk\n', DATA(3));
fprintf(data_ngf1,'R4 4 0 %.11fk\n', DATA(4));
fprintf(data_ngf1,'R5 4 5 %.11fk\n', DATA(5));
fprintf(data_ngf1,'R6 0 8 %.11fk\n', DATA(6));
fprintf(data_ngf1,'R7 6 7 %.11fk\n', DATA(7));
fprintf(data_ngf1,'C 5 7 %.11fu\n', DATA(9));
fprintf(data_ngf1,'Vs 1 0 DC %.11f\n', DATA(8));
fprintf(data_ngf1,'Vaux 8 6 DC 0\n');
fprintf(data_ngf1,'Hd 4 7 vaux %.11fk\n', DATA(11));
fprintf(data_ngf1,'Gb 5 3 (2,4) %.11fm\n', DATA(10));
fclose(data_ngf1);

%%%%%%%%%%%%%%%%%%%SYMBOLIC COMPUTATIONS%%%%%%%%%%%%%%%%%%%

R1 = sym (sprintf('%.11f',DATA(1)*1000));
R2 = sym (sprintf('%.11f',DATA(2)*1000));
R3 = sym (sprintf('%.11f',DATA(3)*1000));
R4 = sym (sprintf('%.11f',DATA(4)*1000));
R5 = sym (sprintf('%.11f',DATA(5)*1000));
R6 = sym (sprintf('%.11f',DATA(6)*1000));
R7 = sym (sprintf('%.11f',DATA(7)*1000));
Vs = sym (sprintf('%.11f',DATA(8)));
C = sym (sprintf('%.11f',DATA(9)*10^-6));
Kb = sym (sprintf('%.11f',DATA(10)/1000));
Kd = sym (sprintf('%.11f',DATA(11)*1000));

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
Icn = 0;

%%%%%Nodal method 1) TABLE FILE%%%%%

diary "Nodal1_tab.tex"
diary on
printf("$I_c$ & %d\n", Icn);
printf("$I_b$ & %.5e\n", Ibn);
printf("$I_R$$_1$ & %e\n",IR1n);
printf("$I_R$$_2$ & %e\n",IR2n);
printf("$I_R$$_3$ & %.5e\n",IR3n);
printf("$I_R$$_4$ & %e\n",IR4n);
printf("$I_R$$_5$ & %.5e\n",IR5n);
printf("$I_R$$_6$ & %e\n",IR6n);
printf("$I_R$$_7$ & %e\n",IR7n);
printf("$V_1$ & %f\n",V1n);
printf("$V_2$ & %f\n",V2n);
printf("$V_3$ & %f\n",V3n);
printf("$V_4$ & %f\n",V4n);
printf("$V_5$ & %f\n",V5n);
printf("$V_6$ & %.5f\n",V6n);
printf("$V_7$ & %.5f\n",V7n);
printf("$V_8$ & %.5f\n",V8n);
diary off

%%%%%Nodal method 2)%%%%%

Vs = sym ('0');
Vx = sym (sprintf('%.11f',V5n-V7n));
syms V1n V2n V3n V4n V5n V6n V7n 

Eqd_Vx = V5n-V7n == Vx;
Eqd_Vs_2 = V1n == Vs;
sn_2 = solve(Eqd_2,Eqd_0_R6,Eqd_0_R7,Eqd_Vs_2,Eqd_Vd,Eqd_Vx,Eqd_3);

%%%%%Nodal method 2) WRITE FILE DATA_2_NG.TXT%%%%%

data_ngf2 = fopen('data_2_ng.txt','w');
fprintf(data_ngf2,'R1 1 2 %.11fk\n', DATA(1));
fprintf(data_ngf2,'R2 2 3 %.11fk\n', DATA(2));
fprintf(data_ngf2,'R3 2 4 %.11fk\n', DATA(3));
fprintf(data_ngf2,'R4 4 0 %.11fk\n', DATA(4));
fprintf(data_ngf2,'R5 4 5 %.11fk\n', DATA(5));
fprintf(data_ngf2,'R6 0 8 %.11fk\n', DATA(6));
fprintf(data_ngf2,'R7 6 7 %.11fk\n', DATA(7));
fprintf(data_ngf2,'Vx 5 7 DC %.6f\n',double(Vx));
fprintf(data_ngf2,'Vs 1 0 DC 0\n');
fprintf(data_ngf2,'Vaux 8 6 DC 0\n');
fprintf(data_ngf2,'Hd 4 7 vaux %.11fk\n', DATA(11));
fprintf(data_ngf2,'Gb 5 3 (2,4) %.11fm\n', DATA(10));
fclose(data_ngf2); 

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
tau = double(C) * Req;

%%%%%Nodal method 2) TABLE FILE%%%%%

diary "Req_tau_tab.tex"
diary on
printf("$I_b$ & %f\n", Ibn);
printf("$I_R$$_1$ & %f\n",IR1n);
printf("$I_R$$_2$ & %f\n",IR2n);
printf("$I_R$$_3$ & %f\n",IR3n);
printf("$I_R$$_4$ & %f\n",IR4n);
printf("$I_R$$_5$ & %.5e\n",IR5n);
printf("$I_R$$_6$ & %f\n",IR6n);
printf("$I_R$$_7$ & %f\n",IR7n);
printf("$V_1$ & %f\n",V1n);
printf("$V_2$ & %f\n",V2n);
printf("$V_3$ & %f\n",V3n);
printf("$V_4$ & %f\n",V4n);
printf("$V_5$ & %f\n",V5n);
printf("$V_6$ & %f\n",V6n);
printf("$V_7$ & %f\n",V7n);
printf("$V_8$ & %f\n",V8n);
printf("$V_X$ & %f\n",double(Vx));
printf("$I_X$ & %.5e\n",Ix);
printf("$R_e$$_q$ & %.6e\n",-Req);
printf("$tau$ & %e\n",-tau);
diary off

%%%%%NATURAL SOLUTION 3)%%%%%

syms vc_n(t)
syms A s

vc_n(t) = A*exp(s*t);

%%%%%NATURAL SOLUTION 3) WRITE FILE DATA_3_NG.TXT%%%%%
data_ngf3 = fopen('data_3_ng.txt','w');
fprintf(data_ngf3,'R1 1 2 %.11fk\n', DATA(1));
fprintf(data_ngf3,'R2 2 3 %.11fk\n', DATA(2));
fprintf(data_ngf3,'R3 2 4 %.11fk\n', DATA(3));
fprintf(data_ngf3,'R4 4 0 %.11fk\n', DATA(4));
fprintf(data_ngf3,'R5 4 5 %.11fk\n', DATA(5));
fprintf(data_ngf3,'R6 0 8 %.11fk\n', DATA(6));
fprintf(data_ngf3,'R7 6 7 %.11fk\n', DATA(7));
fprintf(data_ngf3,'C 5 7 %.11fu\n', DATA(9));
fprintf(data_ngf3,'Vs 1 0 DC 0\n');
fprintf(data_ngf3,'Vaux 8 6 DC 0\n');
fprintf(data_ngf3,'Hd 4 7 vaux %.11fk\n', DATA(11));
fprintf(data_ngf3,'Gb 5 3 (2,4) %.11fm\n', DATA(10));
fprintf(data_ngf3,'.ic v(5)=%.6f v(7)=%f', V5n,V7n);
fclose(data_ngf3);

%%%%%NATURAL SOLUTION 3) NUMERIC%%%%%

A = double(Vx); %initial condition V5(0) = Vx
s = 1/tau; %pois Req e negativo logo tau e negativo

t=0:1e-6:20e-3;
vc_n = A*exp(s*t);
v5_n = vc_n;

%%%%%NATURAL SOLUTION 3) PLOT%%%%%

hf = figure (1);
plot (t*1000, v5_n, "g");
xlabel ("t [ms]");
ylabel ("v5_n [V]");
legend('v5_n(t)','Location','northeast');
print (hf, "v5_n.eps", "-depsc");

%%%%%FORCED SOLUTION 4)%%%%%

Vs_p = sym (0-j);
f = 1000;
w = sym('2000*pi');
Zc = sym (0-1/(w*C)*j);
syms V1n_p V2n_p V3n_p V4n_p V5n_p V6n_p V7n_p 

Eqd4_2 = (V2n_p-V1n_p)/R1 + (V2n_p-V4n_p)/R3 + (V2n_p-V3n_p)/R2 == 0;
Eqd4_0_R6 = (V1n_p-V2n_p)/R1 + (-V4n_p)/R4 + (-V6n_p)/R6 == 0;
Eqd4_0_R7 = (V1n_p-V2n_p)/R1 + (-V4n_p)/R4 + (V6n_p-V7n_p)/R7 == 0;
Eqd4_Vs = V1n_p == Vs_p;
Eqd4_Vd = V4n_p-V7n_p == Kd*(-V6n_p)/R6;
Eqd4_5 = (V5n_p-V7n_p)/Zc + (V5n_p-V4n_p)/R5 + Kb*(V2n_p-V4n_p) == 0;
Eqd4_3 = (V3n_p-V2n_p)/R2 == Kb*(V2n_p-V4n_p);
sn_4 = solve(Eqd4_2,Eqd4_0_R6,Eqd4_0_R7,Eqd4_Vs,Eqd4_Vd,Eqd4_5,Eqd4_3);

%%%%%FORCED SOLUTION 4) DATA FILE AND PLOT%%%%%

diary "phasors_tab.tex"
diary on 
printf("$V_1$ & %f\n",double(abs(sn_4.V1n_p)));
printf("$V_2$ & %f\n",double(abs(sn_4.V2n_p)));
printf("$V_3$ & %f\n",double(abs(sn_4.V3n_p)));
printf("$V_4$ & %f\n",double(abs(sn_4.V4n_p)));
printf("$V_5$ & %f\n",double(abs(sn_4.V5n_p)));
printf("$V_6$ & %f\n",double(abs(sn_4.V6n_p)));
printf("$V_7$ & %f\n",double(abs(sn_4.V7n_p)));
printf("$V_8$ & %f\n",double(abs(sn_4.V6n_p)));
diary off


M_V5n_p = double(abs(sn_4.V5n_p));
A_V5n_p = double(angle(sn_4.V5n_p));

t=0:1e-6:20e-3;
v5_f = M_V5n_p*cos(double(w)*t+A_V5n_p);

hf = figure (2);
plot (t*1000, v5_f, "r");
xlabel ("t [ms]");
ylabel ("v5_f [V]");
legend('v5_f(t)','Location','northeast');
print (hf, "v5_f.eps", "-depsc");

%%%%%SOLUTION 5)%%%%%

%%%%%SOLUTION 5) WRITE FILE DATA_4_NG.TXT%%%%%

data_ngf4 = fopen('data_4_ng.txt','w');
fprintf(data_ngf4,'R1 1 2 %.11fk\n', DATA(1));
fprintf(data_ngf4,'R2 2 3 %.11fk\n', DATA(2));
fprintf(data_ngf4,'R3 2 4 %.11fk\n', DATA(3));
fprintf(data_ngf4,'R4 4 0 %.11fk\n', DATA(4));
fprintf(data_ngf4,'R5 4 5 %.11fk\n', DATA(5));
fprintf(data_ngf4,'R6 0 8 %.11fk\n', DATA(6));
fprintf(data_ngf4,'R7 6 7 %.11fk\n', DATA(7));
fprintf(data_ngf4,'C 5 7 %.11fu\n', DATA(9));
fprintf(data_ngf4,'Vs 1 0 ac 1.0 -90 sin(0 1 1k)\n');
fprintf(data_ngf4,'Vaux 8 6 DC 0\n');
fprintf(data_ngf4,'Hd 4 7 vaux %.11fk\n', DATA(11));
fprintf(data_ngf4,'Gb 5 3 (2,4) %.11fm\n', DATA(10));
fprintf(data_ngf4,'.ic v(5)=%.6f v(7)=%f', V5n,V7n);
fclose(data_ngf4);

%%%%%SOLUTION 5) PLOT%%%%%

t=0:1e-6:20e-3;
v5 = v5_n + v5_f;
vs = sin(double(w)*t);

ti=-5e-3:1e-6:0;
tt = cat(2,ti,t);

v5i = double(sn_1.V5n)*ones(1,size(ti,2));
v5 = cat(2,v5i,v5);

vsi = DATA(8)*ones(1,size(ti,2));
vs = cat(2,vsi,vs);

hf = figure (3);
plot (tt*1000, v5, "r",tt*1000, vs, "b");
xlabel ("t [ms]");
ylabel ("v [V]");
legend('v5(t)','vs(t)','Location','northeast');
print (hf, "v5_vs.eps", "-depsc");

%%%%%FREQUENCY ANALYSIS 6)%%%%%

syms f;
Zc = sym (0-1/(2*sym(pi)*f*C)*j);
syms V1n_p V2n_p V3n_p V4n_p V5n_p V6n_p V7n_p 

Eqd5_2 = (V2n_p-V1n_p)/R1 + (V2n_p-V4n_p)/R3 + (V2n_p-V3n_p)/R2 == 0;
Eqd5_0_R6 = (V1n_p-V2n_p)/R1 + (-V4n_p)/R4 + (-V6n_p)/R6 == 0;
Eqd5_0_R7 = (V1n_p-V2n_p)/R1 + (-V4n_p)/R4 + (V6n_p-V7n_p)/R7 == 0;
Eqd5_Vs = V1n_p == Vs_p;
Eqd5_Vd = V4n_p-V7n_p == Kd*(-V6n_p)/R6;
Eqd5_5 = (V5n_p-V7n_p)/Zc + (V5n_p-V4n_p)/R5 + Kb*(V2n_p-V4n_p) == 0;
Eqd5_3 = (V3n_p-V2n_p)/R2 == Kb*(V2n_p-V4n_p);
sn_5 = solve(Eqd5_2,Eqd5_0_R6,Eqd5_0_R7,Eqd5_Vs,Eqd5_Vd,Eqd5_5,Eqd5_3);

%%%%%FREQUENCY ANALYSIS 6) NUMERIC%%%%%

freq = logspace(-1,6,200);

fh = function_handle(sn_5.V5n_p);
v5_freq = fh(freq); 

fh = function_handle(Vs_p);
vs_freq = fh(freq); 

Vc_p = sn_5.V5n_p - sn_5.V7n_p ;
fh = function_handle(Vc_p);
vc_freq = fh(freq);

%%%%%FREQUENCY ANALYSIS 6) PLOT%%%%%

hf = figure (4);
semilogx(freq,20*log10(abs(v5_freq)), "r",freq,20*log10(abs(vc_freq)), "g",freq,20*log10(abs(vs_freq)),"b");
xlabel ("f [Hz]");
ylabel ("|v| dB");
legend('v5(f)','vc(f)','vs(f)','Location','southwest');
print (hf, "freqresp.eps", "-depsc");

hf = figure (5);
hold on
semilogx(freq,180/pi*angle(v5_freq), "r",freq,180/pi*angle(vc_freq), "g",freq,180/pi*angle(vs_freq),"b");
xlabel ("f [Hz]");
ylabel ("phase [degrees]");
legend('phase.v5(f)','phase.vc(f)','phase.vs(f)','location','northwest');
print (hf, "phase_oct.eps", "-depsc");




