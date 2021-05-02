close all
clear all

format long;
A_vin = 230;
f = 50;
w=2*pi*f;
n=230/1500;


A_vi = A_vin/n;

t=linspace(4e-1, 6e-1, 50000);

vIN = A_vin*cos(w*t);
vI = A_vi*cos(w*t);

hf=figure(1);
plot(t*1000, vI);
title("Input voltage");
xlabel ("t[ms]");
legend("vI");
print (hf,"vI.eps", "-depsc");

%%%%%%%%%%%%%%ENVELOPE DETECTOR%%%%%%%%%%%%%%
Renv=1000e3;
Cenv=1020e-6;

vOenv = zeros(1, length(t));
vOr = zeros(1, length(t));
tOFF = 1/w * atan(1/w/Renv/Cenv) + 0.4;
vOnexp = A_vi*cos(w*tOFF)*exp(-(t-tOFF)/Renv/Cenv);

for i=1:length(t)
  %if (vI(i) > 0) %half rectifier
    %vOr(i) = vI(i); %half rectifier
    vOr(i) = abs(vI(i)); %full rectifier
  %endif
endfor

for i=1:length(t)
  if t(i) < tOFF
    vOenv(i) = vOr(i);
  elseif vOnexp(i) > vOr(i) %qd vOnexp(i) = vOr(i) é t=tON
    vOenv(i) = vOnexp(i);
  else
    %tOFF = tOFF + 1/f; %half rectifier
    tOFF = tOFF + 1/f/2; %full rectifier
    vOnexp = A_vi*abs(cos(w*tOFF))*exp(-(t-tOFF)/Renv/Cenv);
    vOenv(i) = vOr(i);
  endif
endfor

vripple_env = max(vOenv)-min(vOenv);


hf=figure(2);
%plot(t*1000, vOenv,t*1000, vOr);
plot(t*1000, vOenv);
title("Output voltage");
xlabel ("t[ms]");
%legend("envelope","rectified");
legend("envelope");
print (hf,"vOenv_vOr.eps", "-depsc");

%%%%%%%%%%%%%%VOLTAGE REGULATOR%%%%%%%%%%%%%%
Rreg = 951.06e3;
Vt = 0.026;
N = 1;
Is = 1*10^(-14);
n_diodos = 18;
VON = 12/n_diodos;

rd = Vt*N/(Is*exp(VON/(N*Vt)))

VIreg = mean(vOenv);
vireg = vOenv - VIreg;

VOreg = n_diodos*VON;
voreg = n_diodos*rd/(n_diodos*rd+Rreg)*vireg;
vOreg = VOreg + voreg;

vripple_reg = max(vOreg)-min(vOreg)
hf=figure(3);
plot(t*1000, vOreg);
ylim([12-0.5*10^-5 12+0.5*10^-5])
yticks = get (gca, "ytick");
ylabels = arrayfun (@(x) sprintf ("%.6f", x), yticks, "uniformoutput", false);
set (gca, "yticklabel", ylabels)
title("Output voltage");
xlabel ("t[ms]");
legend("regulator");
print (hf,"vOreg.eps", "-depsc");

hf=figure(4);
plot(t*1000, vOreg-12);
title("Output voltage");
xlabel ("t[ms]");
legend("vO-12");
print (hf,"vO-12.eps", "-depsc");


Cost = 1*(Rreg+Renv)/1000 + Cenv*10^6 + 0.1*(4+n_diodos);
Merit = 1/(Cost * (vripple_reg + abs(mean(vOreg-12)) + 10^(-6)));

voavg = mean(vOreg);
vomax = max(vOreg);
vomin = min(vOreg);
vripplereg = vripple_reg;
voenvavg = VIreg;
voenvmax = max(vOenv);
voenvmin = min(vOenv);
vrippleenv = vripple_env;


diary "Envelope_tab.tex"
diary on
printf("$V_I$ & %.3f\n", A_vi);
printf("$V_{max}$ & %.3f\n", voenvmax);
printf("$V_{min}$ & %.3f\n", voenvmin);
printf("$V_{ripple}$ & %.6e\n",vrippleenv);
printf("$V_{avg}$ & %.3f\n",voenvavg);
diary off


diary "Regulator_tab.tex"
diary on
printf("$V_{max}$ & %.5f\n", vomax);
printf("$V_{min}$ & %.5f\n", vomin);
printf("$V_{ripple}$ & %.6e\n",vripplereg);
printf("$V_{avg}$ & %.5f\n",voavg);
printf("$V_{avg}$-12 & %.5f\n",voavg-12);
diary off

diary "Merit_tab.tex"
diary on
printf("Cost & %.3f\n", Cost);
printf("Merit & %.4f\n", Merit);
diary off




