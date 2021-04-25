close all
clear all


A_vin = 230;
f = 50;
w=2*pi*f;
n=10;


A_vi = A_vin/n;

t=linspace(5e-1, 7e-1, 1000);

vIN = A_vin*cos(w*t);
vI = A_vi*cos(w*t);

hf=figure(1);
plot(t*1000, vIN,t*1000, vI);
title("Input voltage");
xlabel ("t[ms]");
legend("vIN","vI");
print (hf,"vIN_vI.eps", "-depsc");

%%%%%%%%%%%%%%ENVELOPE DETECTOR%%%%%%%%%%%%%%
Renv=10e3;
Cenv=10e-6;

vOenv = zeros(1, length(t));
vOr = zeros(1, length(t));
tOFF = 1/w * atan(1/w/Renv/Cenv) + 0.5;
vOnexp = A_vi*cos(w*tOFF)*exp(-(t-tOFF)/Renv/Cenv);

for i=1:length(t)
  if (vI(i) > 0) %half rectifier
    vOr(i) = vI(i); %half rectifier
    %vOr(i) = abs(vI(i)); %full rectifier
  endif
endfor

for i=1:length(t)
  if t(i) < tOFF
    vOenv(i) = vOr(i);
  elseif vOnexp(i) > vOr(i) %qd vOnexp(i) = vOr(i) é t=tON
    vOenv(i) = vOnexp(i);
  else
    tOFF = tOFF + 1/f; %half rectifier
    %tOFF = tOFF + 1/f/2; %full rectifier
    vOnexp = A_vi*abs(cos(w*tOFF))*exp(-(t-tOFF)/Renv/Cenv);
    vOenv(i) = vOr(i);
  endif
endfor

vripple_env = max(vOenv)-min(vOenv)


hf=figure(2);
plot(t*1000, vOenv,t*1000, vOr);
title("Output voltage");
xlabel ("t[ms]");
legend("envelope","rectified");
print (hf,"vOenv_vOr.eps", "-depsc");

%%%%%%%%%%%%%%VOLTAGE REGULATOR%%%%%%%%%%%%%%
Rreg = 10e3;
Vt = 0.026;
N = 1;
Is = 1*10^(-14);
n_diodos = 17;
VON = 0.7;

rd = Vt*N/(Is*exp(VON/(N*Vt)))

VIreg = max(vOenv) - vripple_env/2
vireg = vOenv - VIreg;

VOreg = n_diodos*VON
voreg = n_diodos*rd/(n_diodos*rd+Rreg)*vireg;
vOreg = VOreg + voreg;

vripple_reg = max(vOreg)-min(vOreg)

hf=figure(3);
plot(t*1000, vOreg);
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

mean(vOreg-12)
Quality = 1/vripple_reg + abs(1/mean(vOreg-12))
Cost = 1*(Rreg+Renv)/1000 + Cenv*10^6 + 0.1*(1+n_diodos)
Merit = Quality/Cost





