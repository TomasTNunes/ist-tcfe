close all
clear all

pkg load symbolic

A_vin = 230
f = 50
w=2*pi*f
n=1


A_vi = A_vin/n

t=linspace(0, 2e-1, 1000);

vIN = A_vin*cos(w*t)
vI = A_vi*cos(w*t)


%%%%%%%%%%%%%%ENVELOPE DETECTOR%%%%%%%%%%%%%%
Renv=1
Cenv=1

vOenv = zeros(1, length(t));
vOr = zeros(1, length(t));
tOFF = 1/w * atan(1/w/R/C);
vOnexp = A*cos(w*tOFF)*exp(-(t-tOFF)/R/C);


