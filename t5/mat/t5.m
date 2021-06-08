clear all

R1 = 1000
R2 = 1000
R3 = 315000
R4 = 1000
C1 = 110e-9
C2 = 220e-9


freq = logspace(1,8,100);

Tf = (R1*C1*2*pi*freq*j)./(1+R1*C1*2*pi*freq*j)*(1+R3/R4).*(1./(1+R2*C2*2*pi*freq*j));

f1 = figure();
semilogx(freq,20*log10(abs(Tf)));
xlabel("Frequency [Hz]");
ylabel("Gain [dB]");
title("Gain");
print(f1, "teo_gain.eps", "-depsc");

f2 = figure();
semilogx(freq,180*arg(Tf)/pi);
xlabel("Frequency [Hz]");
ylabel("Phase [Deg]");
title("Phase");
print(f2, "teo_phase.eps", "-depsc");


wL = 1/(R2*C2)
fL = wL/(2*pi)
wH = 1/(R1*C1)
fH = wH/(2*pi)
wO = sqrt(wL*wH)
fO = wO/(2*pi)

gain = abs((R1*C1*2000*pi*j)/(1+R1*C1*2000*pi*j)*(1+R3/R4)*(1/(1+R2*C2*2000*pi*j)))
gain_db = 20*log10(abs(gain))

Z_in = abs(R1 + 1/(j*2000*pi*C1))
Z_out = abs(1/(j*2000*pi*C2+1/R2))

Cost = (R1+R2+R3+R4+15000)/1000 + (4*C1+C2)*1000000 + 13323
gain_deviation = abs(100-gain)
frequency_deviation = abs(fO-1000)
Merit = 1/(Cost*(gain_deviation+frequency_deviation+10^(-6)))



diary "data1_teo.tex"
diary on
printf("$f_{L}$ & %.2f\n", fL);
printf("$f_{H}$ & %.2f\n", fH);
printf("$f_{O}$ & %.2f\n", fO);
diary off


diary "data2_teo.tex"
diary on
printf("$gain$ & %.2f\n", gain);
printf("$gain_{dB}$ & %.3f\n", gain_db);
printf("$Z_{IN}$ & %.2f\n", Z_in);
printf("$Z_{OUT}$ & %.2f\n", Z_out)
diary off



diary "data3_teo.tex"
diary on
printf("$gain_{dev}$ & %.6f\n", gain_deviation);
printf("$freq_{dev}$ & %.4f\n", frequency_deviation);
printf("$Cost$ & %.2f\n", Cost);
printf("$Merit$ & %.6e\n", Merit);
diary off



