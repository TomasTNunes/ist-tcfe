clear all

R1 = 1000
R2 = 1000
R3 = 100000
R4 = 1000
C1 = 220e-9
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


wL = 1/(R1*C1)
wH = 1/(R2*C2)
wO = sqrt(wL*wH)
f = wO/(2*pi)

gain = abs((R1*C1*wO*j)/(1+R1*C1*wO*j)*(1+R3/R4)*(1/(1+R2*C2*wO*j)))
gain_db = 20*log10(abs(gain))

Z_in = abs(R1 + 1/(j*wO*C1))
Z_out = abs(1/(j*wO*C2+1/R2))

Cost = (R1+R2+R3+R4)/1000 + (C1+C2)*1000000
gain_deviation = abs(100-gain)
frequency_deviation = abs(f-1000)
Merit = 1/(Cost*gain_deviation*frequency_deviation)
