%gain stage

VT=25e-3
BFN=178.7
VAFN=69.7
RE1=100
RC1=500
RB1=100000
RB2=20000
VBEON=0.7
VCC=12
RS=100

RB=1/(1/RB1+1/RB2)
VEQ=RB2/(RB1+RB2)*VCC
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1


gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=RB*RS/(RB+RS)

AV1_RE = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB_RE = 20*log10(abs(AV1_RE))
AV1simple_RE = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB_RE = 20*log10(abs(AV1simple_RE))

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=100
ZI1_RE = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX_RE = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX_RE = ro1*(1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB))/(1/RE1+1/(rpi1+RSB) ) 
ZO1_RE = 1/(1/ZX_RE+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

%ouput stage
RE1=100
RL = 8
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

 
Vin=0.01
rpi2=1/gpi2
ro2=1/go2

freq = logspace(1, 8, 700);
gain = zeros(1,700);


for a = 1:1:700
w=2*pi*freq(a);
Zci=1/(j*w*1000*10^(-6));
Zcb=1/(j*w*5000*10^(-6));
Zco=1/(j*w*1995*10^(-6));

Y =[1,0,0,0,0,0,0;
-1/RS, 1/RS+1/Zci, -1/Zci, 0,0,0,0;
0, -1/Zci, 1/Zci+1/RB+1/rpi1, -1/rpi1,0,0,0;
0,0,-gm1-1/rpi1, gm1+1/rpi1 + 1/RE1 + 1/Zcb + 1/ro1, -1/ro1, 0,0;
0,0, gm1, -gm1-1/ro1, 1/ro1 + 1/RC1 + 1/rpi2, -1/rpi2, 0;
0,0,0,0,0, -1/Zco, 1/Zco + 1/RL;
0,0,0,0, -gm2-1/rpi2, gm2 + 1/rpi2 + 1/RE2 + 1/Zco + 1/ro2, -1/Zco;
];

B=[Vin; 0; 0; 0; 0; 0; 0];
X = Y\B;
	
gain_db(a)=20*log10(abs(X(7)/X(1)));

endfor

max_gain_db_3 = max(gain_db) - 3
minn = 1

for a = 1:1:700
if abs(gain_db(a)-max_gain_db_3) < minn
minn = abs(gain_db(a)-max_gain_db_3);
lco=freq(a);
bbb=gain_db(a);
endif
endfor
lco


hf1 = figure (1);
semilogx(freq, gain_db);
xlabel ("f [Hz]");
ylabel ("gain [dB]");
%legend('gain(f) dB','Location','southeast');
print(hf1, "gain_db_teo.eps", "-depsc");














