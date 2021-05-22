%gain stage

VT=25e-3
BFN=178.7
VAFN=69.7
RE1=100
RC1=1000
RB1=80000
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


%DC analysis gain stage

RB=1/(1/RB1+1/RB2) %RB paralelo
VEQ=RB2/(RB1+RB2)*VCC #Voltage thevanin
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1 #IC1 = beta*IB1
IE1=(1+BFN)*IB1 #IE1 = IC1 + IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1

BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
load = 8

#AC analysis gain stage consts
gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7

%DC analysis output stage
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/(RE2);

RSB=RB*RS/(RB+RS)

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))


ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)

file = fopen("../doc/GainStage_OP.txt", "w");
fprintf(file, "VEQ & %e\\\\ \\\hline\n", VEQ);
fprintf(file, "IB1 & %e\\\\ \\\hline\n", IB1);
fprintf(file, "IC1 & %e\\\\ \\\hline\n", IC1);
fprintf(file, "VEQ & %e\\\\ \\\hline\n", IE1);
fprintf(file, "VE1 & %e\\\\ \\\hline\n", VE1);
fprintf(file, "V01 & %e\\\\ \\\hline\n", VO1);
fprintf(file, "VCE & %e\\\\ \\\hline\n", VCE);
fclose(file);

file = fopen("../doc/GainStage_AC.txt", "w");
fprintf(file, "gm1 & %e\\\\ \\\hline\n", gm1);
fprintf(file, "$r \pi 1 & %e\\\\ \\\hline\n", rpi1);
fprintf(file, "r01 & %e\\\\ \\\hline\n", ro1);
fprintf(file, "AV1 & %e\\\\ \\\hline\n", AV1);
fprintf(file, "AV1_DB & %e\\\\ \\\hline\n", AVI_DB);
fprintf(file, "ZI1& %e\\\\ \\\hline\n", ZI1);
fprintf(file, "ZO1 & %e\\\\ \\\hline\n", ZO1);
fclose(file);

file = fopen("../doc/OutputStage_OP.txt", "w");
fprintf(file, "VI2 & %e\\\\ \\\hline\n", VO1);
fprintf(file, "IC2 & %e\\\\ \\\hline\n", IC2);
fprintf(file, "IE2 & %e\\\\ \\\hline\n", IE2);
fprintf(file, "VO2 & %e\\\\ \\\hline\n", VO2);
fclose(file);

file = fopen("../doc/OutputStage_AC.txt", "w");
fprintf(file, "gm2 & %e\\\\ \\\hline\n", gm2);
fprintf(file, "$r \pi 2$ & %e\\\\ \\\hline\n", 1/gpi2);
fprintf(file, "r02 & %e\\\\ \\\hline\n", 1/go2);
fprintf(file, "AV2 & %e\\\\ \\\hline\n", AV2);
fprintf(file, "AV_DB & %e\\\\ \\\hline\n", 20*log10(AV2));
fprintf(file, "ZI2 & %e\\\\ \\\hline\n", ZI2);
fprintf(file, "ZO2 & %e\\\\ \\\hline\n", ZO2);
fclose(file);


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

file = fopen("../doc/Final.txt", "w");
fprintf(file, "AV & %e\\\\ \\\hline\n", AV);
fprintf(file, "AV_BD & %e\\\\ \\\hline\n", AV_DB);
fprintf(file, "ZI & %e\\\\ \\\hline\n", ZI);
fprintf(file, "ZO & %e\\\\ \\\hline\n", ZO);
fclose(file);

jj = 1;

RE1 = 100

TGain = zeros(1, 10);
TGain_DB = zeros(1, 10);
logfreq = zeros(1, 10);

C1 = 1e-3;
C2 = 10e-3;
C3 =2e-3;

vin = 10e-3;

for t = 1:0.1:7

omega = 2*pi*power(10, t);
ZC1 = 1./(I*omega*C1);
ZC2 = 1./(I*omega*C2);
ZC3 = 1./(I*omega*C3);

ZReC1 = 1/(1/RE1 + 1/ZC2);
Zeq = 1/(1/RE2 + 1/(load + ZC3));

eq1 = [RS + ZC1 + RB, -RB, 0, 0, 0, 0, 0];
eq2 = [-RB, RB+rpi1 + ZReC1, 0, -ZReC1, 0, 0, 0];
eq3 = [0, rpi1*gm1, 1, 0, 0, 0, 0];
eq4 = [0, -ZReC1, -ro1, ZReC1+ro1+RC1, -RC1, 0, 0];
eq5 = [0, 0, 0, -RC1, 1/gpi2+RC1+Zeq, 0, -Zeq];
eq6 = [0, 0, 0, 0, 1/gpi2*gm2, 1, 0];
eq7 = [0, 0, 0, 0, -Zeq, -1/go2, Zeq+1/go2];

A = [eq1; eq2; eq3; eq4; eq5; eq6; eq7];

B = [vin; 0; 0; 0; 0; 0; 0];

X = A\B;
Vout = (X(7) - X(5))*Zeq;
TGain(jj) = Vout*load/(load + ZC3)/vin;
TGain_DB(jj) = 20*log10(abs(TGain(jj)));
logfreq(jj) = t;

jj = jj +1;
endfor


hGain = figure();
hACtotal = figure();
plot(logfreq, TGain_DB, "b");
print(hACtotal, "gainACtotal.eps", "-depsc");



