%gain stage

VT=25e-3
BFN=178.7#beta
VAFN=69.7#alpha ??
RE1=100
RC1=1000
RB1=80000
RB2=20000
VBEON=0.7
VCC=12
RS=100


%DC analysis sem approx grosseira
RB=1/(1/RB1+1/RB2) %RB paralelo
VEQ=RB2/(RB1+RB2)*VCC #Voltage thevanin
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1 #IC1 = beta*IB1
IE1=(1+BFN)*IB1 #IE1 = IC1 + IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1 #lei dos nodos
VCE=VO1-VE1

%Capacitor in paralellel with RE
CE1 = 1e-3;

%Capacitor in series with RB
CB1 = 1e-3;

#AC analysis 
gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=RB*RS/(RB+RS)

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))


ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)


AVI_DB_vec = zeros(1, 100);
AVI_vec = zeros(1, 100);
logfreq = zeros(1, 100);
ZO1_vec = zeros(1, 100);
ZI1_vec = zeros(1, 100);

jj= 1;
RE1 = 100;
for freq = 10:100:10e6
omega = 2*pi*freq;

ZCE1 = 1/(I*omega*CE1); %bypass capacitor RE
RE1 = 1/((1/RE1 + 1/ZCE1));

RB = RB + 1/(I*omega*1e-3);
RSB=RB*RS/(RB+RS);


AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2); #gain
AVI_DB = 20*log10(abs(AV1));
AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1);
AVIsimple_DB = 20*log10(abs(AV1simple));

AVI_vec(jj) = AV1;
AVI_DB_vec(jj) = AVI_DB;
logfreq(jj) = log(freq);

ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)));
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)));
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) );
ZO1 = 1/(1/ZX+1/RC1);

ZI1_vec(jj) = ZI1;
ZO1_vec(jj) = ZO1;

jj = jj +1;
endfor

hAC1 = figure();
plot(logfreq, AVI_DB_vec, "b");
print(hAC1, "gainAC1.eps", "-depsc");



%ouput stage

BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7

%DC analysis2
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/(RE2);

%AC analysis2
AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)


%coupling capacitor R2
C2 = 1e-6;

AV2_vec = zeros(1, 100);
ZO2_vec = zeros(1, 100);
ZI2_vec = zeros(1, 100);

jj = 1;

for freq = 10:100:10e6

omega = 2*pi*freq;
ZC2 = 1/(I*omega*C2);

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/(RE2);

%AC analysis2
AV2_vec(jj) = gm2/(gm2+gpi2+go2+ge2);
ZI2_vec(jj) = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
ZO2_vec(jj) = 1/(gm2+gpi2+go2+ge2);

jj = jj +1;
endfor


AV_DBvec = zeros(1, 100);
ZO_vec = zeros(1, 100);

%total
jj= 1;

for freq = 10:100:10e6

omega = 2*pi*freq;
ZC2 = 1/(I*omega*C2);

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/(RE2 + ZC2);

gB = 1/(1/gpi2+ZO1_vec(jj));
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AVI_vec(jj);
AV_DBvec(jj) = 20*log10(abs(AV));
ZI=ZI1_vec(jj);
ZO_vec(jj)=1/(go2+gm2/gpi2*gB+ge2+gB);

jj = jj +1;

endfor

hACtotal = figure();
plot(logfreq, AV_DBvec, "b");
print(hACtotal, "gainACtotal.eps", "-depsc");