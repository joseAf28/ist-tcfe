

% file = fopen("../doc/GainStage_OP.tex", "w");
% fprintf(file, "VEQ & %e\\\\ \\\hline\n", VEQ);
% fprintf(file, "IB1 & %e\\\\ \\\hline\n", IB1);
% fprintf(file, "IC1 & %e\\\\ \\\hline\n", IC1);
% fprintf(file, "VEQ & %e\\\\ \\\hline\n", IE1);
% fprintf(file, "VE1 & %e\\\\ \\\hline\n", VE1);
% fprintf(file, "V01 & %e\\\\ \\\hline\n", VO1);
% fprintf(file, "VCE & %e\\\\ \\\hline\n", VCE);
% fclose(file);

% file = fopen("../doc/GainStage_AC.tex", "w");
% fprintf(file, "gm1 & %e\\\\ \\\hline\n", gm1);
% fprintf(file, "$r \pi 1$ & %e\\\\ \\\hline\n", rpi1);
% fprintf(file, "r01 & %e\\\\ \\\hline\n", ro1);
% fprintf(file, "AV1 & %e\\\\ \\\hline\n", AV1);
% fprintf(file, "$AV1_{DB}$ & %e\\\\ \\\hline\n", AVI_DB);
% fprintf(file, "ZI1& %e\\\\ \\\hline\n", ZI1);
% fprintf(file, "ZO1 & %e\\\\ \\\hline\n", ZO1);
% fclose(file);

% file = fopen("../doc/OutputStage_OP.tex", "w");
% fprintf(file, "VI2 & %e\\\\ \\\hline\n", VO1);
% fprintf(file, "IC2 & %e\\\\ \\\hline\n", IC2);
% fprintf(file, "IE2 & %e\\\\ \\\hline\n", IE2);
% fprintf(file, "VO2 & %e\\\\ \\\hline\n", VO2);
% fclose(file);

% file = fopen("../doc/OutputStage_AC.tex", "w");
% fprintf(file, "gm2 & %e\\\\ \\\hline\n", gm2);
% fprintf(file, "$r \pi 2$ & %e\\\\ \\\hline\n", 1/gpi2);
% fprintf(file, "r02 & %e\\\\ \\\hline\n", 1/go2);
% fprintf(file, "AV2 & %e\\\\ \\\hline\n", AV2);
% fprintf(file, "$AV_{DB}$ & %e\\\\ \\\hline\n", 20*log10(AV2));
% fprintf(file, "ZI2 & %e\\\\ \\\hline\n", ZI2);
% fprintf(file, "ZO2 & %e\\\\ \\\hline\n", ZO2);
% fclose(file);


% %total
% gB = 1/(1/gpi2+ZO1)
% AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
% AV_DB = 20*log10(abs(AV))
% ZI=ZI1
% ZO=1/(go2+gm2/gpi2*gB+ge2+gB)


% file = fopen("../doc/Final.tex", "w");
% fprintf(file, "AV & %e\\\\ \\\hline\n", AV);
% fprintf(file, "$AV_{DB}$ & %e\\\\ \\\hline\n", AV_DB);
% fprintf(file, "ZI & %e\\\\ \\\hline\n", ZI);
% fprintf(file, "ZO & %e\\\\ \\\hline\n", ZO);
% fclose(file);

jj = 1;

Cmax = 1e-6;
Cmin = 220e-9;

R1k = 1e3;
R10k = 10e3;
R100k = 100e3;

invR1 = 3./R1k + 3./R10k;
invR2 = 3./R100k;

invC1 = 3./Cmax;
invC2 = 3./Cmin;


TGain = zeros(1, 10);
TGain_DB = zeros(1, 10);
phase = zeros(1, 10);
logfreq = zeros(1, 10);

C1  = 1./invC1;
C2 = 1./invC2;

R1 = 1./invR1;
R2 = 1./invR2;   

vin = 1.;

for t = 1:0.005:6

omega = 2*pi*power(10, t);
ZC1 = 1./(I*omega*C1);
ZC2 = 1./(I*omega*C2);

eq1 = [1./R1 + 1./ZC2 + 1./ZC1, -1./ZC2];
eq2 = [-1./ZC1, -1./R2];
A = [eq1; eq2];

B = [vin/R1; 0];

X = A\B;
Vout = X(2);
TGain(jj) = Vout*vin;
TGain_DB(jj) = 20*log10(abs(TGain(jj)));
phase(jj) = atan(TGain(jj));
logfreq(jj) = t;

jj = jj +1;
endfor

figure
plot(logfreq, TGain_DB, "b");
title("Gain");
ylabel ("dB");
xlabel ("log10 frequency [Hz]");
print ("gain.eps", "-depsc");

figure
plot(logfreq, phase, "b");
title("Phase");
ylabel ("rad");
xlabel ("log10 frequency [Hz]");
print ("phase.eps", "-depsc");


[Gmax, idxfreqmax] = max(TGain_DB)

freqmax = power(10, logfreq(idxfreqmax))

%impedance input

omega = 2*pi*freqmax;
ZC1 = 1./(I*omega*C1);
ZC2 = 1./(I*omega*C2);

eq1 = [1./R1 + 1./ZC2 + 1./ZC1, -1./ZC2];
eq2 = [-1./ZC1, -1./R2];
A = [eq1; eq2];

B = [vin/R1; 0];

X = A\B;

current_in = -(vin -X(1))/R1
zinput = vin/current_in
