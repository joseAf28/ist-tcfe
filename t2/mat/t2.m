close all
clear all

fid = fopen('../data.txt','r');
for i = 1:8
fskipl(fid);
end

vec = 1:1:10;
#format long

R1 = strsplit(fgetl(fid));
vec(1) = R1 = str2double(R1(4));
R2 = strsplit(fgetl(fid))
vec(2) = R2 = str2double(R2(3))
R3 = strsplit(fgetl(fid))
vec(3) = R3 = str2double(R3(3))
R4 = strsplit(fgetl(fid))
vec(4) = R4 = str2double(R4(3))
R5 = strsplit(fgetl(fid))
vec(5) = R5 = str2double(R5(3))
R6 = strsplit(fgetl(fid))
vec(6) = R6 = str2double(R6(3))
R7 = strsplit(fgetl(fid))
vec(7) = R7 = str2double(R7(3))
Vs = strsplit(fgetl(fid))
vec(8) = Vs = str2double(Vs(3))
Kb = strsplit(fgetl(fid))
vec(9) = Kb = str2double(Kb(3));
Kd = strsplit(fgetl(fid))
vec(10) = Kd = str2double(Kd(3))

fclose(fid);


file = fopen("../sim/data2spice.txt", "w");

for i = 1:10
fprintf(file, "%d \n", vec(i));
end

fclose(file);


#{

format long


printf("\n alinea 1)------------------------------------------\n\n")

line1 = [1, 0, 0, 0, 0, 0, 0];
line2 = [0, 1/R3+1/R2+1/R1, -1/R2, -1/R3, 0, 0, 0];
line3 = [0, Kb+1/R2, -1/R2, -Kb, 0, 0, 0];
line58 = [0, -1/R3, 0, 1/R5+1/R4+1/R3, -1/R5, -1/R7, 1/R7];
line6 = [0, Kb, 0, -Kb-1/R5, 1/R5, 0, 0];
line7 = [0, 0, 0, 0, 0, -1/R6-1/R7, 1/R7];
line8 = [0, 0, 0, 1, 0, Kd/R6, -1];

An = [line1; line2; line3; line58; line6; line7; line8]

bn = [Vs; Vs/R1; 0; 0; 0; 0; 0]

#det(An)

Xn = An\bn


Vx = Xn(5) - Xn(7)


printf("\n alinea 2)-------------------------------------------\n\n")

printf("Nodal Analysis---------------")

A = [1/R1+1/R4+1/R6, -1/R1, 0, -1/R4, 0, -1/R6, 0; -1/R1, 1/R1+1/R2+1/R3, -1/R2, -1/R3, 0, 0, 0; 0, -1/R2-Kb, 1/R2, Kb, 0, 0, 0; -1/R6, 0, 0, 0, 0, 1/R6+1/R7, -1/R7; 0, 0, 0, 0, 1, 0, -1; -Kd/R6, 0, 0, 1, 0, Kd/R6, -1; 1, 0, 0, 0, 0, 0, 0]

det(A)

b = [0; 0; 0; 0; Vx; 0; 0]

Xn2 = A\b

Ix = -(Xn2(4)-Xn2(5))/R5 + Kb*(Xn2(2)-Xn2(4));
Req = Vx/Ix



printf("\n\n alinea 3)------------------------------------------\n\n")

time = 0:1e-6:20e-3;

v6n = Xn2(5)*exp(-time/(Req*C));

hf = figure();
plot(time, v6n, "b");
xlabel("time (s)");
ylabel("v6n(t) (V)");

print(hf, "v6n.eps", "-depsc");



printf("\n\nalinea 4)-------------------------------------------------\n\n")

omega = 2*pi*1000 #rad/s

dlineVs = [1, 0, 0, 0, 0, 0, 0];
dline2 = [-1/R1, 1/R1 + 1/R3+1/R2, -1/R2, -1/R3, 0, 0, 0];
dline58 = [0, -1/R3, 0, 1/R3+1/R4+1/R5, -i*omega*C-1/R5, -1/R7, 1/R7+i*omega*C];
dline3 = [0, Kb+1/R2, -1/R2, -Kb, 0, 0, 0];
dline6 = [0, Kb, 0, -1/R5-Kb, 1/R5+i*omega*C, 0, -i*omega*C];
dline7 = [0, 0, 0, 0, 0, -1/R6-1/R7, 1/R7];
dlineVd = [0, 0, 0, 1, 0, Kd/R6, -1];

dAn = [dline2; dline3; dline6; dline7; dlineVs; dlineVd; dline58]

dbn = [0; 0; 0; 0; -i; 0; 0] #parte real

determinante = det (dAn)

dXn = dAn\dbn

fase = arg(dXn(5))
amplitude = abs(dXn(5))

printf("\n plot solução forçada \n")

time = 0:1e-6:20e-3;

v6force = amplitude*cos(omega*time - fase);

hf2 = figure();
plot(time, v6force, "b;v6(t);");

legend();

xlabel("time (s)");
ylabel("potencial(t) (V)");

print(hf2, "v6_force.eps", "-depsc");


printf("\n\nalinea 5)-------------------------------------------------\n\n")

stepDt = fase/omega

initial_time  = -5e-3:1e-6:0;
final_time = 0:1e-6:20e-3;

v6i = 0*initial_time+Xn(5);
vsi = 0*initial_time + Vs;

v6f = v6n + amplitude*cos(omega*final_time + fase);
vsf = sin(final_time*omega);

hf2 = figure();
plot(initial_time, v6i, "b;v6(t);");
hold on
plot(final_time, v6f, "b");
plot(initial_time, vsi, "g;vs(t);");
plot(final_time, vsf, "g");

legend();

xlabel("time (s)");
ylabel("potencial(t) (V)");

print(hf2, "v6_vs.eps", "-depsc");


printf("\n\nalinea 6)-----------------------------------------------------\n\n")

solution61 = 0.1:2:1e4;
solution62 = 1e4:100:1e6;

solutionVc1 = 0.1:2:1e4;
solutionVc2 = 1e4:100:1e6;
kk = 1

for freq = 0.1:2:1e4
omega = 2*pi*freq; #rad/s

dlineVs = [1, 0, 0, 0, 0, 0, 0];
dline2 = [-1/R1, 1/R1 + 1/R3+1/R2, -1/R2, -1/R3, 0, 0, 0];
dline58 = [0, -1/R3, 0, 1/R3+1/R4+1/R5, -i*omega*C-1/R5, -1/R7, 1/R7+i*omega*C];
dline3 = [0, Kb+1/R2, -1/R2, -Kb, 0, 0, 0];
dline6 = [0, Kb, 0, -1/R5-Kb, 1/R5+i*omega*C, 0, -i*omega*C];
dline7 = [0, 0, 0, 0, 0, -1/R6-1/R7, 1/R7];
dlineVd = [0, 0, 0, 1, 0, Kd/R6, -1];

dAn = [dline2; dline3; dline6; dline7; dlineVs; dlineVd; dline58];

dbn = [0; 0; 0; 0; -i; 0; 0]; #parte real

dXn = dAn\dbn;
solution61(kk) = dXn(5);
solutionVc1(kk) = dXn(5) - dXn(7);
kk = kk +1;

end 

kk = 1;

for freq = 1e4:100:1e6
omega = 2*pi*freq; #rad/s

dlineVs = [1, 0, 0, 0, 0, 0, 0];
dline2 = [-1/R1, 1/R1 + 1/R3+1/R2, -1/R2, -1/R3, 0, 0, 0];
dline58 = [0, -1/R3, 0, 1/R3+1/R4+1/R5, -i*omega*C-1/R5, -1/R7, 1/R7+i*omega*C];
dline3 = [0, Kb+1/R2, -1/R2, -Kb, 0, 0, 0];
dline6 = [0, Kb, 0, -1/R5-Kb, 1/R5+i*omega*C, 0, -i*omega*C];
dline7 = [0, 0, 0, 0, 0, -1/R6-1/R7, 1/R7];
dlineVd = [0, 0, 0, 1, 0, Kd/R6, -1];

dAn = [dline2; dline3; dline6; dline7; dlineVs; dlineVd; dline58];

dbn = [0; 0; 0; 0; -i; 0; 0]; #parte real

dXn = dAn\dbn;
solution62(kk) = dXn(5);
solutionVc2(kk) = dXn(5) - dXn(7);
kk = kk +1;

end 

solution61 = solution61/(-i);
solution62 = solution62/(-i);
argsSol61 = arg(solution61)*(180/pi);
argsSol62 = arg(solution62)*(180/pi);

solutionVc1 = solutionVc1/(-i);
solutionVc2 = solutionVc2/(-i);
argsSolVc1 = arg(solutionVc1)*(180/pi);
argsSolVc2 = arg(solutionVc2)*(180/pi);

vecfreq1 = 20*log10(0.1:2:1e4);
vecfreq2 = 20*log10(1e4:100:1e6);

solutionVs1 = 0*vecfreq1 - i;
solutionVs2 = 0*vecfreq2 -i;

absSol61 = 20*log10(abs(solution61));
absSol62 = 20*log10(abs(solution62));

absSolVs1 = 20*log10(abs(solutionVs1));
absSolVs2 = 20*log10(abs(solutionVs2));

argsSolVs1 = arg(solutionVs1)*(180/pi);
argsSolVs2 = arg(solutionVs2)*(180/pi);

absSolVc1 = 20*log10(abs(solutionVc1));
absSolVc2 = 20*log10(abs(solutionVc2));

hf3 = figure();
plot(vecfreq1, absSol61, "b;abs(T_{v6});");
hold on
plot(vecfreq2, absSol62, "b");
plot(vecfreq1, absSolVc1, "g;abs(T_{vc});");
plot(vecfreq2, absSolVc2, "g");
plot(vecfreq1, absSolVs1, "r;abs(Vs);");
plot(vecfreq2, absSolVs2, "r");

xlabel("20*log(w) [rad/s]")
ylabel("|T| dB");
print(hf3, "absT.eps", "-depsc");

hf4 = figure();
plot(vecfreq1, argsSol61, "b;arg(T_{v6});");
hold on
plot(vecfreq2, argsSol62, "b");
plot(vecfreq1, argsSolVc1, "g;arg(T_{vc});");
plot(vecfreq2, argsSolVc2, "g");
plot(vecfreq1, argsSolVs1, "r;arg(Vs);");
plot(vecfreq2, argsSolVs2, "r");

xlabel("20*log(w) [rad/s]");
ylabel("Phase (degrees)")
print(hf4, "argT.eps", "-depsc");

hf5 = figure();
plot(vecfreq1, argsSolV61, "g");
hold on
plot(vecfreq2, argsSol62, "g");

print(hf5, "argV6.eps", "-depsc");

hf6 = figure()
plot(vecfreq1, absSol61, "b");
hold on
plot(vecfreq2, absSol62, "b");

print(hf6, "absV6.eps", "-depsc");


#}







