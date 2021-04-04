close all
clear all

R1 = 1.04408633697e3
R2 = 2.04051610808e3
R3 = 3.07566747417e3
R4 = 4.05936218175e3 
R5 = 3.05878343538e3
R6 = 2.0603640429e3
R7 = 1.04299566201e3
Vs = 5.18382634375 
C = 1.02590436129e-6
Kb = 7.28659513293e-3
Kd = 8.22752594192e3

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
xlabel("time (ms)");
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

printf("\n\nalinea 5)-------------------------------------------------\n\n")

stepDt = fase/omega

initial_time  = -5e-3:1e-6:0;
final_time = 0:1e-6:20e-3;

v6i = 0*initial_time+Xn(5);
vsi = 0*initial_time + Vs;

v6f = v6n + amplitude*cos(omega*final_time - fase);
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

#{
printf("\n\nalinea 6)-----------------------------------------------------\n\n")

solution6 = 0.1:1:1e4;
kk = 1

for freq = 0.1:1:1e4
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
solution6(kk) = dXn(5) - dXn(7);
kk = kk +1;

end 

solution6 = solution6/(-i);
argsSol6 = arg(solution6)*(180/pi);

vecfreq = 10*log10(0.1:1:1e4);
absSol6 = 20*log10(abs(solution6));

hf3 = figure();
plot(vecfreq, absSol6, "b;absFrq(t);");
print(hf3, "absv6.eps", "-depsc");

hf4 = figure();
plot(vecfreq, argsSol6, "b;argV6;");
print(hf4, "argvc.eps", "-depsc");


#}







