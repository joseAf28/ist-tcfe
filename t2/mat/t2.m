close all
clear all

R1 = 1.04408633697 
R2 = 2.04051610808 
R3 = 3.07566747417 
R4 = 4.05936218175 
R5 = 3.05878343538 
R6 = 2.0603640429 
R7 = 1.04299566201 
Vs = 5.18382634375 
C = 1.02590436129 
Kb = 7.2865951329 
Kd = 8.22752594192 

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

time = 0:1e-2:20;

v6n = Xn(5)*exp(-time/(Req*C));

hf = figure();
plot(time, v6n, "b");
xlabel("time (ms)");
ylabel("v6n(t) (V)");

print(hf, "v6n.eps", "-depsc");



printf("\n\nalinea 4)-------------------------------------------------\n\n")

omega = 2*pi*10e3 #rad/s

dlineVs = [1, 0, 0, 0, 0, 0, 0];
dline2 = [-1/R1, 1/R1 + 1/R3+1/R2, -1/R2, -1/R3, 0, 0, 0];
dline58 = [0, -1/R3, 0, 1/R3+1/R4+1/R5, -i*omega*C-1/R5, -1/R7, 1/R7+i*omega*C];
dline3 = [0, Kb+1/R2, -1/R2, -Kb, 0, 0, 0];
dline6 = [0, Kb, 0, -1/R5-Kb, 1/R5+i*omega*C, 0, -i*omega*C];
dline7 = [0, 0, 0, 0, 0, -1/R6-1/R7, 1/R7];
dlineVd = [0, 0, 0, 1, 0, Kd/R6, -1];

dAn = [dline2; dline3; dline6; dline7; dlineVs; dlineVd; dline58]

dbn = [1; 0; 0; 0; 0; 0; 0] #parte imaginaria

determinante = det (dAn)

dXn = dAn\dbn

fase = arg(dXn(5))
amplitude = abs(dXn(5))

printf("\n\nalinea 5)-------------------------------------------------\n\n") #corrigir unidades de tempo
initial_time = time = -5:1e-2:0-1e-2;
final_time = 0:1e-2:20;

v6 = 0*initial_time+Xn(5);
