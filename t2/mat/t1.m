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


printf("\n alinea a)------------------------------------------\n\n")

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


bVa = Xn(5) - Xn(7)


printf("\n alinea b)-------------------------------------------\n\n")

#{

bline10 = [-1/R1, 0, -1/R4, 0, -1/R6, 0];
bline2 = [1/R1+1/R2+1/R3, -1/R2, -1/R3, 0, 0, 0];
#bline2 = [0, 0, -1, 0, Kd/R6, 1]
bline3 = [-Kb-1/R2, 1/R2, Kb, 0, 0, 0];
#bline7 = [d, 0, 1/R7, -1/Kd - 1/R7]
bline7 = [0, 0, 0, 0, Kd/R6+1/R7, -1/R8];
bline586 = [-1/R3 + Kb, 0, 1/R4 + 1/R3 - Kb, 0, -1/R7, 1/R7];
#bline586 = [0, 0, -1, 0, Kd/R6, 1]
blineV = [0, 0, 0, 1, 0, -1]

bAn = [bline10; bline2; bline3; bline7; bline586; blineV]

det(bAn)

bbn = [0; 0; 0; 0; 0; bVa]

bXn = bAn\bbn


mline1 = [1 - Kb*R3, Kb*R3, 0, 0]
mline2 = [0, -R4, R6 + R7 - Kd + R4, 0]
mline3 = [-R3, R3 + R1 + R4,  -R4, 0]
mline4 = [-R5, 0, Kb, 1]

meshA = [mline1; mline2; mline3; mline4]

det(meshA)

meshb = [0; 0; 0; -R5]

meshIn = meshA\meshb

#}

#evaluate equivalen resistance equals R5

print("\n\n alinea c)------------------------------------------")

time = 0:1e-2:20;

v6n = bVa*exp(-t/(R5*C));

hf = figure();
plot(time, v6n, "b");
xlabel ("t");
ylabel("v6n(t)");

print(hf, "v6n.eps", "-depsc");




hf = figure()



printf("\n\nalinea d)-------------------------------------------------")

omega = 2*pi*2*10e3 #rad/s

dline1 = [1, 0, 0, 0, 0, 0, 0];
dline2 = [0, 1/R2 + 1/R3+1/R2, -1/R2, -1/R3, 0, 0, 0];
dline58 = [0, -1/R3, 0, 1/R3+1/R4+1/R5, -i*omega*C-1/R5, -1/R7, 1/R7+i*omega*C];
dline3 = [0, Kb+1/R2, -1/R2, -Kb, 0, 0, 0];
dline6 = [0, Kb, 0, -1/R5-Kb, 1/R5+i*omega*C, 0, -i*omega*C];
dline7 = [0, 0, 0, 0, 0, -1/R6-1/R7, 1/R7];
dlineVd = [0, 0, 0, 1, 0, Kd/R6, -1];

dAn = [dline1; dline2; dline58; dline3; dline6; dline7; dlineVd]

dbn = [1; 1/R1; 0; 0; 0; 0; 0] #parte imaginaria

determinante = det (dAn)

dXn = dAn\dbn


#{

filename = "../doc/op_nodal_tab.tex";

file = fopen(filename, 'w');

exampleI = "$I_{R%d}$ & %e\\\\ \\\hline\n";
exampleV = "$V_%d$ & %e\\\\ \\\hline\n";

veci = [0, 0, 0, 0, 0, 0, 0];
vecv = [0, 0, 0, 0, 0, 0, 0];

vecv(1) = Xn(1);
vecv(2) = Xn(2);
vecv(3) = Xn(3);
vecv(4) = Xn(4);
vecv(5) = Xn(4)-R6*(Xn(4)-Xn(1))/(R6+R7);
vecv(6) = Xn(3) + R2*(Kb*Xn(3));
vecv(7) = Xn(4) + Va;

veci(1) = (vecv(3) - vecv(7))/R1;
veci(2) = (vecv(6) - vecv(3))/R2;
veci(3) = -vecv(3)/R3;
veci(4) = -vecv(4)/R4;
veci(5) = -vecv(2)/R5;
veci(6) = (vecv(4) - vecv(5))/R6;
veci(7) = (vecv(5) - vecv(1))/R7;

IVc = Id - veci(7);


for i = 1:7
fprintf(file, exampleI, i, veci(i))
end

fprintf(file, "$I_{Vc}$ & %e\\\\ \\\hline\n", IVc)
fprintf(file, "$I_d$ & %e\\\\ \\\hline\n", Id)

fprintf(file, "$V_0$ & 0.0\\\\ \\\hline\n")

for i= 1:7
fprintf(file, exampleV, i, vecv(i))
end

#}