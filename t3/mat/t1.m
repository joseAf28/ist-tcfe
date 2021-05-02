close all
clear all

R1 = 1.04408633697e3 
R2 = 2.04051610808e3
R3 = 3.07566747417e3
R4 = 4.05936218175e3
R5 = 3.05878343538e3
R6 = 2.0603640429e3
R7 = 1.04299566201e3
Va = 5.18382634375 
Id = 1.02590436129e-3
Kb = 7.2865951329e-3
Kc = 8.22752594192e3

format long

printf("\nNodal Analysis------------------------------------------\n\n")

An = [0, 1/R5, Kb, 0; 0, 0, Kb-1/R1-1/R3, 1/R1; 1/(R6+R7), 0, 1/R1, -1/R4-1/R1-1/(R6+R7); -1/Kc, 1/R5, 1/R3, 1/R4]

bn = [Id; -Va/R1; Va/R1; Id]

Xn = An\bn

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

printf("\nMesh Analysis------------------------------------------\n\n")

Am = [(R1+R3+R4)*(Kb*R3 - 1)/(Kb*R3) - R3, -R4; -R4*(Kb*R3 - 1)/(Kb*R3), R4 + R6 + R7 - Kc]
bm = [-Va; 0]

Ym = inv(Am)*bm
%[Ib Ic]


Ia = (Kb*R3 - 1)*Ym(1)/(Kb*R3)
Ib = Ym(1)
Ic = Ym(2)
Id = Id;

printf("\n\n")

veci = [0, 0, 0, 0, 0, 0, 0];

veci(1) = Ia;%R1
veci(2) = Ib;%R2
veci(3) = Ia - Ib;%R3
veci(4) = Ic - Ia;%R4
veci(5) = Ib - Id;%R5
veci(6) = Ic;%R6
veci(7) = Ic;%R7

IVc = Id - Ic;

printf("\n\n")

vecv = [0, 0, 0, 0, 0, 0, 0];

vecv(1) = -Kc*Ic;
vecv(2) = -veci(5)*R5;
vecv(3) = -R3*veci(3);
vecv(4) = (Ia - Ic)*R4;
vecv(5) = vecv(4) - Ic*R6;
vecv(6) = vecv(3) + R2*Ib;
vecv(7) = vecv(3) - Ia*R1;


filename = "../doc/op_mesh_tab.tex";

file = fopen(filename, 'w');

for i = 1:7
fprintf(file, exampleI, i, veci(i))
end

fprintf(file, "$I_{Vc}$ & %e\\\\ \\\hline\n", IVc)
fprintf(file, "$I_d$ & %e\\\\ \\\hline\n", Id)

fprintf(file, "$V_0$ & 0.0\\\\ \\\hline\n")

for i= 1:7
fprintf(file, exampleV, i, vecv(i))
end










