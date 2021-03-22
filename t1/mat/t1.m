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


printf("\nMesh Analysis------------------------------------------\n\n")

Am = [(R1+R3+R4)*(Kb*R3 - 1)/(Kb*R3) - R3, -R4; -R4*(Kb*R3 - 1)/(Kb*R3), R4 + R6 + R7 - Kc]
bm = [-Va; 0]

Ym = inv(Am)*bm
%[Ib Ic Vk]


Ia = (Kb*R3 - 1)*Ym(1)/(Kb*R3)
Ib = Ym(1)
Ic = Ym(2)
Id = Id;

printf("\n\n")

veci = [0, 0, 0, 0, 0];

veci(1) = Ym(1) - Ia;%3
veci(2) = Ym(2) - Ia;%4
veci(3) = Ym(1) - Id;%5
veci(4) = Ym(2);%6
veci(5) = Ym(2);%7

Vector = veci(1)

printf("\n\n")

vecv = [0, 0, 0, 0, 0, 0, 0, 0];

vecv(1) = Kc*Ic;
vecv(2) = -veci(3)*R5;
vecv(3) = -R3*veci(1);
vecv(4) = (Ia - Ic)*R4;
vecv(5) = vecv(4) - Ic*R6;
vecv(6) = vecv(3) - R2*Ib;
vecv(7) = vecv(3) + Ia*R1;
vecv(8) = vecv(4);


filename = "../doc/op_octave.tex";

file = fopen(filename, 'w');

exampleI = "I%d & %e\\\\ \\\hline\n";
exampleV = "V%d & %e\\\\ \\\hline\n";

for i = 1:5
fprintf(file, exampleI, i, veci(i))
end

for i= 1:8
fprintf(file, exampleV, i, vecv(i))
end










