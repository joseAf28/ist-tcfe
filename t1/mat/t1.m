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

I3 = Ym(1) - Ia
I4 = Ym(2) - Ia
I5 = Ym(1) - Id
I6 = Ym(2)
I7 = Ym(2)

printf("\n\n")

V1 = Kc*Ic
V2 = -I5*R5
V3 = -R3*I3
V4 = (Ia - Ic)*R4
V5 = V4 - Ic*R6
V6 = V3 - R2*Ib
V7 = V3 + Ia*R1
V8 = V4




