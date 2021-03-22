% script para experimentar o octave

R1 = 1.04408633697 
R2 = 2.04051610808 
R3 = 3.07566747417 
R4 = 4.05936218175 
R5 = 3.05878343538 
R6 = 2.0603640429 
R7 = 1.04299566201 
Va = 5.18382634375 
Id = 1.02590436129 
Kb = 7.2865951329 
Kc = 8.22752594192 

printf("Nodal Analysis---------------")

A = [0, 1/R5, Kb, 0; 0, 0, Kb-1/R1-1/R3, 1/R1; 1/(R6+R7), 0, 1/R1, -1/R4-1/R1-1/(R6+R7); -1/Kc, 1/R5, 1/R3, 1/R4]

b = [Id; -Va/R1; Va/R1; Id]

Xn = A\b

