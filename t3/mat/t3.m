

%Envelope Detector

Von = 0.75;













%Voltage regulator


%ripple que vem do voltage regulator
%Vs average voltage source of Envelope detector
Vs = 1;

rd = 40.1743; %ohm determinada a partir do voltage regulator do spice

Rvr = 3000; %ohm

vri = 16*rd/(16*rd + Rvr)*Vs; %amplitude: uses 16 diodes in series

%meu plot faço depois..........
%cos(wt) somente com a parte positiva??
