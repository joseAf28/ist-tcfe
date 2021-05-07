

%Envelope Detector

% f=50; %Hz
% w=2*pi*f;
% R=1e3 %ohm
% C=8e-6 %Farads


% t=linspace(0, 5e-2, 10000);
% A = 34.5; %Volt
% vS=A*cos(w*t);
% vO = abs(vS);

% VON=0.75
% vlim =3*VON

% for i=1:length(t)
%   if vO(i) > vlim
%     vO(i) = vO(i)-vlim;
%   else
%     vO(i) = 0;
%   endif
% endfor

% figure
% plot(t*1000, vO)
% title("Output voltage v_o(t)")
% xlabel ("t[ms]")
% ylabel ("v_o[V]")
% print ("vo.eps", "-depsc");

%envelope detector
t=linspace(0, 10e-2, 10000);
f=50; %Hz
w=2*pi*f;
A = 34.5; %Volts
R = 1e3; %Ohm
C = 8e-6; %Farad

VON=0.75;
vlim =3*VON;

B = A - vlim;


vS = A * cos(w*t);
vS0 = abs(vS);
vOhr = zeros(1, length(t));
vO = zeros(1, length(t));

figure
for i=1:length(t)
  if (vS0(i) > vlim)
    vOhr(i) = vS0(i)-vlim;
  else
    vOhr(i) = 0;
  endif
endfor

plot(t*1000, vOhr)
hold

counter = 0;
tOFF = 1/w * atan(1/w/R/C);
texp = linspace(0, 1e-2, 1000);


function f = f(B, w, tOFF, tstart, R, C)
f = B*cos(w*tOFF)*exp(-(tstart - tOFF)/R/C);
endfunction

for(j = 0:9)
for i=1:length(texp)
k = i + j*length(texp);
    valueexp = f(B, w, tOFF, texp(i), R, C);
  if texp(i) < tOFF
    vO(k) = vOhr(i);
  elseif valueexp > vOhr(i)
    vO(k) = valueexp;
  else 
    vO(k) = vOhr(i);
  endif
endfor
endfor

plot(t*1000, vO)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc");



Vaverage = mean(vO)

vIncrement = vO-Vaverage;

plot(t*1000, vIncrement)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc");



%Voltage regulator


%ripple que vem do voltage regulator
%Vs average voltage source of Envelope detector

rd = 40.1743; %ohm determinada a partir do voltage regulator do spice

Rvr = 3000; %ohm


vri = 16*rd/(16*rd + Rvr)*vIncrement + 0.75*16; %amplitude: uses 16 diodes in series

ripple = max(vri) - min(vri)

vriaverage = mean(vri)

plot(t*1000, vri)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("voltageRegulator.eps", "-depsc");

%meu plot faço depois..........
%cos(wt) somente com a parte positiva??


