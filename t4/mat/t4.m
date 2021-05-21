%Envelope detector--------------

t=linspace(2e-2, 12e-2, 10000);
f=50; %Hz
w=2*pi*f;
A = 20; %Volts


R = 7e3; %Ohm
C = 11e-6; %Farad


VON=12.021610/18
vlim =2*VON;

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

hold on

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
title("Output voltage Envelope Detector v(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc");

hold off

Vaverage = mean(vO)

%Voltage regulator--------------
% vii = 6.475163e+00;
% iir = vii/3000;


vii = 6.475163e+00;
iir = vii/3000; %Id

rd = (1.610751/iir)/18 %ohm determinada a partir do voltage regulator do spice
Rvr = 3000; %ohm


vIncrement = vO-Vaverage;
vri = 18*rd/(18*rd + Rvr)*vIncrement + 18*VON; %amplitude: uses 18 diodes in series

ripple = max(vri) - min(vri)

vriaverage = mean(vri)

cost = 18 + 22 * 0.1
vdev = mean(vri-12)

merit = 1/ cost*(ripple + vdev + 1e-6)

figure
plot(t*1000, vri)
title("Output voltage of Voltage Regulator v(t)")
xlabel ("t[ms]")
print ("voltageRegulator.eps", "-depsc");

vri = vri - 18*VON;
plot(t*1000, vri)
title("Output voltage Deviation v(t)")
xlabel ("t[ms]")
print ("Deviation.eps", "-depsc");

