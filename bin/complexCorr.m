
%% "complexCorr" does a complex correlation betwen two vector
%      time series.
%
function [cor, dir]=complexCorr(n,u1,v1,u2,v2)



u1bar = mean(u1);
v1bar = mean(v1);
u2bar = mean(u2);
v2bar = mean(v2);


c12 = complex(0,0);
speed1 = complex(0,0);
speed2 = complex(0,0);

for i=1:n
   w1 = complex(u1(i)-u1bar,v1(i)-v1bar);
   w2 = complex(u2(i)-u2bar,v2(i)-v2bar);
   c12 = c12+w1*conj(w2);
   speed1 = speed1 + w1*conj(w1);
   speed2 = speed2 + w2*conj(w2);
end

c12 = c12/sqrt(speed1*speed2);
creal = (c12+conj(c12))*0.5;
cimag = imag(c12);
cor = sqrt(creal*creal+cimag*cimag);
dir = atan2(cimag,creal)*180./pi;

