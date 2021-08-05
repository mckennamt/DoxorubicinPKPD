function Ct = modelCt2C2P(kpars,Cp,t)

% function to calculate the concentration in the extravascular,
% intracelluar space as a function of the kinetic parameters, 
% the concentration in the BLOOD and the time 

% kpars = [K1 k2]

% 3/14/11 - JUF

k2 = kpars(2);

f1 = kpars(1)*(exp(-k2*t));

Ct = conv(f1*mean(diff(t)),Cp);

Ct = Ct(1:length(t));