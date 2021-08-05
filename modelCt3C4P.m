function Ct = modelCt3C4P(kpars,Cp,t)

% function to calculate the concentration in the extravascular,
% intracelluar space as a function of the kinetic parameters, 
% the concentration in the BLOOD and the time 

% kpars = [K1 Vd k3 k4]

% 3/14/11 - JUF

k2 = kpars(1)/kpars(2); % 6/24/11 - JGW

alpha1=(k2+kpars(3)+kpars(4)-sqrt(((k2+kpars(3)+kpars(4)).^2)-(4*k2*kpars(4))))/2;
alpha2=(k2+kpars(3)+kpars(4)+sqrt(((k2+kpars(3)+kpars(4)).^2)-(4*k2*kpars(4))))/2;

f1 = (kpars(1)/(alpha2-alpha1))*(((kpars(3)+kpars(4)-alpha1)*exp(-alpha1*t))+((alpha2-kpars(3)-kpars(4))*exp(-alpha2*t)));
Ct = conv(f1*mean(diff(t)),Cp);

Ct = Ct(1:length(t));