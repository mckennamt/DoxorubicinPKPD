function dSdt = chemoCompartmentModel_drugOff(t,S,parms)


% extract parameters
% compartment model parameters
k12 = parms(1);
k21 = parms(2);
k23 = parms(3);

% drug compartments
dSdt(1) = k21*S(2) - k12*S(1);
dSdt(2) = k12*S(1) - k21*S(2) - k23*S(2);
dSdt(3) = k23 * S(2);%-.025*S(3);

dSdt = dSdt';