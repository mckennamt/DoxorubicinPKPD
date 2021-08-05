function dSdt = chemoCompartmentModel_drugOn(t,S,parms,AIF,AIFt)


% extract parameters
% compartment model parameters
k12 = parms(1);
k21 = parms(2);
k23 = parms(3);

% drug compartments
AIFc = interp1(AIFt,AIF,t);
dSdt(1) = k12*AIFc - k21*S(1) - k23*S(1);
dSdt(2) = k23*S(1);

dSdt = dSdt';

