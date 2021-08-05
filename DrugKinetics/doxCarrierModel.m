

function dSdt = doxCarrierModel(t, S, parms, inputFcn, Tspan)

inputVal = interp1(Tspan, inputFcn, t);

% extract parameters
k12 = parms(1);
k21 = parms(2);
k23 = parms(3);
%kc = parms(2);
K_half = parms(4);


% diffusion into cell, carrier into cell, diffusion out
%dSdt(1) = k12*inputVal + kc*inputVal./(K_half + inputVal) - k21*S(1) - k23*S(1);
dSdt(1) = k12*inputVal./(K_half + inputVal) - k21*S(1) - k23*S(1);
dSdt(2) = k23*S(1);

dSdt = dSdt';
