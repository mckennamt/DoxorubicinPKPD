function dCidt = twoPeakDoxModel(t,C, parms, inputFcn, tSpan)

cVal = interp1(tSpan, inputFcn, t);

% extract Parameters
k1 = parms(1);
k2 = parms(2);
Ki = parms(3);
k3 = parms(4);

dCidt = (k1*cVal + k2*cVal/(Ki+cVal) - k3*C(1));

end