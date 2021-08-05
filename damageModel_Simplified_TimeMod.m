

function dSdt = damageModel_Simplified_TimeMod(t,S,parms,...
    fixedParms,maxDrugConc, drugVec, drugTime)


% extract parameters

% fixed parameters
kp = fixedParms(1);
theta = fixedParms(2);


if fixedParms(3) == 1
    % stretch Cb
    % normalize drug vector to maximum if we want to follow Cb curve
    %drugConc = interp1(drugTime,drugVec,t/delay);
    %
    %k_tot = kp + (-kd-kp) * (drugConc./maxDrugConc);
    %
    %if k_tot<-kd
    %    k_tot = -kd;
    %end
    
    kd = parms(1);
    theta = parms(2);
    k_tot = (kp - kd).*(t>=0) + ...
        kp.*(t<0);
    
elseif fixedParms(3) == 2
    % shift Cb
    % normalize drug vector to maximum if we want to follow Cb curve
    %if t>=delay
    %    drugConc = interp1(drugTime,drugVec,(t-delay));
    %else
    %    drugConc = 0;
    %end
    %
    %k_tot = kp + (-kd-kp) * (drugConc./maxDrugConc);
    %
    %if k_tot<-kd
    %    k_tot = -kd;
    %end
    
    % ktot parameters
    kd = parms(1);
    delay = parms(2);
    delay2 = parms(3);
    
    % linear transition
    k_tot = (kp + (-kd-kp)/delay * t) * (t>=0 & t<=delay) + ...
        kp  * (t<0) + ...
        (-kd + (kd+kp)/(delay2-delay) * (t-delay)) * (t>delay & t<=delay2) + ...
        kp * (t>delay2);
    %-kd * (t>delay & t<delay2) + ...
    
elseif fixedParms(3) == 3
    
    % ktot parameters
    kd = parms(1);
    delay = parms(2);
    delay2 = parms(3);
    
    % 'exponential' transitions
    rec = delay2*(t);%.*(t>d2);
    
    k_tot = (kp - ((kp+kd).*(1-exp(-delay.*t))).*(exp(-rec))) .* (t>=0) + ...
        kp .* (t<0);
    
elseif fixedParms(3) == 4
    % biexponential model
    kd = parms(1);
    r = parms(2);
    r2 = r/parms(3);
    theta = parms(4);
    
    %%biexponential function
    %k_tot = (kp - kd.*(exp(-ttr.*t) - exp(-ttd.*t))).* (t>=0) + ...
    %    kp .* (t<0);
    
    maxTau = (r-r2)^-1*log(r/r2);
    normMe = exp(-r*maxTau)-exp(-r2*maxTau);
    
    if r~=r2 && abs(normMe) >1e-5
        k_tot = (kp - kd/normMe.*(exp(-r.*t) - exp(-r2.*t))).* (t>=0) + ...
            kp .* (t<0);
    else
        maxTau = 1/r;
        normMe = r^2*maxTau*exp(-r*maxTau);
        k_tot = (kp - kd/normMe.*r^2.*t.*exp(-r.*t)).* (t>=0) + ...
            kp .* (t<0);
    end
    
    
elseif fixedParms(3) == 5
    % monoexponential model
    kd = parms(1);
    r = parms(2);
    theta = parms(3);
    
    %k_tot = (kp - (kd*r^2.*t.*exp(-r.*t))).* (t>0) + ...
    %    kp .* (t<=0);
    %k_tot = (kp - (kd.*t.*exp(-r.*t))).* (t>0) + ...
    %    kp .* (t<=0);
    
    maxTau = 1/r;
	normMe = maxTau*exp(-r*maxTau);
    k_tot = (kp - (kd/normMe*t.*exp(-r.*t))).* (t>0) + ...
        kp .* (t<=0);
    
elseif fixedParms(3) == 6
    % exponential transitions
    kd = parms(1);
    r = parms(2);
    theta = parms(3);
    
    tau = log(2)/r;
    k_tot = (kp - (kd*(exp(r.*t)-1))).* (t>0 & t<=tau) + ...
        (kp - (kd*exp(-r.*(t-tau)))).* (t>tau) + ...
        kp .* (t<=0);
    
elseif fixedParms(3) == 7
    % log transition
    kd = parms(1);
    r = parms(2);
    tau = parms(3);
    theta = parms(4);
    
    c = (1-exp(-r*tau));
    k_tot = (kp - (kd*(1-exp(-r.*t)))).* (t>0 & t<=tau) + ...
        (kp - (kd*c*exp(-r.*(t-tau)))).* (t>tau) + ...
        kp .* (t<=0);
    
elseif fixedParms(3) == 8
    % exponential transitions
    kd = parms(1);
    r = parms(2);
    theta = parms(3);
    
    %tau = log(2)/r;
    %k_tot = (kp - (kd*(exp(r.*t)-1))).* (t>0 & t<=tau) + ...
    %    (kp - (kd*exp(-r2.*(t-tau)))).* (t>tau) + ...
    %    kp .* (t<=0);
    k_tot = (kp - (kd*exp(-r.*(t)))).* (t>0) + ...
        kp .* (t<=0);
    
end

dSdt(1) = k_tot*S(1)*(1-S(1)/theta);

dSdt = dSdt';

