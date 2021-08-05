

function [N,J] = damageModel_analytic(N0, t, parms, fixedParms)


% extract parameters

% fixed parameters
kp = fixedParms(1);
theta = fixedParms(2);


if fixedParms(3) == 1
    
    kd = parms(1);
    theta = parms(2);
    
    k_tot = kp - kd;
    
    if any(t<=0)
        t_pre = t-min(t);
        N_preTreat = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t_pre));
        
        % get number of cells at t=0 for post-treatment initialization
        t_eval = 0-min(t);
        N0 = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t_eval));
        t_trans = t;
        
    else
        N_preTreat = zeros(size(t));
        
        t_eval = min(t);
        N0 = theta*N0 ./ (N0 + (theta-N0).*exp(-k_tot.*0));
        t_trans = t - t_eval;
        
    end
    
    N_treat = theta*N0 ./ (N0 + (theta-N0).*exp(-k_tot.*t_trans));
    
    N = N_treat.*(t>0) + N_preTreat.*(t<=0);
    
    if nargout > 1
        J = zeros(length(t), length(parms));
    end
    
elseif fixedParms(3) == 5
    % monoexponential model
    kd = parms(1);
    r = parms(2);
    theta = parms(3);
    
    if any(t<=0)
        t_pre = t-min(t);
        N_preTreat = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t_pre));
        
        % get number of cells at t=0 for post-treatment initialization
        t_eval = 0-min(t);
        N0 = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t_eval));
        t_eval = 0;
        
    else
        N_preTreat = zeros(size(t));
        t_eval = min(t);
        %d0 = t(1).*(kp+kd.*exp(1-r.*t(1))) + kd/r.*exp(1-r.*t(1));
        
    end
    
    d0 = t_eval.*(kp+kd.*exp(1-r.*t_eval)) + kd/r.*exp(1-r.*t_eval);
    C = (theta*N0 / (theta-N0))*exp(-d0);
    
    d = t.*(kp+kd.*exp(1-r.*t)) + kd/r.*exp(1-r.*t);
    N_treat = C.*exp(d)./(1+C.*exp(d)./theta);
    
    N = N_treat.*(t>0) + N_preTreat.*(t<=0);
    
    if nargout > 1
        J = zeros(length(t), length(parms));
        
        if kd<1e-8
            kd = 1e-4;
        end
        
        % Jacobian
        alpha = (1+C.*exp(d)./theta).^-1;
        
        dCdkd = exp(1-exp(1)*kd/r)*N0*theta./(N0*r-r*theta);
        dCdr = exp(1-exp(1)*kd/r)*kd*N0*theta./(r^2*(theta-N0));
        dCdtheta = -exp(-kd*exp(1)/r) * N0.^2 ./((N0-theta)^2);
        
        dddkd = exp(1-r.*t).*(1+r.*t)./r;
        dddr = -exp(1-r.*t).*kd.*(1+r.*t.*(1+r.*t))./r.^2;
        dddtheta = 0;
        
        %dNdkd
        J(:,1) = dCdkd.*exp(d).*alpha + ...
            C.*alpha.*exp(d).*dddkd +...
            -C.*exp(d).*alpha.^2.*(dCdkd.*exp(d) + C.*exp(d).*dddkd).*theta^-1;
        
        %dNdr
        J(:,2) = dCdr.*exp(d).*alpha + ...
            C.*alpha.*exp(d).*dddr +...
            -C.*exp(d).*alpha.^2.*(dCdr.*exp(d) + C.*exp(d).*dddr).*theta^-1;
        
        %dNdtheta
        J(:,3) = dCdtheta.*exp(d).*alpha + ...
            C.*alpha.*exp(d).*dddtheta +...
            -C.*exp(d).*alpha.^2.*(dCdtheta.*exp(d)./theta + C.*exp(d).*dddtheta./theta - ...
            C.*exp(d)./(theta^2));
        
    end
    
    
elseif fixedParms(3) == 2
    % shift Cb
    % normalize drug vector to maximum if we want to follow Cb curve
    % No analytic solution (yet)
    
    % ktot parameters
    kd = parms(1);
    delay = parms(2);
    delay2 = parms(3);
    
    % linear transition
    k_tot = (kp + (-kd-kp)/delay * t) * (t>=0 & t<=delay) + ...
        kp  * (t<0) + ...
        (-kd + (kd+kp)/(delay2-delay) * (t-delay)) * (t>delay & t<=delay2) + ...
        kp * (t>delay2);
    
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

