function M = mhgOptimalM(kappa)
% compute the minimum M to get the correct
% mhg estimation

% for kappa <=8
forLowKappa=[1 16 21 25 28 31 34 37 40];

if kappa <=8
    M = forLowKappa(ceil(kappa)+1);
elseif kappa <=32
    M = ceil(-0.015*kappa^2 + 2.5*kappa + 22);
else
    M = ceil(-0.0041*kappa^2 + 2*kappa + 28);
    if M > 140
        M = 140;
    end 
end

