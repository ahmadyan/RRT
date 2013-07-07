function c =  pnCurrent ( Upn,  Iss,  Ute)
% The function computes the exponential pn-junction current.
    c = Iss * (exp (Upn / Ute) - 1);
end
