function c = pnConductance ( Upn, Iss, Ute)
% The function computes the exponential pn-junction current's derivative.
c = Iss * exp (Upn / Ute) / Ute;
end

