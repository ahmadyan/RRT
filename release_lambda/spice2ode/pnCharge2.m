function q = pnCharge2(Uj, Cj, Vj, Mj)
% This function computes the pn-junction depletion charge with no
%   linearization factor given.
  if (Uj <= 0)
    q = Cj * Vj / (1 - Mj) * (1 - exp ((1 - Mj) * log (1 - Uj / Vj)));
  else
    q = Cj * Uj * (1 + Mj * Uj / 2 / Vj);
  end
end