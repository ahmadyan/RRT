function c = pnCapacitance(  Uj,  Cj,  Vj, Mj,  Fc )
%pnCapacitance Computes pn-junction depletion capacitance.
  if (Uj <= Fc * Vj)
    c = Cj * exp (-Mj * log (1 - Uj / Vj));
  else
    c = Cj * exp (-Mj * log (1 - Fc)) * (1 + Mj * (Uj - Fc * Vj) / Vj / (1 - Fc));
  end
end

