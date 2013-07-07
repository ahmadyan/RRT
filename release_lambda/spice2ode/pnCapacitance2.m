function c = pnCapacitance2 ( Uj, Cj, Vj, Mj)
% This function computes the pn-junction depletion capacitance with
% no linearization factor given.
    if (Uj <= 0)
        c = Cj * exp (-Mj * log (1 - Uj / Vj));
    else
        c = Cj * (1 + Mj * Uj / Vj);
    end
end