function q = pnCharge ( Uj,  Cj,  Vj, Mj,  Fc)
%Computes pn-junction depletion charge.

  if (Uj <= Fc * Vj),
    a = 1 - Uj / Vj;
    b = exp ((1 - Mj) * log (a));
    q = Cj * Vj / (1 - Mj) * (1 - b);
  else
    % Approach 1
    % a = 1 - Fc;
    % b = exp ((1 - Mj) * log (a));
    % a = exp ((1 + Mj) * log (a));
    % c = 1 - Fc * (1 + Mj);
    % d = Fc * Vj;
    % e = Vj * (1 - b) / (1 - Mj);
    % q = Cj * (e + (c * (Uj - d) + Mj / 2 / Vj * (sqr (Uj) - sqr (d))) / a);

    %Approach 2: this variant is numerically more stable
    a = 1 - Fc;
    b = exp (-Mj * log (a));
    f = Fc * Vj;
    c = Cj * (1 - Fc * (1 + Mj)) * b / a;
    d = Cj * Mj * b / a / Vj;
    e = Cj * Vj * (1 - a * b) / (1 - Mj) - d / 2 * f * f - f * c;
    q = e + Uj * (c + Uj * d / 2);
  end
end
