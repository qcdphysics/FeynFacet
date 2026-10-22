<|"Description" -> 
  "NLO hard coefficients of l q -> q + X (tagged quark), \
Processes/LeptonQuark", "Basis" -> "{regular, D0, D1, delta} with D0 = \
[1/(1-w)]_+, D1 = [log(1-w)/(1-w)]_+, delta = delta(1-w)", 
 "Scales" -> "LR = log(muR^2/s) (one common scale for renormalization, the \
incoming PDF and the fragmentation)", "Cross section" -> "sigma = s E \
dsigma/d^3p = Born delta(1-w) + Prefactor (c_regular + c_D0 D0 + c_D1 D1 + \
c_delta delta), distributions' coefficients with the prefactor at w = 1", 
 "Variables" -> "s = (l+k)^2, v = 1 + t/s, w = -u/(s+t); Unpolarized = spin \
average, Polarized = helicity difference", 
 "Results" -> <|"Unpolarized" -> 
    <|"Born" -> (2*alphaEM^2*eq^2*(1 + v^2))/(s*(-1 + v)^2*v), 
     "Prefactor" -> (alphaEM^2*alphaS*CF*eq^2)/(Pi*s*v*w), 
     "Coefficients" -> {-((2 - 4*v + 4*v^2 - 2*v^3 - w + 2*v*w - v^2*w + 
           2*v^4*w - 2*v^2*w^2 - 2*v^5*w^2 + v^2*w^3 - v^4*w^3 + 4*v^5*w^3 - 
           4*v^5*w^4 + 2*v^5*w^5)/((-1 + v)^2*(-1 + v*w)*(1 - v + v*w))) + 
        (LR*(2 - 2*v + 2*v^2 - 2*v^3 - w - 4*v*w + 5*v^2*w + 2*v^3*w + 
           2*v^4*w - 6*v^2*w^2 - 2*v^3*w^2 - 2*v^4*w^2 - 2*v^5*w^2 + 
           v^2*w^3 + 4*v^3*w^3 - 3*v^4*w^3 + 6*v^5*w^3 + 4*v^4*w^4 - 
           8*v^5*w^4 + 4*v^5*w^5))/((-1 + v)^2*(-1 + v*w)*(1 - v + v*w)) - 
        (2*(-1 - v^2 + 5*v*w - 2*v^2*w + v^3*w + 2*v^2*w^2 - 2*v^3*w^2 + 
           2*v^3*w^3)*Log[1 - v])/(-1 + v)^2 - 
        ((2 - 2*v + 2*v^2 - 2*v^3 - w - 4*v*w + 5*v^2*w + 2*v^3*w + 2*v^4*w - 
           6*v^2*w^2 - 2*v^3*w^2 - 2*v^4*w^2 - 2*v^5*w^2 + v^2*w^3 + 
           4*v^3*w^3 - 3*v^4*w^3 + 6*v^5*w^3 + 4*v^4*w^4 - 8*v^5*w^4 + 
           4*v^5*w^5)*Log[v])/((-1 + v)^2*(-1 + v*w)*(1 - v + v*w)) - 
        ((2 - 2*v + 2*v^2 - 2*v^3 - w - 8*v*w + 9*v^2*w + 2*v^3*w + 2*v^4*w - 
           10*v^2*w^2 - 2*v^3*w^2 - 2*v^4*w^2 - 2*v^5*w^2 + v^2*w^3 + 
           8*v^3*w^3 - 7*v^4*w^3 + 6*v^5*w^3 + 8*v^4*w^4 - 8*v^5*w^4 + 
           4*v^5*w^5)*Log[1 - w])/((-1 + v)^2*(-1 + v*w)*(1 - v + v*w)) + 
        (2*(1 + v^2 + 2*v^2*w^2)*Log[w])/(-1 + v)^2 + 
        (2*(-1 - v^2 + 3*v*w - 2*v^2*w + v^3*w - 2*v^3*w^2 + 2*v^3*w^3)*
          Log[1 - v*w])/(-1 + v)^2 + (2*v*w*Log[1 - v + v*w])/(1 - v + v*w), 
       (-4*LR*(1 + v^2))/(-1 + v)^2 + (4*(1 + v^2)*Log[1 - v])/(-1 + v)^2 + 
        (4*(1 + v^2)*Log[v])/(-1 + v)^2, (8*(1 + v^2))/(-1 + v)^2, 
       (-8*(1 + v^2))/(-1 + v)^2 + ((1 + v^2)*Log[v]^2)/(-1 + v)^2 + 
        LR*((-3*(1 + v^2))/(-1 + v)^2 - (2*(1 + v^2)*Log[v])/(-1 + v)^2) + 
        Log[1 - v]*((3*(1 + v^2))/(-1 + v)^2 + (2*(1 + v^2)*Log[v])/
           (-1 + v)^2)}|>, "Polarized" -> 
    <|"Born" -> (-2*alphaEM^2*eq^2*(1 + v))/(s*(-1 + v)*v), 
     "Prefactor" -> (alphaEM^2*alphaS*CF*eq^2)/(Pi*s*v*w), 
     "Coefficients" -> {((-1 + v)*(-1 + w)*(1 + v*w))/(1 - v + v*w) - 
        (LR*(-2 + 2*v^2 + w + 3*v*w - 10*v^2*w + 2*v^3*w + v*w^2 + 
           9*v^2*w^2 - 6*v^3*w^2 + 4*v^3*w^3))/((-1 + v)*(1 - v + v*w)) + 
        (2*(-1 - v + 3*v*w - v^2*w + 2*v^2*w^2)*Log[1 - v])/(-1 + v) + 
        ((-2 + 2*v^2 + w + 3*v*w - 10*v^2*w + 2*v^3*w + v*w^2 + 9*v^2*w^2 - 
           6*v^3*w^2 + 4*v^3*w^3)*Log[v])/((-1 + v)*(1 - v + v*w)) + 
        ((-2 + 2*v^2 + w + 7*v*w - 14*v^2*w + 2*v^3*w + v*w^2 + 13*v^2*w^2 - 
           6*v^3*w^2 + 4*v^3*w^3)*Log[1 - w])/((-1 + v)*(1 - v + v*w)) - 
        (2*(1 + v)*Log[w])/(-1 + v) - (2*(1 + v*w)*(-1 - v + 2*v*w)*
          Log[1 - v*w])/(-1 + v) - (2*v*w*Log[1 - v + v*w])/(1 - v + v*w), 
       (4*LR*(1 + v))/(-1 + v) - (4*(1 + v)*Log[1 - v])/(-1 + v) - 
        (4*(1 + v)*Log[v])/(-1 + v), (-8*(1 + v))/(-1 + v), 
       (8*(1 + v))/(-1 + v) - ((1 + v)*Log[v]^2)/(-1 + v) + 
        Log[1 - v]*((-3*(1 + v))/(-1 + v) - (2*(1 + v)*Log[v])/(-1 + v)) + 
        LR*((3*(1 + v))/(-1 + v) + (2*(1 + v)*Log[v])/(-1 + v))}|>|>|>
