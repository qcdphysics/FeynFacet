<|"Description" -> "NLO hard coefficients of l q -> q + X with transverse \
spin transfer from the incoming to the tagged quark (Transverse), and the \
spin average (Unpolarized), Processes/Collins", 
 "Basis" -> "{regular, D0, D1, delta} with D0 = [1/(1-w)]_+, D1 = \
[log(1-w)/(1-w)]_+, delta = delta(1-w)", "Scales" -> "LR = log(muR^2/s); LI = \
log(muI^2/s), LO = log(muO^2/s), LL = log(muL^2/s) for the incoming quark, \
the outgoing quark and the photon-in-lepton subtraction", 
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
           (-1 + v)^2)}|>, "Transverse" -> 
    <|"Born" -> (4*alphaEM^2*eq^2)/(s*(-1 + v)^2), 
     "Prefactor" -> (alphaEM^2*alphaS*CF*eq^2)/(Pi*s*v*w), 
     "Coefficients" -> {(-2*v*w)/(-1 + v*w) + (4*LO*v*w*(1 + v*w))/
         (-1 + v)^2 + (2*LL*v*w*(1 + v^2 - 2*v*w - 2*v^2*w + 2*v^2*w^2))/
         ((-1 + v)^2*(-1 + v*w)) + (2*(-1 + 3*v - 3*v^2 + v^3 - 2*v*w + 
           2*v^2*w - 4*v^2*w^2)*Log[1 - v])/(-1 + v)^2 - 
        (2*v*w*(-1 + v^2 - 2*v*w - 2*v^2*w + 4*v^2*w^2)*Log[v])/
         ((-1 + v)^2*(-1 + v*w)) - (2*v*(1 - 2*v + v^2 - 4*w - v*w + 
           2*v^2*w - v^3*w + v*w^2 - 2*v^2*w^2 + v^3*w^2 + 4*v^2*w^3)*
          Log[1 - w])/((-1 + v)^2*(-1 + v*w)) + (2*v*(3 + v^2)*w*Log[w])/
         (-1 + v)^2 - (2*(-1 + 2*v - v^2 + v*w + 2*v^2*w + v^3*w - 4*v^2*w^2)*
          Log[1 - v*w])/(-1 + v)^2 + 2*(1 - v + v*w)*Log[1 - v + v*w], 
       (-4*LI*v)/(-1 + v)^2 - (4*LO*v)/(-1 + v)^2 + (8*v*Log[1 - v])/
         (-1 + v)^2 + (8*v*Log[v])/(-1 + v)^2, (16*v)/(-1 + v)^2, 
       (-16*v)/(-1 + v)^2 - (3*LI*v)/(-1 + v)^2 + (2*v*Log[v]^2)/(-1 + v)^2 + 
        LO*((-3*v)/(-1 + v)^2 - (4*v*Log[v])/(-1 + v)^2) + 
        Log[1 - v]*((6*v)/(-1 + v)^2 + (4*v*Log[v])/(-1 + v)^2)}|>|>, 
 "ResolvedBorn" -> 
  <|"Transfer" -> <|"q qp -> q" -> 
      <|"NN" -> -1/4*((-1 + Nc)*(1 + Nc)*(-4*s^2 - 4*s*t - 4*t^2 + D*t^2))/
          (Nc^2*t^2), "XX" -> ((-1 + Nc)*(1 + Nc)*(s + t)*
          (4*s^2 + 4*s*t - 4*t^2 + D*t^2)*(t + u))/(4*Nc^2*s*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*s^3*t + 8*s^2*t^2 + D*s*t^3 - 
           4*t^4 + D*t^4 + 8*s^3*u + 12*s^2*t*u + 4*s*t^2*u - 4*t^3*u + 
           D*t^3*u))/(8*Nc^2*s*t^2*u)|>, "q qpb -> q" -> 
      <|"NN" -> -1/4*((-1 + Nc)*(1 + Nc)*(-4*s^2 - 4*s*t - 4*t^2 + D*t^2))/
          (Nc^2*t^2), "XX" -> ((-1 + Nc)*(1 + Nc)*(s + t)*
          (4*s^2 + 4*s*t - 4*t^2 + D*t^2)*(t + u))/(4*Nc^2*s*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*s^3*t + 8*s^2*t^2 + D*s*t^3 - 
           4*t^4 + D*t^4 + 8*s^3*u + 12*s^2*t*u + 4*s*t^2*u - 4*t^3*u + 
           D*t^3*u))/(8*Nc^2*s*t^2*u)|>, "q q -> q" -> 
      <|"NN" -> -1/4*((-1 + Nc)*(1 + Nc)*(-4*Nc*s^3 + 4*s^2*t - 2*D*s^2*t - 
            8*Nc*s^2*t + 28*s*t^2 - 12*D*s*t^2 + D^2*s*t^2 - 8*Nc*s*t^2 + 
            D*Nc*s*t^2 + 16*t^3 - 8*D*t^3 + D^2*t^3 - 4*Nc*t^3 + D*Nc*t^3))/
          (Nc^3*t^2*(s + t)), "XX" -> ((-1 + Nc)*(1 + Nc)*
          (4*Nc*s^3 - 4*s^2*t + 2*D*s^2*t + 8*Nc*s^2*t + 20*s*t^2 - 
           8*D*s*t^2 + D^2*s*t^2 + D*Nc*s*t^2 + 16*t^3 - 8*D*t^3 + D^2*t^3 - 
           4*Nc*t^3 + D*Nc*t^3)*(t + u))/(4*Nc^3*s*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*Nc*s^4*t - 4*s^3*t^2 + 
           2*D*s^3*t^2 + 12*Nc*s^3*t^2 + 16*s^2*t^3 - 6*D*s^2*t^3 + 
           D^2*s^2*t^3 + 8*Nc*s^2*t^3 + D*Nc*s^2*t^3 + 36*s*t^4 - 
           16*D*s*t^4 + 2*D^2*s*t^4 - 4*Nc*s*t^4 + 2*D*Nc*s*t^4 + 16*t^5 - 
           8*D*t^5 + D^2*t^5 - 4*Nc*t^5 + D*Nc*t^5 + 8*Nc*s^4*u - 8*s^3*t*u + 
           4*D*s^3*t*u + 20*Nc*s^3*t*u - 12*s^2*t^2*u + 6*D*s^2*t^2*u + 
           16*Nc*s^2*t^2*u + 20*s*t^3*u - 8*D*s*t^3*u + D^2*s*t^3*u + 
           D*Nc*s*t^3*u + 16*t^4*u - 8*D*t^4*u + D^2*t^4*u - 4*Nc*t^4*u + 
           D*Nc*t^4*u))/(8*Nc^3*s*t^2*(s + t)*u)|>, 
     "q qb -> q" -> <|"NN" -> ((-1 + Nc)*(1 + Nc)*(4*Nc*s^3 + 4*s^2*t - 
           2*D*s^2*t + 4*Nc*s^2*t - 20*s*t^2 + 8*D*s*t^2 - D^2*s*t^2 + 
           4*Nc*s*t^2 - D*Nc*s*t^2 - 8*t^3 + 2*D*t^3))/(4*Nc^3*s*t^2), 
       "XX" -> ((-1 + Nc)*(1 + Nc)*(s + t)*(4*Nc*s^3 + 4*s^2*t - 2*D*s^2*t + 
           4*Nc*s^2*t + 28*s*t^2 - 12*D*s*t^2 + D^2*s*t^2 - 4*Nc*s*t^2 + 
           D*Nc*s*t^2 + 8*t^3 - 2*D*t^3)*(t + u))/(4*Nc^3*s^2*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*Nc*s^4*t + 4*s^3*t^2 - 
           2*D*s^3*t^2 + 8*Nc*s^3*t^2 + 32*s^2*t^3 - 14*D*s^2*t^3 + 
           D^2*s^2*t^3 + D*Nc*s^2*t^3 + 36*s*t^4 - 14*D*s*t^4 + D^2*s*t^4 - 
           4*Nc*s*t^4 + D*Nc*s*t^4 + 8*t^5 - 2*D*t^5 + 8*Nc*s^4*u + 
           8*s^3*t*u - 4*D*s^3*t*u + 12*Nc*s^3*t*u + 12*s^2*t^2*u - 
           6*D*s^2*t^2*u + 4*Nc*s^2*t^2*u + 28*s*t^3*u - 12*D*s*t^3*u + 
           D^2*s*t^3*u - 4*Nc*s*t^3*u + D*Nc*s*t^3*u + 8*t^4*u - 2*D*t^4*u))/
         (8*Nc^3*s^2*t^2*u)|>, "q qb -> qb" -> 
      <|"NN" -> ((-1 + Nc)*(1 + Nc)*(20*s^2 - 10*D*s^2 + D^2*s^2 + 20*s*t - 
           10*D*s*t + D^2*s*t + 8*t^2 - 2*D*t^2))/(4*Nc^3*s*(s + t)), 
       "XX" -> ((-1 + Nc)*(1 + Nc)*(20*s^2 - 10*D*s^2 + D^2*s^2 + 20*s*t - 
           10*D*s*t + D^2*s*t - 8*t^2 + 2*D*t^2)*(t + u))/(4*Nc^3*s^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(20*s^3*t - 10*D*s^3*t + D^2*s^3*t + 
           40*s^2*t^2 - 20*D*s^2*t^2 + 2*D^2*s^2*t^2 + 12*s*t^3 - 8*D*s*t^3 + 
           D^2*s*t^3 - 8*t^4 + 2*D*t^4 + 40*s^3*u - 20*D*s^3*u + 
           2*D^2*s^3*u + 60*s^2*t*u - 30*D*s^2*t*u + 3*D^2*s^2*t*u + 
           20*s*t^2*u - 10*D*s*t^2*u + D^2*s*t^2*u - 8*t^3*u + 2*D*t^3*u))/
         (8*Nc^3*s^2*(s + t)*u)|>, "q g -> q" -> 
      <|"NN" -> ((-8*s^2 + 4*D*s^2 - 8*s*t + 4*D*s*t + 16*t^2 - 8*D*t^2 + 
           D^2*t^2)*(2*Nc^2*s^2 + 2*Nc^2*s*t - t^2 + Nc^2*t^2))/
         (4*(-2 + D)*Nc^2*s*t^2*(s + t)), 
       "XX" -> ((-8*s^2 + 4*D*s^2 - 8*s*t + 4*D*s*t + 16*t^2 - 8*D*t^2 + 
           D^2*t^2)*(2*Nc^2*s^2 + 2*Nc^2*s*t - t^2 + Nc^2*t^2)*(t + u))/
         (4*(-2 + D)*Nc^2*s^2*t^2*u), "Collins" -> 
        ((-8*s^2 + 4*D*s^2 - 8*s*t + 4*D*s*t + 16*t^2 - 8*D*t^2 + D^2*t^2)*
          (2*Nc^2*s^2 + 2*Nc^2*s*t - t^2 + Nc^2*t^2)*(s*t + t^2 + 2*s*u + 
           t*u))/(8*(-2 + D)*Nc^2*s^2*t^2*(s + t)*u)|>, 
     "q qp -> qp" -> <|"NN" -> 0, "XX" -> 0, "Collins" -> 0|>, 
     "q qb -> qp" -> <|"NN" -> 0, "XX" -> 0, "Collins" -> 0|>, 
     "q qb -> qpb" -> <|"NN" -> 0, "XX" -> 0, "Collins" -> 0|>, 
     "qb qpb -> qb" -> <|"NN" -> -1/4*((-1 + Nc)*(1 + Nc)*(-4*s^2 - 4*s*t - 
            4*t^2 + D*t^2))/(Nc^2*t^2), 
       "XX" -> ((-1 + Nc)*(1 + Nc)*(s + t)*(4*s^2 + 4*s*t - 4*t^2 + D*t^2)*
          (t + u))/(4*Nc^2*s*t^2*u), "Collins" -> 
        ((-1 + Nc)*(1 + Nc)*(4*s^3*t + 8*s^2*t^2 + D*s*t^3 - 4*t^4 + D*t^4 + 
           8*s^3*u + 12*s^2*t*u + 4*s*t^2*u - 4*t^3*u + D*t^3*u))/
         (8*Nc^2*s*t^2*u)|>, "qb qp -> qb" -> 
      <|"NN" -> -1/4*((-1 + Nc)*(1 + Nc)*(-4*s^2 - 4*s*t - 4*t^2 + D*t^2))/
          (Nc^2*t^2), "XX" -> ((-1 + Nc)*(1 + Nc)*(s + t)*
          (4*s^2 + 4*s*t - 4*t^2 + D*t^2)*(t + u))/(4*Nc^2*s*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*s^3*t + 8*s^2*t^2 + D*s*t^3 - 
           4*t^4 + D*t^4 + 8*s^3*u + 12*s^2*t*u + 4*s*t^2*u - 4*t^3*u + 
           D*t^3*u))/(8*Nc^2*s*t^2*u)|>, "qb qb -> qb" -> 
      <|"NN" -> -1/4*((-1 + Nc)*(1 + Nc)*(-4*Nc*s^3 + 4*s^2*t - 2*D*s^2*t - 
            8*Nc*s^2*t + 28*s*t^2 - 12*D*s*t^2 + D^2*s*t^2 - 8*Nc*s*t^2 + 
            D*Nc*s*t^2 + 16*t^3 - 8*D*t^3 + D^2*t^3 - 4*Nc*t^3 + D*Nc*t^3))/
          (Nc^3*t^2*(s + t)), "XX" -> ((-1 + Nc)*(1 + Nc)*
          (4*Nc*s^3 - 4*s^2*t + 2*D*s^2*t + 8*Nc*s^2*t + 20*s*t^2 - 
           8*D*s*t^2 + D^2*s*t^2 + D*Nc*s*t^2 + 16*t^3 - 8*D*t^3 + D^2*t^3 - 
           4*Nc*t^3 + D*Nc*t^3)*(t + u))/(4*Nc^3*s*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*Nc*s^4*t - 4*s^3*t^2 + 
           2*D*s^3*t^2 + 12*Nc*s^3*t^2 + 16*s^2*t^3 - 6*D*s^2*t^3 + 
           D^2*s^2*t^3 + 8*Nc*s^2*t^3 + D*Nc*s^2*t^3 + 36*s*t^4 - 
           16*D*s*t^4 + 2*D^2*s*t^4 - 4*Nc*s*t^4 + 2*D*Nc*s*t^4 + 16*t^5 - 
           8*D*t^5 + D^2*t^5 - 4*Nc*t^5 + D*Nc*t^5 + 8*Nc*s^4*u - 8*s^3*t*u + 
           4*D*s^3*t*u + 20*Nc*s^3*t*u - 12*s^2*t^2*u + 6*D*s^2*t^2*u + 
           16*Nc*s^2*t^2*u + 20*s*t^3*u - 8*D*s*t^3*u + D^2*s*t^3*u + 
           D*Nc*s*t^3*u + 16*t^4*u - 8*D*t^4*u + D^2*t^4*u - 4*Nc*t^4*u + 
           D*Nc*t^4*u))/(8*Nc^3*s*t^2*(s + t)*u)|>, 
     "qb q -> qb" -> <|"NN" -> ((-1 + Nc)*(1 + Nc)*(4*Nc*s^3 + 4*s^2*t - 
           2*D*s^2*t + 4*Nc*s^2*t - 20*s*t^2 + 8*D*s*t^2 - D^2*s*t^2 + 
           4*Nc*s*t^2 - D*Nc*s*t^2 - 8*t^3 + 2*D*t^3))/(4*Nc^3*s*t^2), 
       "XX" -> ((-1 + Nc)*(1 + Nc)*(s + t)*(4*Nc*s^3 + 4*s^2*t - 2*D*s^2*t + 
           4*Nc*s^2*t + 28*s*t^2 - 12*D*s*t^2 + D^2*s*t^2 - 4*Nc*s*t^2 + 
           D*Nc*s*t^2 + 8*t^3 - 2*D*t^3)*(t + u))/(4*Nc^3*s^2*t^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(4*Nc*s^4*t + 4*s^3*t^2 - 
           2*D*s^3*t^2 + 8*Nc*s^3*t^2 + 32*s^2*t^3 - 14*D*s^2*t^3 + 
           D^2*s^2*t^3 + D*Nc*s^2*t^3 + 36*s*t^4 - 14*D*s*t^4 + D^2*s*t^4 - 
           4*Nc*s*t^4 + D*Nc*s*t^4 + 8*t^5 - 2*D*t^5 + 8*Nc*s^4*u + 
           8*s^3*t*u - 4*D*s^3*t*u + 12*Nc*s^3*t*u + 12*s^2*t^2*u - 
           6*D*s^2*t^2*u + 4*Nc*s^2*t^2*u + 28*s*t^3*u - 12*D*s*t^3*u + 
           D^2*s*t^3*u - 4*Nc*s*t^3*u + D*Nc*s*t^3*u + 8*t^4*u - 2*D*t^4*u))/
         (8*Nc^3*s^2*t^2*u)|>, "qb q -> q" -> 
      <|"NN" -> ((-1 + Nc)*(1 + Nc)*(20*s^2 - 10*D*s^2 + D^2*s^2 + 20*s*t - 
           10*D*s*t + D^2*s*t + 8*t^2 - 2*D*t^2))/(4*Nc^3*s*(s + t)), 
       "XX" -> ((-1 + Nc)*(1 + Nc)*(20*s^2 - 10*D*s^2 + D^2*s^2 + 20*s*t - 
           10*D*s*t + D^2*s*t - 8*t^2 + 2*D*t^2)*(t + u))/(4*Nc^3*s^2*u), 
       "Collins" -> ((-1 + Nc)*(1 + Nc)*(20*s^3*t - 10*D*s^3*t + D^2*s^3*t + 
           40*s^2*t^2 - 20*D*s^2*t^2 + 2*D^2*s^2*t^2 + 12*s*t^3 - 8*D*s*t^3 + 
           D^2*s*t^3 - 8*t^4 + 2*D*t^4 + 40*s^3*u - 20*D*s^3*u + 
           2*D^2*s^3*u + 60*s^2*t*u - 30*D*s^2*t*u + 3*D^2*s^2*t*u + 
           20*s*t^2*u - 10*D*s*t^2*u + D^2*s*t^2*u - 8*t^3*u + 2*D*t^3*u))/
         (8*Nc^3*s^2*(s + t)*u)|>, "qb g -> qb" -> 
      <|"NN" -> ((-8*s^2 + 4*D*s^2 - 8*s*t + 4*D*s*t + 16*t^2 - 8*D*t^2 + 
           D^2*t^2)*(2*Nc^2*s^2 + 2*Nc^2*s*t - t^2 + Nc^2*t^2))/
         (4*(-2 + D)*Nc^2*s*t^2*(s + t)), 
       "XX" -> ((-8*s^2 + 4*D*s^2 - 8*s*t + 4*D*s*t + 16*t^2 - 8*D*t^2 + 
           D^2*t^2)*(2*Nc^2*s^2 + 2*Nc^2*s*t - t^2 + Nc^2*t^2)*(t + u))/
         (4*(-2 + D)*Nc^2*s^2*t^2*u), "Collins" -> 
        ((-8*s^2 + 4*D*s^2 - 8*s*t + 4*D*s*t + 16*t^2 - 8*D*t^2 + D^2*t^2)*
          (2*Nc^2*s^2 + 2*Nc^2*s*t - t^2 + Nc^2*t^2)*(s*t + t^2 + 2*s*u + 
           t*u))/(8*(-2 + D)*Nc^2*s^2*t^2*(s + t)*u)|>, 
     "qb qpb -> qpb" -> <|"NN" -> 0, "XX" -> 0, "Collins" -> 0|>, 
     "qb q -> qpb" -> <|"NN" -> 0, "XX" -> 0, "Collins" -> 0|>, 
     "qb q -> qp" -> <|"NN" -> 0, "XX" -> 0, "Collins" -> 0|>|>, 
   "Channels" -> {"q qp -> q", "q qpb -> q", "q q -> q", "q qb -> q", 
     "q qb -> qb", "q g -> q", "q qp -> qp", "q qb -> qp", "q qb -> qpb", 
     "qb qpb -> qb", "qb qp -> qb", "qb qb -> qb", "qb q -> qb", "qb q -> q", 
     "qb g -> qb", "qb qpb -> qpb", "qb q -> qpb", "qb q -> qp"}, 
   "Variables" -> "s, t (u = -s-t), D, Nc; divided by gs^4; first incoming \
parton polarized, tagged parton c"|>|>
