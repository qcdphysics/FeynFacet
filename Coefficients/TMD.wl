<|"Description" -> "NLO TMD beam functions, TMD PDFs and soft functions, \
Processes/TMD (notation in its README: as = alphaS/(4 Pi), Lp, Lv, LQ, Lh)", 
 "Beam" -> <|"qq" -> 2*as*CF*(1 + Lp - x + Lp*x) - 
     2*as*CF*Lp*LQ*delta[1 - x] - 4*as*CF*Lp*plusD[0, 1 - x], 
   "qg" -> -(as*(2*(-1 + x)*x + Lp*(1 - 2*x + 2*x^2))), 
   "gq" -> 2*as*CF*x - (2*as*CF*Lp*(2 - 2*x + x^2))/x, 
   "gg" -> (4*as*CA*Lp*(-1 + 2*x - x^2 + x^3))/x - 
     2*as*CA*Lp*LQ*delta[1 - x] - 4*as*CA*Lp*plusD[0, 1 - x]|>, 
 "TMDPDF" -> <|"qq" -> 2*as*CF*(1 + Lp - x + Lp*x) - 
     (as*CF*(12*Lh*Lp + 6*Lp^2 + Pi^2)*delta[1 - x])/6 - 
     4*as*CF*Lp*plusD[0, 1 - x], 
   "qg" -> -(as*(2*(-1 + x)*x + Lp*(1 - 2*x + 2*x^2))), 
   "gq" -> 2*as*CF*x - (2*as*CF*Lp*(2 - 2*x + x^2))/x, 
   "gg" -> (4*as*CA*Lp*(-1 + 2*x - x^2 + x^3))/x - 
     (as*CA*(12*Lh*Lp + 6*Lp^2 + Pi^2)*delta[1 - x])/6 - 
     4*as*CA*Lp*plusD[0, 1 - x]|>, 
 "Soft" -> <|"q" -> -1/3*(as*CF*(6*Lp^2 + 12*Lp*Lv + Pi^2)), 
   "g" -> -1/3*(as*CA*(6*Lp^2 + 12*Lp*Lv + Pi^2))|>|>
