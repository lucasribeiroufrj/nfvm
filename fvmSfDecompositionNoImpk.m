function [Ef, Tf] = fvmSfDecompositionNoImpk(Sf, Sfnorm, nf, e_CN)
%fvmSfDecompositionForGhost Over-relaxed approach.
%   See [p. 244, Moukalled - 2015].

cosTheta = e_CN ^ nf;
Ef = (Sfnorm/cosTheta) * e_CN;
Tf = Sf - Ef;

end

