function [Ef, Tf] = fvmSfDecomposition(Sf, nf, impK, eCN)
%fvmSfDecomposition Over-relaxed approach [p. 244, Moukalled - 2015]

    sf = impK.' * Sf;
    cosTheta = eCN ^ nf; % Inner product
    Ef = (sf.norm()/cosTheta) * eCN;
    Tf = sf - Ef;
end
