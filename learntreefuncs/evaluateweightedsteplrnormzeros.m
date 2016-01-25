function val = evaluateweightedsteplrnormzeros(alpha, egam, egbv, egfmat, egtmat, degam, degbv, params, egweights)

% this might be confined to ceninds later, don't need this now
egam = egam + alpha*degam;
egam = orthonormalize(egam, params.proj_dimensions, params.reconstruction_size);
egbv = egbv + alpha*degbv;
val = evallrfun(egfmat, egtmat, egam, egbv, params, egweights);
