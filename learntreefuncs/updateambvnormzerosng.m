function [egam, egbv, alpha, gam, gbv] = updateambvnormzerosng(egam, egbv, egfmat, egtmat, weights, params, egweights)

weps = params.weps;
nfeats = params.nfeats;
wal = params.maxgradientstep;
nzeros = params.nzeros;

proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;

tfn = size(egfmat, 1);

% all globals
    %
    % our functional update is A_n+1=(A_n+w D)/norm(A_n+w D)
    % so that A is always normalized; this means we need to be
    % quite careful about gradients.
    %

    % Computes gradient wrt A numerically.
    egam = orthonormalize(egam, proj_dimensions, reconstruction_size);

    % Compute the gradient on a subset of the data (for speed)
    % Chose All
    % subsample = 1:numel(weights);
    % Chose Highest Weight
    %[~,idx] = sort(weights, 'descend');
    %subsample = idx(1:min(numel(idx),params.maxnodegradient));
    % Chose Random Sample
    subsample = randsample(numel(weights), min(numel(weights),params.maxnodegradient));

    [sgv, gstars, egbs, ~] = evallrfun(egfmat(subsample,:), egtmat(subsample,:), egam(subsample,:), egbv(subsample,:), params, egweights(subsample));

    wgbv=-gstars+(1/(1-weps))*(egbs-weps*gstars);
    % We force the last nzeros rows of the gradient to be zero
    wgbv(:, proj_dimensions-nzeros+1:end) = 0;
    wgbv = wgbv .* repmat(egweights(subsample), [1,size(wgbv,2)]);

    wgbm = zeros(numel(subsample), nfeats*proj_dimensions);

    % We force the last nzeros rows of the gradient to be zero
    for i=1:proj_dimensions-nzeros
        wgbm(:, nfeats*(i-1)+1:nfeats*i)=((wgbv(:, i))*ones(1, nfeats)).*(egfmat(subsample,:));
    end
    wgbm = wgbm .* repmat(egweights(subsample), [1,size(wgbm,2)]);

    gam=zeros(1, proj_dimensions*reconstruction_size);
    gradeps=1e-3;

    % We don't need to touch the setrows.
    for lwgi=params.setrows*reconstruction_size+1:proj_dimensions*reconstruction_size
        tempegam = egam(subsample,:);
        tempegam(:, lwgi) = egam(subsample, lwgi) + weights(subsample) * gradeps;
        tempegam = orthonormalize(tempegam, proj_dimensions, reconstruction_size);
        wgv = evallrfunaonly(egtmat(subsample,:), tempegam, egbs, params, egweights(subsample,:));
        gam(lwgi)=(wgv-sgv)/gradeps;
    end

    gbv=sum((weights(subsample)*ones(1, proj_dimensions*nfeats+proj_dimensions)).*[wgbm, wgbv], 1);

    degam = -(weights*ones(1, proj_dimensions*reconstruction_size)).*(ones(tfn, 1)*gam);
    % but we need two properties.  first, degam.am=0 so that the magnitude
    % of am doesn't grow, and second, that the whole gradient is
    degbv = -weights*gbv;
    %
    % now take a step - don't need howfaralpha, psd is covered
    %
    % phase 1: blank search for best value over
    % we look for descent -
    % makesvals % don't need to do this, svals are up to date
    % now an alpha step
    alpha=0;
    cs=evaluateweightedsteplrnormzeros(alpha, egam, egbv, egfmat, egtmat, degam, degbv, params, egweights);

    alpha=wal;
    ce=evaluateweightedsteplrnormzeros(alpha, egam, egbv, egfmat, egtmat, degam, degbv, params, egweights);
    stepflag=0;
    fail=0;
    while (stepflag==0)&&(fail==0)
        alpha=wal/2;
        ci=evaluateweightedsteplrnormzeros(alpha, egam, egbv, egfmat, egtmat, degam, degbv, params, egweights);

        tvw=(3*cs-4*ci+ce)/(4*cs-8*ci+4*ce);
        tv=max(min(tvw, 1), 0);
        if tv==0
            cmin=cs;
        elseif tv==1
            cmin=ce;
            alpha=tv*wal;
        else
            alpha=tv*wal;
            cmin=evaluateweightedsteplrnormzeros(alpha, egam, egbv, egfmat, egtmat, degam, degbv, params, egweights);
        end
        if ~(cmin<cs)
            wal=wal*(1/2);
            ce=ci;
        else
            stepflag=1;
           % fprintf(1, 'step %d\n', nnode.alpha);
        end
        if wal<1e-6
            fail=1;
            disp('* no step *');
            alpha=0;
        end
    end
    if fail==0
        egam = egam + alpha*degam;
        egbv = egbv + alpha*degbv;

        % now this needs a renormalization
        egam = orthonormalize(egam, proj_dimensions, reconstruction_size);
    end
