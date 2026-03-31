function varargout = beta1_vec(hnd_denoise_w1, w_star, nu_W1_vec, nu_W1_mism_vec, nu_vec, iterInf, Np2)

    % Run output denoiser for (vector) inputs "nu_W1_vec"

    Vz2 = var(w_star);
    muZ2 = mean(w_star);

    if isempty(nu_W1_mism_vec)
        nu_W1_mism_vec = nu_W1_vec; 
    end 

    beta1_avg_vec = zeros(1, length(nu_W1_vec));

    for ll = 1:length(nu_W1_vec)

        % ------- Adapted from GAMPMATLAB package -------
        % https://sourceforge.net/projects/gampmatlab/ 
        % Authors: Sundeep Rangan, Philip Schniter and others
        % Generate w_star = p1 + e (e is the error)
        % with E[e] = 0, var(e) = nu_W1_vec, and Cov[p1,e] = 0 (uncorrelated)
        % Construct e = a * w_star + b * Np2 + c
        % and find p1 = w_star - e
        a = min(1, nu_W1_mism_vec(ll) / Vz2);
        b = sqrt(nu_W1_mism_vec(ll) * (1 - a));
        c = -a * muZ2;

        p1 = w_star - (a * w_star + b * Np2 +c);

        [w1hat, beta1] = hnd_denoise_w1(p1.', nu_W1_vec(ll), nu_vec, iterInf);
        beta1_avg_vec(ll) = mean(beta1);
    end

    %% Return different numbr of output depending if the function is used with vector inputs or not
    varargout{1} = beta1_avg_vec; 
    if isscalar(nu_W1_vec)
        varargout{2} = w1hat; 
        varargout{3} = p1; 
    end 

end
