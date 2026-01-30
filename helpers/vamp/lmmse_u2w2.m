function [u2d_hat, alpha2_hat, w2_hat, beta2_hat, u2d_hat_part] = lmmse_u2w2( ...
        Pcfg, ...
        r2, nu_U2, ...
        p2, nu_W2, ...
        w_intf, n, np, ...
        idx_data, ...
        Hcfft, ... %Eigenvalues of circulant channel A
        Hpfft, ... %Eigenvalues of circulant precoder P
        n_it_thresh, ... %Number of iterations until annealing is turned off
        iter, ...
        nu_N1sim, ... %Pre-intensity noise on whole simulation bandwidth [-B,+B]
        Lambda_D_col, ...
        gam_Bd_over_Bsam, ...
        N_os, ...
        ind_sic, ...
        PtxdB, ...
        annealcfg, ...
        Adata, ...
        rGVAMP, ...
        flag_trig, ...
        u2d_hat_part ... %May be supplied if only the messages (w2hat, beta2) need to be calculated
    )

    %% Perform LMMSE denoising for model W = Adata*U + N1 + s, where N1 is colored CSCG and s is known interference.
    % Note: nu_N1sim := B' * nu1, i.e., the optical noise power before filtering, i.e., over the whole "simulation bandwidth [-B, +B]"
    % Thus, nuN1/gamma = nu_N1sim

    u2d_hat = []; %init
    w2_hat = []; %init

    %% Configuration for variance annealing
    if iter > 1E15 %Annealing OFF; enter here for state evolution / EXIT charts
        L_D = Lambda_D_col;
        gamma = gam_Bd_over_Bsam; %Colored noise, relative bandwidth: Bd/Bsam = Bd/(Nos*B)

    elseif iter > n_it_thresh || flag_trig %Annealing OFF, if more than n_it_thresh iterations have passed, or trigger flag is active
        nu_N1sim = nu_N1sim + (ind_sic <= 1) * annealcfg{3}(1); %Allow for small additional (white) annealing
        L_D = Lambda_D_col + (ind_sic <= 1) * annealcfg{3}(1); %Allow for small additional (white) annealing
        gamma = gam_Bd_over_Bsam; %Colored noise, relative bandwidth: Bd/Bsam = Bd/(Nos*B)

    else %Annealing ON
        minAnneal = annealcfg{2}(3);
        nu_N1sim = max(minAnneal, annealcfg{2}(1) * exp(annealcfg{2}(2) * PtxdB)); %Optical noise annealing
        gamma = 1; %White noise
    end

    if any(contains(["cDC", "cDR", "cU", "cO", "_"], Pcfg)) %Supports only FFT-based precoder

        %% Annealing
        if iter > n_it_thresh || flag_trig %Use actual LMMSE denoiser
            nu_N1simU = nu_N1sim;
            LD_U = L_D;
            nu_N1simW = nu_N1sim;
            LD_W = L_D;
        else %Mismatch LMMSE denoiser, rely more on denoised u2hat
            nu_N1simU = nu_N1sim;
            LD_U = nu_N1sim;
            nu_N1simW = 1/4 * nu_N1sim; %1/4 of noise of U
            LD_W = 1/4 * nu_N1sim; %1/4 of noise of U
        end

        %% Conditional variances
        alpha2_hat = nu_U2 * (nu_N1simU + nu_W2) / (2^(rGVAMP) * N_os * nu_U2 + nu_N1simU + nu_W2); %rGVAMP = 1 for real-valued prior

        nu_N1_W = (gamma * nu_N1simW); %Variance of filtered noise (effective variance)
        beta2_hat = nu_W2 / (nu_W2 + nu_N1simW) * (nu_N1_W + alpha2_hat * np / n * nu_W2 / (nu_W2 + nu_N1simW));

        if ~isempty(r2) && ~isempty(p2) % Compute only of r2 and p2 are supplied
            if isempty(u2d_hat_part) % Compute u2d_hat_part if it is not supplied as a function argument
                L_N1tilde_inv_U = 1 ./ (nu_W2 + LD_U);
                tmp1 = conj(Hcfft) .* L_N1tilde_inv_U .* ((1 / sqrt(2 * n)) * fft((p2 - w_intf), 2 * n)); %Note: unitary 2n-FFT
                tmp2 = sum([tmp1(1:n), tmp1(n + 1:end)], 2);

                interm_hat = 1 / sqrt(N_os) * (sqrt(n) * ifft(conj(Hpfft) .* tmp2)); %Note: unitary n-inverse-FFT
                u2d_hat_part = interm_hat(idx_data); %Remove know, i.e., apply T matrix
            end

            if rGVAMP %Real-valued prior
                u2d_hat = alpha2_hat * (2^(rGVAMP) * real(u2d_hat_part) + 1 / nu_U2 * r2);
            else
                u2d_hat = alpha2_hat * (u2d_hat_part + 1 / nu_U2 * r2);
            end

            % -------- For debugging --------------
            % Fm = 1/sqrt(2* n )*dftmtx(2*n);
            % if isscalar(LnU)
            %     L_N1til_U_D = diag((nu_W2 + LnU)*ones(2*n,1));
            % else
            %     L_N1til_U_D = diag((nu_W2 + LnU));
            % end
            % u2d_hat_L = alpha2_hat * (real(Adata'*(((Fm'*L_N1til_U_D*Fm)/2^rGVAMP)\(p2 - w_intf))) + 1/nu_U2 * r2);
            % if norm(u2d_hat_L - u2d_hat) > 1E-9; error('calculation wrong'); end
            % ------------------------------------

            if nargout == 4 %If 4 output arguments are required, calculate also messages (w2_hat, beta2_hat)
                u2d_hat_padded = zeros(n, 1);
                u2d_hat_padded(idx_data) = u2d_hat; %T*x
                L_N1tilde_inv_W = 1 ./ (nu_W2 + LD_W);

                x2d_t1 = Hpfft .* (1 / sqrt(n) * fft(u2d_hat_padded, n)); %Note: unitary n-FFT
                u2d_t3tmp = sqrt(2 * n) / sqrt(N_os) * ifft(L_N1tilde_inv_W .* Hcfft .* [x2d_t1; x2d_t1]); %Unitary 2n-inverse-FFT

                w2_hat = nu_W2 * (u2d_t3tmp + ifft(L_N1tilde_inv_W .* fft(w_intf - p2, 2 * n), 2 * n)) + p2;

                % -------- For debugging --------------
                % Fm = 1/sqrt(2* n )*dftmtx(2*n);
                % if isscalar(LnW)
                %     L_N1til_W_D = diag((nu_W2 + LnW)*ones(2*n,1));
                % else
                %     L_N1til_W_D = diag((nu_W2 + LnW));
                % end
                % w2_hat_L = nu_W2 * ( ((Fm'*L_N1til_W_D*Fm)) \ (Adata * u2d_hat + w_intf - p2) ) + p2;
                % if norm(w2_hat_L - w2_hat) > 1E-9;  error('calculation wrong'); end
                %  ------------------------------------

            end
        end

    else
        error('Precoder not supported.');
    end
end
