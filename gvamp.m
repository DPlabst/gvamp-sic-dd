% <header>
% /*********************************************************************************************************
%  * Information Rates of Approximate Message Passing for Bandlimited Direct-Detection Channels
%  * Daniel Plabst; Institute for Communications Engineering (ICE), Technical University of Munich, Germany
%  * Mohamed Akrout; University of Manitoba, Winnipeg, Canada
%  * Version r2: Date: 2026-01-28
%  *********************************************************************************************************/
% </header>

function [filename, vver] = gvamp(nvArgs)
    arguments

        nvArgs.m string % Real modulation: 'M-ASK-o' or 'M-MB-o-nu' (order M, offset o \ in [0, 1] and Maxwell - Boltzmann shaping parameter nu)
        nvArgs.ch double = 0 % 0=Optical Fiber
        nvArgs.S double % Total number of SIC stages S >=1
        nvArgs.st (1, :) double = [] % Simulate subset of SIC stages, e.g., S = 4 and st = [2,3,4]

        nvArgs.P string % Precoder: {cO, cU, cDC, cDR, _}; cO=circul. orthogonal (real), cU=circul. unitary (complex), cDC=chromatic dispersion (complex), cDR = chromatic dispersion with hermitian symmetry (real),"_" = none

        nvArgs.runGVAMP logical = 1 % Run GVAMP
        nvArgs.runEXIT logical = 1 % Run EXIT

        nvArgs.rGVAMP logical = 1 % If real-valued prior, the LMMSE and input denoiser uses a real-value constraint

        nvArgs.ps string % Pulseshaping filter {RC, RRC}
        nvArgs.alp double % Roll-off factor \in [0.01, 1]
        nvArgs.Nsp double % Pulseshaping filter span (choose large enough of roll-off factor is small)

        nvArgs.n double % Block length
        nvArgs.L double % Fiber length in [m], for L > 0, chromatic dispersion will occur
        nvArgs.Rs double % Symbol rate in [Bd] or [Sym/s]

        nvArgs.Nr double % Number of GVAMP recovery tries per block with different initialization
        nvArgs.Nb double % Number of blocks to be transmitted
        nvArgs.Npi double = 0 % Optional Npi x N_span symbols within a block are pilots

        nvArgs.Ptx string % Transmit power vector in dB, e.g.,'0:1:15'
        nvArgs.nu (1, :) double = [1, 0] %Noise variances of [optical white CSCG, electrical real WGN] over bandwidth [-B,B]

        nvArgs.d1 double % Initial damping of mean
        nvArgs.d2 double % Initial damping of variance
        nvArgs.aD double = 1 % [0/1] automatic damping
        nvArgs.it double % Maximum number of GVAMP iterations
        nvArgs.itEXIT double % Maximum number of EXIT iterations

        nvArgs.RxBW logical = false % [0/1]: Apply PD filter to filter out-of-band noise (Filter bandwidth is adjusted to signal bandwidth) 
        nvArgs.ADCfilt double = 0 % Apply ADC filter with bandwidth [-B, +B]

        nvArgs.outdencfg (2, 1) cell; %Cell-array: first cell use (mismatched) [nu_1, nu_1] as variance in output denoiser, second cell: generalized Xi2 integration method {'SPA', 'Davies'} when electrical noise > 0
        nvArgs.dampcfgStepW (1, 3) double; %[minWindow, maxWindow, [0/1]:reduce step window with increasing SIC stage]
        nvArgs.annealcfg (3, 1) cell; %Cell-array {[n_it_thresh, reduce threshold vs. SIC stage], [alpha, beta, minAnneal], [anneal if iter > switchIterThresh]}

        nvArgs.VAMP_dbg_mode double = 0 %[0,1,2], 0=none, 1=plot rate trace, 2=plot internal variance vs. actual variances
        nvArgs.EXIT_dbg_mode double = 0 %[0,1], 0=none, 1=plot exit charts
        nvArgs.run_simp double = 0 % [0,1], 0=off (default), 1=on (accelerates GVAMP iterations by switch debiasing and s-gmi optimization off)

        nvArgs.cons double = 0 % [0/1] Console output (iterations, rate, error, ...)

        nvArgs.save_file double = 0 % Save results to file
        nvArgs.save_trace double = 0 % Save iteration traces also in MAT files (increases filesize)
        nvArgs.save_logAPP double = 0 % Save logAPPs (increases filesize)

    end

    addpath('helpers/vamp/'); %Packages
    addpath('helpers/comm'); %Helper files from SP rates
    addpath('helpers/other'); %Mainly for printing and organization
    addpath('helpers/other/gampmatlab'); %For GAMPMATLAB project
    addpath('helpers/other/gx2'); %For generalized Xi2 PDFs

    %% -------------------------------- Settings  --------------------------------
    vver = "r2"; %VAMP version

    % --- Settings ------
    overwrite_flag = 1; % Overwrite saved simulation files
    VAMP_dbg_mode = nvArgs.VAMP_dbg_mode;
    EXIT_dbg_mode = nvArgs.EXIT_dbg_mode;

    % Setups that can accerlerate GVAMP
    if nvArgs.run_simp
        debias_nmse = 0;
        gvamp_sgmi = 0;
    else
        debias_nmse = 1;
        gvamp_sgmi = 1;
    end

    % - Choices: {'RC','RRC'}
    ps_filter = struct('name', char(nvArgs.ps), 'alpha', nvArgs.alp);

    filename = strcat( ...
        vver, ',', ... %Version
        nvArgs.m, ... %Modulation format
        ',S=', num2str(nvArgs.S), ... %Number of SIC stages
        ',C=', num2str(nvArgs.ch), ... %Used channel matrix A
        ',Rb=', num2str(nvArgs.RxBW), ... %Bool of filter before PD is active
        ',n1=', num2str(nvArgs.nu(1)), ... %Optical noise on whole bandwidth [-B,B]
        ',n2=', num2str(nvArgs.nu(2)), ... %Electrical noise on whole bandwidth [-B,B]
        ',n1m=', num2str(nvArgs.outdencfg{1}(1)), ... %Denoiser config for mismatch
        ',n2m=', num2str(nvArgs.outdencfg{1}(2)), ... %Denoiser config for mismatch
        ',P=', nvArgs.P, ... %Precoder
        ',ps=', ps_filter.name, ... %Pulseshaping filter
        ',a=', num2str(round(ps_filter.alpha, 2)), ... %Roll-off factor
        ',Nsp=', num2str(nvArgs.Nsp), ... %Filter span
        ',Npi=', num2str(nvArgs.Npi), ... %Append Npi x Nsp pilots at block end
        ',n=', num2str(nvArgs.n), ... %Block length
        ',L=', num2str(round(nvArgs.L / 1E3, 2)), ... %Fiber length [km]
        ',Rs=', num2str(round(nvArgs.Rs / 1E9, 2)), ... %Symbol rate [GBd]
        ',Nr=', num2str(nvArgs.Nr), ... %Number of independent GVAMP tries per block
        ',Nb=', num2str(nvArgs.Nb), ... %Number of transmitted blocks
        ',I=', num2str(nvArgs.it) ... %Number of maximum iterations
    );

    % ------------ Check if file exists already ----------
    if overwrite_flag == 0
        folderPath = strcat('results/', vver, '/');
        fullFilePath = fullfile(folderPath, strcat(filename, '.txt'));

        if exist(fullFilePath, 'file') == 2 % Check if the file exists
            fprintf('Simulation file already exists... skipping\n');
            return;
        end
    end

    % ---------- System Parameters  ---------
    N_os = 2; %Oversampling factor
    N = nvArgs.n;
    M = N_os * N;
    L_fiber = nvArgs.L;
    R_sym = nvArgs.Rs;
    N_span = nvArgs.Nsp;
    N_block = nvArgs.Nb; %Number of TX blocks
    N_recov = nvArgs.Nr; %Number of tries per block
    Pcfg = nvArgs.P;
    N_pil = round(nvArgs.Npi * N_span); %Append Npi * N_span pilots and block end
    S_SIC = nvArgs.S; % Number of SIC stages
    maxIt = nvArgs.it; %Maximum amount of iterations
    damp = nvArgs.d1; %Damping for mean
    dampGam = nvArgs.d2; %Damping for variance
    vamptol = 1E-6; %Stopping tolerance for GVAMP

    if isempty(nvArgs.st)
        s_vec = 1:S_SIC; % Simulate all S stages or select a subset of [1, ... , S_SIC]
    else
        if any(nvArgs.st > S_SIC); error('Wrong SIC levels supplied.'); end
        s_vec = nvArgs.st;
    end

    % ------------- Noise ----------------
    TXnoise = struct('set', 'on', 'var_n_pre', nvArgs.nu(1)); %Pre-intensity noise (Optical noise)
    RXnoise = struct('set', 'on', 'var_n_post', nvArgs.nu(2)); %Post-intensity noise (Electrical noise)
    P_tx_vec_dB = str2num(nvArgs.Ptx); % TX power range
    P_tx_vec = 10.^(P_tx_vec_dB / 10);
    L_SNR = length(P_tx_vec);

    if TXnoise.var_n_pre ~= 0
        x_dB_TX_PW = 10 * log10(P_tx_vec / TXnoise.var_n_pre); %Define SNR by pre-intensity noise
    elseif TXnoise.var_n_pre == 0 && RXnoise.var_n_post > 0
        x_dB_TX_PW = 10 * log10(P_tx_vec / RXnoise.var_n_post); %Define SNR by post-intensity noise
    end

    % ---------------- Pre-allocate tensors  ----------------
    IEXIT_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);

    IqYX_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);
    REQITER_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);
    APP_cost_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);
    SER_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);
    HDMI_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);
    MMSESE_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);
    MMSEVAMP_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov);

    P_tx_wave_measured_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov); %TX
    P_tx_sym_measured_ten = zeros(S_SIC, length(P_tx_vec), N_block, N_recov); %TX

    % ---- Pre-allocate long vectors for parfor per SIC stage -----
    IqYX_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);
    APP_cost_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);
    SER_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);
    HDMI_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);
    rateEXIT_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);

    mmseEXIT_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);
    mmse_VAMP_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);
    reqiter_VAMP_lvec = zeros(length(P_tx_vec) * N_block * N_recov, 1);

    rate_VAMP_iter_lvec = zeros(length(P_tx_vec) * N_block * N_recov, nvArgs.it);
    rate_EXIT_iter_lvec = zeros(length(P_tx_vec) * N_block * N_recov, nvArgs.itEXIT);

    nmse_VAMP_iter_lvec = zeros(length(P_tx_vec) * N_block * N_recov, nvArgs.it);
    nmse_EXIT_iter_lvec = zeros(length(P_tx_vec) * N_block * N_recov, nvArgs.itEXIT);

    P_tx_wave_measured = zeros(length(P_tx_vec) * N_block * N_recov, 1); %TX
    P_tx_sym_measured = zeros(length(P_tx_vec) * N_block * N_recov, 1); %TX

    % ---------- Generate filters and matrices -----------
    h_tx = gen_tx_filter(ps_filter.name, ps_filter.alpha, N_os, N_span); %TX filter
    [a_comb] = gen_fiber_response(0, N_os, h_tx, L_fiber, R_sym); %Combined filter of TX filter and CD channel
    L_A = fft(a_comb, M).'; %Eigenvalues of channel

    % Optical PD input filter
    if nvArgs.RxBW
        gD = fir1(100, (1 + ps_filter.alpha) / N_os);
    else
        gD = 1;
    end

    [S, pU, mod_alph] = gen_mod_alphabet(nvArgs.m); %Generate modulation alphabet and PMF

    % Tensors for log APPs
    log_APP_lvec = zeros(length(P_tx_vec) * N_block * N_recov, mod_alph.Q, N / S_SIC); %tensor for logAPPs per SIC stage

    %A = gen_chan(nvArgs.ch, a_comb, N_span, N_os, N, M); % Generate GVAMP channel

    rng(1, "twister"); % set rng
    [L_P] = gen_prec(nvArgs.ch, Pcfg, N); % Generate precoder (can also output a P matrix)
    %A = A * P; % Apply precoder to channel

    for ind_sic = s_vec %Loop over SIC stages

        % Variance annealing
        n_it_thresh = round(nvArgs.annealcfg{1}(1)) * exp(- nvArgs.annealcfg{1}(2) * (ind_sic - 1)); %Switch OFF annealing after n_it_thresh iterations
        n_it_min = 0; %Run GVAMP for at least minIt iterations
        autoDampCFG = configure_autodamp(nvArgs, ind_sic); %Configure automatic damping

        % ---------- Compute Adata ------------
        [idx_data, ~, ~, ~, Nprime] = gen_sic_mat(N_pil, N_span, N, [], ind_sic, S_SIC); 
        %[idx_data, idx_pil, Adata, Apil, Nprime] = gen_sic_mat(N_pil, N_span, N, A, ind_sic, S_SIC); %Legacy

        trAAh = 1/2 * sum(sum(reshape(abs(L_A).^2, N, 2), 2) .* abs(L_P).^2); 
        trAdAdh = Nprime / N * trAAh; 
        %%rAAh_leg = sum(conj(A) .* A, 'all'); %Legacy: Energy of (full) A
        %%trAdAdh_leg = sum(conj(Adata) .* Adata, 'all'); %Legacy: Energy of Adata

        print_sys_params(mod_alph, ps_filter, TXnoise, RXnoise, L_fiber, R_sym, nvArgs); % Print settings
        fprintf('Measured (normalized) power per symbol: %.2f\n', 1 / N * trAAh / N_os); %Print power
        if abs(10 * log10(N_os) - 10 * log10(1 / N * trAAh)) > 0.2; error('Power constraint violation'); end %Sanity check

        fprintf('\n')
        fprintf('Starting iterations ...\n')

        for ind_n_x_rb = 1:(L_SNR * N_block * N_recov) %can also do parfor over all SNRs x number of blocks x number of trials per block

            w_star = [];
            [ind_n, rb, rr] = ind2sub([L_SNR, N_block, N_recov], ind_n_x_rb);
            rng(rb + 123456, "twister"); % set seed

            % Print progress
            fprintf('SIC %d/%d | SNR %d/%d | Block %d/%d | Try %d/%d\n', ...
                ind_sic, S_SIC, ...
                ind_n, length(P_tx_vec), ...
                rb, N_block, ...
                rr, N_recov);

            %% Run transmitter
            [Ualph_sc, u0, n_tx, u_ind, P_tx_sym, mean_U_prior, var_U_prior] = comm_transmitter(P_tx_vec(ind_n), TXnoise.var_n_pre, N_os, S, mod_alph.name, N, mod_alph, pU);

            % ------- FFT: Apply channel to calculate known interference --------
            u0_pil_padded = u0.';
            u0_pil_padded(idx_data) = 0; %zero unknown indices
            x2d_t1_pil = L_P .* (1 / sqrt(N) * fft(u0_pil_padded, N)); %Note: unitary n-FFT
            w_intf = sqrt(2 * N) / sqrt(N_os) * ifft(L_A .* [x2d_t1_pil; x2d_t1_pil]); %Unitary 2n-inverse-FFT
            %%s_intf = u0(idx_pil).'; %Legacy: Known interference or pilots
            %%w_intf = Apil * s_intf; %Legacy: Known interference via pilots, SIC stages

            P_tx_sym_measured(ind_n_x_rb) = P_tx_sym; %Measure power

            %% -------- VAMP ---------

            % ------- FFT: Apply channel: W_nonoise = A*U --------
            %u0_padded = zeros(N, 1);
            %u0_padded(idx_data) = u0; %T*x
            x2d_t1 = L_P .* (1 / sqrt(N) * fft(u0.', N)); %Note: unitary n-FFT
            w_nonoise = sqrt(2 * N) / sqrt(N_os) * ifft(L_A .* [x2d_t1; x2d_t1]); %Unitary 2n-inverse-FFT
            %w_nonoise_leg = A * (u0).'; %Legacy

            P_tx_wave_measured(ind_n_x_rb) = mean(abs(w_nonoise).^2); %Measure power

            n_tx = conv(n_tx, gD, 'same'); %Filter optical noise to bandwidth [-(1+alpha)*B/2, +(1+alpha)*B/2]
            w_noisy = w_nonoise + n_tx.'; %Add filtered noise

            SNRmeasdB(ind_n_x_rb) = 10 * log10(mean(abs(w_nonoise).^2) / mean(abs(n_tx).^2)); %Measure (some) SNR in dB

            w_noisy_up = interp(w_noisy, 2); %Interpolate by x2 before squaring
            y_inter = abs(w_noisy_up).^2; %Squaring causes bandwidth doubling

            %% --------- DAC Filtering --------
            if nvArgs.ADCfilt == 1

                gRX = fir1(1000, 0.5); %DAC filter to bandwidth [-B, +B]

                % Extend (pad) to avoid transients at border edges
                w_noisy_ext = [w_noisy(end - N_os * N_span + 1:end); w_noisy; w_noisy(1:N_os * N_span)];
                w_noisy_up_ext = interp(w_noisy_ext, 2); %Interpolate
                y_inter_ext = abs(w_noisy_up_ext).^2; %Apply square-law detector: squaring causes bandwidth doubling

                y_interfilt = conv(y_inter_ext, gRX, 'same');
                y_interfilt_ds = downsample(y_interfilt, 2);

                n_rx = sqrt(RXnoise.var_n_post) * randn(N_os * N, 1); %Generate real white Gaussian noise
                y_noisy = y_interfilt_ds(N_os * (N_span) + 1:N_os * N_span + M) + n_rx; %Select middle part

                %% Take absolute value of measurements
                % When nu1 > 0 and nu2 = 0, y_inter_ext \in R_0^+.
                % When filtering intensities with the DAC filter, y can have negative values at the block edges 
                % The surrogate q(y|x) assigns zero likelihood to negative y, which causes problems in damping. 
                % For simplificty, we take the absolute value, another
                % approach models small, but nonzero, electrical noise in
                % q(y|x).
                if TXnoise.var_n_pre > 0 && RXnoise.var_n_post == 0
                    y_noisy = abs(y_noisy); 
                end 

            elseif nvArgs.ADCfilt == 0
                y_inter_ds = downsample(y_inter, 2);
                n_rx = sqrt(RXnoise.var_n_post) * randn(N_os * N, 1); %Generate real white Gaussian noise
                y_noisy = y_inter_ds + n_rx; %Don't filter; just sample to Nos = 2 (may result in aliasing)
            end

            % -------------- GVAMP Initial messages ---------
            rng(rr + 2 * 123456, "twister"); %Different initializations with offset

            pvar = 1e1 * P_tx_vec(ind_n) * 1 / M * trAdAdh; %for p1
            rvar = 2 * P_tx_vec(ind_n); %for r1

            phat = sqrt(pvar / 2) .* (randn(M, 1) + 1j * randn(M, 1)); %initialize mean for p1 randomly; always complex

            if nvArgs.rGVAMP %If prior is real-valued
                rhat = sqrt((rvar)) .* (randn(Nprime, 1)); %real-valued init
            else
                rhat = sqrt((rvar) / 2) .* (randn(Nprime, 1) + 1j * randn(Nprime, 1)); %complex-valued init
            end

            % Scale noise variances, e.g., for mismatched denoising
            nu_p_vec = [nvArgs.outdencfg{1}(1) * TXnoise.var_n_pre; nvArgs.outdencfg{1}(2) * RXnoise.var_n_post];

            %% ------------ Handles ---------
            hnd_denoise_w1 = @(p1, nu_Z1, nu_p_vec, niter) denoise_w1(y_noisy.', p1, nu_Z1, nu_p_vec, nvArgs.outdencfg{2});
            hnd_denoise_u1 = @(r1, nu_U1, rGVAMP) denoise_u1(pU, Ualph_sc, r1, nu_U1, rGVAMP);

            hnd_s = @(x) x(ind_sic:S_SIC:end); %Returns symbols of current stage among all symbols
            hnd_s_pil = @(x0) x0(idx_data(ind_sic:S_SIC:end)); %Returns symbols of current stage among all data symbols

            if debias_nmse
            hnd_debias = @(in) disambig1Dfft(in, hnd_s_pil(u0)); %Remove phase ambiguities; From GAMPMATLAB package
            else
                hnd_debias = @(in) in.'; %note: transpose
            end
            hnd_nmse = @(x1) mean(abs(hnd_debias(hnd_s(x1)).' - hnd_s_pil(u0)).^2) / P_tx_vec(ind_n); %Normalized MSE

            hnd_h2 = @(pe) - ((1 - (pe + eps)) .* log2(1 - (pe + eps)) + (pe + eps) .* log2((pe + eps))); %Binary entropy function
            hnd_hdmi = @(pe) max(zeros(length(pe), 1), log2(length(Ualph_sc)) - (hnd_h2(pe) + pe .* (log2(length(Ualph_sc) - 1)))); %Hard decision mismatched rate (cf. Fano's inequality)
            hnd_sgmi = @(r1, nu_U1, rGVAMP) comm_gmi(r1, nu_U1, u_ind, 1, pU, idx_data, ind_sic, S_SIC, hnd_denoise_u1, rGVAMP); %Optimize s-parameter
            hnd_gmi = @(r1, nu_U1, rGVAMP) comm_gmi(r1, nu_U1, u_ind, 0, pU, idx_data, ind_sic, S_SIC, hnd_denoise_u1, rGVAMP); %Do not optimize s-parameter
            hnd_ser = @(u1) comm_calc_ser(u1, Ualph_sc, u_ind, idx_data, ind_sic, S_SIC);

            if gvamp_sgmi
                hnd_gvamp_rate = hnd_sgmi; %with optimization over s
            else
                hnd_gvamp_rate = hnd_gmi; %without optimization over s
            end

            %% Damping configuration
            DAMPcfg = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % Damp only the extrinsic coming from output denoiser (p2, nuW2)

            %% ----------- LMMSE denoiser -------------------
            if nvArgs.RxBW %Check if bandwidth-limited noise
                id = 0:1 / M:(1 - 1 / M);
                ps_filt_marg = 0.01; %Slightly model the PD filter larger (avoids divisions by zero in eigenvalues)
                gam_BdBsam = (1 + (ps_filter.alpha + ps_filt_marg)) / 2;
                idx2 = (id > (gam_BdBsam * N / M) & id < (1 - gam_BdBsam * N / M)); %Out of band
                Lambda_D = (TXnoise.var_n_pre) * ones(1, M).'; %var_n_pre is the variance of the noise on the whole simulation bandwidth [-B, +B]
                Lambda_D(idx2) = 0;
            else %Otherwise flat
                gam_BdBsam = 1; %flat
                Lambda_D = (TXnoise.var_n_pre) * ones(1, M).'; %white
            end

            hnd_lmmse_u2w2 = @(r2, nu_U2, p2, nu_W2, iter, flag_trig, rGVAMP, u2d_hat_part) ...
                lmmse_u2w2(Pcfg, r2, nu_U2, p2, nu_W2, w_intf, N, Nprime, idx_data, L_A, L_P, n_it_thresh, ...
                iter, TXnoise.var_n_pre, Lambda_D, gam_BdBsam, N_os, ind_sic, P_tx_vec_dB(ind_n), nvArgs.annealcfg, [], rGVAMP, flag_trig, u2d_hat_part);
            % ----------------

            %Handle for damp cost function (Note: Must take Adata as an argumenet, if channel is not FFT-based)
            hnd_damp_cost = @(r1, nu_X1, x1, alpha1) calc_damp_cost(y_noisy.', N, r1, nu_X1, x1, alpha1, TXnoise.var_n_pre, gam_BdBsam, nu_p_vec, pU, [], w_intf, trAdAdh, idx_data, L_A, L_P, Pcfg, N_os, nvArgs.rGVAMP, hnd_denoise_u1);

            % ------------ Define ground truth signals -----------
            x_star = u0.';
            w_star = w_noisy; %Optical noise included in w_star

            %% Run state evolution / EXIT
            if nvArgs.runEXIT
                [rate_EXIT, mmse_EXIT, rate_EXIT_iter, mmse_EXIT_iter] = exit.run_exit( ...
                    nu_p_vec, ...
                    w_star, ...
                    x_star, ...
                    hnd_denoise_w1, ...
                    hnd_denoise_u1, ...
                    hnd_lmmse_u2w2, ...
                    hnd_gmi, ...
                    idx_data, ...
                    nvArgs.itEXIT, ...
                    P_tx_vec(ind_n), ...
                    nvArgs.rGVAMP, ...
                    ind_n_x_rb, ... %used as seed
                    EXIT_dbg_mode, ...
                    Nprime, ...
                    M ...
                );

                mmseEXIT_lvec(ind_n_x_rb) = mmse_EXIT;
                rateEXIT_lvec(ind_n_x_rb) = rate_EXIT;

                %% Iteratation vectors
                rate_EXIT_iter_lvec(ind_n_x_rb, :) = rate_EXIT_iter;
                nmse_EXIT_iter_lvec(ind_n_x_rb, :) = mmse_EXIT_iter;

            end

            %% Run GVAMP
            if nvArgs.runGVAMP

                [IqYX, SER, NMSE, IqYX_vamp_iter, NMSE_vamp_iter, APP_cost, reqIter, u1_best, log_APP] = run_gvamp_asym( ...
                    nu_p_vec, ...
                    phat.', ...
                    pvar.', ...
                    rhat.', ...
                    rvar.', ...
                    DAMPcfg, ...
                    hnd_denoise_w1, ...
                    hnd_denoise_u1, ...
                    hnd_lmmse_u2w2, ...
                    hnd_gvamp_rate, ...
                    hnd_ser, ...
                    hnd_nmse, ...
                    hnd_damp_cost, ...
                    damp, ...
                    dampGam, ...
                    autoDampCFG, ...
                    n_it_min, ...
                    maxIt, ...
                    vamptol, ...
                    nvArgs.cons, ...
                    x_star, ...
                    w_star, ...
                    idx_data, ...
                    nvArgs.rGVAMP, ...
                    VAMP_dbg_mode, ...
                    EXIT_dbg_mode, ...
                    ind_sic ...
                );

                %% Scalar results
                IqYX_lvec(ind_n_x_rb) = IqYX;
                SER_lvec(ind_n_x_rb) = SER;
                HDMI_lvec(ind_n_x_rb) = hnd_hdmi(SER);
                APP_cost_lvec(ind_n_x_rb) = APP_cost;
                mmse_VAMP_lvec(ind_n_x_rb) = NMSE;
                reqiter_VAMP_lvec(ind_n_x_rb) = reqIter;

                rate_VAMP_iter_lvec(ind_n_x_rb, :) = IqYX_vamp_iter;
                nmse_VAMP_iter_lvec(ind_n_x_rb, :) = NMSE_vamp_iter;
                log_APP_lvec(ind_n_x_rb, :, :) = log_APP;

            end

        end

        % Reshape to tensors:
        IqYX_ten(ind_sic, :, :, :) = reshape(IqYX_lvec, L_SNR, N_block, N_recov);
        SER_ten(ind_sic, :, :, :) = reshape(SER_lvec, L_SNR, N_block, N_recov);
        HDMI_ten(ind_sic, :, :, :) = reshape(HDMI_lvec, L_SNR, N_block, N_recov);
        APP_cost_ten(ind_sic, :, :, :) = reshape(APP_cost_lvec, L_SNR, N_block, N_recov);
        REQITER_ten(ind_sic, :, :, :) = reshape(reqiter_VAMP_lvec, L_SNR, N_block, N_recov);
        IEXIT_ten(ind_sic, :, :, :) = reshape(rateEXIT_lvec, L_SNR, N_block, N_recov);

        MMSESE_ten(ind_sic, :, :, :) = reshape(mmseEXIT_lvec, L_SNR, N_block, N_recov);
        MMSEVAMP_ten(ind_sic, :, :, :) = reshape(mmse_VAMP_lvec, L_SNR, N_block, N_recov);
        P_tx_wave_measured_ten(ind_sic, :, :, :) = reshape(P_tx_wave_measured, L_SNR, N_block, N_recov);
        P_tx_sym_measured_ten(ind_sic, :, :, :) = reshape(P_tx_sym_measured, L_SNR, N_block, N_recov);

        %% Then tensors save performance over iterations
        RATEVAMP_iter_ten(ind_sic, :, :, :, :) = reshape(rate_VAMP_iter_lvec, L_SNR, N_block, N_recov, nvArgs.it);
        RATESE_iter_ten(ind_sic, :, :, :, :) = reshape(rate_EXIT_iter_lvec, L_SNR, N_block, N_recov, nvArgs.itEXIT);

        MMSEVAMP_iter_ten(ind_sic, :, :, :, :) = reshape(nmse_VAMP_iter_lvec, L_SNR, N_block, N_recov, nvArgs.it);
        MMSESE_iter_ten(ind_sic, :, :, :, :) = reshape(nmse_EXIT_iter_lvec, L_SNR, N_block, N_recov, nvArgs.itEXIT);

        %% For APPs:
        log_APP_ten(ind_sic, :, :, :, :, :) = reshape(log_APP_lvec, L_SNR, N_block, N_recov, mod_alph.Q, N / S_SIC);

    end

    % ----------------------------------------------------
    % --------------------- SAVE ------------------------
    % ----------------------------------------------------

    %% Process data - Choose rate with maximum LR
    [~, minRR] = min(APP_cost_ten, [], 4, 'linear'); %Block with *minimum* cost among different trials

    IqsYX_vec = mean(IqYX_ten(minRR), 3); %Apply best true and average over blocks
    IqYX = mean(IqsYX_vec, 1); % Average rate across all SIC stage: 1 x SNR matrix;

    SERs_vec = mean(SER_ten(minRR), 3); %Apply best true and average over blocks
    SER = mean(SERs_vec, 1); % Average SER across all SIC stage: 1 x SNR matrix;

    HDMIs_vec = mean(HDMI_ten(minRR), 3); %Apply best true and average over blocks
    HDMI = mean(HDMIs_vec, 1); % Average HDMI rate across all SIC stage: 1 x SNR matrix;

    if nvArgs.save_logAPP
    APP_cost_vec = mean(APP_cost_ten(minRR), 3); %Apply best true and average over blocks
    AvgAPP_cost = mean(APP_cost_vec, 1); % Average log likelihood across all SIC stage: 1 x SNR matrix;
    end

    REQITER_vec = mean(REQITER_ten(minRR), 3); %Apply best true and average over blocks
    REQITER = mean(REQITER_vec, 1); % Average log likelihood across all SIC stage: 1 x SNR matrix;

    MMSEVAMP_vec = mean(MMSEVAMP_ten(minRR), 3); %Apply best true and average over blocks
    MMSEVAMP = mean(MMSEVAMP_vec, 1); %Average over SIC stages:  1 x SNR matrix;

    %% Just average
    MMSE_SE_vec = mean(MMSESE_ten, [3, 4]); % Average SER across all SIC stage: 1 x SNR matrix;
    MMSE_SE = mean(MMSE_SE_vec, 1); %Average over SIC stages:  1 x SNR matrix;

    IEXIT_vec = mean(IEXIT_ten, [3, 4]); % Average SER across all SIC stage: 1 x SNR matrix;
    IEXIT = mean(IEXIT_vec, 1); %Average over SIC stages:  1 x SNR matrix;

    %% log-APPs: among the Nr trials, this selects the best according to the internal metric APP_cost_ten
    if nvArgs.save_logAPP
        log_APP_ten_temp1 = reshape(log_APP_ten, [], mod_alph.Q, N / S_SIC);
        log_APP_ten_temp2 = log_APP_ten_temp1(minRR, :, :);
        log_APP_ten_sel = reshape(log_APP_ten_temp2, S_SIC, L_SNR, N_block, mod_alph.Q, N / S_SIC); %Dim: S_SIC x L_SNR x N_block x Alph_Size (Q) x Blocklength (N)
    end

    %% Sanity-check transmit powers
    h_steep_convd = conv(conj(fliplr(h_tx)), h_tx); %Correlated TX filters
    P_tx_corr = sum(pU .* abs(S - sum(pU .* S)).^2) + abs(sum(pU .* S))^2 * (sum(h_steep_convd(1:N_os:end)) / N_os);
    SNRdB = x_dB_TX_PW + 10 * log10(P_tx_corr);

    P_tx_measured_dB = 10 * log10(mean(P_tx_wave_measured_ten, [1, 3, 4]));

    % For debugging purposes, compare with P_tx_vec_measured
    dev_tx_pow = abs(SNRdB - P_tx_measured_dB);

    if any(dev_tx_pow > 0.3) %Allow at most to be 0.3 dB apart
        warning('Potential problem with average power computation at TX.')
    end

    %% Print back:
    header_str = strcat("IqYXs_", string(1:S_SIC));
    [header_str(:), IqsYX_vec]

    if nvArgs.save_file

        %% CSV:
        %% Build results matrix:
        % SNR x Results
        HEADER = ['SNR', ...
                      'IqYX', ... %Overall rate
                      strcat("IqYXs_", string(1:S_SIC)), ... %Individual SIC stage rates
                      'EXIT', ...
                      strcat("EXITs_", string(1:S_SIC)), ...
                      'SER', ...
                      strcat("SERs_", string(1:S_SIC)), ...
                  ];

        RES = real([SNRdB; ...
                        IqYX; ...
                        IqsYX_vec; ...
                        IEXIT; ...
                        IEXIT_vec; ...
                        SER; ...
                        SERs_vec; ...
                    ])';

        writematrix([HEADER; RES], 'results/' + vver + '/'+filename + '.txt');

        %% Reduce filesize before storing
        if ~nvArgs.save_trace %Clear traces 
            clear *_iter*
        end

        if ~nvArgs.save_logAPP %Clear log_APP
            clear log_APP_ten*
        end

        clear rem_pil hnd_*
        clear L_A id idx2 gD Lambda_D arg_rnd_* idx_data idx_pil A Aup A GmU GmD GmV AD GdU AV h3 h4 Adata Apil Acomb_k Acomb_k2 phat u_ind w_intf P w_noisy w_noisy_up u0 w_star w_nonoise w_intf y_inter y_noisy n_tx
        clear s_intf rhat F arg_rnd simWDM nvArgs L_P h_steep_convd h_tx h_tx_lim  h_comb_pad a_comb
        clear alpha1_SEs alpha2_SEs beta1_SEs beta2_SEs rate_SEs z_star x_star n_rx u0_pil_padded x2d_t1 x2d_t1_pil u1_best
        clear D1 D2 Fh Fh_A Fh_B Fh_C Fh_D Fn
        clear P_tx_sym_measured P_tx_wave_measured
        clear APP_cost_ten  HDMI_ten reqiter_VAMP_lvec REQITER_ten init_cfg GZF y_inter_ds 
        clear log_APP_ten_temp1 log_APP_ten_temp2 log_APP 
        clear RXnoise minRR ps_filter cluster P_tx_sym_measured_ten P_tx_wave_measured_ten autoDampCFG *_lvec*

        save('results/' + vver + '/' + filename + '.mat'); % Save as .MAT
    end

end
