function [IqYX, SER, NMSE, IqYX_itvec, NMSE_itvec, ...
              damp_cost, reqIter, u1best, log_APP] ...
    = run_gvamp_asym(nu_init, ...
    p1, nu_W1, r1, nu_U1, ...
    damp_cfg, ...
    hnd_denoise_w1, ...
    hnd_denoise_u1, ...
    hnd_lmmse_u2w2, ...
    hnd_sgmi, ...
    hnd_ser, ...
    hnd_nmse, ...
    hnd_damp_cost, ...
    damp, ... %mean damping 
    dampV, ... %variance damping 
    autoDampCFG, ...
    nitMin, ...
    nitMax, ...
    tol, ...
    flagConsOutput, ...
    u_star, ...
    w_star, ...
    idx_data, ...
    rGVAMP, ...
    VAMP_dbg_mode, ...
    EXIT_dbg_mode, ...
    ind_sic ...
)

%% Note
% This version of GVAMP is assymetric and simplified to one Turbo loop. It works when the extrinsic (r1, nuU1) is independent
% of the exstrinsic (r2, nuU2). This is for example the case for allpass channels.
% Otherwise, one should use the more general GVAMP [https://arxiv.org/abs/1612.01186] which runs two Turbo loops.

dbg_nu_U1_mc = [];
dbg_nu_W1_mc = [];

% VAMP
i = 0; % iteration counter (for damping)
k = 0; % true iteration counter (without damping)
fail_count = 0;

stop = false; % stop flag

r1_prev = [];
r2_prev = [];
p1_prev = [];
p2_prev = [];
u1_prev = [];
u2_prev = [];
w1_prev = [];
w2_prev = [];
alpha1_prev = [];
nu_U1_prev = [];
nu_W1_prev = [];
nu_U2_prev = [];
nu_W2_prev = [];

dbg_alpha1_alg = [];
dbg_alpha1_mc  = [];
dbg_alpha2_alg = [];
dbg_alpha2_mc = [];
dbg_nu_U1_alg = [];
dbg_nu_U1_mc = [];
dbg_nu_U2_alg = [];
dbg_nu_U2_mc = [];
dbg_beta1_alg = [];
dbg_beta1_mc = [];
dbg_beta2_alg = [];
dbg_beta2_mc = [];
dbg_nu_W1_alg = [];
dbg_nu_W1_mc = [];
dbg_nu_W2_alg = [];
dbg_nu_W2_mc = [];

NMSE_itvec = inf * ones(1, nitMax);
SER_itervec = inf * ones(1, nitMax);
IqYX_itvec = 0 * ones(1, nitMax);
IqYX_itvec_k = 0 * ones(1, nitMax);
damp_cost_itvec = inf * ones(1, nitMax);

% ------ Taken from GAMPMATLAB package ------
gamMin = 1E-8; %1E-8
gamMax = 1E+11; %1E+11
clipVar = @(nu0) 1 / min(max(1 / nu0, gamMin), gamMax); %Clip variances

dampMean = @(in, old, damp) damp * in + (1 - damp) * old; %Damp mean
dampVar = @(in, old, damp) 1 / (damp * 1 / in + (1 - damp) * 1 / old); %Damp precision
% -------------------------------------------

nu(1) = nu_init(1); %Optical noise
nu(2) = nu_init(2); %Electrical noise

trig = 0;
n_trig = 0;

while ~stop

    i = i + 1; %Iteration counter for damping
    k = k + 1; %Iteration counter for complexity

    %---------------------------------------------------------------------------
    %   ___        _               _     ____                   _
    %  / _ \ _   _| |_ _ __  _   _| |_  |  _ \  ___ _ __   ___ (_)___  ___ _ __
    % | | | | | | | __| '_ \| | | | __| | | | |/ _ \ '_ \ / _ \| / __|/ _ \ '__|
    % | |_| | |_| | |_| |_) | |_| | |_  | |_| |  __/ | | | (_) | \__ \  __/ |
    %  \___/ \__,_|\__| .__/ \__,_|\__| |____/ \___|_| |_|\___/|_|___/\___|_|
    %                 |_|
    % --------------------------------------------------------------------------
    if VAMP_dbg_mode >= 2 || EXIT_dbg_mode
        dbg_nu_W1_alg(i) = nu_W1;
        dbg_nu_W1_mc(i) = var(p1 - w_star.');
    end

    [w1hat, beta1_vec] = hnd_denoise_w1(p1, nu_W1, nu, i); %Calculate (w1hat, beta1)
    beta1 = mean(beta1_vec);

    if VAMP_dbg_mode >= 2
        dbg_beta1_alg(i) = beta1;
        dbg_beta1_mc(i) = var(w1hat - w_star.');
    end

    if damp_cfg(1) && ~isempty(w1_prev)
        w1hat = dampMean(w1hat, w1_prev, damp); % Damp
    end

    % Convert to extrinsics
    p2 = (nu_W1 * w1hat - beta1 * p1) / (nu_W1 - beta1);
    nu_W2 = clipVar(beta1 * nu_W1 / (nu_W1 - beta1)); %clip variances if required

    if VAMP_dbg_mode >= 2 || EXIT_dbg_mode
        dbg_nu_W2_alg(i) = nu_W2;
        dbg_nu_W2_mc(i) = var(p2 - w_star.');
    end

    if damp_cfg(2) && ~isempty(nu_W2_prev)
        nu_W2 = dampVar(nu_W2, nu_W2_prev, dampV); % Damp
    end

    if damp_cfg(3) && ~isempty(p2_prev)
        p2 = dampMean(p2, p2_prev, damp); %Damp
    end

    % -------------------------------------------
    %  _     __  __ __  __ ____  _____   ___
    % | |   |  \/  |  \/  / ___|| ____| |_ _|
    % | |   | |\/| | |\/| \___ \|  _|    | |
    % | |___| |  | | |  | |___) | |___   | |
    % |_____|_|  |_|_|  |_|____/|_____| |___|
    %
    % --------------------------------------------
    % Part I of LMMSE denoising
    % Dummy inputs, because extrinsic (nuU1, r1) does not depend on (nuU2, r2)
    dummy_r2 = 1; %Scalar works
    dummy_nu_U2 = 1;
    [u2hat, alpha2, ~, ~, u2d_hat_part] = hnd_lmmse_u2w2(dummy_r2.', dummy_nu_U2, p2.', nu_W2, i, trig, rGVAMP, []); % Compute (u2hat, alpha2)

    if VAMP_dbg_mode >= 2
        dbg_alpha2_alg(i) = alpha2;
        dbg_alpha2_mc(i) = var(u2hat - u_star(idx_data)); %all column
    end

    u2hat = u2hat.'; %Transpose
    alpha2 = alpha2.';

    if damp_cfg(10) && ~isempty(u2_prev)
        u2hat = dampMean(u2hat, u2_prev, damp); %damp
    end

    % Calculate extrinsics
    r1 = (dummy_nu_U2 * u2hat - alpha2 * dummy_r2) / (dummy_nu_U2 - alpha2);
    nu_U1 = clipVar((alpha2 * dummy_nu_U2) / (dummy_nu_U2 - alpha2));

    if damp_cfg(11) && ~isempty(nu_U1_prev)
        nu_U1 = dampVar(nu_U1, nu_U1_prev, dampV); %damp
    end

    if damp_cfg(12) && ~isempty(r1_prev)
        r1 = dampMean(r1, r1_prev, damp); %damp
    end

    if VAMP_dbg_mode >= 2
        dbg_nu_U1_alg(i) = nu_U1;
        dbg_nu_U1_mc(i) = var(r1 - u_star(idx_data).');
    end

    % -----------------------------------------------------------------------
    %  ___                   _     ____                   _
    % |_ _|_ __  _ __  _   _| |_  |  _ \  ___ _ __   ___ (_)___  ___ _ __
    %  | || '_ \| '_ \| | | | __| | | | |/ _ \ '_ \ / _ \| / __|/ _ \ '__|
    %  | || | | | |_) | |_| | |_  | |_| |  __/ | | | (_) | \__ \  __/ |
    % |___|_| |_| .__/ \__,_|\__| |____/ \___|_| |_|\___/|_|___/\___|_|
    %           |_|
    % -----------------------------------------------------------------------
    [u1hat, alpha1_vec] = hnd_denoise_u1(r1, nu_U1, rGVAMP); %Calculate (u1hat, alpha1)
    alpha1 = mean(alpha1_vec); %Scalar

    if VAMP_dbg_mode >= 2
        dbg_alpha1_alg(i) = alpha1;
        dbg_alpha1_mc(i) = var(u1hat - u_star(idx_data).'); %row vectors
    end

    if i > 1 && k > 1 %Calculate and save some metrics 
        NMSE_itvec(i) = hnd_nmse(u1hat); %NMSE
        SER_itervec(i) = hnd_ser(u1hat.'); %SER
        [IqYX_tmp, M_APP_tmp] = hnd_sgmi(r1, nu_U1, rGVAMP); %GMI
        IqYX_itvec(i) = IqYX_tmp;
        IqYX_itvec_k(k) = IqYX_tmp;
    end

    damp_cost_itvec(i) = hnd_damp_cost(r1, nu_U1, u1hat, alpha1); % Cost for automatic damping

    if damp_cfg(4) && ~isempty(u1_prev)
        u1hat = dampMean(u1hat, u1_prev, damp); %damp
    end

    % Calculate extrinsics
    r2 = (nu_U1 * u1hat - alpha1 * r1) / (nu_U1 - alpha1);
    nu_U2 = clipVar((nu_U1 * alpha1) / (nu_U1 - alpha1));

    if VAMP_dbg_mode >= 2
        dbg_nu_U2_alg(i) = nu_U2;
        dbg_nu_U2_mc(i) = var(r2 - u_star(idx_data).');
    end

    if damp_cfg(5) && ~isempty(nu_U2_prev)
        nu_U2 = dampVar(nu_U2, nu_U2_prev, dampV); %damp 
    end

    if damp_cfg(6) && ~isempty(r2_prev)
        r2 = dampMean(r2, r2_prev, damp); %damp
    end

    % --------------------------------------------------------
    %  _     __  __ __  __ ____  _____   ___ ___ 
    % | |   |  \/  |  \/  / ___|| ____| |_ _|_ _|
    % | |   | |\/| | |\/| \___ \|  _|    | | | | 
    % | |___| |  | | |  | |___) | |___   | | | | 
    % |_____|_|  |_|_|  |_|____/|_____| |___|___|                                    
    % -------------------------------------------------------
    % Part II of LMMSE denoising
    [~, ~, w2hat, beta2] = hnd_lmmse_u2w2(r2.', nu_U2, p2.', nu_W2, i, trig, rGVAMP, u2d_hat_part); %Compute (w2hat, beta2) and reused partially computed u2d_hat_part

    if VAMP_dbg_mode >= 2
        dbg_beta2_alg(i) = beta2;
        dbg_beta2_mc(i) = var(w2hat - w_star); %all column
    end

    w2hat = w2hat.';
    beta2 = beta2.';

    if damp_cfg(7) && ~isempty(w2_prev)
        w2hat = dampMean(w2hat, w2_prev, damp); %damp 
    end

    %Calculate extrinsics 
    p1 = (nu_W2 * w2hat - beta2 * p2) / (nu_W2 - beta2);
    nu_W1 = clipVar((beta2 * nu_W2) / (nu_W2 - beta2));

    if damp_cfg(8) && ~isempty(nu_W1_prev)
        nu_W1 = dampVar(nu_W1, nu_W1_prev, dampV); %damp 
    end

    if damp_cfg(9) && ~isempty(p1_prev)
        p1 = dampMean(p1, p1_prev, damp); %damp 
    end

    % --------------------------------------------------------------
    % --------------------------------------------------------------
    % --------------------------------------------------------------
    % ------------------------ DONE --------------------------------
    % --------------------------------------------------------------
    % --------------------------------------------------------------
    % --------------------------------------------------------------

    %% Plot progress
    if i > 1
        change = (abs(damp_cost_itvec(i) - damp_cost_itvec(i - 1)) ./ abs(damp_cost_itvec(i)));
    else
        change = inf;
    end
    
    if flagConsOutput == 1 %Show progress on console
        if mod(i, 1) == 0 && i ~= 1
            fprintf("Iter:%d " + ... % Report true iteration count k 
                "Rate:%.3f " + ...
                "SER:%.3f " + ...
                "NMSE:%.3f " + ...
                "DampCost:% .3f " + ...
                "Step: %.2f " + ...
                "\n", k, IqYX_itvec(i), SER_itervec(i), 10 * log10(NMSE_itvec(i)), damp_cost_itvec(i), damp);
        end
    end

    % Stopping criterion 
    if i > nitMin && (k >= nitMax || damp <= autoDampCFG.stepTol || dampV <= autoDampCFG.stepTol || (change < tol)) 

        stop = true;
        u1best = u1hat; %Do not return
        IqYX = IqYX_itvec(i); %rate
        [~, log_APP] = hnd_sgmi(r1, nu_U1, rGVAMP); %GMI
        SER = SER_itervec(i); %ser
        damp_cost = damp_cost_itvec(i); %damping cost 
        NMSE = NMSE_itvec(i); %nmse 
        reqIter = find(IqYX_itvec_k >= 0.995 * IqYX, 1); %Number of iterations to reach 99.5 % of ultimate rate

        % Some plots
        plot_gvamp_dbg(VAMP_dbg_mode, EXIT_dbg_mode, n_trig, IqYX_itvec, IqYX, i, ...
                       dbg_alpha1_alg, dbg_alpha1_mc, ...
                       dbg_alpha2_alg, dbg_alpha2_mc, ...
                       dbg_nu_U1_alg, dbg_nu_U1_mc, ...
                       dbg_nu_U2_alg, dbg_nu_U2_mc, ...
                       dbg_beta1_alg, dbg_beta1_mc, ...
                       dbg_beta2_alg, dbg_beta2_mc, ...
                       dbg_nu_W1_alg, dbg_nu_W1_mc, ...
                       dbg_nu_W2_alg, dbg_nu_W2_mc); 

    else 
        % Fail: Decrease damping and repeat step if...
        if any(isnan(r1)) ... %Message are NaN OR
                || (((i - 1) - autoDampCFG.stepWindow > 1) ... % iteration count it larger than step window
                && ~(damp <= autoDampCFG.stepMin) ... %mean damping is smaller than threshold 
                && ~(dampV <= autoDampCFG.stepMin) ... %variance damping is smaller than threshold 
                && ~(damp_cost_itvec(i) < max(damp_cost_itvec((i - 1) - autoDampCFG.stepWindow:i - 1))) ... % damp cost at iteration i is not smaller than maximum cost of previous stepWindow iterations
                && i >= autoDampCFG.dampoffs) %if damping 

            %Count failures; from GAMPMATLAB package
            fail_count = fail_count + 1;
            if fail_count > autoDampCFG.maxBadSteps %If number of failures > max. number of bad steps
                fail_count = 0;
                autoDampCFG.stepMax = max(autoDampCFG.stepMin, autoDampCFG.maxstepDecr * autoDampCFG.stepMax); %decrease maximum step size
            end

            tmp_damp = max(autoDampCFG.stepMin, autoDampCFG.stepDec * damp); %adjust damping 
            tmp_dampGam = max(autoDampCFG.stepMin, autoDampCFG.stepDec * dampV); %adjust damping 

            % Activate one-time trigger, when iteration count larger than 10 * exp(- (ind_sic - 1))
            % -> Deactivates annealing in LMMSE denoiser
            % -> Do not repeat step when trigger gets activated 
            if tmp_damp < 1 && trig == 0 && i > 10 * exp(- (ind_sic - 1)) 
                fail_count = 0;
                autoDampCFG.stepMax = 1; %Reset maximum step size 
                trig = 1;
                n_trig = i; %Iteration number where trigger got activated (and annealing got turned OFF)
            else %Otherwise, repeat step 

                damp = tmp_damp;
                dampV = tmp_dampGam;
                i = i - 1; % Unsuccessful step (go back) 

                %Go step back
                r1 = r1_prev;
                r2 = r2_prev;
                p1 = p1_prev;
                p2 = p2_prev;
                u1hat = u1_prev;
                u2hat = u2_prev;
                w1hat = w1_prev;
                w2hat = w2_prev;
                alpha1 = alpha1_prev;
                nu_U1 = nu_U1_prev;
                nu_W1 = nu_W1_prev;
                nu_U2 = nu_U2_prev;
                nu_W2 = nu_W2_prev;

                r1_prev = r1_prev2;
                r2_prev = r2_prev2;
                p1_prev = p1_prev2;
                p2_prev = p2_prev2;
                u1_prev = u1_prev2;
                u2_prev = u2_prev2;
                w1_prev = w1_prev2;
                w2_prev = w2_prev2;
                alpha1_prev = alpha1_prev2;
                nu_U1_prev = nu_U1_prev2;
                nu_W1_prev = nu_W1_prev2;
                nu_U2_prev = nu_U2_prev2;
                nu_W2_prev = nu_W2_prev2;
            end

        else
            % If successful
            if i >= autoDampCFG.dampoffs
                damp = min(autoDampCFG.stepMax, autoDampCFG.stepInc * damp); %increase damp 
                dampV = min(autoDampCFG.stepMax, autoDampCFG.stepInc * dampV); %increase damp V
            end

            %Update older variables
            r1_prev2 = r1_prev;
            r2_prev2 = r2_prev;
            p1_prev2 = p1_prev;
            p2_prev2 = p2_prev;
            u1_prev2 = u1_prev;
            u2_prev2 = u2_prev;
            w1_prev2 = w1_prev;
            w2_prev2 = w2_prev;
            nu_U1_prev2 = nu_U1_prev;
            nu_W1_prev2 = nu_W1_prev;
            nu_U2_prev2 = nu_U2_prev;
            nu_W2_prev2 = nu_W2_prev;
            alpha1_prev2 = alpha1_prev;

            %Update old variables
            r1_prev = r1;
            r2_prev = r2;
            p1_prev = p1;
            p2_prev = p2;
            u1_prev = u1hat;
            u2_prev = u2hat;
            w1_prev = w1hat;
            w2_prev = w2hat;
            alpha1_prev = alpha1;
            nu_U1_prev = nu_U1;
            nu_W1_prev = nu_W1;
            nu_U2_prev = nu_U2;
            nu_W2_prev = nu_W2;

        end

    end

end
