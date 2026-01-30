% <header>
% /*********************************************************************************************************
%  * Information Rates of Approximate Message Passing for Bandlimited Direct-Detection Channels
%  * Daniel Plabst; Institute for Communications Engineering (ICE), Technical University of Munich, Germany
%  * Mohamed Akrout; University of Manitoba, Winnipeg, Canada
%  * Version r2: Date: 2026-01-29
%  *********************************************************************************************************/
% </header>
clear; tic;

for modstr = ["ASK-0.2"]
    for M = [16]
        [filename, vver] = gvamp( ...
            m = num2str(M) + "-" + modstr, ... % Modulation
            ch = 0, ... % Optical channel
            P = 'cO', ... % Precoder: cO=circul. orthogonal (real)
            runGVAMP = 1, ... % Run GVAMP
            runEXIT = 1, ... % Run EXIT
            run_simp = 0, ... % Run simplified GVAMP without s-GMI optimization 
            S = 1, ... % Number of SIC stages
            ...%st = [1,2,4], ... % Compute subset of all stages
            Ptx = '10', ... % Transmit powers in dB, e.g.,'-5:1:20' (mind the SNR definition)
            n = 4 * 512, ... % Block length
            Nb = 4, ... % Number of blocks to be transmitted
            Nr = 1, ... % Number of recovery trials per block
            it = 250, ... % Maximum number of GVAMP iterations
            itEXIT = 100, ... % Maximum number of EXIT iterations
            d1 = 1, ... % Initial damping mean
            d2 = 1, ... % Initial damping mean
            aD = 1, ... % Autodamp
            ps = 'RRC', ... % Pulseshaping filter
            alp = 0.01, ... % Roll-off factor
            Nsp = 250, ... % Pulseshaping filter span
            Npi = 0, ... % Use Npi * Nsp pilots at block edge
            L = 4E3, ... % Fiber length in [m]
            Rs = 300E9, ... % Symbol rate in [Bd] or [Sym/s]
            nu = [0, 1], ... % Noise variances of [optical white CSCG, electrical real WGN].
            RxBW = 0, ... % [0/1]: Apply PD filter to filter out-of-band noise 
            ADCfilt = 1, ... %[0/1]: Apply ADC filter with bandwidth [-B, +B]
            outdencfg = {[0, 1], 'SPA'}, ... % Use scaled [nu(1), nu(2)] in output denoiser, calc. generalized Xi2 via method {SPA, Davies} when electrical noise > 0
            dampcfgStepW = [3, 3, 1], ... [minWindow, maxWindow, [0/1]:reduce step window with increasing SIC stage]
            annealcfg = {[100, 1], [0.47, 0.20, eps], [0.015]}, ... % Cell-array {[n_it_thresh, reduce threshold vs. SIC stage], [alpha, beta, minAnneal], [anneal if iter > switchIterThresh]}
            cons = 1, ... % Console output
            save_file = 1, ... % Save results in MAT and TXT file
            VAMP_dbg_mode = 1, ... % [0,1,2]; 0=none, 1=plot rate trace, 2=plot internal variance vs. actual variances
            EXIT_dbg_mode = 0, ... % [0,1]; 0=none, 1=plot exit charts
            save_trace = 0 ...
        );
    end
end

toc;
load('results/' + vver + '/' + filename + '.mat'); %load results
