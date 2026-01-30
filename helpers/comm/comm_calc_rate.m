function [MI, log_b_x1_norm] = comm_calc_rate(opt_var, r1, nu_U1, u_ind, pX, data_ind, ind_sic, S_SIC, hnd_denoise_u1, rGVAMP)

    % Optimized variance
    nu_U1 = max(1E-15, nu_U1 + opt_var); %Lower bound variance
    idx_data = 1:length(r1);
    idx_stage_s = idx_data(1:(S_SIC - (ind_sic - 1)):end);

    %% Stable implementation in log domain
    [~, ~, log_b_x1_norm] = hnd_denoise_u1(r1(idx_stage_s), nu_U1, rGVAMP);
    
    log_b_x1_norm_res = log_b_x1_norm.';

    u_fin = u_ind(data_ind(idx_stage_s)).'; %Data indices corresponding to current SIC stage
    q_XgY_sel_log2 = +log2(exp(1)) * log_b_x1_norm_res(sub2ind(size(log_b_x1_norm_res), 1:size(log_b_x1_norm_res, 1), u_fin)).'; %This does the same as a for loop

    hqXgY = -mean(q_XgY_sel_log2);
    MI = sum(-pX .* log2(pX)) - hqXgY; %Estimate mismatched rate

end
