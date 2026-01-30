function [Iq_XY, log_APP] = comm_gmi(r1, nu_U1, u_ind, opt_flag, pU, data_ind, ind_sic, S_SIC, hnd_denoise_x1, rGVAMP)
    %Returns (s)_GMI and APP metrics

    if isnan(nu_U1) || isinf(nu_U1) %Handle NaN or inf
        nu_U1 = 1;
    end

    if opt_flag == 1 % Optimize s-parameter
        fun = @(opt_var) -comm_calc_rate(opt_var, r1, nu_U1, u_ind, pU, data_ind, ind_sic, S_SIC, hnd_denoise_x1, rGVAMP);

        fun = @(p) fun(p);
        options = optimset('Display', 'off', 'FunValCheck', 'off', 'MaxFunEvals', 30, 'MaxIter', 30);
        [param_opt] = fminbnd(fun, -5 * abs(nu_U1), 5 * abs(nu_U1), options); %Interval from nu_U1 +/- [5*nu_U1]

        %Iq_XY = -fun(param_opt); % [TODO: Updated Nov. 19]
        [Iq_XY, log_APP] = comm_calc_rate(param_opt, r1, nu_U1, u_ind, pU, data_ind, ind_sic, S_SIC, hnd_denoise_x1, rGVAMP);
        Iq_XY = max([0, Iq_XY]); %Clip

    elseif opt_flag == 0 % Do not optimize
        opt_var = 0;
        [Iq_XY, log_APP] = comm_calc_rate(opt_var, r1, nu_U1, u_ind, pU, data_ind, ind_sic, S_SIC, hnd_denoise_x1, rGVAMP);
        if length(Iq_XY) > 1
            error('Wrong dimensions.')
        end
        Iq_XY = max([0, Iq_XY]); %Clip

    end

end
