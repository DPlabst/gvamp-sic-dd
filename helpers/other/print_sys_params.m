function print_sys_params(mod_alph, ps_filter, TXnoise, RXnoise, L_fiber, R_sym, nvArgs)
    % Print system configuration in a nicely formatted table

    fprintf('%s\n', '*** C-band Fiber-Optic Communication System ***');
    fprintf('%s\n', repmat('-', 1, 60));

    fprintf('%-25s | %-30s\n', 'Parameter', 'Value');
    fprintf('%s\n', repmat('-', 1, 60));
    
    % Modulation info
    fprintf('%-25s | %d - %s\n', 'Modulation', mod_alph.Q, mod_alph.name);
    fprintf('%-25s | %.2f\n', 'Constellation offset', mod_alph.off);
    fprintf('%-25s | %.2f\n', 'Gaussian Shaping', mod_alph.nu);
    
    % Transmitter pulse shape and precoder 
    fprintf('%-25s | %s\n', 'Pulse shape', ps_filter.name);
    fprintf('%-25s | %.3f\n', 'Roll-off factor', ps_filter.alpha);

    if strcmp(nvArgs.P, 'cO')
        fprintf('%-25s | %s\n', 'Precoder', 'Real orth. (FFT) + rand. phase of eigs.');
    elseif strcmp(nvArgs.P, '_')
        fprintf('%-25s | %s\n', 'Precoder', 'None');
    end 

    % TX noise
    if TXnoise.var_n_pre > 0
        fprintf('%-25s | %.3f\n', 'TX AWGN variance', TXnoise.var_n_pre);
    else
        fprintf('%-25s | %s\n', 'TX AWGN variance', 'None');
    end

    % RX noise
    if RXnoise.var_n_post > 0
        fprintf('%-25s | %.3f\n', 'RX AWGN variance', RXnoise.var_n_post);
    else
        fprintf('%-25s | %s\n', 'RX AWGN variance', 'None');
    end

    % Fiber link
    if L_fiber > 0
        fprintf('%-25s | %.1f km\n', 'Fiber length', L_fiber / 1e3);
        fprintf('%-25s | %.1f GBd\n', 'Symbol rate', R_sym / 1e9);
    else
        fprintf('%-25s | %s\n', 'Fiber', 'No fiber');
    end

    % PD RX filter
    if nvArgs.RxBW
        fprintf('%-25s | %s\n', 'PD Filter', 'Active');
    else
        fprintf('%-25s | %s\n', 'PD Filter', 'None');
    end

    % RX filter
    if nvArgs.ADCfilt
        fprintf('%-25s | %s\n', 'ADC Filter', 'Active');
    else
        fprintf('%-25s | %s\n', 'ADC Filter', 'None');
    end

    % RX detector
    fprintf('%-25s | %s\n', 'Detector', 'Square-law detector');

    fprintf('%s\n\n', repmat('-', 1, 60));

    % RX filter
    if nvArgs.run_simp
        fprintf('%s\n', 'Note: Running simplified GVAMP without s-GMI optimization and debiasing');
    end

    pause(0.5)

end
