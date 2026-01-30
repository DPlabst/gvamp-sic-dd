function [p_vec, P] = gen_prec(ch, Pcfg, N)

    p_vec = []; %Initialize

    if ch == 0 ... % For fiber channels
            && any(contains(["cDC", "cDR", "cU", "cO", "_"], Pcfg)) % Some all-pass precoders

        p_abs = 1; %Overall channel needs to be band-edge symmetric
        arg_rnd_full = zeros(N, 1); %Phase

        if strcmp(Pcfg, 'cU') %Choose uniformly in [0,2*pi]
            arg_rnd_full = rand(N, 1) * 2 * pi;

        elseif strcmp(Pcfg, 'cO') %Choose uniformly in [0,2*pi], but with real-value constraint
            arg_rnd_halve = 2 * pi * rand(N / 2, 1);
            phi0 = pi * (randi(2, 2, 1) - 1);
            arg_rnd_full = [phi0(1); arg_rnd_halve(2:end); phi0(2); -flipud(arg_rnd_halve(2:end))]; %Hermitian symmetric

        elseif strcmp(Pcfg, 'cDC') %Mimic CD
            fvec = fftshift(1 / N * (-N / 2:1:(N / 2 - 1))).';
            CDgain = 50; %Adjust artificial CD
            arg_rnd_full = (fvec * CDgain).^2;

        elseif strcmp(Pcfg, 'cDR') %Mimic CD, but with real-value constraint
            fvec = 1 / N * [(0:1:(N / 2 - 1))].';
            CDgain = 6; %Adjust artificial CD
            arg_rnd_halve = (fvec * CDgain).^2;
            phi0 = pi * (randi(2, 2, 1) - 1);
            arg_rnd_full = [phi0(1); arg_rnd_halve(2:end); phi0(2); -flipud(arg_rnd_halve(2:end))]; %Hermitian symmetric

        elseif strcmp(Pcfg, '_')
            %No precoding
        end

        %% Construct precoder
        p_vec = p_abs .* exp(1j * arg_rnd_full); %
        Ptx = sum(abs(p_vec).^2); %Normalize power
        p_vec = sqrt(N / Ptx) .* p_vec;

        % Save time: Don't do explicit
        F = 1 / sqrt(N) * dftmtx(N); %Unitary DFT matrices
        P = F' * (p_vec .* F);

        P = diag(1 ./ sqrt(diag(P * P'))) * P; %Fine for all circulant matrices!

    else
        error('Combination of precoder and channel is not supported.')
    end

end
