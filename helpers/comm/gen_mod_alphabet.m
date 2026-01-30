function [S, pX, mod_alph] = gen_mod_alphabet(mod_alph_str)

    %% Generates constellation symbols and PMF
    modf_str_arr = split(mod_alph_str, '-');

    if length(modf_str_arr) == 4
        mod_alph = struct('name', char(modf_str_arr(2)), 'Q', str2double(modf_str_arr(1)), 'off', str2double(modf_str_arr(3)), 'nu', str2double(modf_str_arr(4))); %For shaping
    elseif length(modf_str_arr) == 3
        mod_alph = struct('name', char(modf_str_arr(2)), 'Q', str2double(modf_str_arr(1)), 'off', str2double(modf_str_arr(3)), 'nu', 0); %For offset
    else
        mod_alph = struct('name', char(modf_str_arr(2)), 'Q', str2double(modf_str_arr(1)), 'off', 0, 'nu', 0);
    end

    %% Available alphabets
    if strcmp(mod_alph.name, 'ASK')
        % ------ Q-ASK ----
        S = linspace(-1, 1, mod_alph.Q) + mod_alph.off;
        pX = 1 / mod_alph.Q * ones(1, mod_alph.Q); %Generate uniform input distribution

    elseif strcmp(mod_alph.name, 'MB') % "Maxwell-Boltzmann"
        % ------ Gaussian ----
        S = linspace(-1, 1, mod_alph.Q);
        gam = mod_alph.nu;
        pX = exp(-S.^2 * (gam));
        pX = pX ./ sum(pX); %Normalize
        pX(pX < 1e-8) = 1E-8;
        pX = pX ./ sum(pX); %Normalize again

        alpha = sqrt(var(S, 1) / sum(pX .* abs(S).^2)); %Scale power
        S = alpha * S + mod_alph.off; %Symmetric with offset

    elseif strcmp(mod_alph.name, 'QAM') % 
        m = sqrt(mod_alph.Q);
        Sb = linspace(- (m), (m), (m));
        S = Sb + 1j * Sb';
        S = S(:).';
        pX = 1 / mod_alph.Q * ones(1, mod_alph.Q); %Generate uniform input distribution
    end

    S = 1 / sqrt(sum(pX .* abs(S).^2)) * S; %Normalize constellation to unit power

    %  ------ Sanity Check power constraint -----
    if (sum(pX .* abs(S).^2) - 1) > 5 * eps
        error('Power constraint violated.');
    end

end
