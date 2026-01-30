function [idx_data, idx_pil, Adata, Apil, Nprime] = gen_SIC_mat(N_pil, N_span, L_sequ, A, ind_sic, S_SIC)
    %Generate pilot indices

    idx_data = 1:L_sequ; %Indices of data
    idx_pil_edge = [round(linspace(L_sequ - N_pil + 1, L_sequ, N_pil))]; %add pilots at edge (if required)
    idx_pil_sic = [];

    for jj = 1:ind_sic - 1
        idx_pil_sic = [idx_pil_sic, 1 + (jj - 1):S_SIC:L_sequ - N_pil]; %
    end

    idx_pil = [idx_pil_edge, idx_pil_sic]; %Stitch together known symbols through pilots and SIC
    idx_data(idx_pil) = []; %Remove indices of known data

    if ~isempty(A)
        Adata = A(:, idx_data); %Columns corresponding to unknown data
        Apil = A(:, idx_pil); %Columns corresponding to known data (pilots or known interference)
    else
        Adata = [];
        Apil = []; 
    end 

    Nprime = length(idx_data); %size(Adata, 2);

    SE = (L_sequ) / (N_span + L_sequ); %Spectral efficiency

end
