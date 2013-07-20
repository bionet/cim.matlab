%IDENTIFY_H_IDEAL_IAF_TRIG Identify channels in cascade with an asynchronous delta-sigma modulator in a MIMO circuit.
%   [PH_REC SPK_NUM] = IDENTIFY_H_ASDM_MIMO(T,T_PH,U,W,L,B,D,K,Z) identifies
%   the projection PH_REC of the linear filters H in cascade with an
%   asynchronous delta-sigma modulator (ASDM) in a multi-input multi-output
%   (MIMO) circuit. PH is the projection of H onto the reproducing kernel
%   Hilbert space (RKHS) of the trigonometric polynomials with bandwidth W and
%   order L, and is of time course T_PH. When L is equal to infinity, PH is the
%   The ASDM has a bias B, a Schmitt trigger threshold D, and an integration
%   capacitance K. The signal U with a time course T is encoded into the ASDM
%   output Z. The identified filter projection PH_REC is returned, as well as
%   the number of total spikes SPK_NUM used for identification.
%
%   Author:               Yevgeniy B. Slutskiy
%
%   Revision Author:      Chung-Heng Yeh
%
%   Copyright 2012-2014   Yevgeniy B. Slutskiy and Chung-Heng Yeh

% This version: using multiple signals only
function [Ph_rec, spk_num] = identify_h_asdm_mimo(t, t_Ph, u, W, L, b, d, k, z)

    dt = t(2)-t(1);                      % set the time step
    T  = 2*pi*L/W;                       % compute the period of the trigonometric space
    em = @(m,t) exp(1j*m*W/L*t)/sqrt(T); % declare function handle for the basis of the space
    M  = size(u,1);                      % get the number of filters
    N  = size(u,3);                      % get the number of signals

    u_l = fft(u,[],2)*sqrt(T)/numel(t);  % reconstruct the signal coefficients
    u_l = u_l(:, [end-L+1:end 1:L+1],:); % get the relevant Fourier coefficients

    % declare function handle for getting spike time
    get_tk = @(i) t( [diff(z(i,:))~=0 false] );
    % get spike time
    tk = arrayfun( get_tk, 1:N,'UniformOutput',false);
    % get the number of total tspikes
    spk_num = sum(cellfun(@numel,tk));
    % function handle for computing q of each stimuli
    com_q = @(i) sign(z(i,1))*(-1).^(2:numel(tk{i})).*(2*k*d - b*diff(tk{i}));
    % compute q of each stimuli and merge them to a single vector
    q = cell2mat(arrayfun( com_q, 1:N, 'UniformOutput',false))';

    % build the matrix B block by block; the matrix B is the block diagonal
    % matrix comprised of Phi^i; please refer to Eq. (13)
    B = zeros(numel(q),(2*L+1)*N);                 % allocate memory for the matrix B
    row_idx = 0;                                   % row index for each Phi^i
    for i=1:N
        % fill in each of Phi^i
        row_idx = row_idx(end)+(1:numel(tk{i})-1); % compute the row index of each Phi^i
        col_off = (i-1)*(2*L+1);                   % compute the column offset of each Phi^i 
        % for basis elements -L:0 (the rest are conjugate symmetric)
        for j=-L:0
            if j == 0
                B(row_idx,col_off+j+L+1) = diff(tk{i})';
            else
                B(row_idx,col_off+j+L+1) = sqrt(T)/(1j*j*W/L)*...
                     (em(j,tk{i}(2:end))-em(j,tk{i}(1:end-1)));
            end
        end
        % use the symmetry of the Phi matrix to fill in the remaining entries
        B(row_idx,col_off+(L+2:2*L+1)) = conj(B(row_idx,col_off+(L:-1:1)));
    end

    % build the matrix U block by block; please refer Eq. (13)
    U = zeros((2*L+1)*N,(2*L+1)*M);
    for i=1:N
        % fill in each of U^i
        col_off = 1+(i-1)*(2*L+1);
        for j=-L:L
            % fill each of u^i_l into U
            U( col_off+L+j, (j+L)*M+(1:M) ) = transpose(u_l(:,j+L+1,i));
        end
    end
    Phi = conj(B*U);

    % recover the Fourier coefficients of all filters
    h = pinv(Phi,1e-9)*q;
    % reshape h so that each column of h represents a single filter
    h = transpose(reshape(h,M,2*L+1));
    % allocate memory for the reconstructed filters
    Ph_rec = zeros(M,length(t_Ph));
    for m=1:M
        h_l_pad = [h(L+1:end,m); zeros(numel(t_Ph)-2*L-1,1); h(1:L,m)]';
        Ph_rec(m,:) = real(ifft(h_l_pad))*numel(t_Ph)/sqrt(T);
    end
    % shift the filters according to the preferred time course
    idx = floor(mod(t_Ph(1)-T/2,T)/dt)+1;
    Ph_rec = Ph_rec(:,[idx:end 1:idx-1]);

    display(['Number of spikes needed to invert the matrix: ' num2str(size(u,1)*(2*L+2))]);
    display(['Number of spikes used: ' num2str(spk_num)]);
end
