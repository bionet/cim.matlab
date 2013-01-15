%IDENTIFY_H_IDEAL_IAF Identify linear dendritic processing filter in cascade with a spiking neuron model.
%   [PH WINDOWS] = IDENTIFY_H_IDEAL_IAF(DT,T_PH,T,U,W,B,D,K,TK,TAU_1,
%   TAU_2,N_MAX) identifies the projection PH of the linear dendritic
%   processing filter H in cascade with an IAF neuron. PH is the projection
%   of H onto the input signal space and is of time course T_PH. The
%   ideal integrate-and-fire (IAF) neuron has a bias B, a firing threshold
%   D, and a capacitance K. The [Filter]-[Ideal IAF] neural circuit encodes
%   a signal U with a time course T_U into a sequence of spike times TK.
%   DT is the time-step of the numeric computation. At most N_MAX temporal
%   windows with length (TAU_2-TAU_1) are used in the identification.
%   IDENTIFY_H_IDEAL_IAF(...) returns the identified filter PH, as well as
%   the picked temporal windows WINDOWS.
%
%   [PH WINDOWS H_REC_N] = IDENTIFY_H_IDEAL_IAF(...,'Calc_MSE_N',true)   
%   identifies filter using the first i windows of WINDOWS, for i from 1 to
%   N_MAX where N is the maximum number of windows.
%
%   Author:               Yevgeniy B. Slutskiy
%
%   Revision Author:      Chung-Heng Yeh
%
%   Copyright 2010-2012   Yevgeniy B. Slutskiy and Chung-Heng Yeh

function [Ph windows varargout] = ...
         identify_h_ideal_iaf(dt, t_Ph, t, u, W, b, d, ... 
         k, tk, tau_1, tau_2, N_max, varargin)

    % Handle the optional input parameters, e.g., Calc_MSE_N     
    p = inputParser;
    addParamValue(p,'Calc_MSE_N',false,@islogical);
    p.KeepUnmatched = true;

    parse(p,varargin{:});
    UnmatchedParam = fieldnames(p.Unmatched);
    if ~isempty(UnmatchedParam)
        error(['"',UnmatchedParam{1},'" is not a valid parameter. ' ...
               'Please use "Calc_MSE_N" instead...']);
    end

    % Check the number of output arguments
    if p.Results.Calc_MSE_N
        error(nargoutchk(1, 3, nargout, 'struct'));
    else
        error(nargoutchk(1, 2, nargout, 'struct'));
    end 

    % Find the non-overlapping spike windows. Keep in mind that for every
    % window we will have to look at most tau_1 seconds into the future and at
    % most tau_2 seconds into the past. Need to guarantee that we have enough 
    % of the input signal to look at.

    % Compute the length of the temporal window
    W_length = dt*round((tau_2-tau_1)/dt); 
    % Get the first window anchor (the spike based on which we pick the window)
    window_anchor = tk(tk>=t(1) + W_length);    % find all feasible tk that can be used as window anchors
    window_start  = window_anchor(2) - dt;  	% get the start time of the first window

    % Initialize temporal window parameters
    windows = {};             % cell array of temporal windows
    tk_in_window = {};        % cell array of spikes in these window
    all_tk = [];              % array of spikes contained in all temporal window

    % While still have enough of the corresponding input signal, keep
    % generating temporal windows.
    while window_start+W_length<t(end)

        % Generate the temporal window
        window_end = window_start+W_length;
        windows{end+1} = [window_start window_end];

        % Find all spikes falling in the current temporal window
        tk_in_window{end+1} = tk( window_start<=tk & tk<=window_end );

        % Find the biggest gap in the combined spike train
        all_tk = sort([all_tk, tk_in_window{end} - window_start]); % sort combined spikes
        gaps = diff(all_tk);                                       % compute gaps
        biggest_gap_idx = find(gaps == max(gaps),1);               % find the first index for the biggest gap

        % Find the next anchor at least a gap_spike_time away from the end of
        % current temporal window
        gap_spike_time = dt*round((all_tk(biggest_gap_idx) + all_tk(biggest_gap_idx+1))/(2*dt));
        window_anchor = tk(tk>window_end + gap_spike_time);

        % If the next anchor exists (may run out of spikes)
        if numel(window_anchor)~=0
            window_start = window_anchor(1) - gap_spike_time;      % get the start time of next temporal window
        else
            break;
        end
    end
    
    % Find and Remove windows that contain fewer than 2 spikes     
    idx = cellfun( @numel, tk_in_window ) < 2;
    tk_in_window = tk_in_window(~idx);
    windows = windows(~idx);

    % Set the maximum number of windows to be used.
    N = min(numel(windows),N_max);                                       
    tk_in_window = tk_in_window(1:N);
    windows = windows(1:N);

    % Compute q of each window, and merge all q to a single vector
    q_N = cell2mat(cellfun( @(x) k*d-b*diff(x), tk_in_window, 'UniformOutput',false))';

    % Compute T (Corollary 3 of [1])
    T = (tau_1+tau_2)/2;

    % Compute all t_k^j-tau+T (Corollary 3 of [1])
    sk_N = cell2mat( cellfun( @(tk_in_win,win) tk_in_win(1:end-1)-mean(win)+T,...
                     tk_in_window,windows,'UniformOutput',false));

    % Initialize ensemble matrices G_N
    G_N = zeros(numel(sk_N));
    
    % Compute the matrix G = [G1; G2; ...; GN] block by block. Each block 
    % Gi = [G1i; G2i; ... Gki] associates with the spikes in the i-th 
    % window. Entries of each row G1i are integrals of different shifted 
    % version of u(t) on the same interval [t^i_k t^i_(k+1)]. To compute
    % every entry efficiently, we calculate the integral of u(t), U(t), one 
    % time, and find the integral of u(t-s) on [a b] by computing 
    % U(b-s)-U(a-s). The above idea can be vectorized for entries of the 
    % same row, since they have the same integral interval.
    
    int_u = cumtrapz(t,u);                  % compute the integral of u
    row = 1;                                % initialize row index
    for i = 1:N  
        % get the lower bound of the first integral. Need to convert the
        % domain from continuous value to discrete index.
        low_idx = time2idx( tk_in_window{i}(1)-sk_N-t(1), dt, 1, numel(u) );
        for k = 2:numel(tk_in_window{i})
            % get the upper bound of the integral interval
            up_idx = time2idx( tk_in_window{i}(k)-sk_N-t(1), dt, 1, numel(u) );
            % compute the explicit integral
            G_N(row,:) = int_u(up_idx)-int_u(low_idx);
            % the upper bound becomes the lower of the next interval
            low_idx = up_idx;              
            row = row + 1;                 % update row index 
        end
    end
    % Calculate ck_N    
    ck_N = pinv(G_N,1e-7)*q_N;   

    % Recover the filter
    Ph = synsig(dt,t_Ph, ck_N, sk_N, W ); 
    
    % If error as a function of the number of windows is requested
    if p.Results.Calc_MSE_N
        h_rec_N = zeros(N,length(t_Ph));
        G_ix_end = cumsum( cellfun( @numel, tk_in_window ) - 1 );
        for i = 1:N-1
            % get the matrix G corresponding to first i windows
            G = G_N(1:G_ix_end(i), 1:G_ix_end(i));      % get the matrix for a population of i neurons
            q = q_N(1:G_ix_end(i));                     % get the q_i
            ck = pinv(G,1e-7)*q;                        % calculate cks
            sk = sk_N(1:G_ix_end(i));                   % get sks
            h_rec_N(i,:) = synsig(dt,t_Ph, ck, sk, W ); % recover the filter
        end
        h_rec_N(N,:) = Ph(1,:);                         % the last entry is the recovery for all N windows     
        varargout = {h_rec_N};
    end
    
    %V = synsig( DT, T_SIG, S, S_T, OMEGA ) synthesizes the band-limited 
    %   signal V of time course T_SIG from samples S at samples time S_T.
    %   The time step is DT, and the frequency support of V is bounded by
    %   [-OMEGA, OMEGA].
    function v = synsig(dt, t_sig, s, s_t, Omega)
        v = zeros(size(t_sig));
        % Naive method 
        if dt >= 5e-6;
            for j = 1:numel(s)
                v = v + s(j)* Omega/pi * sinc( Omega/pi*(t_sig- s_t(j) ));
            end
        % FFT-based method. Faster and more accurate only for small dt. 
        else
            v( round((s_t-t_sig(1))/dt)+1 ) = s;
            t_sinc = dt*(-numel(t_sig)+1:numel(t_sig));
            v = fftfilt(v, sinc(t_sinc*Omega/pi)*Omega/pi);
            v = v(end-length(t_sig)+1:end);
        end
    end

    %IDX = time2idx( TIME, DT, LOWERBOUND, UPPERBOUND ) discretizes the 
    %   continuous-valued TIME to IDX with step DT. IDX is hard limited by 
    %   the upper bound UPPERBOUND and lower bound LOWERBOUND.
    function idx = time2idx( time, dt, lowerBound, upperBound )
        idx = round( time/dt ) + 1;
        idx(idx < lowerBound) = lowerBound;
        idx(idx > upperBound) = upperBound;
    end
end % end of function
