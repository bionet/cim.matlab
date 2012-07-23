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
        error(['"',UnmatchedParam{1},'" is not a valid paprameter. ' ...
               'Please use "Calc_MSE_N" instead...']);
    end

    % Check the number of ouput arguments
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


    % Compute total number of spikes that are contained in temopral windows
    num_tk = sum( cellfun(@numel,tk_in_window) );


    % Comupte the matrix G_N

    G_N = zeros(num_tk-numel(windows));   % initialize ensemble matrices G_N
    q_N = zeros(num_tk-numel(windows),1); % initialize ensemble matrices q_N

    shifted_u = zeros(size(t));
    % build the matrix G_N column by column hence the integral of the shifted
    % stimulus will only be computed once.
    col_idx = 1;


    for i = 1:N                 % column block index
        % get the shift imposed by spikes from the window i
        shift = round((tk_in_window{i} - windows{i}(1) + tau_1)/dt);
        % for every column
        for k_cntr=1:length(tk_in_window{i})-1              
            if shift(k_cntr)>0
                shifted_u = [zeros(1,shift(k_cntr)) u(1:end-shift(k_cntr))];
            else
                shifted_u = [u(abs(shift(k_cntr))+1:end) zeros(1,abs(shift(k_cntr)))];
            end
            % compute the integral of shifted stimulus
            shifted_u_int = dt*cumtrapz(shifted_u);

            % compute G_N, row block by row block
            row_offset = 0;
            for j=1:N
               idx = round( ( tk_in_window{j} - t(1) ) / dt );
               idx( idx<1 ) = 1;
               idx( idx>length(t) ) = length(t);

               % Read out the value from the integral of the shifted stimulus
               G_N( row_offset+1:row_offset+numel(tk_in_window{j})-1, col_idx ) ...
                   = diff( shifted_u_int(idx) );
               row_offset = row_offset + numel(tk_in_window{j})-1;
            end 
            q_N(col_idx) = k*d - b*diff( tk_in_window{i}(k_cntr:k_cntr+1) );
            col_idx = col_idx + 1;
        end            
    end

    % Calculate ck_N    
    ck_N = pinv(G_N)*q_N;                      


    % Recover the filter
    Ph = zeros(size(t_Ph));                
    ck_counter = 1;
    for j=1:N
        for m=1:length(tk_in_window{j})-1
            Ph = Ph + ck_N(ck_counter) .* W/ pi * ...
                    sinc( W*(t_Ph - tk_in_window{j}(m) + windows{j}(1) - tau_1 )/pi);
            ck_counter = ck_counter + 1;
        end
    end


    % If error as a function of the number of windows is requested
    if p.Results.Calc_MSE_N
        h_rec_N = zeros(N,length(t_Ph));
        G_ix_end = length(tk_in_window{1}) - 1;     % get the end index of the population matrix for i=1
        for i = 1:N-1
            % get the matrix G corresponding to first i windows
            G = G_N(1:G_ix_end, 1:G_ix_end);    % get the matrix for a population of i neurons
            q = q_N(1:G_ix_end);                % get the q_i
            Ginv = pinv(G);                     % compute the pseudoinverse of G
            ck = Ginv*q;                        % calculate cks

            ck_counter = 1;
            for j=1:i
                for m=1:length(tk_in_window{j})-1
                    h_rec_N(i,:) = h_rec_N(i,:) + ck(ck_counter)*W/pi*...
                                   sinc(W*(t_Ph - tk_in_window{j}(m) + windows{j}(1) - tau_1 )/pi);
                    ck_counter = ck_counter + 1;
                end
             end
             G_ix_end = G_ix_end + length(tk_in_window{i+1}) - 1;	% update the end index of the population matrix
        end
        h_rec_N(N,:) = Ph(1,:);              % the last entry is the recovery for all N windows     
        varargout = {h_rec_N};
    end

end % end of function
