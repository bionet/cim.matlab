%IDENTIFY_H_IDEAL_IAF_TRIG Identify linear dendritic processing filter in cascade with a spiking neuron model.
%   [PH_REC SPK_NUM] = IDENTIFY_H_IDEAL_IAF_TRIG(T,T_PH,U,W,L,B,D,K,S) identifies
%   the projection PH_REC of the linear filter H in cascade with an
%   integrate-and fire (IAF) neuron. PH is the projection of H onto the
%   reproducing kernel Hilbert space of the trigonometric polynomials with
%   bandwidth W and order L, and is of time course T_PH. When L is equal to
%   infinity, PH is the projection of H onto the Paley-Wiener space. The ideal
%   IAF neuron has a bias B, a firing threshold D, and a capacitance K. The
%   signal U with a time course T is encoded into a sequence of spike times S.
%   The identified filter projection PH_REC is returned, as well as the number
%   of total spikes used for identification.
%
%   IDENTIFY_H_IDEAL_IAF_TRIG('RandomThreshold',true) identifies the filter
%   using the algorithm for the random threshold case.
%
%   IDENTIFY_H_IDEAL_IAF_TRIG('SMOOTHFACTOR', NUMBER) uses NUMBER as the smooth
%   factor for the random threshold case. The defacult value of NUMBER is 1e-11.
%
%   IDENTIFY_H_IDEAL_IAF_TRIG('TAU', TAU) takes only the spikes within the
%   window TAU for identification. By default, TAU is set to cover the entire
%   time course T of the signal U.
%
%   Author:               Yevgeniy B. Slutskiy
%
%   Revision Author:      Chung-Heng Yeh
%
%   Copyright 2012-2014   Yevgeniy B. Slutskiy and Chung-Heng Yeh

function [Ph_rec, spk_num] = identify_h_ideal_iaf_trig(...
          t, t_Ph, u, W, L, b, d, k, s, varargin)

    % Handle the optional input parameters
    p = inputParser;
    addParamValue(p,'SmoothFactor',1e-11,@isnumeric);
    addParamValue(p,'RandomThreshold',false,@islogical);
    addParamValue(p,'Tau',[Inf Inf],@(x )isnumeric(x) & numel(x)==2 );
    parse(p,varargin{:});
    % compute the window for faithful identificaiton
    if all( p.Results.Tau == [Inf Inf] )
        tau = [t(1) t(end)];
    else
        tau = p.Results.Tau;
    end

    % Main body of the function starts here
    dt = t(2)-t(1);               % compute the time set

    % function handle for getting spike time of each stimuli
    get_tk = @(i) t( s(i,:)==1 & t >= tau(1) & t <= tau(2) );
    % get spike time
    tk = arrayfun( get_tk, 1:size(s,1),'UniformOutput',false);
    % get the number of total spikes
    spk_num = sum(sum(s));
    % function handle for computing q of each stimuli
    com_q  = @(i) k*d - b*diff(tk{i});
    % compute q of each stimuli and merge them to a single vector
    q = cell2mat(arrayfun( com_q, 1:size(u,1), 'UniformOutput',false))';
    % If L is infinity, it implies that the Paley-Wiener space is considered.
    if isinf(L),
        Phi = zeros(numel(q));                                       % allocate memory for the Phi matrix
        sk = cell2mat(arrayfun( @(i) tk{i}(1:end-1),...              % merge all 'tk-alpha' to a single vector sk
                 1:numel(tk),'UniformOutput', false ) )';
        % Compute the Phi matrix row by row
        row = 1;
        for i = 1:size(u,1)                                          % i: index of the stimuli
            int_u = dt*cumtrapz(u(i,:));                             % compute the integral of the i-th stimulus
            for k = 1:numel(tk{i})-1                                 % k: index of spike of the i-th stimulus
                for m = 1:numel(sk)                                  % m: index of representation function
                    idx = time2idx( [tk{i}(k) tk{i}(k+1)] -sk(m)-t(1), dt, 1, numel(t) );
                    % integral of f from a to b is equal to F(b) - F(a)
                    Phi(row,m) = int_u( idx(2) ) - int_u( idx(1) );
                end
                row = row + 1;
            end
        end
        h_l = pinv(Phi,1e-9)*q;                                      % compute the h_l vector
        Ph_rec = zeros(size(t_Ph));                                  % allocate memory for reconstructed filter
        % Reconstruct the filter
        for m=1:numel(sk)
            Ph_rec = Ph_rec + h_l(m) .* W/pi*sinc(W*(t_Ph - sk(m))/pi);
        end
    else
        T = 2*pi*L/W;                                                % find the period of the RKHS
        e_l = @(m,t) exp(1j*m*W/L*t)/sqrt(T);                        % declare function handle for the basis of the RKHS
        u_l = fft(u,[],2)*sqrt(T)/numel(t);                          % reconstruct Fourier coefficient of the stimuli
        Phi = zeros(spk_num-size(u,1),2*L+1);                        % allocate memory for the Phi matrix
        % Compute the Phi matrix row block by row block. Each block associates
        % with the spikes of a single stimulus.
        row_idx = 0;
        for i = 1:size(u,1)
            % compute the row index of the present block
            row_idx = row_idx(end)+(1:numel(tk{i})-1);
            % for basis elements -L:0 (the rest are conjugate symmetric)
            for m=-L:0
                if m == 0
                    Phi(row_idx,L+1) = u_l(i,1)*diff(tk{i});
                else
                    Phi(row_idx,m+L+1) = u_l(i,end+m+1)*sqrt(T)/(1j*m*W/L)...
                                   *(e_l(m,tk{i}(2:end))-e_l(m,tk{i}(1:end-1)));
                end
            end
            % use the symmetry of the Phi matrix to fill in the remaining entries
            Phi(row_idx,L+2:end) = conj(Phi(row_idx,L:-1:1));
        end
        Phi = conj(Phi);
        % Recover h_l
        if p.Results.RandomThreshold  % if solving the random thresholds problem
            PhiPhi = Phi'*Phi;
            h_l = (PhiPhi + numel(q)*p.Results.SmoothFactor*eye(size(PhiPhi)))\(Phi'*q);
        else
            h_l = pinv(Phi,1e-9)*q;
        end

        % Reconstruct the filter using the inverse Fourier transform
        h_l_pad = [h_l(L+1:end); zeros(numel(t_Ph)-2*L-1,1); h_l(1:L)]';
        Ph_rec = real(ifft(h_l_pad))*numel(t_Ph)/sqrt(T);
        % Shift the reconstructed filter according to the preferred time course
        idx = floor(mod(t_Ph(1)+T/2,T)/dt)+1;
        Ph_rec = Ph_rec([idx:end 1:idx-1]);

        display(['Number of spikes needed: ' num2str(2*L+2)]);
    end
    display(['Number of spikes used: ' num2str(spk_num)]);

    function idx = time2idx( time, dt, lowerBound, upperBound )
    %IDX = time2idx( TIME, DT, LOWERBOUND, UPPERBOUND ) discretizes the
    %   continuous-valued TIME to IDX with step DT. IDX is hard limited by
    %   the upper bound UPPERBOUND and lower bound LOWERBOUND.
        idx = round( time/dt ) + 1;
        idx(idx < lowerBound) = lowerBound;
        idx(idx > upperBound) = upperBound;
    end
end
