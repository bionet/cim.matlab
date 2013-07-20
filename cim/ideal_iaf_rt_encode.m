%IDEAL_IAF_RT_ENCODE Encode a temporal signal using an ideal integrate-and-fire neuron with random thresholds
%   [S V DELTA_V] = IDEAL_IAF_RT_ENCODE(U,T,B,D,K) encodes a temporal signal U
%   with a time course T using an ideal integrate-and-fire neuron with a bias B,
%   a capacitance K, and a random threshold varying from spike to spike. The
%   membrane voltage of the neuron has a random initial value, and  is computed
%   using the trapezoidal rule for integration. The output of encoder contains
%   an array of spike time S, the voltage trace V, as well as an array of
%   random thresholds associated with each spikes.
%
%   IDEAL_IAF_RT_ENCODE(...,'DeltaVar',DELTAVAR) resets the threshold to an
%   sample of a Gaussian distribution of mean D and variance DETAVAR after each
%   spike generation. The default vale of DELTAVAR is zero.
%
%   IDEAL_IAF_RT_ENCODE(...,'RandomInit',false) encodes the signal with the
%   membrane voltage initialized to zero.
%
%   IDEAL_IAF_RT_ENCODE is a simplified version of IAF_ENCODE in the
%   ted.matlab toolbox. Interested users may wish to have a look at the
%   Time Encoding and Decoding Machine toolbox in the Bionet repository.
%
%   Authors: Lev Givon, Yevgeniy B. Slutskiy and Chung-Heng Yeh
%
%   Copyright 2012-2014 Lev Givon, Yevgeniy B. Slutskiy and Chung-Heng Yeh

function [s, v, delta_v] = ideal_iaf_rt_encode(u, t, b, d, k, varargin)
    % Handle the optional input arguments
    p = inputParser;
    addParamValue(p,'DeltaVar',0,@isnumeric);
    addParamValue(p,'RandomInit',true,@islogical);
    parse(p,varargin{:});

    % Main body of the function starts here
    sigma = p.Results.DeltaVar;

    dt = t(2)-t(1);                                                  % get the time step
    s  = zeros(size(t));                                             % allocate a vector for the spike train of IAF neuron
    v  = zeros(size(t));                                             % allocate a vector for the membrane voltage of IAF neuron
    delta_v = d + sigma*randn();                                     % initialize the random threshold
    if p.Results.RandomInit,
        v(1) = 0.9*delta_v*rand();                                   % initialize the membrane voltage with random value
    end

    % Perform the numerical integration using the trapzoidal method
    for i = 2:numel(t)
       v(i) = v(i-1) + dt/k*(b+0.5*(u(i)+u(i-1)));                   % update the membrane voltage
        if v(i) >= delta_v(end)                                      % spike detection
            v(i) = v(i)-delta_v(end);                                % reset the membrane voltage
            s(i)  = 1;                                               % record spike
            delta_v = [delta_v d+sigma*randn()];                     % get a new random threshold
        end
    end
end
