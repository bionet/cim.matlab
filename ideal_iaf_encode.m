%IDEAL_IAF_ENCODE Encode a temporal signal using an ideal integrate-and-fire neuron
%   [S V] = IDEAL_IAF_ENCODE(U,T,B,D,K) encodes a temporal signal U with a
%   time course T using an ideal integrate-and-fire neuron with a bias B,
%   a firing threshold D, and a capacitance K. The membrane voltage of the
%   neuron has an initial value equal to zero and is computed using the
%   trapezoidal rule for integration. The output of the neuron is returned
%   as an array of spike time S, as well as the voltage trace V.
%
%   IDEAL_IAF_ENCODE is a simplified version of IAF_ENCODE in the
%   ted.matlab toolbox. Interested users may wish to have a look at the
%   Time Encoding and Decoding Machine toolbox in the Bionet repository.
%
%   Authors: Lev Givon, Yevgeniy B. Slutskiy and Chung-Heng Yeh
%
%   Copyright 2010-2012 Lev Givon, Yevgeniy B. Slutskiy and Chung-Heng Yeh

function [s, v]  = ideal_iaf_encode(u, t, b, d, k)

    dt = t(2)-t(1);                 % get the time-step dt
    s = [];                         % initialize the spikes list
    v = zeros(1, length(t));        % initialize the voltage trace
    e = 0;                          % initialize the integration error term

    for i=2:length(t)

        % compute the membrane voltage
        v(i) = e + v(i-1) + ( b + 0.5*(u(i)+u(i-1)) )*dt/k;     % 6/17/2012

        % if the voltage is above the threshold d, generate a spike
        if v(i) >= d                
            s(end+1) = t(i);        % generate a spike
            v(i) = v(i)-d;          % reset the voltage
            e = v(i);               % keep track of the integration error           
        else
            e = 0;
        end
    end

end % end of function
