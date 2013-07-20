%ASDM_ENCODE Encode a temporal signal using an asynchronous sigma/delta modulator
%   [Z, V] = IDEAL_IAF_RT_ENCODE(U,T,B,D,K) encodes a temporal signal U
%   with a time course T using an asynchronous sigma/delta modulator (ASDM)
%   with bias B, Schmitt trigger threshold D, and integration capacitance K.
%   The integrator voltage of the ASDM is initialized to -D, and is computed
%   using the rectangular rule for integration. The output of this function
%   contains the ASDM output Z, as well as the integrator voltage trace V.
%
%   ASDM_ENCODE is a simplified version of ASDM_ENCODE in the ted.matlab
%   toolbox. Interested users may wish to have a look at the Time Encoding and
%   Decoding Machine toolbox in the Bionet repository.
%
%   Authors: Lev Givon, Yevgeniy B. Slutskiy and Chung-Heng Yeh
%
%   Copyright 2012-2014 Lev Givon, Yevgeniy B. Slutskiy and Chung-Heng Yeh

function [z, v]  = asdm_encode(u, t, b, d, k)

    dt = t(2)-t(1);                                       % get the time step
    v = zeros(1, length(t));                              % allocate memory for the integrator output
    z = zeros(1, length(t));                              % allocate memory for the ASDM output
    v(1) = -d;                                            % initialize the integrator output
    z(1) = -b;                                            % initialize the ASDM output

    % Perform the numerical integration using the rectangular method
    for i=2:length(t)
        v(i) = v(i-1) + dt/k*(u(i) - z(i-1));             % update the integrator output
        if abs(v(i)) >= d && sign(v(i)) ~= sign(z(i-1))   % if the integrator output reaches the threshold
            z(i)  = -z(i-1);                              % switch the phase of ASDM
        else
            z(i) = z(i-1);                                % remain at the same phase
        end
    end
end
