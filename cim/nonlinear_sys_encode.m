%NONLINEAR_SYS_ENCODE Encode a temporal signal using a nonlinear dynamic system which has a stable limit cycle
%   [S VARARGOUT] = NONLINEAR_SYS_ENCODE(U,T,SYS_ODE,INIT_STATE) encodes a
%   temporal signal U with a time course T using a nonlinear dynamic system
%   which has a stable limit cycle and an analytic system equation. The state
%   variable of the dynamic system is initialized to INIT_STATE, and the system
%   equation is implemented in the function handle SYS_ODE, which satisfies the
%   following equation,
%
%       dy/dt = SYS_ODE( u(t), y(t) ).
%
%   The output of the neuron is returned as an array of spike time S, as well
%   as the trace of each state variables. The spike is defined as the times
%   that local maximums of the first state variable occurs.
%
%   Authors: Yevgeniy B. Slutskiy and Chung-Heng Yeh
%
%   Copyright 2012-2014 Yevgeniy B. Slutskiy and Chung-Heng Yeh

function [s varargout] = nonlinear_sys_encode( u, t, sys_ode, init_state)
    dt = t(2) - t(1);                         % set the time step
    y = zeros(numel(init_state),numel(t));    % allocate memory for the state variable of the dynamic system
    y(:,1) = init_state;                      % initialize the state variable of the dynamic system
    for i = 2:numel(t)
        y(:,i) = y(:,i-1) + dt*sys_ode(u(i),y(:,i-1));
    end
    % spike detection
    s = [false y(1,1:end-2)<y(1,2:end-1)&y(1,2:end-1)>y(1,3:end) false];
    % pack the trace of each state variable into output argument list
    varargout = cell(1,numel(init_state));
    for i = 1:numel(init_state)
        varargout{i} = y(i,:);
    end
end
