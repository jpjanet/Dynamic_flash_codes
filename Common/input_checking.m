function [flag, status ] = input_checking (testcase,gain,tspan,tol,clock,...
    clock_reps,jtest,integrator,stats,...
    plotting, fn, print_figs)
% Test if input is correct
%   This script returns 0 if the inputs are all of the required form and
%   have intrepratable types, and 1 otherwise. A string is return as status
%   to help indentify any problems.

flag = 0;
status = [];

% Check testcase
if ~ismember(testcase,[1,2,3])
    flag = 1;
    status = 'Testcase must be 1, 2 or 3.';
end

% Check gain
if any([(gain <= 0),~isnumeric(gain), ~isscalar(gain)])
    flag = 1;
    if isempty(status)
        status = 'Gain must be positive a double.';
    else
        status = [status, ' Gain must be positive a double.'];
    end
end

% Check tspan
if any([all(tspan < 0),~isnumeric(tspan),~(all(diff(tspan)>0)),~isvector(tspan)])
    flag = 1;
    if isempty(status)
        status = 'Tspan must be positive, monotonic vector.';
    else
        status = [status,' Tspan must be positive, monotonic vector.'];
    end
end

% Check tol
if any([max(tol <= 0),~isnumeric(tol),~(isvector(tol)||isscalar(tol))])
    flag = 1;
    if isempty(status)
        status = 'Tolerance must be positive scalar or vector.';
    else
        status = [status, ' Tolerance must be positive scalar or vector.'];
    end
end


% Check clock
if ~ismember(clock,[0,1])
    flag = 1;
    if isempty(status)
        status = 'Clock must be either 0 or 1.';
    else
        status = [status, 'Clock must be either 0 or 1.'];
    end
end

if (~clock && ~isscalar(tol))
    flag = 1;
    if isempty(status)
        status = 'Vector tolerance only supported when clock is on.';
    else
        status = [status, ' Vector tolerance only supported when clock is on.'];
    end    
    
end


% Check Clock_reps
if any([(clock_reps <= 0),~isnumeric(clock_reps), ~isscalar(clock_reps),mod(clock_reps,1)])
    flag = 1;
    if isempty(status)
        status = 'clock_reps must be positive a integer.';
    else
        status = [status, ' clock_reps must be positive a integer.'];
    end
end

% Check jtest
if ~ismember(jtest,[0,1])
    flag = 1;
    if isempty(status)
        status = 'jtest must be either 0 or 1.';
    else
        status = [status, ' jtest must be either 0 or 1.'];
    end
end

% Check integrator
if  ((~strcmp(integrator,'ode23t')) && (~strcmp(integrator,'ode15s'))) 
    flag = 1;
    if isempty(status)
        status = 'integrator must be either ode23t or ode15s.';
    else
        status = [status, ' integrator must be either ode23t or ode15s.'];
    end
end

% Check stats
if  ((~strcmp(stats,'on')) && (~strcmp(stats,'off'))) 
    flag = 1;
    if isempty(status)
        status = 'Stats must be either on or off.';
    else
        status = [status, ' Stats must be either on or off.'];
    end
end

% Check plotting
if ~ismember(plotting,[0,1])
    flag = 1;
    if isempty(status)
        status = 'plotting must be either 0 or 1.';
    else
        status = [status, ' plotting must be either 0 or 1.'];
    end
end

% Check fn
if any([(fn <= 0),~isnumeric(fn), ~isscalar(fn),mod(fn,1)])
    flag = 1;
    if isempty(status)
        status = 'fn must be either a positive integer.';
    else
        status = [status,' fn must be either a positive integer.'];
    end
end


% Check print_figs
if ~ismember(print_figs,[0,1])
    flag = 1;
    if isempty(status)
        status = 'print_figs must be either 0 or 1.';
    else
        status = [status, 'print_figs must be either 0 or 1.'];
    end
end



end

