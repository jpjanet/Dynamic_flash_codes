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
if any([(tol <= 0),~isnumeric(tol),~(isvector(tol)||isscalar(tol))])
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
if any([(gain <= 0),~isnumeric(gain), ~isscalar(gain)])
    flag = 1;
    if isempty(status)
        status = 'Gain must be positive a double.';
    else
        status = [status, ' Gain must be positive a double.'];
    end
end


end

