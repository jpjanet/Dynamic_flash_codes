%% Dynamic Flash Equations 
% Jon Paul Janet
% 07/07/2015

clc, clear all; close all;
%% Set up and solve dynamic flash equations 
% This is the top level script that calls various functions needed to
% simulate dynamic flash vessels using the two models in [1]. Each call
% will solve both systems (the full and simplified models) and return the
% results from both for comparison.
% See the readme.txt file, or Appendix C in the accompanying report.  

%% User-specified parameters
% This section is used to control which  solver and test case are being
% used. This is the only section that is intended to be edited.

% 1) test case and general setup options
testcase = 1; % Switch between varying temp (=1) and varying feed (=2)
              % and varying feed with control (=3) 
integrator = 'ode15s'; % Switch to select integrator, either 'ode15s' 
                       % or 'ode23t', string.
stats = 'off'; % Switch to provide solver output details 
               % (number of evals etc), either 'on' or 'off', string.              

gain = 3; % Set proportional controller gain, double. Should be positive,
          % setting a higher gain incurs serious computational time
          % penalty at high tolerances. Recommended range = [0.5,3]. This
          % parameter has no effect if testcase != 3.
              
tspan = [0,48]';% Integration time span, either a vector of two doubles
                % [begin, end]', or a vector of desired sample times.

% 2) Solver and debug options
tol = [5E-5]; % Provide abs error tolerance for integration, this can 
              % be a double or a vector of doubles. If a single tolerance 
              % is given, both systems are simulated for the given test
              % case. In the case that a vector is provided, the option for
              % timing (clock, below) must be switched on. This will 
              % cause the code to repeat for each tolerance, time the run,
              % and provide a comparison across ALL test cases and solvers.
% tol = [5E-2,5E-3,5E-4,5E-5,5E-6,5E-7,5E-8,5E-9]';             
                        
clock = 1;  % Toggle if the runs should be timed, binary.  Options are =1 
            % for yes, and =0 for no. If yes (=1) is selected, 
            % you may optionally  use a vector of tolerances. This enables
            % automatic for time-error comparisons. 
            
clock_reps = 5; % Number of times to repeat the timer (natural #). This  
                 % has no effect if the clock is set to zero.
                  
include3 = 0; % Switch to include the third test case in the timed runs. 
              % The options are =1 for yes and =0 for no. 
              % This only applies in the case that tol is a vector and
              % clock = 1, such that the time-error comparison code is
              % being run. Test case 3 takes significantly longer to solve
              % compared to cases 1 and 2, and hence leaving this disabled
              % might decrease computational time significantly.

jtest  = 0; % Switch to activate numerical checking of provided 
            % Jacobian functions, binary. On = 1, off = 0.
            


% 3) plotting and postprocessing control

plotting = 1; % Post-process results and generate plots, binary.  =1  
              % for yes, =0 for no.

fn = 1; % Figure number, integer. This will start drawing at figure = fn. 
        % No effect if plotting = 0. 

print_figs = 0; % Option to automatically save plots, binary. =1 for yes,
                % =0 for no. Will write to current directory.

imageF = '-depsc'; % Image format, string. Only used if print_figs = 1.
                   % See the Matlab function help for options. -depsc
                   % produces vector eps output (recommended)

comp_string = ['prop  ';'n-but ';'n-pent';'n-hex ']; % Names of the 
                                                     % components to 
                                                     % display.  


                                                     
                                                     
%% Input Checking

% Check the input parameters
[flag, status ] = input_checking (testcase,gain,tspan,tol,clock,...
                                  clock_reps,jtest,integrator,stats,...
                                  plotting, fn, print_figs);
                              
if flag ~= 0
    disp(['Parameter setting failue .',status])
    return
else
    disp('Input passed verification.')
end

% Print details to terminal 
disp(['Intialisation, solver = ',integrator,', verbose = ',stats])
fprintf('\n')
if ~clock
    disp('No clock -  run will not be timed')
    fprintf('\n')
end
disp('Calling System 1')
fprintf('\n')

%% Call soliving routines
if (clock) && (~isscalar(tol))
     [ fn, TIME_INFO_15s_1,TIME_INFO_15s_2,TIME_INFO_15s_3,...
           TIME_INFO_23t_1,TIME_INFO_23t_2,TIME_INFO_23t_3]=timed_runs(fn,tspan,jtest,tol,clock,stats,clock_reps,gain,include3);
    
else
    [tt1,Mi,x1,y1,Ml,Mv,M1,L1,V1,T1,P1,Pout1,F1,w1,rho_R,sol_t1] = s1_explicit_function(tspan,jtest,stats, integrator,testcase,tol(1),clock,clock_reps,gain);
    disp('~~~~ System 1 complete ~~~~')
    fprintf('\n')
    fprintf('\n')
    disp('Calling System 2')
    fprintf('\n')
    [tt2,x2,y2,M2,L2,V2,T2,P2,F2,w2,sol_t2] = s2_explicit_function(tspan,jtest,stats, integrator,testcase,tol(1),clock,clock_reps,gain);
    disp('~~~~ System 2 complete ~~~~')
    fprintf('\n')
    fprintf('\n')
    
end

% Post-Process
if plotting && (isscalar(tol))
    fprintf('\n')
    disp('Preparing postprocessing')
    [fn] = flash_postprocessing( tt1,Mi,x1,y1,Ml,Mv,M1,L1,V1,T1,P1,Pout1,F1,w1,rho_R,tt2,x2,y2,M2,L2,V2,T2,P2,F2,w2,testcase,comp_string ,integrator,imageF, print_figs ,fn);
else
    disp('Postprocessing off')
end








%% References
% [1]: Biegler, L. T. (2010). 
% "Nonlinear programming: concepts, algorithms, and applications to
%  chemical processes, volume 10. SIAM, Philadelphia, PA.




