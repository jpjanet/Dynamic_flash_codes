%% Dynamic Flash Equations 
% Jon Paul Janet
% 07/07/2015

clc, clear all; close all;
%% Set up and solve dynamic flash equations 
% This is the top level script that calls various functions need to
% simulate dynamic flash vessels using the two models in [1]. Each call
% will solve both systems (the full and simplified models) and return the
% results from both for comparison

%% User-specified parameters
% This section is used to contol which  solver  and test case is being
% used. This is the only section that is intended to be edited.

% 1) test case options
testcase = 2; % Swtich between varying temp (=1) and varying feed (=2)
              % and varying feed with control (=3) 
              

gain = 3; % Set proportional controller gain, double. Should be positive,
          % setting a higher gain inccurs serious computational time
          % penalty at high tolerances. Reccomnended range = [0.5,3]. This
          % parameter has no effect if testcase != 3.
              
tspan = [0,48]';% Integration time span, either a vector of two doubles
                % [begin, end]', or a vector of desired sample times

% 2) Solver and debug options
tol = [5E-5];
%tol = [5E-2,5E-7,5E-7]; % Provide abs error tolerance for integration
                        % this can either be a double or a vector of
                        % doubles. If a vector is provided, the option for
                        % timing (" clock, below ") must be switched on.                    
                        
clock = 0;  % Toggle if the runs should be timed, binary.  Options are =1 
            % for yes, and =0 for no. If yes (=1) is selected, 
            % you may optionally  use a vector of tolerances.
            % This provides allows for time-error comparisons.
            
clock_reps = 10; % Number of times to repeat the timer (integer). This  
                  % has no effect if the clock is set to zero

jtest  = 0; % Switch to activate numerical checking of provided 
            % Jacobian functions, binary. On = 1, off = 0.
            
integrator = 'ode15s'; % Switch to select integrator, either 'ode15s' 
                       % or 'ode23t', string.
stats = 'off'; % Switch to provide solver output details 
               % (number of evals etc), either 'on' or 'off', string.

% 3) plotting and postprocessing control

plotting = 1; % Post-process results and generate plots, binary.  =1 yes 
              % for yes, =0 no.

fn = 1; % Figure number, integer. This will start drawing at figure = fn.

print_figs = 0; % Option to automatically save plots, binary. =1 for yes,
                % =0 for no. Will write to current directory.

imageF = '-depsc'; % Image format, string. Only used if print_figs = 1.
                   % See the Matlab function help for options. -depsc
                   % produces vector eps output (reccomended)

comp_string = ['prop  ';'n-but ';'n-pent';'n-hex ']; % Names of the 
                                                     % components to 
                                                     % display.  


                                                     
                                                     
%% Input Checking

% Check the input parameters
[flag, status ] = input_checking (testcase,gain,tspan,tol,clock,...
                                  clock_reps,jtest,integrator,stats,...
                                  plotting, fn, print_figs);

% Print details to terminal 
disp(['Intialisation, solver = ',integrator,', verbose = ',stats])
fprintf('\n')
disp(['Test case = ',num2str(testcase)])
fprintf('\n')
if ~clock
    disp('No clock -  run will not be timed')
    fprintf('\n')
end
disp('Calling System 1')
fprintf('\n')

%% Call soliving routines
if clock
     [ fn, TIME_INFO_15s_1,TIME_INFO_15s_2,TIME_INFO_15s_3,...
           TIME_INFO_23t_1,TIME_INFO_23t_2,TIME_INFO_23t_3]=timed_runs(fn,tspan,jtest,tol,clock,stats,clock_reps,gain);
    
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
if plotting
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





% 
%  [tt2,Mi2,x2,y2,Ml2,Mv2,M2,L2,V2,T2,P2,Pout2,F2,w2,rho_R,sol_t2] = s1_explicit_function(tspan,jtest,stats, integrator,2,tol,clock,clock_reps,gain);
%  
%  [tt3,Mi3,x3,y3,Ml3,Mv3,M3,L3,V3,T3,P3,Pout3,F3,w3,rho_R,sol_t3] = s1_explicit_function(tspan,jtest,stats, integrator,3,tol,clock,clock_reps,gain);
%  
%  
%  
%  set(0, 'DefaultAxesFontSize', 24);
% set(0, 'DefaultLineLineWidth', 3);
% global Nc
% 
% FR3 = V3./F3; % Flash Ratio 1
% FR2 = V2./F2; % Flash Ratio 2
% SR3 = (y3.*(V3*ones(1,Nc)))./(x3.*(L3*ones(1,Nc))); % Split ratios by species
% SR2= (y2.*(V2*ones(1,Nc)))./(x2.*(L2*ones(1,Nc))); %  Split Ratios by species
% SP3 = SR3(:,2)./SR3(:,3); % Seperation power butane/pentane 1
% SP2 = SR2(:,2)./SR2(:,3); % Seperation power butane/pentane 2
% figure(fn)
% clf(fn)
% set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
%     0.8500,0.3250,0.0980]);
% plot(tt2,L2,'--^',tt3,L3,'-o',[0,48],[L2(1),L2(1)],':k')
% axis tight
% title(['Liquid outflow, control compared to free system'])
% xlabel('Time, t [h]')
% ylabel('Liquid phase outflow, [kmol h^{-1} ]')
% legend('Control off (test case 2)','Control on (test case 3)','Target')
% grid on
% fn=fn+1;
% 
% figure(fn)
% clf(fn)
% set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
%     0.8500,0.3250,0.0980]);
% plot(tt2,T2,'--^',tt3,T3,'-o')
% axis tight
% title(['Temperature, control compared to free system'])
% xlabel('Time, t [h]')
% ylabel('Temperature, [k]')
% legend('Control off (test case 2)','Control on (test case 3)')
% grid on
% fn=fn+1;
% 
% 
% 
% figure(fn)
% clf(fn)
% set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
%     0.8500,0.3250,0.0980]);
% plot(tt2,FR2,'--^',tt3,FR3,'-o')
% axis tight
% title('Flash ratio, control compared to free system')
% xlabel('Time, t [h]')
% ylabel('Flash ratio, [ ]')
% grid on
% legend('Control off (test case 2)','Control on (test case 3)')
% fn=fn+1;
% figure(fn)
% clf(fn)
% set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
%     0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560]);
% plot(tt2,SP2,'--^',tt3,SP3,'-o')
% axis tight
% title('Separation power, control compared to free system')
% xlabel('Time, t [h]')
% ylabel('Separation power n-butane to n-pentane, [ ]')
% grid on
% legend('Control off (test case 2)','Control on (test case 3)')
% fn=fn+1;

