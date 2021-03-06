function [ fn, TIME_INFO_15s_1,TIME_INFO_15s_2,TIME_INFO_15s_3,...
    TIME_INFO_23t_1,TIME_INFO_23t_2,TIME_INFO_23t_3] = timed_runs(fn,tspan,jtest,tol,clock,stats,clock_reps,gain,inlcude3)
%Runs the simulation for both the vector of tolerances in tol, for both
%systems, both solvers and all test cases.
% This code repeatedly runs the simulations in s1_explicit_function and
% s2_explicit_function for test cases 1, 2 and 3 (iff include=1), for
% each of the tolerances in tol, and a number of repeats given in
% clock_reps. The average time for solving for each case is returned as
% an output for each case.
% Due to long simulation times, Test Case 3 will not simulate for
% tolerances tighter than 5e-7

integrator = 'ode15s';
TIME_INFO_15s_1= [];
TIME_INFO_15s_2= [];
TIME_INFO_15s_3= [];
for l = 1:(2+inlcude3)
    testcase = l;
    if l ~= 3
        for j = 1:length(tol)
            disp('Calling System 1, ode15s')
            fprintf('\n')
            [tt1,~,~,~,~,~,~,~,~,~,~,~,~,~,~,sol_t1] = s1_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 1 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            disp('Calling System 2, ode15s')
            fprintf('\n')
            [tt2,~,~,~,~,~,~,~,~,~,sol_t2] = s2_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 2 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            if l ==1
                TIME_INFO_15s_1 = [TIME_INFO_15s_1;sol_t1,length(tt1),sol_t2,length(tt2)];
            elseif l==2
                TIME_INFO_15s_2 = [TIME_INFO_15s_2;sol_t1,length(tt1),sol_t2,length(tt2)];
            end
        end
    elseif l==3
        for j = 1:sum(tol>= 5e-7)
            disp('Calling System 1, ode15s')
            fprintf('\n')
            [tt1,~,~,~,~,~,~,~,~,~,~,~,~,~,~,sol_t1] = s1_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 1 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            disp('Calling System 2, ode15s')
            fprintf('\n')
            [tt2,~,~,~,~,~,~,~,~,~,sol_t2] = s2_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 2 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            TIME_INFO_15s_3 = [TIME_INFO_15s_3;sol_t1,length(tt1),sol_t2,length(tt2)];
        end
    end
end

integrator = 'ode23t';
TIME_INFO_23t_1= [];
TIME_INFO_23t_2= [];
TIME_INFO_23t_3= [];
for l=1:(2+inlcude3)
    testcase =l;
    if l ~= 3
        for j = 1:length(tol)
            disp('Calling System 1, ode23t')
            fprintf('\n')
            [tt1,~,~,~,~,~,~,~,~,~,~,~,~,~,~,sol_t1] = s1_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 1 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            disp('Calling System 2, ode23t')
            fprintf('\n')
            [tt2,~,~,~,~,~,~,~,~,~,sol_t2] = s2_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 2 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            if l==1
                TIME_INFO_23t_1 = [TIME_INFO_23t_1;sol_t1,length(tt1),sol_t2,length(tt2)];
            elseif l==2
                TIME_INFO_23t_2 = [TIME_INFO_23t_2;sol_t1,length(tt1),sol_t2,length(tt2)];
            end
        end
    elseif l==3
        for j = 1:sum(tol>= 5e-7)
            disp('Calling System 1, ode23t')
            fprintf('\n')
            [tt1,~,~,~,~,~,~,~,~,~,~,~,~,~,~,sol_t1] = s1_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 1 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            disp('Calling System 2, ode23t')
            fprintf('\n')
            [tt2,~,~,~,~,~,~,~,~,~,sol_t2] = s2_explicit_function(tspan,jtest,stats, integrator,testcase,tol(j),clock,clock_reps,gain);
            disp('~~~~ System 2 complete ~~~~')
            fprintf('\n')
            fprintf('\n')
            TIME_INFO_23t_3 = [TIME_INFO_23t_3;sol_t1,length(tt1),sol_t2,length(tt2)];
        end
    end
end

figure(fn)
clf(fn)
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultLineLineWidth', 3);
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;]);
loglog(tol,TIME_INFO_15s_1(:,1),'--sq',tol,TIME_INFO_15s_1(:,3),'--o',...
    tol,TIME_INFO_23t_1(:,1),':^',tol,TIME_INFO_23t_1(:,3),':x','Markersize',15)
legend('ode15s - System 1 ','ode15s - System 2 ',...
    'ode23t - System 1 ','ode23t - System 2 ')
xlabel('Absolute error tolerance')
ylabel('Solution time, [s]')
title(' Abs. error tolerance and solution time, Test Case 1')
fn=fn+1;

figure(fn)
clf(fn)
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultLineLineWidth', 3);
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;]);
loglog(tol,TIME_INFO_15s_2(:,1),'--d',tol,TIME_INFO_15s_2(:,3),'--+',...
    tol,TIME_INFO_23t_2(:,1),':p',tol,TIME_INFO_23t_2(:,3),':h','Markersize',15)
legend('ode15s - System 1 ','ode15s - System 2 ',...
    'ode23t - System 1 ','ode23t - System 2 ')
xlabel('Absolute error tolerance')
ylabel('Solution time, [s]')
title(' Abs. error tolerance and solution time, Test Case 2')
fn=fn+1;

figure(fn)
clf(fn)
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultLineLineWidth', 3);
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;]);
subplot(121)
loglog(tol,TIME_INFO_15s_1(:,4),'--sq',tol,TIME_INFO_15s_1(:,2),'--o',...
    tol,TIME_INFO_23t_1(:,4),':^',tol,TIME_INFO_23t_1(:,2),':x','Markersize',15)
legend('ode15s - System 1 ','ode15s - System 2 ',...
    'ode23t - System 1 ','ode23t - System 2 ')
xlabel('Absolute error tolerance')
ylabel('Number of timesteps')
title(' Test Case 1')
subplot(122)
loglog(tol,TIME_INFO_15s_2(:,2),'--sq',tol,TIME_INFO_15s_2(:,4),'--o',...
    tol,TIME_INFO_23t_2(:,2),':^',tol,TIME_INFO_23t_2(:,4),':x','Markersize',15)
legend('ode15s - System 1 ','ode15s - System 2 ',...
    'ode23t - System 1 ','ode23t - System 2 ')
xlabel('Absolute error tolerance')
ylabel('Number of timesteps')
title(' Test Case 2')
xlabel('Absolute error tolerance')
fn=fn+1;


if inlcude3
    figure(fn)
    clf(fn)
    set(0, 'DefaultAxesFontSize', 24);
    set(0, 'DefaultLineLineWidth', 3);
    set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
        0.8500,0.3250,0.0980;]);
    
    loglog(tol(1:sum(tol>= 5e-7)),TIME_INFO_15s_3(:,2),'--sq',tol(1:sum(tol>= 5e-7)),TIME_INFO_15s_3(:,4),'--o',...
        tol(1:sum(tol>= 5e-7)),TIME_INFO_23t_3(:,2),':^',tol(1:sum(tol>= 5e-7)),TIME_INFO_23t_3(2,4),':x',...
        tol(2:end),TIME_INFO_15s_2(2:end,1),'k--d','Markersize',15)
    legend('test case 3 - ode15s - System 1 ','test case 3 - ode15s - System 2 ',...
        'test case 3 - ode23t - System 1 ','test case 3 - ode23t - System 2 ', 'test case 2 - ode15s - System 1 ')
    xlabel('Absolute error tolerance')
    ylabel('Solution time, [s]')
    title(' Solution time and error tolerance for test Case 3')
    fn=fn+1;
end


end

