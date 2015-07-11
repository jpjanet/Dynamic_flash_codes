
function [tt,Mi,x,y,Ml,Mv,M,L,V,T,P,Pout,F,w,rho_R,sol_t] = s1_explicit_function(tspan,jtest,Stats, integrator, testcase, tol,clock,clock_reps,gain)
%This function runs the simulation specified in the inputs using the 
% system 1 equations.
%   The input is defined in the dynamic_flash_mainfile.m script
%   The outputs are the various state varaibles as vectors at the times
%   specified in the vector tt.

global Nc C_l C_v
Nc = 4;   % Number of speices

%% Print the problem details to terminal
disp(['System 1. Solving with tol = ',num2str(tol),'. Test case = ',num2str(testcase),'. Integrator = ',integrator])
if clock
    disp(['Timer is on with reps = ',num2str(clock_reps)])
end



%% Species Data
% The Antoine eq coefficients are given in a matrix, one row per
% component
Ant_coeff = [3.928, 803.99, 247.04;  ...
    3.932, 935.77, 238.00;  ...
    3.978, 1064.84, 232.00; ...
    4.000, 1170.88, 224.32];
Molar_weights = [44.10; 58.12; 72.15; 86.20];

% Pure component critical properties stored as rows,
% with Tc [K], Pc [bar], Vc [cm^3/mol], Zc []
Critical_props = [369.80, 42.48, 200.00, 0.276; ...
    425.12, 37.96, 255.00, 0.274; ...
    469.80, 33.70, 311.00, 0.268; ...
    507.60, 30.25, 368.00, 0.264];



%% Initial conditions
% First, consisent initial conditions are determined using the
% Rachford-Rice routine, as described in Section 4. The base-case
% conditions are used:
w0 = [0.1;0.2;0.3;0.4]; % Inital molar fraction
F0 = 100; % Feedrate in [kmol/h]
Tref = 366.5; % Reference T in K
P0 = 6.89; % Initial pressure in  Bar(a)
Pout0 =  1; % Initial external pressue in Bar (a)
PSI0 = 0.1219; % Initial guess from source
t = 0; % Initial time
CT = 1.5; % volume, m3
phi0 = [F0*w0;CT;Pout0];

% Call the Rachford-Rice functions
if clock
    tic
    for i = 1: clock_reps % Repeat the timer if desired
        [z0, dotz0,CT, C_l, C_v, rho_R, ierr ] = s1_initial_condition_function(Tref,w0,F0,P0,Pout0,Critical_props,PSI0);
    end
    init_t = toc;
    init_t = init_t/clock_reps; % Repeat the timer if desired
    disp(['Initial conditions found in ',num2str(init_t),' seconds'])
else
    [z0, dotz0,CT, C_l, C_v, rho_R, ierr ] = s1_initial_condition_function(Tref,w0,F0,P0,Pout0,Critical_props,PSI0);
end


if (ierr == 0) % check inital condition set-up
    disp('Successful steady-state conditions calculated by Rachford-Rice routine')
    normEr = norm(SYSTEM1(t,z0,dotz0,phi0, rho_R,C_l,C_v));
    disp(['Error at t = 0 is ',num2str(normEr)]);
else
    disp('Unsuccesful steady-state condition calculation')
    normEr = norm(SYSTEM1(t,z0,dotz0,phi0, rho_R,C_l,C_v));
    disp(['Error at t = 0 is ',num2str(normEr)]);
    disp('Please ensure operation in two-phase region at t = 0')
end



%% External function wrappers

% reform as function of X
x0 = z0([(1:3*Nc + 5),(3*Nc + 7)]);


% Apply controls and Inhomogenieties
if testcase ==1
    Tcon = @(t,X) (TcontrolX_varying(t,X,Tref));
    phifun = @(t) (s1_phi(t,F0,w0,CT,Pout0,@FeedFucntion));
elseif testcase == 2
    Tcon = @(t,X) (TcontrolX(t,X,Tref));
    phifun = @(t) (s1_phi(t,F0,w0,CT,Pout0,@FeedFucntion_varying));
elseif testcase == 3
    LTarget = x0(3*Nc + 4);
    Tcon = @(t,X) Tcontrol_v(t,X,Tref,LTarget,gain,1);
    phifun = @(t) (s1_phi(t,F0,w0,CT,Pout0,@FeedFucntion_varying));
else
    disp('Error: testcase undefined')
    return
end


% system equations
System1_Explicit_RHS = @(t,x) (SYSTEM1_E_RHS(t,x,Tcon,phifun,rho_R,C_l,C_v));

% Jacobian
System1_Explicit_J = @(t,x) (SYSTEM1_E_J(t,x,Tcon,phifun,rho_R));


% Test explicit system Jacobian:
if (jtest)
    F_handle = @(x)(System1_Explicit_RHS(0,x));
    J = System1_Explicit_J(0,x0);
    delta = 0.00001;
    nJ = numjac(F_handle,x0,delta);
    G = J - nJ;
    disp(['Difference between numerical (central dif) Jacobian and Analytic =  ',num2str(norm(G))])
    disp(['With delta   = ',num2str(delta)])
    figure(2)
    clf()
    imagesc(abs(nJ-J))
    colorbar
end

% Mass matrix
System1_Explicit_M =  @(t,x)(SYSTEM1_E_M(t,x,Tcon));

% Mass matrrix is constant
M = System1_Explicit_M(0,x0);



% Set up ODE options
options = odeset('Mass',M,'Jacobian',System1_Explicit_J,'Stats',Stats,'AbsTol',tol,'RelTol',1E-10);

if ~clock
    clock_reps = 1;
    sol_t = 'clock_off';
end
% Solve DAE
disp(['Hand-off to solver sucessful. Repetitions = ', num2str(clock_reps),'. Sysem size  = ',num2str(length(x0)), '. Solver output = '])
fprintf('\n')

if strcmp(integrator,'ode15s');
    if clock
        tic
        for i = 1:clock_reps
            [tt,X] = ode15s(System1_Explicit_RHS,tspan,x0,options);
        end
        sol_t = toc;
        sol_t=sol_t/clock_reps;
    else
        [tt,X] = ode15s(System1_Explicit_RHS,tspan,x0,options);
    end
elseif strcmp(integrator,'ode23t');
    if clock
        tic
        for i = 1:clock_reps
            [tt,X] = ode23t(System1_Explicit_RHS,tspan,x0,options);
        end
        sol_t = toc;
        sol_t=sol_t/clock_reps;
    else
        [tt,X] = ode23t(System1_Explicit_RHS,tspan,x0,options);
    end
else
    disp('Invalid integrator')
    return
end

if clock
    disp(['Solved in (average of) ',num2str(sol_t),' seconds'])
else
    disp(['Solver complete'])
end
fprintf('\n')

% Post-Process
[Mi,x,y,Ml,Mv,M,L,V,T,P,Pout,F,w] = s1_X_to_variable(tt,X,Tcon,phifun);

end
