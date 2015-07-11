clc, clear all;
% Species are 1)propane 2)butane 3)pentate 4)hexane

global Nc C_l C_v
Nc = 4;   % Number of speices
jtest = 1;
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

              

%% Problem Input
w0 = [0.1;0.2;0.3;0.4]; % inital molar fraction
F0 = 100; % feedrate in [kmol/h]
Tref = 366.5; % reference T in K
P0 = 6.89; % Bar(a)
Pout0 =  1; % Bar (a)
PSI0 = 0.1219; % Initial guess from source
t = 0; % Initial time
CT = 1.5; % volume, m3
phi0 = [F0*w0;CT;Pout0];


[z0, dotz0,CT, C_l, C_v, rho_R, ierr ] = s1_initial_condition_function(Nc,Tref,w0,F0,P0,Pout0,Ant_coeff,Molar_weights, Critical_props,PSI0);





if (ierr == 0) % check inital condition set-up
    disp('Successful steady-state conditions calculated by Rachford-Rice Routine')
    normEr = norm(SYSTEM1(t,z0,dotz0,phi0, rho_R,C_l,C_v));
    disp(['Error at t = 0 is ',num2str(normEr)]);
else
    disp('Unsuccesful steady-state condition calculation')
    normEr = norm(SYSTEM1(t,z0,dotz0,phi0, rho_R,C_l,C_v));
    disp(['Error at t = 0 is ',num2str(normEr)]);
    disp('Please ensure operation in two-phase region at t = 0')
end







%% External function wrappers
% Inhomogenieties
phifun = @(t) (s1_phi(t,F0,w0,CT,Pout0,@FeedFucntion));
 
% Apply controls
Tcon = @(t,X) (TcontrolX(t,X,Tref));


% reform as function of X
x0 = z0([(1:3*Nc + 5),(3*Nc + 7)]);
dotx0=dotz0([(1:3*Nc + 5),(3*Nc + 7)]);





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
[maxval, maxloc] = max(abs(G(:)));
[maxloc_row, maxloc_col] = ind2sub(size(G), maxloc);
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


% Set up ode15s options
options = odeset('Mass',M,'Jacobian',System1_Explicit_J,'AbsTol',1E-10,'MaxStep',0.5,'MassSingular','yes','BDF','on');
tspan = [0,20];

% Solve DAE
[tt,X] = ode15s(System1_Explicit_RHS,tspan,x0,options);
%  if jtest
%      for i = 1:length(tt)
%          disp(norm(System1_Explicit_RHS(tt(i),X(i,:)')))
%      end
%  end



% Post-Process
[Mi,x,y,Ml,Mv,M,L,V,T,P,Pout,F,w] = s1_X_to_variable(tt,X,Tcon,phifun);
K = VLEFUN(T(end),P(end));




return






