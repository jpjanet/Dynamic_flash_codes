function  [fn] = flash_postprocessing( tt1,Mi,x1,y1,Ml,Mv,M1,L1,V1,T1,P1,...
                                       Pout1,F1,w1,rho_R,tt2,x2,y2,M2,L2,...
                                       V2,T2,P2,F2,w2,testcase,...
                                       comp_string,integrator,imageF,...
                                       print_figs,fn)
%Used to draw graphs based on the simulation results from
%dynamic_flash_mainfile.m
% The inputs are the simulation timessteps tt1 and t22, for System 1 and
% Sytem 2 respectively, and the physical variables for each system. The
% function also takes arguments that dictate if the figures are saved or
% only displayed, and what format to save the figures as (print_figs,
% imageF respectively).

% Set plot defaults
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultLineLineWidth', 3);
global Nc

FR1 = V1./F1; % Flash Ratio 1
FR2 = V2./F2; % Flash Ratio 2
SR1 = (y1.*(V1*ones(1,Nc)))./(x1.*(L1*ones(1,Nc))); % Split ratios by species
SR2= (y2.*(V2*ones(1,Nc)))./(x2.*(L2*ones(1,Nc))); %  Split Ratios by species
SP1 = SR1(:,2)./SR1(:,3); % Seperation power butane/pentane 1
SP2 = SR2(:,2)./SR2(:,3); % Seperation power butane/pentane 2
VOL_v = Mv./GASDENSE(y1, T1, P1); % System 1, vapour inventory volume
VOL_l = zeros(length(tt1),1);
for j=1:length(tt1)
    VOL_l(j) =  Ml(j)/LIQDENSE(x1(j,:)', T1(j), P1(j),rho_R);
end

%% Save figures:
if print_figs
disp(['Saving figures as ',imageF, ' in CD'])
else
disp('Not saving figures')
end

%% Plot temperature or inflow profile
set(0, 'DefaultLineLineWidth', 1.5);
if testcase==1
figure(fn)
clf(fn)
set(gcf,'color','white');
z = zeros(size(tt1'));
col = T1';  % Vary color with temp
surface([tt1';tt1'],[T1';T1'],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',3);
colormap jet
title('Temperature input function')
xlabel('Time, t [h]')
ylabel('Temperature, T [K]')
axis([0 48 310 420 ] )
line([0,48],[T1(1), T1(1)],'color','k','LineStyle','--','LineWidth',1.5)
line([15,15],[T1(1),T1(1)+50],'color','r','LineStyle','-','LineWidth',1)
text(15 + 1 ,T1(1) + 25,'T^{Ref} + 50 K','color','r','FontSize',24)
text(45 - 5.5 ,T1(1) - 25,'T^{Ref} - 50 K','color','b','FontSize',24)
text(24 -2.5,T1(1) + 5,['T^{Ref} = ',num2str(T1(1)),' K ',],'color','k','FontSize',24)
line([45,45],[T1(1),T1(1)-50],'color','b','LineStyle','-','LineWidth',1)
fn=fn+1;
elseif testcase ==2
figure(fn)
clf(fn)
set(gcf,'color','white');
plot(tt1,w1)
grid on
%axis([0 48 0 0.5])
axis tight
title('Composition input function')
xlabel('Time, t [h]')
ylabel('Inflow molar fraction, w [ ]')
legend(comp_string)
fn =fn+1;
end
    
%% Comparison Plots
% Number of timesteps plot
figure(fn)
clf(fn)
set(gcf,'color','white');
plot(tt1,[1:1:length(tt1)],'r--^',tt2,[1:1:length(tt2)],'b--sq')
title(['Number of timesteps, Test Case ',num2str(testcase),', integrator = ',integrator])
xlabel('Time, t [h]')
ylabel('Cumulative number of timesteps')
axis([0 48 0 max([length(tt2);length(tt1)])] )
legend('System 1','System 2')
if print_figs 
print(['Timesteps_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF,'-painters')
end
fn=fn+1;

% Plot composition-time space for each species
figure(fn)
clf(fn)
subplot(121)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560]);
plot(tt1,x1,'--^',tt2,x2,':sq')
title(['Liquid phase, Test Case ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Molar fraction, [ ]')
grid on
%legend([repmat('S1 ', 4,1),comp_string,;repmat('S2', 4,1),comp_string,],'Location','northoutside','Orientation','horizontal')
%legend([repmat('S1 ', 4,1),comp_string,;repmat('S2 ', 4,1),comp_string])
subplot(122)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560]);
plot(tt1,y1,'--^',tt2,y2,':sq')
title(['Gas phase, Test Case ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Molar fraction, [ ]')
grid on
legend([repmat('S1 ', 4,1),comp_string,;repmat('S2 ', 4,1),comp_string],'Location','best')
if print_figs
print(['Composition_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
end
fn=fn+1;

% Plot x-y space for each species
figure(fn)
clf(fn)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560]);
plot(x1,y1,'--^',x2,y2,':sq',[0:0.1:0.7],[0:0.1:0.7],'--k')
hold on
plot(ones(4,1)*x1(1,:),ones(4,1)*y1(1,:),'o','MarkerSize',16)
axis([0 0.7 0 0.7])
title(['x-y plot, Test Case ',num2str(testcase)])
xlabel('Liquid molar fraction, x [ ]')
ylabel('Vapour molar fraction, y [ ]')
legend([repmat('S1 ', 4,1),comp_string,;repmat('S2 ', 4,1),comp_string],'Location','best')
grid on
if print_figs
print(['XY_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
end
fn=fn+1;


% Phase outflows
figure(fn)
clf(fn)
subplot(121)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980]);
plot(tt1,V1/V1(1),'--^',tt1,L1/L1(1),'--x',tt2,V2/V2(1),':sq',tt2,L2/L2(1),':o')
axis tight
title(['Dimensionless outflows, Test Case  ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Dimensionless phase outflow, [ ]')
legend('S1 V/V(0)','S1 L/L(0)','S2 V/V(0)','S2 L/L(0)')
grid on
subplot(122)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;]);
plot(tt1,M1/M1(1),'--^',tt2,M2/M2(1),':sq')
axis tight
title(['Dimensionless hold-ups, Test Case  ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Dimensionless hold-up, [ ]')
legend('S1 M/M(0)','S2 M/M(0)')
grid on
if print_figs
print(['Outflows_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
end
fn =fn+1;


% hold-up and volume split for S1
figure(fn)
clf(fn)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;0.9290,0.6940,0.1250]);
plot(tt1,M1,'--^',tt1,Mv,'--d',tt1,Ml,'--x')
axis tight
title(['Hold-up by phase for System 1, Test Case ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Phase hold-up, [kmol]')
legend('Total','Vapour','Liquid')
if print_figs
print(['S1comp_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
end
fn =fn+1;


% figure(fn)
% clf(fn)
% set(gcf,'color','white')
% if testcase ==3
% ha = area(tt1(1:55:end),[VOL_l(1:55:end),VOL_v(1:55:end)]/1.5);
% else
%     ha = area(tt1,[VOL_l,VOL_v]/1.5);
% end
% 
% ha(1).FaceColor = [0,0.4470, 0.7410];
% ha(2).FaceColor = [0.8500,0.3250,0.0980];
% axis tight
% legend('Liquid','Vapour')
% title(['Hold-up by volume for System 1, Test Case  ',num2str(testcase)])
% xlabel('Time, t [h]')
% ylabel('Fraction volume occupied, [ ]')
% if print_figs
% print(['S1volfrac_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
% end
% fn =fn+1;






% Flash ratio and separtion ratio
figure(fn)
clf(fn)
subplot(121)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980]);
plot(tt1,FR1,'--^',tt2,FR2,':sq')
axis tight
title(['Flash ratio, Test Case ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Flash ratio, [ ]')
grid on
legend('System 1','System 2')
subplot(122)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560]);
plot(tt1,SR1,'--^',tt2,SR2,':sq')
axis tight
title(['Separation ratio, Test Case ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Separation ratio, [ ]')
grid on
legend([repmat('S1 ', 4,1),comp_string,;repmat('S2 ', 4,1),comp_string],'Location','best')
if print_figs
print(['FRSR_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
end
fn=fn+1;

% Separation power
figure(fn)
clf(fn)
set(gcf,'color','white','defaultAxesColorOrder',[0,0.4470, 0.7410;...
    0.8500,0.3250,0.0980]);
plot(tt1,SP1,'--^',tt2,SP2,':sq')
if testcase ==2 
axis([0 48 0.9*max(SP2) 1.1*max(SP2)])
end 
title(['Separation power, Test Case ',num2str(testcase)])
xlabel('Time, t [h]')
ylabel('Separation power n-butane to n-pentane, [ ]')
grid on
legend('System 1','System 2')
if print_figs
print(['SP_Testcase_',num2str(testcase),'_Integrator_',integrator],imageF)
end
fn =fn+1;





end

