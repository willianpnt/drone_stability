clear all
close all
clc

%% Aircraft and Engine data
% Reference system pos: 1/4 of root chord from the leading edge

tailsit.general.m = 2.340;
tailsit.general.cbar = 0.172;
tailsit.general.density = 1.050;
tailsit.general.Sref = 0.102; 
tailsit.general.cg = [0.001 0]; % CG Position
tailsit.general.speed = 18; % True airspeed
tailsit.general.ME_min = 5;
tailsit.general.Iyy = 0.016;
tailsit.general.pdin = 0.5*tailsit.general.density*(tailsit.general.speed)^2;

%% Wing data

tailsit.wing.b = 0.6; % Wingspan
tailsit.wing.S = 0.102; % Wing's Area
tailsit.wing.AR = tailsit.wing.b^2/tailsit.wing.S; %Aspect Ratio
tailsit.wing.cr = 0.2; % Root chord
tailsit.wing.ct = 0.14; % Tip chord
tailsit.wing.cbar = 0.17; % Mean aerodynammic chord
tailsit.wing.position = [0 0]; %Wing's CA position
tailsit.wing.incidence = -2; % Wing's incidence
tailsit.wing.max_deflexion = 15;
tailsit.wing.min_deflexion = -15;
tailsit.wing.deflexion = tailsit.wing.min_deflexion:5:tailsit.wing.max_deflexion; % Elevon's deflexion

tailsit.wing.delta_0 = readtable('Polars/A1_0.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f'); % Polar aerodymanic: no elevon's deflexion
tailsit.wing.delta_5 = readtable('Polars/A1_5.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f'); % Polar aerodynamic: 5o elevon's deflexion
tailsit.wing.delta_10 = readtable('Polars/A1_10.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f');  % Polar aerodynamic: 10o elevon's deflexion
tailsit.wing.delta_15 = readtable('Polars/A1_15.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f');  % Polar aerodynamic: 15o elevon's deflexion
tailsit.wing.delta_5neg = readtable('Polars/A1_5neg.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f');  % Polar aerodynamic: -5o elevon's deflexion
tailsit.wing.delta_10neg = readtable('Polars/A1_10neg.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f');  % Polar aerodynamic: -10o elevon's deflexion
tailsit.wing.delta_15neg = readtable('Polars/A1_15neg.txt','HeaderLines',5,'Format','%f%f%f%f%f%f%f%f%f%f%f%f%f');
tailsit.wing.alpha = tailsit.wing.delta_0.Variables;
tailsit.wing.alpha = tailsit.wing.alpha(:,1);
tailsit.general.alpha = tailsit.wing.alpha - tailsit.wing.incidence;

%% Static Stability

% Longitudinal

% Aerodynamic matrix
for delta_e = tailsit.wing.min_deflexion:tailsit.wing.max_deflexion
    [L(:,delta_e+16),D(:,delta_e+16),M(:,delta_e+16)] = compute_aero(tailsit, delta_e);
end

for i = 1:size(M,1)
    deltae_trim(i,:) = interp1(M(i,:),tailsit.wing.min_deflexion:tailsit.wing.max_deflexion,0);
    CL_trim(i,:) = interp1(tailsit.wing.min_deflexion:tailsit.wing.max_deflexion, L(i,:),deltae_trim(i));
    CD_trim(i,:) = interp1(tailsit.wing.min_deflexion:tailsit.wing.max_deflexion, D(i,:),deltae_trim(i));
    ME(i,:) = polyfit(L(:,i),M(:,i),1);
end 
ME = -ME(:,1)*100;
xPN = ME*tailsit.general.cbar/100 + tailsit.general.cg(1);

tailsit.static.LMatrix = L;
tailsit.static.DMatrix = D;
tailsit.static.MMatrix = M;
tailsit.static.deltae_trim = deltae_trim;
tailsit.trim = table(tailsit.wing.alpha, deltae_trim, CL_trim, CD_trim, (CL_trim./CD_trim), ME, xPN);
tailsit.trim.Properties.VariableNames = {'alpha', 'delta_e', 'C_L trim', 'C_D trim', 'Eff. Aero', 'Static Margin','Neutral Point Pos.'};
writetable(tailsit.trim,'A1_trim_polar.csv');
%% Dynamic Stability

% Flight Condition definition
alpha = 0;
trim_polar = tailsit.trim.Variables;
delta_e = interp1(trim_polar(:,1),trim_polar(:,2),alpha);
Cl = interp1(trim_polar(:,1),trim_polar(:,3),alpha);
Cd = interp1(trim_polar(:,1),trim_polar(:,4),alpha);
Cw0 = tailsit.general.m*9.81/(tailsit.general.pdin*tailsit.general.Sref);
h = tailsit.general.cg(1)/tailsit.general.cbar;
u0 = tailsit.general.speed;
w0 = 0;
Cz = -Cl*cos(alpha*pi/180)-Cd*sin(alpha*pi/180);
Cx = Cl*sin(alpha*pi/180)-Cd*cos(alpha*pi/180);
theta0 = 0; %acos(-Cz/Cw0);
tailsit.general.thrust = (Cw0*sin(theta0)-Cx)*tailsit.general.pdin*tailsit.general.Sref;


% Aerodynamic derivate
[Cl_fl, Cd_fl, Cm_fl] = compute_aero(tailsit,delta_e);
p_Cl = polyfit(tailsit.general.alpha*pi/180, Cl_fl,1);
Cl0 = p_Cl(2);
Cla = p_Cl(1);
Clde = polyfit((tailsit.wing.min_deflexion:tailsit.wing.max_deflexion)*pi/180,L(floor(alpha)+6,:),1);
Clde = Clde(1);

p_Cd = fit2(Cl_fl,Cd_fl);
k = p_Cd(2);
Cd0 = p_Cd(1);
Cda = 2*k*Cl*Cla;
Cdde = polyfit((tailsit.wing.min_deflexion:tailsit.wing.max_deflexion)*pi/180,D(floor(alpha)+6,:),1);
Cdde = Clde(1);

p_Cm = polyfit(tailsit.general.alpha*pi/180, Cm_fl,1);
Cm0 = p_Cm(2);
Cma = p_Cm(1);
Cmde = polyfit((tailsit.wing.min_deflexion:tailsit.wing.max_deflexion)*pi/180,M(floor(alpha)+6,:),1);
Cmde = Clde(1);

% Stability derivates
% Dimensioneless
F_target = Cla/(2*pi);
func = @(k)theodorsen(k)-F_target;
k = secante(func,0,1,1e-6,100);
[Fk,Gk] = theodorsen(k);    

Cladot = pi + 2*pi*Gk/k;
Clq = -2*Cla*(h-3/4);
Cmadot = pi*(h-1/2) + 2*pi*Gk/k*(h-1/4);
Cmq = -2*Cla*(h-1/2)^2;

% Dimensional Derivate

% u derivate
Cxu = 0;
Czu = 0;
Cmu = 0;

Xu = tailsit.general.density*u0*tailsit.general.Sref*Cw0*sin(theta0)+...
    0.5*tailsit.general.density*u0*tailsit.general.Sref*Cxu;
Zu = -tailsit.general.density*u0*tailsit.general.Sref*Cw0*cos(theta0)+...
    0.5*tailsit.general.density*u0*tailsit.general.Sref*Czu;
Mu = 0.5*tailsit.general.density*u0*tailsit.general.Sref*   Cmu*tailsit.general.cbar;
    
% w derivate
Cxa = Cl*cos(alpha*pi/180) + Cla*sin(alpha*pi/180) + Cd*sin(alpha*pi/180)-Cda*cos(alpha*pi/180);
Cza = -(Cla*cos(alpha*pi/180) - Cl*sin(alpha*pi/180)) + Cda*sin(alpha*pi/180) + Cd*cos(alpha*pi/180);

Xw = 0.5*tailsit.general.density*u0*tailsit.general.Sref*Cxa;
Zw = 0.5*tailsit.general.density*u0*tailsit.general.Sref*Cza;
Mw = 0.5*tailsit.general.density*u0*tailsit.general.Sref*Cma*tailsit.general.cbar;

%wdot derivate
Cxadot = Cladot*sin(alpha*pi/180);
Czadot = -(Cladot*cos(alpha*pi/180));

Xwdot = 0.25*tailsit.general.density*tailsit.general.cbar*tailsit.general.Sref*Cxadot;
Zwdot = 0.5*tailsit.general.density*tailsit.general.cbar*tailsit.general.Sref*Czadot;
Mwdot = 0.5*tailsit.general.density*tailsit.general.Sref*Cmadot*tailsit.general.cbar^2;

%q derivate
Cxq = Clq*sin(alpha*pi/180);
Czq = -(Clq*cos(alpha*pi/180));

Xq = 0.25*tailsit.general.density*u0*tailsit.general.Sref*Cxq*tailsit.general.cbar;
Zq = 0.25*tailsit.general.density*u0*tailsit.general.Sref*Czq*tailsit.general.cbar;
Mq = 0.25*tailsit.general.density*u0*tailsit.general.Sref*Cmq*tailsit.general.cbar^2;

%delta_e derivate
Cxde = Clde*sin(alpha*pi/180) - Cdde*cos(alpha*pi/180);
Czde = -Clde*cos(alpha*pi/180) - Cdde*sin(alpha*pi/180);

Xde = tailsit.general.pdin*tailsit.general.Sref*Cxde;
Zde = tailsit.general.pdin*tailsit.general.Sref*Czde;
Mde = tailsit.general.pdin*tailsit.general.Sref*tailsit.general.cbar*Cmde;

%Matrix Longitudinal

ALong = [Xu/tailsit.general.m Xw/tailsit.general.m 0 -9.81*cos(theta0);...
    Zu/(tailsit.general.m-Zwdot) Zw/(tailsit.general.m-Zwdot) (Zq+tailsit.general.m*u0)/(tailsit.general.m-Zwdot)...
    -tailsit.general.m*9.81*sin(theta0)/(tailsit.general.m-Zwdot);...
    1/tailsit.general.Iyy*(Mu + (Mwdot*Zu)/(tailsit.general.m-Zwdot)) 1/tailsit.general.Iyy*(Mw + (Mwdot*Zw)/(tailsit.general.m-Zwdot))...
    1/tailsit.general.Iyy*(Mq + (Mwdot*(Zq + tailsit.general.m*u0))/(tailsit.general.m-Zwdot)) ...
    -Mw*tailsit.general.m*9.81*sin(theta0)/(tailsit.general.Iyy*(tailsit.general.m-Zwdot));...
    0 0 1 0];

BLong = [Xde/tailsit.general.m; Zde/(tailsit.general.m-Zwdot); Mde/tailsit.general.Iyy+Mwdot/tailsit.general.Iyy*(Zde)/(tailsit.general.m-Zwdot); 0];

Mode_Long = eig(ALong);
tailsit.dynamic.MLong = ALong;
[~,index] = min(abs(imag(Mode_Long)));
tailsit.dynamic.phugoid.eigenvalue = [Mode_Long(index); conj(Mode_Long(index))];
tailsit.dynamic.phugoid.Frequency = imag(Mode_Long(index));
tailsit.dynamic.phugoid.Damping = -cos(atan2(abs(imag(Mode_Long(index))),real(Mode_Long(index))));
tailsit.dynamic.phugoid.wn = abs(Mode_Long(index));
[~,index] = max(abs(imag(Mode_Long)));
tailsit.dynamic.short_period.eigenvalue = [Mode_Long(index); conj(Mode_Long(index))];
tailsit.dynamic.short_period.Frequency = imag(Mode_Long(index));
tailsit.dynamic.short_period.Damping = -cos(atan2(abs(imag(Mode_Long(index))),real(Mode_Long(index))));
tailsit.dynamic.short_period.wn = abs(Mode_Long(index));
[tailsit.dynamic.tfnum, tailsit.dynamic.tfden] = ss2tf(ALong, BLong, eye(4), zeros(4,1));
tailsit.dynamic.Gude = tf(tailsit.dynamic.tfnum(1,:),tailsit.dynamic.tfden);
tailsit.dynamic.Gwde = tf(tailsit.dynamic.tfnum(2,:),tailsit.dynamic.tfden);
tailsit.dynamic.Guwdote = tf(tailsit.dynamic.tfnum(3,:),tailsit.dynamic.tfden);
tailsit.dynamic.Gqde = tf(tailsit.dynamic.tfnum(4,:),tailsit.dynamic.tfden);
%% Simulink Simulation

t_sim = 100;
sim('long.slx');

%% Plots

        figure()
        plot(tailsit.general.alpha', deltae_trim, 'LineWidth', 1.5);
        xlabel('\alpha [deg.]'); ylabel('\delta_e [deg.]');
        title('Elevon''s deflexion to trim the drone');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        plot(tailsit.general.alpha', CL_trim, 'LineWidth', 1.5);
        xlabel('\alpha [deg.]'); ylabel('C_{Ltrim} [-]');
        title('C_L trimmed');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        plot(CD_trim,CL_trim,'LineWidth', 1.5);
        ylabel('C_{Ltrim} [-]'); xlabel('C_{Dtrim} [-]');
        title('Drag Polar');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        plot(tailsit.general.alpha', CL_trim./CD_trim,'LineWidth', 1.5);
        xlabel('\alpha [-]'); ylabel('\eta [-]');
        title('Aerodynamic Efficiency');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        hold on
        plot(CL_trim, ME,'LineWidth', 1.5);
        plot(CL_trim, ones(size(CL_trim,1),1)*tailsit.general.ME_min,'--','LineWidth', 1.5)
        xlabel('C_{Ltrim} [-]'); ylabel('ME [%]');
        title('Static Margin');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);
        hold off

        figure() % Longitudinal Abacus
        hold on
        legend('on');
        set(legend,'FontSize',10);
        set(legend,'Location','Northeast Outside');
        for i = 1:11
            leg = strcat('\delta_e = ',num2str(-12+2*i),'°');
             plot(L(:,2*(i+2)),M(:,2*(i+2)),'LineWidth',1.5,'DisplayName',leg);
        end
        title('Longitudinal Trim Abacus')
        xlabel('C_L [-]'); ylabel('C_M [-]');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);
        hold off

        figure()
        plot(AoA_Sim.time,AoA_Sim.signals.values,'LineWidth',1.5);
        xlabel('t [s]'); ylabel('\alpha [deg.]');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        plot(u_Sim.time,u_Sim.signals.values,'LineWidth',1.5);
        xlabel('t [s]'); ylabel('u [m/s]');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        plot(w_Sim.time,w_Sim.signals.values,'LineWidth',1.5);
        xlabel('t [s]'); ylabel('w [m/s]');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);
        
        figure()
        plot(wdot.time,wdot.signals.values,'LineWidth',1.5);
        xlabel('t [s]','Interpreter','latex'); ylabel('$\dot{w}$ $[m/s^2]$','Interpreter','latex');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        plot(theta_Sim.time,theta_Sim.signals.values,'LineWidth',1.5);
        xlabel('t [s]'); ylabel('\theta [deg.]');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);

        figure()
        t = linspace(0,t_sim,1000);
        y = step(tailsit.dynamic.Gude,t);
        plot(t,y,'LineWidth',1.5);
        xlabel('t [s]'); ylabel('\Deltau [m/s]');
        grid on
        grid minor
        box on
        set(gcf,'Color', [1 1 1]);
% 
%         figure()
%         t = linspace(0,t_sim,1000);
%         y = step(tailsit.dynamic.Gwde,t);
%         plot(t,y,'LineWidth',1.5);
%         xlabel('t [s]'); ylabel('\Deltaw [m/s]');
%         grid on
%         grid minor
%         box on
%         set(gcf,'Color', [1 1 1]);
       
%     end
% end