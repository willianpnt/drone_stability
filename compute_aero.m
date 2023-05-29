function [Cl,Cd,Cm] = compute_aero(tailsit,delta_e)

alpha = tailsit.general.alpha;
xcg = tailsit.general.cg(1);
zcg = tailsit.general.cg(2);
cbar = tailsit.general.cbar;

delta_sup = 5*ceil(delta_e/5);
delta_inf = 5*floor(delta_e/5);

if abs(delta_sup) > 15 || abs(delta_inf) > 15
    Cl = zeros(15,1);
    Cd = zeros(15,1);
    Cm = zeros(15,1);
    return
elseif isnan(delta_sup) || isnan(delta_inf)
    Cl = zeros(15,1)*NaN;
    Cd = zeros(15,1)*NaN;
    Cm = zeros(15,1)*NaN;
    return
end

if delta_sup < 0
    name_sup = strcat('delta_',num2str(-delta_sup),'neg');
    polar_sup = tailsit.wing.(name_sup).Variables;
    name_inf = strcat('delta_',num2str(-delta_inf),'neg');
    polar_inf = tailsit.wing.(name_inf).Variables;
elseif delta_inf < 0
    name_sup = strcat('delta_',num2str(delta_sup));
    polar_sup = tailsit.wing.(name_sup).Variables;
    name_inf = strcat('delta_',num2str(-delta_inf),'neg');
    polar_inf = tailsit.wing.(name_inf).Variables;
else    
    name_sup = strcat('delta_',num2str(delta_sup));
    polar_sup = tailsit.wing.(name_sup).Variables;
    name_inf = strcat('delta_',num2str(delta_inf));
    polar_inf = tailsit.wing.(name_inf).Variables;
end

if delta_sup == delta_inf
    polar_w = polar_sup;
else
    for i = 1:min(size(polar_sup,1),size(polar_inf,1))
        for k=1:13 
       polar_w(i,k) = interp1([delta_inf, delta_sup],[polar_inf(i,k), polar_sup(i,k)],delta_e,'pchip'); 
    end
    end
end

Cl = interp1(polar_w(:,1), polar_w(:,3), alpha,'pchip');
Cd = interp1(polar_w(:,1), polar_w(:,6), alpha, 'pchip'); 
Cmw = interp1(polar_w(:,1), polar_w(:,9), alpha, 'pchip'); 
alpha = pi/180 * alpha;

Cm = (Cl.*cos(alpha) +Cd.*sin(alpha))*xcg/cbar + (Cl.*sin(alpha)-Cd.*cos(alpha))*zcg/cbar + Cmw;

return
end