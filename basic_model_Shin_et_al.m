
x = 0:1:72; % nm did 600 to compare to Shin et al (still very high e- temperatures, still think due to bad C_e calculations (multiplicity of bands?) and non-3D no-flux condition constraining high quantities of charge carriers to a localized space
t = 0:0.01:500; % ps

m = 0;
%options = odeset('MaxOrder',5);
sol = pdepe_opt(m,@pdefun,@pdeic,@pdebc,x,t); % ,options

figure; surf(x,t,log(10.^(sol(:,:,1))*1e21)/log(10),'edgecolor','none');
xlabel('x (nm)','FontSize',20);
ylabel('t (ps)','FontSize',20);
zlabel('log(carrier concentration)');
figure; surf(x,t,sol(:,:,2),'edgecolor','none');
xlabel('x (nm)','FontSize',20);
ylabel('t (ps)','FontSize',20);
zlabel('Lattice Temperature (K)');
figure; surf(x,t,sol(:,:,3),'edgecolor','none');
xlabel('x (nm)','FontSize',20);
ylabel('t (ps)','FontSize',20);
zlabel('Electron Temperature (K)');

%% RUNNING INTO STABILITY ISSUE --> CHANGE EVERYTHING TO nm, eV, ps

function [c,f,s] = pdefun(x,t,u,dudx)

JpeV = 1.60218e-19; % J/eV
R = 0.508; % dimensionless
A = pi*(84.7458*1000)^2; % nm^2
F = 1.3741/1000*(100^2)/JpeV/(1e9)^2; % eV/nm^2
temporal_stdev = 0.700/2/sqrt(2*log(2)); % ps
h = (4.135667696e-15)*1e12; % eV*ps
freq = 2.9979e8/515e-9/1e12; % 1/ps
sigma = (1/5.97e-2); % nm

E_c = 0.66; % eV
E_v = 0; % eV
E_g = E_c-E_v;

m_o = 9.10938356e-31; % kg
m_e = 0.22*m_o; % kg

n_e = 10^(u(1)); % 1/nm^3

n_g = 2.4e13*(100^3)/(1e9)^3; % 1/nm^3
n_auger_piecewise = 1e18*(100^3)/(1e9)^3; % 1/nm^3
n_lim = (1e20)*(100^3)/(1e9)^3; % 1/nm^3
n_lim2 = (5e19)*(100^3)/(1e9)^3; % 1/nm^3
n_lim3 = (1e13)*(100^3)/(1e9)^3; % 1/nm^3
n_lim4 = (5e13)*(100^3)/(1e9)^3; % 1/nm^3

n_e_norm = n_e*(1e7)^3; % 1/cm^3
n_e_norm_lim1 = n_lim*(1e7)^3; % 1/cm^3
n_e_norm_lim2 = n_lim2*(1e7)^3; % 1/cm^3
n_e_norm_lim3 = n_lim3*(1e7)^3; % 1/cm^3
n_e_norm_lim4 = n_lim4*(1e7)^3; % 1/cm^3
n_g_norm = n_g*(1e7)^3; % 1/cm^3

tempCe_lim1 = (log(20000)/log(10) - 3.876)/0.4073;
tempCe_lim2 = (log(19950)/log(10) - 3.876)/0.4073;
tempCe_lim3 = (log(78)/log(10) - 3.876)/0.4073;
tempCe_lim4 = (log(100)/log(10) - 3.876)/0.4073;
tempCe = (log(u(3))/log(10) - 3.876)/0.4073;
p1Ce =   3.765e-05;
p2Ce =   0.0006056;
p3Ce =    0.004413;
p4Ce =     0.01654;
p5Ce =     0.02823;
p6Ce =    -0.03744;
p7Ce =      -0.147;
p8Ce =     -0.2213;
p9Ce =      0.4762;
p10Ce =      4.234;
if (u(3) > 20000)
    C_e_lim1 = p1Ce*tempCe_lim1^9 + p2Ce*tempCe_lim1^8 + p3Ce*tempCe_lim1^7 + p4Ce*tempCe_lim1^6 + p5Ce*tempCe_lim1^5 + p6Ce*tempCe_lim1^4 + p7Ce*tempCe_lim1^3 + p8Ce*tempCe_lim1^2 + p9Ce*tempCe_lim1^1 + p10Ce;
    C_e_lim2 = p1Ce*tempCe_lim2^9 + p2Ce*tempCe_lim2^8 + p3Ce*tempCe_lim2^7 + p4Ce*tempCe_lim2^6 + p5Ce*tempCe_lim2^5 + p6Ce*tempCe_lim2^4 + p7Ce*tempCe_lim2^3 + p8Ce*tempCe_lim2^2 + p9Ce*tempCe_lim2^1 + p10Ce;
    C_e = 10^((C_e_lim1 - C_e_lim2)/(tempCe_lim1 - tempCe_lim2)*(tempCe - tempCe_lim2) + C_e_lim2)/JpeV/(1e9)^3/n_e; % eV/K
elseif (u(3) < 78)
    C_e_lim3 = p1Ce*tempCe_lim3^9 + p2Ce*tempCe_lim3^8 + p3Ce*tempCe_lim3^7 + p4Ce*tempCe_lim3^6 + p5Ce*tempCe_lim3^5 + p6Ce*tempCe_lim3^4 + p7Ce*tempCe_lim3^3 + p8Ce*tempCe_lim3^2 + p9Ce*tempCe_lim3^1 + p10Ce;
    C_e_lim4 = p1Ce*tempCe_lim4^9 + p2Ce*tempCe_lim4^8 + p3Ce*tempCe_lim4^7 + p4Ce*tempCe_lim4^6 + p5Ce*tempCe_lim4^5 + p6Ce*tempCe_lim4^4 + p7Ce*tempCe_lim4^3 + p8Ce*tempCe_lim4^2 + p9Ce*tempCe_lim4^1 + p10Ce;
    C_e = 10^((C_e_lim3 - C_e_lim4)/(tempCe_lim3 - tempCe_lim4)*(tempCe - tempCe_lim4) + C_e_lim4)/JpeV/(1e9)^3/n_e; % eV/K
else
    C_e = 10^(p1Ce*tempCe^9 + p2Ce*tempCe^8 + p3Ce*tempCe^7 + p4Ce*tempCe^6 + p5Ce*tempCe^5 + p6Ce*tempCe^4 + p7Ce*tempCe^3 + p8Ce*tempCe^2 + p9Ce*tempCe^1 + p10Ce)/JpeV/(1e9)^3/n_e; % eV/K
end

if (n_e > n_auger_piecewise)
    C_eeh = (10^190.9124*n_e_norm^(0.6698*log(n_e_norm)/log(10) - 21.5429 - 1))/(n_e_norm^2-n_g_norm^2); % cm^6/s
    C_hhe = (10^120.1814*n_e_norm^(0.4574*log(n_e_norm)/log(10) - 13.6881 - 1))/(n_e_norm^2-n_g_norm^2); % cm^6/s
    k = ((C_eeh + C_hhe)*n_e_norm^2)/1e12; % function of n_e, Dominici et. al., results in 1/ps units
elseif (n_e <= n_auger_piecewise && n_e >= n_g)
    C_eeh = (10^(-30.9373)*n_e_norm^(0.0110*log(n_e_norm)/log(10) + 2.6402 - 1))/(n_e_norm^2-n_g_norm^2);
    C_hhe = (10^(-28.9705)*n_e_norm^(0.0114*log(n_e_norm)/log(10) + 2.6257 - 1))/(n_e_norm^2-n_g_norm^2);
    k = ((C_eeh + C_hhe)*n_e^2)/1e12;
else
    k = 0;
end

n_e_diff_100 = (log(n_e_norm)/log(10)-17.27)/2.215; % (n-17.27)/2.215
n_e_diff_100_lim1 = (log(n_e_norm_lim1)/log(10)-17.27)/2.215;
n_e_diff_100_lim2 = (log(n_e_norm_lim2)/log(10)-17.27)/2.215;
n_e_diff_100_lim3 = (log(n_e_norm_lim3)/log(10)-17.27)/2.215;
n_e_diff_100_lim4 = (log(n_e_norm_lim4)/log(10)-17.27)/2.215;
p1_100 = 21.34;
p2_100 = 113;
p3_100 = 218.2;
p4_100 = 164.1;
p5_100 = -20.31;
p6_100 = 246.6;
log_diff_100k_lim1 = p1_100*n_e_diff_100_lim1^5 + p2_100*n_e_diff_100_lim1^4 + p3_100*n_e_diff_100_lim1^3 + p4_100*n_e_diff_100_lim1^2 + p5_100*n_e_diff_100_lim1^1 + p6_100;
log_diff_100k_lim2 = p1_100*n_e_diff_100_lim2^5 + p2_100*n_e_diff_100_lim2^4 + p3_100*n_e_diff_100_lim2^3 + p4_100*n_e_diff_100_lim2^2 + p5_100*n_e_diff_100_lim2^1 + p6_100;
log_diff_100k_lim3 = p1_100*n_e_diff_100_lim3^5 + p2_100*n_e_diff_100_lim3^4 + p3_100*n_e_diff_100_lim3^3 + p4_100*n_e_diff_100_lim3^2 + p5_100*n_e_diff_100_lim3^1 + p6_100;
log_diff_100k_lim4 = p1_100*n_e_diff_100_lim4^5 + p2_100*n_e_diff_100_lim4^4 + p3_100*n_e_diff_100_lim4^3 + p4_100*n_e_diff_100_lim4^2 + p5_100*n_e_diff_100_lim4^1 + p6_100;

if (n_e > n_lim)
    log_diff_100k = (log_diff_100k_lim1 - log_diff_100k_lim2)/(n_e_diff_100_lim1 - n_e_diff_100_lim2)*(n_e_diff_100-n_e_diff_100_lim2) + log_diff_100k_lim2;
elseif (n_e < n_lim3)
    log_diff_100k = (log_diff_100k_lim4 - log_diff_100k_lim3)/(n_e_diff_100_lim4 - n_e_diff_100_lim3)*(n_e_diff_100-n_e_diff_100_lim3) + log_diff_100k_lim3;
else
    log_diff_100k = p1_100*n_e_diff_100^5 + p2_100*n_e_diff_100^4 + p3_100*n_e_diff_100^3 + p4_100*n_e_diff_100^2 + p5_100*n_e_diff_100^1 + p6_100;
end 

n_e_diff_300 = (log(n_e_norm)/log(10)-17.73)/2.478; % (n-17.73)/2.478
n_e_diff_300_lim1 = (log(n_e_norm_lim1)/log(10)-17.73)/2.478;
n_e_diff_300_lim2 = (log(n_e_norm_lim2)/log(10)-17.73)/2.478;
n_e_diff_300_lim3 = (log(n_e_norm_lim3)/log(10)-17.73)/2.478;
n_e_diff_300_lim4 = (log(n_e_norm_lim4)/log(10)-17.73)/2.478;
p1_300 = 7.81;
p2_300 = 33.92;
p3_300 = 49.47;
p4_300 = 25.18;
p5_300 = -3.337;
p6_300 = 62.31;
log_diff_300k_lim1 = p1_300*n_e_diff_300_lim1^5 + p2_300*n_e_diff_300_lim1^4 + p3_300*n_e_diff_300_lim1^3 + p4_300*n_e_diff_300_lim1^2 + p5_300*n_e_diff_300_lim1^1 + p6_300;
log_diff_300k_lim2 = p1_300*n_e_diff_300_lim2^5 + p2_300*n_e_diff_300_lim2^4 + p3_300*n_e_diff_300_lim2^3 + p4_300*n_e_diff_300_lim2^2 + p5_300*n_e_diff_300_lim2^1 + p6_300;
log_diff_300k_lim3 = p1_300*n_e_diff_300_lim3^5 + p2_300*n_e_diff_300_lim3^4 + p3_300*n_e_diff_300_lim3^3 + p4_300*n_e_diff_300_lim3^2 + p5_300*n_e_diff_300_lim3^1 + p6_300;
log_diff_300k_lim4 = p1_300*n_e_diff_300_lim4^5 + p2_300*n_e_diff_300_lim4^4 + p3_300*n_e_diff_300_lim4^3 + p4_300*n_e_diff_300_lim4^2 + p5_300*n_e_diff_300_lim4^1 + p6_300;

if (n_e > n_lim)
    log_diff_300k = (log_diff_300k_lim1 - log_diff_300k_lim2)/(n_e_diff_300_lim1 - n_e_diff_300_lim2)*(n_e_diff_300 - n_e_diff_300_lim2) + log_diff_300k_lim2;
elseif (n_e < n_lim3)
    log_diff_300k = (log_diff_300k_lim4 - log_diff_300k_lim3)/(n_e_diff_300_lim4 - n_e_diff_300_lim3)*(n_e_diff_300 - n_e_diff_300_lim3) + log_diff_300k_lim3;
else
    log_diff_300k = p1_300*n_e_diff_300^5 + p2_300*n_e_diff_300^4 + p3_300*n_e_diff_300^3 + p4_300*n_e_diff_300^2 + p5_300*n_e_diff_300^1 + p6_300;
end 

%D_e_old = ((log_diff_300k - log_diff_100k)/(300 - 100)*(u(3)-100) + log_diff_100k)/(100^2)*(1e9)^2/1e12
D_e = (10^((log(log_diff_300k) - log(log_diff_100k))/(300 - 100)*(u(2)-100) + log(log_diff_100k))/(100^2)*(1e9)^2)/1e12; % nm^2/ps Young et. al.

n_i = 5.323/72.630*6.022e23*(100^3)/(1e9)^3; % 1/nm^3

logLatTemp = log(u(2))/log(10);
logRoomTemp = log(293.7)/log(10);

K_i_temp_adj = (logLatTemp - 1.976)/0.7905;
K_i_temp_adj_RT = (logRoomTemp - 1.976)/0.7905;
K_i_temp_lim_1 = (log(1290)/log(10) - 1.976)/0.7905;
K_i_temp_lim_2 = (log(1300)/log(10) - 1.976)/0.7905;
p1Ki =     0.01207;
p2Ki =     0.04158;
p3Ki =    -0.00384;
p4Ki =     -0.1541;
p5Ki =      0.1266;
p6Ki =      0.0513;
p7Ki =      -1.026;
p8Ki =      0.3893;
log_therm_cond_fit_RT = p1Ki*K_i_temp_adj_RT^7 + p2Ki*K_i_temp_adj_RT^6 + p3Ki*K_i_temp_adj_RT^5 + p4Ki*K_i_temp_adj_RT^4 + p5Ki*K_i_temp_adj_RT^3 + p6Ki*K_i_temp_adj_RT^2 + p7Ki*K_i_temp_adj_RT^1 + p8Ki;
log_therm_cond_fit_lim_1 = p1Ki*K_i_temp_lim_1^7 + p2Ki*K_i_temp_lim_1^6 + p3Ki*K_i_temp_lim_1^5 + p4Ki*K_i_temp_lim_1^4 + p5Ki*K_i_temp_lim_1^3 + p6Ki*K_i_temp_lim_1^2 + p7Ki*K_i_temp_lim_1^1 + p8Ki;
log_therm_cond_fit_lim_2 = p1Ki*K_i_temp_lim_2^7 + p2Ki*K_i_temp_lim_2^6 + p3Ki*K_i_temp_lim_2^5 + p4Ki*K_i_temp_lim_2^4 + p5Ki*K_i_temp_lim_2^3 + p6Ki*K_i_temp_lim_2^2 + p7Ki*K_i_temp_lim_2^1 + p8Ki;
if (u(2) > 1300)
    log_therm_cond_fit = (log_therm_cond_fit_lim_2 - log_therm_cond_fit_lim_1)/(K_i_temp_lim_2 - K_i_temp_lim_1)*(K_i_temp_adj - K_i_temp_lim_1) + log_therm_cond_fit_lim_1;
else
    log_therm_cond_fit = p1Ki*K_i_temp_adj^7 + p2Ki*K_i_temp_adj^6 + p3Ki*K_i_temp_adj^5 + p4Ki*K_i_temp_adj^4 + p5Ki*K_i_temp_adj^3 + p6Ki*K_i_temp_adj^2 + p7Ki*K_i_temp_adj^1 + p8Ki; % (LT - 1.976)/0.7905
end
K_i = (10^(log_therm_cond_fit))/JpeV*100/1e9/1e12; % eV/nm/K/ps
K_i_RT = (10^(log_therm_cond_fit_RT))/JpeV*100/1e9/1e12; % eV/nm/K/ps

xi = 1; % assume. indirect bandgap

C_i_temp_adj = (logLatTemp - 2.08)/0.6879;
C_i_temp_adj_RT = (logRoomTemp - 2.08)/0.6879;
C_i_temp_lim_1 = (log(1190)/log(10) - 2.08)/0.6879;
C_i_temp_lim_2 = (log(1200)/log(10) - 2.08)/0.6879;
p1Ci =    -0.03258;
p2Ci =     -0.1327;
p3Ci =      0.1028;
p4Ci =       0.948;
p5Ci =       0.326;
p6Ci =       -2.35;
p7Ci =      -1.014;
p8Ci =       4.331;
p9Ci =       3.812;
C_i_RT = (4.184*(p1Ci*C_i_temp_adj_RT^8 + p2Ci*C_i_temp_adj_RT^7 + p3Ci*C_i_temp_adj_RT^6 + p4Ci*C_i_temp_adj_RT^5 + p5Ci*C_i_temp_adj_RT^4 + p6Ci*C_i_temp_adj_RT^3 + p7Ci*C_i_temp_adj_RT^2 + p8Ci*C_i_temp_adj_RT^1 + p9Ci))/6.022e23/JpeV;
C_i_lim_1 = (4.184*(p1Ci*C_i_temp_lim_1^8 + p2Ci*C_i_temp_lim_1^7 + p3Ci*C_i_temp_lim_1^6 + p4Ci*C_i_temp_lim_1^5 + p5Ci*C_i_temp_lim_1^4 + p6Ci*C_i_temp_lim_1^3 + p7Ci*C_i_temp_lim_1^2 + p8Ci*C_i_temp_lim_1^1 + p9Ci))/6.022e23/JpeV;
C_i_lim_2 = (4.184*(p1Ci*C_i_temp_lim_2^8 + p2Ci*C_i_temp_lim_2^7 + p3Ci*C_i_temp_lim_2^6 + p4Ci*C_i_temp_lim_2^5 + p5Ci*C_i_temp_lim_2^4 + p6Ci*C_i_temp_lim_2^3 + p7Ci*C_i_temp_lim_2^2 + p8Ci*C_i_temp_lim_2^1 + p9Ci))/6.022e23/JpeV;
if (u(2) > 1200)
    C_i = (C_i_lim_2 - C_i_lim_1)/(C_i_temp_lim_2 - C_i_temp_lim_1)*(C_i_temp_adj - C_i_temp_lim_1) + C_i_lim_1;
else
    C_i_fit = p1Ci*C_i_temp_adj^8 + p2Ci*C_i_temp_adj^7 + p3Ci*C_i_temp_adj^6 + p4Ci*C_i_temp_adj^5 + p5Ci*C_i_temp_adj^4 + p6Ci*C_i_temp_adj^3 + p7Ci*C_i_temp_adj^2 + p8Ci*C_i_temp_adj^1 + p9Ci; % cal/K/mol, (LT - 2.08)/0.6879
    C_i = (4.184*C_i_fit)/6.022e23/JpeV; % eV/K
end

tauTemp = (u(3)-650)/260.2;
tauTemp_lim_1 = (1090-650)/260.2;
tauTemp_lim_2 = (1100-650)/260.2;
p1tau = 0.001878;
p2tau = -0.007429;
p3tau = 0.03456;
p4tau = -0.2661;
p5tau = -13.51;
tau_fit_lim_1 = (p1tau*tauTemp_lim_1^4 + p2tau*tauTemp_lim_1^3 + p3tau*tauTemp_lim_1^2 + p4tau*tauTemp_lim_1^1 + p5tau);
tau_fit_lim_2 = (p1tau*tauTemp_lim_2^4 + p2tau*tauTemp_lim_2^3 + p3tau*tauTemp_lim_2^2 + p4tau*tauTemp_lim_2^1 + p5tau);
if (u(3) > 1100)
    tau_fit = (10^((tau_fit_lim_2 - tau_fit_lim_1)/(tauTemp_lim_2 - tauTemp_lim_1)*(tauTemp - tauTemp_lim_1) + tau_fit_lim_1))*1e12; % ps
else
    tau_fit = (10^(p1tau*tauTemp^4 + p2tau*tauTemp^3 + p3tau*tauTemp^2 + p4tau*tauTemp^1 + p5tau))*1e12; % ps
end

v_s = 5400*1e9/1e12; % nm/ps

g = (pi^2*m_e*n_e*v_s^2/6/tau_fit/u(2))*(1e3)^2/JpeV; % kg/nm/ps^3/K*(1000)^2/JpeV for correction factor

n_e_dot = (1-R)*F/h/freq/sigma/temporal_stdev/sqrt(2*pi)*exp(-(t-10*temporal_stdev)^2/2/temporal_stdev^2)*exp(-x/sigma); % 1/nm^3/ps

c = [log(10)*n_e;n_i*C_i;n_e*C_e];
f = [D_e*dudx(1)*log(10)*(log(10)+1)*n_e; K_i*dudx(2); 1e-300*dudx(3)];
s = [-k*n_e + n_e_dot;... %  
    g*(u(3)-u(2)) + xi*k*n_e*(C_e*u(3)+E_g);...
    D_e*C_e*log(10)*n_e*dudx(1)*dudx(3) - g*(u(3)-u(2)) + n_e_dot*(h*freq-E_g)];
% 
% if C_e*u(3) > h*freq
% C_e*u(3) - h*freq
% end

% t
% x
% u
% s(1)/c(1)
% s(2)/c(2)
% if (u(3) < 293.7)
% stuff = [t;x;C_i;C_e;log(D_e)/log(10);K_i;k;tau_fit;g;10^u(1);u(2);u(3);log(10)*n_e*dudx(1);dudx(2);dudx(3)]
% assignin('base','stuff',stuff);
% 
% s(3)/c(3) % suddenly spikes. why?
% end

end

function u0 = pdeic(x)
n_g = 2.4e13*(100^3)/(1e9)^3; % log_10 1/nm^3
u0 = [log(n_g)/log(10); 293.7; 293.7]; % 1/nm^3, K, K
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
n_g = log(2.4e13*(100^3)/(1e9)^3)/log(10); % log_10 1/nm^3
pl = [0; 0; 0]; % 
pr = [0; ur(2) - 293.7; 0]; %  ur(1) - n_g; ur(2) - 293.7; ur(3) - 293.7

ql = [1;1;1];
qr = [1;1;1];
end