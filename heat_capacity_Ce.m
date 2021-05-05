JpeV = 1.60218e-19; % J/eV

E_c = 0.66*JpeV; % J
E_v = 0*JpeV; % J
E_g = E_c-E_v;

m_o = 9.10938356e-31; % kg
m_h = 0.34*m_o; % kg
m_e = 0.22*m_o; % kg
k_B = 8.617333262145e-5*JpeV; % J/K
hbar = 4.135667696e-15/2/pi*JpeV; % J*s

N_V = @(Te) 2.*(m_h.*k_B.*Te./2./pi./hbar.^2).^(3/2);
N_C = @(Te) 2.*(m_e.*k_B.*Te./2./pi./hbar.^2).^(3/2);
chem_pot = @(Te) (E_c + E_v)/2 + k_B.*Te./2.*log(N_V(Te)./N_C(Te)); % J
fermi_dist = @(E,Te) (exp((E-chem_pot(Te))./k_B./Te) + 1).^(-1); % unitless pop/pop
dens_states = @(E,Te) (E > E_c).*m_e./hbar.^3./pi.^2.*sqrt((2.*m_e.*(E-E_c))) + N_C(Te); % 1/J/m^3 ||| 
dFdTe = @(E,Te) ((E-chem_pot(Te)).*exp((E-chem_pot(Te))./k_B./Te))./(k_B.*Te.^2.*(exp((E-chem_pot(Te))./k_B./Te) + 1).^2); % 1/K
d2FdTe2 = @(E,Te) -((E-chem_pot(Te))^2*exp((E-chem_pot(Te))/(k_B.*Te)))/(k_B^2*Te.^4*(exp((E-chem_pot(Te))/(k_B.*Te)) + 1)^2) + (2*(E-chem_pot(Te))^2*exp(2*(E-chem_pot(Te))/(k_B.*Te)))/(k_B^2*Te.^4*(exp((E-chem_pot(Te))/(k_B.*Te)) + 1)^3) - (2*(E-chem_pot(Te))*exp((E-chem_pot(Te))/(k_B.*Te)))/(k_B*Te.^3*(exp((E-chem_pot(Te))/(k_B.*Te)) + 1)^2);

C_e = @(Te) integral(@(E)dFdTe(E,Te).*dens_states(E,Te).*E,-1*JpeV,5*JpeV); % J/K/m^3

tempRange = 78:20000;
energyRange = (0.66:0.01:2)*JpeV; % J
C_e_res = zeros(1,length(tempRange));
dens_states_res = zeros(length(tempRange),length(energyRange));

for i = 1:length(energyRange)
    for j = 1:length(tempRange)
        dens_states_res(j,i) = dens_states(energyRange(i),tempRange(j));
    end
end

for i = 1:length(tempRange)
    C_e_res(i) = C_e(tempRange(i));
end

log_temp_range = log(tempRange)/log(10);
log_C_e_res = log(C_e_res)/log(10);

C_e_res_d1 = diff(C_e_res)./diff(tempRange);
C_e_res_d2 = diff(C_e_res_d1)./diff(tempRange(1:end-1))./diff(tempRange(1:end-1));

logtempCe = (log_temp_range - 3.876)/0.4073;
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
C_e_res_fit = p1Ce*tempCe.^9 + p2Ce*tempCe.^8 + p3Ce*tempCe.^7 + p4Ce*tempCe.^6 + p5Ce*tempCe.^5 + p6Ce*tempCe.^4 + p7Ce*tempCe.^3 + p8Ce*tempCe.^2 + p9Ce*tempCe.^1 + p10Ce;

figure;
plot(tempRange,C_e_res/JpeV);
% hold on;
% plot(tempRange,C_e_res_fit/JpeV);

h = (4.135667696e-15)*1e12; % eV*ps
freq = 2.9979e8/515e-9/1e12; % 1/ps
n_g = 2.4e13*(100^3)/(1e9)^3; % 1/nm^3

figure;
plot(tempRange(2:end), h*freq./(C_e_res(2:end)/JpeV/(1e9)^3/n_g));