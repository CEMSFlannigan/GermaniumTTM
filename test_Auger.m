
[t,y] = ode15s(@vdp1000,[1e-9 1],1e20);
plot(t,log(y(:,1))/log(10),'-o');


function dydt = vdp1000(t,y)

n_i = 2.4e13;
n_auger_piecewise = 1e18;

if (y > n_auger_piecewise)
    C_eeh = (10^190.9124*y^(0.6698*log(y)/log(10) - 21.5429 - 1))/(y^2-n_i^2); % cm^6/s
    C_hhe = (10^120.1814*y^(0.6698*log(y)/log(10) - 21.5429 - 1))/(y^2-n_i^2); % cm^6/s
    k = (((C_eeh + C_hhe))*y^2); % function of n_e, Dominici et. al., results in 1/ps units
elseif (y <= n_auger_piecewise && y >= n_i)
    C_eeh = (10^(-30.9373)*y^(0.0110*log(y)/log(10) + 2.6402 - 1))/(y^2-n_i^2);
    C_hhe = (10^(-28.9705)*y^(0.0114*log(y)/log(10) + 2.6257 - 1))/(y^2-n_i^2);
    k = (((C_eeh + C_hhe))*y^2);
else
    k = 0;
end

ndot = 1e22*abs(sin((1/1e-5)*t));

dydt = -k*y;

end