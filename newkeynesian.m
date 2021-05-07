% New keynesian model
% impusive response






eps_tau=randn(100,1);
tau_bar=;

tau=[];
tau=[tau sigma_tau*eps_tau(1)];
for i=2:T
    tau=[tau tho_tau^(i-1)*sigma_tau*eps_tau(i)];
    L_taget=(beta*mkt_pwr)^(1/(1+gamma))*1./(1+tau_bar+rho_tau*tau(i-1)).^(1/(1+gamma));
end
tau=tau+tau_bar;


L_t=(1+sigma_tau*eps_tau(i)):