% Phillips curves with different money growth rate means

clear all
clc

% parameterization
beta=0.98;
% time preference (discounting factor)

gamma=0.5;
% inverse of the FRISCH

M1=1;
% initial level of money supply (normalised to 1)

Z1=1;
% initial level of productivity(normalised to 1)

T=5000;
% number of periods

rho_mktpwr=0.9;
% 1/rho is the mark up term (Monopolistic competition)

init=T/10;
% starting point for the analysis

% set the parameters for the Markov process (tau and g are AR(1))
g_bar = 0.02;
rho_tau = 0.99;
rho_g = 0.95;
sigma_tau = 0.001;
sigma_g = 0.001;

tau_bar1 = 0.02;
tau_bar2=0.06;
tau_bar3=0.1;
tau_bar=[tau_bar1, tau_bar2, tau_bar3];


% draw random numbers from a std normal distibution to simulate shocks
eps_tau=randn(T,1);
eps_g=randn(T,1);


% create the results matrix
taus=zeros(T,3);
gs=zeros(T,1);
Money=zeros(T,3);
Z=zeros(T,1);
% assign to the result matrix the initial values
taus(1,:)=tau_bar;
gs(1,1)=g_bar;
Money(1,:)=M1;
Z(1,1)=Z1;

% draw shocks from the std normal distribution
eps_tau=randn(T,1);
eps_g=randn(T,1);

for t=1:3
    for i=2:T
        taus(i,t)=rho_tau*(taus(i-1,t)-tau_bar(t))+sigma_tau*eps_tau(i,1); 
        % deviation from the mean
   
        taus(i,t)=taus(i,t)+tau_bar(t); 
        % tau is computed as deviation from the mean + the mean
    
        taus(i,t)=max(taus(i,t),-0.99); 
        % 1+tau cannot be negative
    
        Money(i,t) = Money(i-1,t)*(1+taus(i,t)); 
        % money supply path
    end
end


for i = 2:T
    gs(i)=rho_g*(gs(i-1)-g_bar)+sigma_g*eps_g(i,1);
    % deviation from the mean
    
    gs(i)=gs(i)+g_bar; 
    % g is computed as deviation from the mean + the mean
    
    gs(i,1)=max(gs(i,1),-0.99); 
    % 1+g cannot be negative
    
    Z(i)= Z(i-1)*(1+gs(i));
    % productivity path
end


% compute the target of labor effort and the effective labor
% create two matrix to results in
L_target=zeros(T,3);
L_t=zeros(T,3);
% assign initial level
for t=1:3

    L_target(1,t)=(rho_mktpwr*beta/(1+tau_bar(t)+rho_tau^2*tau_bar(t)))^(1/(1+gamma));
L_t(1,t)= ((1+tau_bar(t)+rho_tau*tau_bar(t)+sigma_tau*eps_tau(1))/(1+tau_bar(t)+rho_tau*tau_bar(t)))*L_target(1,t);
end


for t=1:3
    for i=2:T
        L_target(i,t)=(rho_mktpwr*beta/(1+tau_bar(t)+rho_tau^2*(taus(i-1,t)-tau_bar(t))))^(1/(1+gamma));
        % target for labor effort
        
        L_t(i,t)=((1+tau_bar(t)+rho_tau*(taus(i-1,t)-tau_bar(t))+sigma_tau*eps_tau(i))/(1+tau_bar(t)+rho_tau*(taus(i-1,t)-tau_bar(t))))*L_target(i,t);    
        % realised labor effort
    end
end

figure(1)
scatter(L_t(init:T,1),taus(init:T,1))
title("Tau vs Labor (tau mean=0.02)")
lsline


figure(2)
scatter(L_t(init:T,2),taus(init:T,2))
title("Tau vs Labor (tau mean=0.06)")
lsline


figure(3)
scatter(L_t(init:T,3),taus(init:T,3))
title("Tau vs Labor (tau mean=0.10)")
lsline

figure(4)
scatter(L_t(init:T,1),taus(init:T,1))
hold on
scatter(L_t(init:T,2),taus(init:T,2))
scatter(L_t(init:T,3),taus(init:T,3))
legend("0.02","0.06","0.10")
hold off



figure(5)
scatter(L_t(init:T,1),eps_tau(init:T,1))
hold on
scatter(L_t(init:T,2),eps_tau(init:T))
scatter(L_t(init:T,3),eps_tau(init:T))
legend("0.02","0.06","0.10")
hold off
