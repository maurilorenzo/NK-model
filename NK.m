% NK.m
% New Keynesian model with Monopolistic competition, sticky prices and
% information frictions



clear all
clc

% parameterization
beta=0.98;
% time preferences

gamma=0.5;
% inverse of the FRISCH

M1=1;
% initial level of money supply (normalised to 1)

Z1=1;
% initial level of productivity(normalised to 1)

T=1000;
% number of periods

rho_mktpwr=0.9;
% 1/rho is the mark up term (Monopolistic competition)

init=T/10;
% starting point for the analysis

% set the parameters for the Markov process (tau and g are AR(1))
tau_bar = 0.02;
g_bar = 0.02;
rho_tau = 0.1;
rho_g = 0.95;
sigma_tau = 0.05;
sigma_g = 0.001;


% create the results matrix
taus=zeros(T,1000);
gs=zeros(T,1000);
Money=zeros(T,1000);
Z=zeros(T,1000);
L_target=zeros(T,1000);
L_t=zeros(T,1000);
Y_t=zeros(T,1000);
L_variation=zeros(T-1,1000);
Ygrowth_10yearsavg=zeros(T-10,1000);
tau_10yearsavg=zeros(T-9,1000);
P_t=zeros(T,1000);
corr_MvsP=zeros(1000,1);
corr_MvsL=zeros(1000,1);
corr_tauvsL=zeros(1000,1);
corr_tauvsLvariation=zeros(1000,1);
corr_tauvsYLR=zeros(1000,1);

% assign to the result matrix the initial values
taus(1,:)=0;
gs(1,:)=0;
Money(1,:)=M1;
L_target(1,:)=(rho_mktpwr*beta/(1+tau_bar+rho_tau^2*tau_bar))^(1/(1+gamma));

% MONTECARLO simulation
for s=1:1000
    % draw random numbers from a std normal distibution to simulate shocks
    eps_tau=randn(T,1);
    eps_g=randn(T,1);

    for i=2:T
        taus(i,s)=rho_tau*taus(i-1,s)+sigma_tau*eps_tau(i,1);
        % deviation from the mean
        
        gs(i,s)=rho_g*gs(i-1,s)+sigma_g*eps_g(i,1);
        % deviation from the mean
    
        gs(i,s)=gs(i,s)+g_bar;    
        % g is computed as deviation from the mean + the mean
    
        taus(i,s)=max(taus(i,s),-0.9);   
        % 1+tau cannot be negative   
         
        gs(i,s)=max(gs(i,s),-0.9);
        % 1+g cannot be negative 
    end
 
    % compute level of Money and productivity
    Z(1,s)=Z1*(1+gs(1,s)); % initial value of g
    for i = 2:T
        Z(i,s)= Z(i-1,s)*(1+gs(i,s));
        Money(i,s)= Money(i-1,s)*(1+tau_bar+taus(i-1,s));
    end
    
    % compute the target of labor effort and the effective labor
    % assign initial level to L_t
    L_t(1,:)=((1+tau_bar+rho_tau*tau_bar+sigma_tau*eps_tau(1))/(1+tau_bar+rho_tau*tau_bar))*L_target(1,s);
    for i=2:T
        L_target(i,s)=(rho_mktpwr*beta/(1+tau_bar+rho_tau^2*taus(i-1,s)))^(1/(1+gamma));
        L_t(i,s)=((1+tau_bar+rho_tau*taus(i-1,s)+sigma_tau*eps_tau(i))/(1+tau_bar+rho_tau*taus(i-1,s)))*L_target(i,s);
  end
    
    % compute production and price level
    Y_t(:,s)=Z(:,s).*L_t(:,s); 
    % production=income=consumption=real balance 
    P_t(:,s)=Money(:,s).*(1+tau_bar+taus(:,s))./Y_t(:,s); 
    % price level assuming CIA holds (transfers are received before the goods market)

    % compute LR averages (10 years rolling window)
    for i=1:T-10
       Ygrowth_10yearsavg(i,s)=10^-1*(log(Y_t(i+10))-log(Y_t(i)));
    end
    
    for i=1:T-9
       tau_10yearsavg(i,s)=sum(taus(i:i+9,s))/10;
    end
    

    % compute labor percentage variation
    for i=1:T-1
        L_variation(i,s)=L_t(i+1,s)/L_t(i,s)-1;
    end
    
    % compute correlations
    MvsP=corrcoef(Money(init:T,s),L_target(init:T));
    corr_MvsP(s)=MvsP(2);
    MvsL=corrcoef(Money(init:T,s),L_t(init:T,s));
    corr_MvsL(s)=MvsL(2);
    tauvsL=corrcoef(taus(init:T,s)+tau_bar,L_t(init:T,s));
    corr_tauvsL(s)=tauvsL(2);
    tauvsLvariation=corrcoef(taus(init+1:T,s)+tau_bar,L_variation(init:T-1,s));
    corr_tauvsLvariation(s)=tauvsLvariation(2);
    tauvsY_LR=corrcoef(tau_10yearsavg(init:T-10,s)+tau_bar, Ygrowth_10yearsavg(init:T-10,s));
    corr_tauvsYLR(s)=tauvsY_LR(2);
    
end


% compute mean of the correlations
meanMvsP=mean(corr_MvsP)
meanMvsL=mean(corr_MvsL)
meantauvsL=mean(corr_tauvsL)
meanstauvsL_variation=mean(corr_tauvsLvariation)
meantauvsY_LR=mean(corr_tauvsYLR)
% plot the histograms of the corrlations
figure(1)
hist(corr_MvsP)
title("Money vs Price")

figure(2)
hist(corr_MvsL)
title("Money vs Labor")

figure(3)
hist(corr_tauvsL)
title("Tau vs Labor")

figure(4)
hist(corr_tauvsLvariation)
title("Tau vs percentage variation of Labor")

figure(5)
hist(corr_tauvsYLR)
title("Tau vs Output (LR averages)")