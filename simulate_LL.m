function [y,mu] = simulate_LL(c,d,mu1,Ht,Qt,T) 
%% ADVANCED ECONOMETRICS
%
%  SIMULATE DATA FROM LINEAR PARAMETER-DRIVEN LOCAL-LEVEL MODEL
%
%  Description: 
%  This code snippet shows how to simulate data from a  
%  linear parameter-driven local-level (LL) model given by:
%
%  y(t) = c + mu(t) + epsilon(t)
%
%  mu(t+1) = d + mu(t) + eta(t)
%
%  where eps ~ NID(0,sigma_eps^2)  and  eta ~ NID(0,sigma_eta^2)
%
%% Parameter Values
%
%    d = 0;      % intercept parameter in update equation
%    c = 0;     % autoregressive parameter in update equation
%    sigma_eta = 15099;  % standard error of innovations in update equation
%                           Qt
%    sigma_eps = 1469.1;  % standard error of innovations in observation equation
%                            Ht
%    mu1 = mu1; % define initial value for time series x
%    T = T length of denerated observations

%%  Generate Innovations

    eta =  Qt * randn(T,1); % generate a vector of T random normal 
                                % variables with variance sigma_v^2

    epsilon = Ht * randn(T,1); % generate a vector of T random normal 
                                % variables with variance sigma_eps^2

%% 4. Define Time Series Vector

    mu = zeros(T,1); % define vector of zeros of length T
    y = zeros(T,1); % define vector of zeros of length T
    
%% 5. Define Initialization for Time Series

    mu(1) = mu1;

%% 6. Generate Time Series

    for t=1:T % start recursion from t=2 to t=T
       
       mu(t+1) = d + mu(t) + eta(t); % update equation
        
       y(t) = c + mu(t) +  epsilon(t); % observation equation
        
    end % end recursion
end








