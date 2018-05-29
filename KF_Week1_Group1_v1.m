%%  TIME SERIES ECONOMETRICS
%
%   ASSIGNMENT 1: LOCAL LEVEL MODEL AND KALMAN FILTER
%   Charlotte Taman, Femke Vedder, Rose Barzilai, Zuzana Leova (Group 1) 

%% 0. Clean Workspace and Command Window

clear all        %clear workspace
clc              %clear command window

%% 0. Read Data
A=importdata('Nile.dat');
A.data;
fid = fopen('Nile.dat','r');
datacell = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 1, 'Collect', 1);
fclose(fid);
A.data = datacell{1};
data=A.data(:,1);

%% 1. Setup

    T=size(data,1);  % sample size

%% MLE of parameters

%% Part C - QMLE
%% 2. Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',100000000); % maximum number of iterations    
        %options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
%% 3. Initial Parameter Values to parameters that has to be estimated

    theta_ini = [0;0]; 
        
%% 4. Parameter Space Bounds
        lb = [10^(-7); 10^(-7)]; 
        ub = [inf; inf];
      
%% 5. Optimize Log Likelihood Criterion

      [theta_hat,llik_val,exitflag,hessian]=...
          fmincon(@(theta) - llik_fun_new(data,theta),theta_ini,[],[],[],[],lb,ub,[],options);
    
display('parameter sigma2_eta')
theta_hat(1)
    
display('parameter omega')
theta_hat(2)
 
%% Parameter Initialisation
    y = data;
    sigma_eps = 15099; 
    sigma_eta = 1469.1;
    
    p = 10^7;
    a = 0; 
    P0 = 10^7;
    a0 = 0; 
    H = 1;
    mT=1;
                      
%% Apply Kalman Filter
     [~,at,alpha_hat] = kf_smooth(y,H,1,0,0,sigma_eps,a0,P0,0);
   

%% PART E: Replicate Slide 65
x =1871:1:1970;

figure(4)
plot(x(2:100),alpha_hat(2:100),'b')  % plot the time-series y in blue 'b'  
hold on
plot(x(2:100),y(2:100),'ko')
hold on
plot(x(2:100),at(2:100),'r')

%% Auxiliary filter
% Regression part
dt = zeros(T,2); %scaled analogue of the distance
for tt = 1:T
    %tt = 43;
    shock = zeros(T,1);
    shock(tt) = 10; % put the shock on the first time index t=1
    Xt = shock; %dummy variables expressing a shock
    shock2 = zeros(T,1);
    shock2(tt:end) = 10;
    Xt2 = shock2;
    
    %KFS procedure applied
    [~,at,alpha_hat,score_lik,ut,ut_star] = kf_smooth_adj(y,H,1,0,0,sigma_eps,a0,P0,Xt); 
    [~,at2,alpha_hat2,score_lik2,ut2,ut_star2] = kf_smooth_adj(y,H,1,0,0,sigma_eps,a0,P0,Xt2); 
    %[llik,alphat,at,score_lik] = kf_smooth(yt,Ht,mT,c,d,Qt,a0,P0,Xt) 
    
    hessian = [0.0101 -0.0263;-0.0263 0.188]; %taken from the paper(in R wrong results for some reason)
    dt(tt,1) = hessian(1,1) *  score_lik(1,2)/sqrt(hessian(1,1));
    
    %plot the comparison plots for ut and ut_star + rt and rt_star
    figure(tt)
    plot(ut,'b')
    hold on
    plot(ut_star,'r')
end

