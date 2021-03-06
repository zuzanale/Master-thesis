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
                     
%% 3. Initial Parameter Values to parameters that has to be estimated

    theta_ini = [0;0]; 
        
%% 4. Parameter Space Bounds
        lb = [10^(-7); 10^(-7)]; 
        ub = [inf; inf];
      
%% 5. Optimize Log Likelihood Criterion
     [theta_hat,llik,exitflag,output,lambda,grad,hessian]=...
          fmincon(@(theta) - llik_fun_new(data,theta),theta_ini,[],[],[],[],lb,ub,[],options);

display('parameter sigma2_eta')
theta_hat(1)
    
display('parameter omega')
theta_hat(2)
 
%% Parameter Initialisation
    y = data;
    sigma_eps = 15098.65; 
    sigma_eta = 1469.163;
    
    p = 10^7;
    %a = 0; 
    a =1120 ;
    %P0 = 10^7;
    P0 = 286379470;
    a0 = 0; 
    H = 1;
    mT=1;
                      
%% Auxiliary filter
% Regression part
dt = zeros(T,2); %scaled analogue of the distance

for tt = 2:T
    
    my_field = strcat('score_',num2str(tt));
    variable.(my_field) =zeros(T,2);
    my_field2 = strcat('score2_',num2str(tt));
    variable.(my_field2) =zeros(T,2);
    
    %tt = 43;
    shock = zeros(T,1);
    shock(tt) = 1; % put the shock on the first time index t=1
    Xt = shock; %dummy variables expressing a shock
    shock2 = zeros(T,1);
    shock2(tt:end) = 1;
    Xt2 = shock2;
    
    %KFS procedure applied
    [~,at,alpha_hat,score_lik,ut,ut_star,rt,rt_star,St,st] = kf_smooth_adj(y,H,1,0,0,sigma_eps,a0,P0,Xt); 
    [~,at2,alpha_hat2,score_lik2,ut2,ut_star2] = kf_smooth_adj(y,H,1,0,0,sigma_eps,a0,P0,Xt2); 
    %[llik,alphat,at,score_lik] = kf_smooth(yt,Ht,mT,c,d,Qt,a0,P0,Xt) 
    
    %hessian = [0.0101 -0.0263;-0.0263 0.188]; %taken from the paper(in R wrong results for some reason)
    %hessian = [1.677492e-07 1.717391e-07;1.717391e-07 1.686288e-06]; %version fro R using KFAS
    dt(tt,1) = (-hessian(1,1))^(-1) *  score_lik(1,2)/((-sqrt(hessian(1,1)))^(-1));
    dt(tt,2) = (-hessian(2,2))^(-1) *  score_lik(2,2)/((-sqrt(hessian(2,2)))^(-1));
    
    dt2(tt,1) = (-hessian(1,1))^(-1) * score_lik(1,3);
    dt2(tt,2) = (-hessian(2,2))^(-1) * score_lik(2,3);
  
    dt3(tt,1) =  (-hessian(1,1))^(-1) * score_lik2(1,3);
    dt3(tt,2) = (-hessian(2,2))^(-1) * score_lik2(2,3);
   
    slope_changes(tt,1) =  score_lik(2,3);
    shocks =  score_lik(1,3);
    
    variable.(my_field) = score_lik;
    variable.(my_field2)= score_lik2;
    %plot the comparison plots for ut and ut_star + rt and rt_star
    %figure(tt)
    %subplot(2,1,1)
    %plot(ut,'b')
    %hold on
   % plot(ut_star,'r')
    
   % subplot(2,1,2)
   % plot(rt,'b')
  %  hold on
  %  plot(rt_star,'r')
  
  %figure(1);subplot(2,1,1);stem(dt(:,2));subplot(2,1,2);stem(dt(:,1))
  %figure(2);subplot(2,1,1);stem(dt2(:,2));subplot(2,1,2);stem(dt2(:,1))
  %figure(3);subplot(2,1,1);stem(dt(:,1),dt(:,2));subplot(2,1,2);stem(dt(:,2),dt(:,1))
end

%%
%  SIMULATE DATA FROM LINEAR PARAMETER-DRIVEN LOCAL-LEVEL MODEL
%
%  Description: 
%  This code snippet shows how to simulate data from a  
%  linear parameter-driven local-level (LL) model given by:
%
%  x(t) = c + mu(t) + epsilon(t)
%
%  mu(t+1) = d + mu(t) + eta(t)
%
%  where eps ~ NID(0,sigma_eps^2)  and  eta ~ NID(0,sigma_eta^2)
%

T = 100; %length of simulated time series
c = 0;
d = 0;
mu1 = 1120; % define initial value for state observation mu
Ht = sigma_eps;
Qt = sigma_eta;

%simulate observations y^k
[x,mu] = simulate_LL(c,d,mu1,Ht,Qt,T); 

%% simulate an additive shocks
x2 = x; %create a copy of simulated time series
x2(15) = min(x2) - 0.3 * mean(x2);
x2(45) = max(x2) + 0.2 * mean(x2);
x2(70:90) = x2(70:90) + min(x2) + 0.15 * mean(x2) - x2(70:90) * 0.12;
  
%% apply our algorithm for outlier detecion
dt = zeros(T,2); %scaled analogue of the distance

for tt = 1:(T-10)
    
    my_field = strcat('score_',num2str(tt));
    variable.(my_field) =zeros(T,2);
    my_field2 = strcat('score2_',num2str(tt));
    variable.(my_field2) =zeros(T,2);
    
    %tt = 43;
    shock = zeros(T,1);
    shock(tt) = 10; % put the shock on the first time index t=1
    Xt = shock; %dummy variables expressing a shock
    shock2 = zeros(T,1);
    shock2(tt:(tt+10)) = shock2(tt:(tt+10)) + 10;
    Xt2 = shock2;
    
    %KFS procedure applied
    [~,at,alpha_hat,score_lik,ut,ut_star,rt,rt_star] = kf_smooth_adj(x2,Ht,1,0,0,Qt,a0,P0,Xt); 
    [~,at2,alpha_hat2,score_lik2,ut2,ut_star2] = kf_smooth_adj(x2,Ht,1,0,0,Qt,a0,P0,Xt2); 
    %[llik,alphat,at,score_lik] = kf_smooth(yt,Ht,mT,c,d,Qt,a0,P0,Xt) 
    
    %hessian = [0.0101 -0.0263;-0.0263 0.188]; %taken from the paper(in R wrong results for some reason)
    hessian = [2.097215  5.352072;5.352072 36.698336]; %version fro R using KFAS
    %dt(tt,1) = hessian(1,1) *  score_lik(1,2)/sqrt(hessian(1,1));
    %dt(tt,2) = hessian(2,2) *  score_lik(2,2)/sqrt(hessian(2,2));
    
    dt2(tt,1) = ((score_lik(1,3))/(hessian(1,1)))/(hessian(1,1)^(-0.5));
    dt2(tt,2) =((score_lik(2,3))/(hessian(2,2)))/(hessian(2,2)^(-0.5));
    
    dt3(tt,1) = ((score_lik2(1,3))/(hessian(1,1)))/(hessian(1,1)^(-0.5));
    dt3(tt,2) =((score_lik2(2,3))/(hessian(2,2)))/(hessian(2,2)^(-0.5));
   
    slope_changes(tt,1) =  score_lik(2,3);
    shocks =  score_lik(1,3);
    
    variable.(my_field) = score_lik;
    variable.(my_field2)= score_lik2;
    %plot the comparison plots for ut and ut_star + rt and rt_star
    %figure(tt)
    %subplot(2,1,1)
    %plot(ut,'b')
    %hold on
   % plot(ut_star,'r')
    
   % subplot(2,1,2)
   % plot(rt,'b')
  %  hold on
  %  plot(rt_star,'r')
  
  %figure(1);subplot(2,1,1);stem(dt(:,2));subplot(2,1,2);stem(dt(:,1))
  %figure(2);subplot(2,1,1);stem(dt2(:,2));subplot(2,1,2);stem(dt2(:,1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% envelope simulation approach
% Monte CaRLO SIMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 99; %number of simulations
T = 100; %length of simulated time series
j = 2; %number of different types of shocks
c = 0;
d = 0;
mu1 = 1120; % define initial value for state observation mu
Ht = sigma_eps;
Qt = sigma_eta;
hessian = [2.097215  5.352072;5.352072 36.698336]; %version fro R using KFAS
dtk = [zeros(T,j)]; % j is a number of different types of shocks

for k=1:K
    
    %simulate observations y^k
    [xk,mux] = simulate_LL(c,d,mu1,Ht,Qt,T); 
    
    %iterate for each shock at time t from T
    for tt = 1:T
        
        % create an additive shock at time tt
        shock = zeros(T,1);
        shock(tt) = 10; % put the shock on the first time index t=1
        Xt = shock; %dummy variables expressing a shock
    
         %KFS procedure applied
        [~,at,alpha_hat,score_lik,ut,ut_star,rt,rt_star] = kf_smooth_adj(xk,Ht,1,c,d,Qt,a0,P0,Xt); 
    
        %calculate delta ll - score
        dtk(tt,1,k) = sqrt(hessian(1,1)) * (score_lik(1,3))/hessian(1,1);
        dtk(tt,2,k) = sqrt(hessian(2,2)) * (score_lik(2,3))/hessian(2,2);
    end
end

%% produce confidence intervals/simulation envelopes
tint99_1=zeros(k,2); % t intervals
tint99_2=zeros(k,2); % t intervals
tint95_1=zeros(k,2); % t intervals
tint95_2=zeros(k,2); % t intervals

for i=1:k
    
    %% 99% interval
    %test dt for first type of shock
    [h,p,ci,stats] = ttest(dtk(:,1,i),0.01);
    tint99_1(i,:) = ci; %assign lower and upper bound to the variable
    
    %test dt for second type of shock 
     [h,p,ci,stats] = ttest(dtk(:,2,i),0.01);
    tint99_2(i,:) = ci; %assign lower and upper bound to the variable
    
     %% 95% interval
    %test dt for first type of shock
    [h,p,ci,stats] = ttest(dtk(:,1,i),0.05);
    tint95_1(i,:) = ci; %assign lower and upper bound to the variable
    
    %test dt for second type of shock 
     [h,p,ci,stats] = ttest(dtk(:,2,i),0.05);
    tint95_2(i,:) = ci; %assign lower and upper bound to the variable

end

figure(1)
stem(dt2(:,1))
hold on
plot(tint95_2,'b')

figure(2)
stem(dt2(:,1))
hold on
plot(tint99_2,'b')