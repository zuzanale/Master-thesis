function [llik,at,alphat,score_lik] = kf_smooth(yt,Ht,mT,c,d,Qt,a0,P0,Xt) 
%kf_smooth(y,H,1,0,0,sigma_eps,a0,P0,Xt); 

s0=0;
S0=0;
%% Kalman filter, matrix notation
vy = yt; %v stands for vector, m for a matrix
T = size(vy,1);
at = zeros(T,1);
Pt = zeros(T,1);
at_star = zeros(T,1);
Pt_star = zeros(T,1);
vt = zeros(T,1);
Kt = zeros(T,1);
Ft = zeros(T,1);
Et = zeros(T,1); % for auxiliary smoother recursion
Rt = zeros(T,1); % for auxiliary smoother recursion
score_lik = zeros(size(Ht,1)+size(Qt,1),2);
%mL = zeros(T,1);
mZ = 1;
mR = 1;


%% filtering
%Calculate all values of the kalman filter 
at(1,1) =  a0;
Pt(1,1) = P0;

for t = 1:T
    if t==T
        vt(t,1) = yt(t,1) - c - mZ * at(t,1);
    	Ft(t,1) = mZ * Pt(t,1) * mZ' + Ht;
        Kt(t,1) = mT * Pt(t,1) * mZ'/Ft(t,1);
    else
        vt(t,1) = yt(t,1) - c - mZ * at(t,1);
    	Ft(t,1) = mZ * Pt(t,1) * mZ' + Ht;
        Kt(t,1) = mT * Pt(t,1) * mZ'/Ft(t,1);
        at(t+1,1) = mT * at(t,1) + Kt(t,1)*vt(t,1) + d; 
        Pt(t+1,1) = mT * Pt(t,1)*mT' + mR*Qt*mR' - Kt(t,1)* Ft(t,1) * Kt(t,1)';
    end
end


  
%% loglikelihood evaluation due to (Q)MLE
l=  -(1/2)*T*log(2*pi) -(1/2)*sum(log(Ft) +((vt.^2)./Ft)); 
llik =mean(l);

%% Kalman smoothing
Lt = mT * ones(T,1) - Kt * mZ;
Ltrans = Lt';

rt = zeros(T,1);
rt_star = zeros(T,1);
alphat= zeros(T,1);
Nt = zeros(T,1);
Nt_star = zeros(T,1);
Vt = zeros(T,1);
ut = zeros(T,1);
Dt = zeros(T,1);
ut_star = zeros(T,1);
Dt_star = zeros(T,1);
epshat = zeros(T,1);
etahat = zeros(T,1);

r = 0;
N = 0;

%initialization
rt(T,1)=r;
Nt(T,1) = N;
alphat(T,1) = at(T,1) + Pt(T,1)*rt(T,1) + d;
Vt(T,1) = Pt(T,1) - Pt(T,1)*N*Pt(T,1) - c;

for t = T:-1:1
   if t==1
       %alphat(t,1) = at(t,1) + Pt(t,1)*rt(t-1,1) + d;
       %Vt(t,1) = Pt(t,1) - Pt(t,1)*Nt(t-1,1)*Pt(t,1) - c; 
       ut(t,1) = Ft(t,1)\vt(t,1) - Kt(t,1)'*rt(t,1);
       Dt(t,1) = inv(Ft(t,1)) + Kt(t,1)^2*Nt(t,1);
       epshat(t,1) = Ht * ut(t,1);
       etahat(t,1) = Qt * mR'* rt(t,1);
   else
        rt(t-1,1) = mZ'/Ft(t,1) * vt(t,1) + Ltrans(1,t) * rt(t,1);
        Nt(t-1,1) = mZ'/Ft(t,1) * mZ + Ltrans(1,t) * Nt(t,1)*Lt(t,1); 
    
        alphat(t,1) = at(t,1) + Pt(t,1)*rt(t-1,1) + d;
        Vt(t,1) = Pt(t,1) - Pt(t,1)*Nt(t-1,1)*Pt(t,1)  - c;
        ut(t,1) = Ft(t,1)\vt(t,1) - Kt(t,1)'*rt(t,1);
        Dt(t,1) = inv(Ft(t,1)) + Kt(t,1)^2*Nt(t,1);
        epshat(t,1) = Ht * ut(t,1);
        etahat(t,1) = Qt * mR'* rt(t,1);
   end
end    

%%  disturbance smoothing % this part should be rewritten 
%for t = T:-1:1
     %ut(t,1) = Ft(t,1)\vt(t,1) - Kt(t,1)'*rt(t,1);
     %Dt(t,1) = inv(Ft(t,1)) + Kt(t,1)^2*Nt(t,1);
     %epshat(t,1) = Ht * ut(t,1);
     %etahat(t,1) = Qt * mR'* rt(t,1);
     %varianceEpsilont = Ht*ones(size(Dt,1),1) - Ht^2*Dt;
     %varianceEtat = Qt*ones(size(Nt,1),1) - Qt^2*Nt;
%end

%% Score function calculation
score_lik(2,1)=1/2*sum(trace(ut*ut'-Dt) * Ht); %derivated by Ht
score_lik(1,1)=1/2*sum(trace(rt*rt'-Nt) * Qt); %derivated by Qt

%% Auxiliary filter
% initial values
st(1,1) = s0;
St(1,1) = S0;
bt = zeros(T,1);
At = zeros(T,1);
At(1,1) = mT * a0;
V_x = zeros(T,1);

%% correction of KF
for t = 1:T
    if t==1
        V_x(t,1) = -Xt(t) - c - mZ * At(t,1); 
        At(t+1,1) = mT * At(t,1) + Kt(t,1)*V_x(t,1) + d;
    else
        st(t,1) = st(t-1,1) + V_x(t,1)'/Ft(t,1) * vt(t,1);
        St(t,1) = St(t-1,1) + V_x(t,1)'/Ft(t,1) * V_x(t,1);
        At(t+1,1) = mT * At(t,1) + Kt(t,1)*V_x(t,1) + d;
        V_x(t,1) = -Xt(t,1) - c - mZ * At(t,1); 
        if st(t,1)==0
            bt(t,1)=0; 
        else
            bt(t,1) = St(t,1)\st(t,1); 
        end
    end
end
for t=1:T-1
    at_star(t+1,1) = at(t+1,1) + At(t+1,1) * bt(t,1);%t+1 index in the paper
    Pt_star(t+1,1) = Pt(t+1,1) - At(t+1,1) * St(t,1) * At(t+1,1)';%t+1 index in the paper
end
%% Auxiliary Smoother

%initialization
Rt(T,1)=0;

for t = T:-1:1
   if t==1
       Et(t,1)=Ft(t,1)\V_x(t,1) - Kt(t,1) * Rt(t,1);   
       ut_star(t,1) = ut(t,1) + Et(t,1) * bt(t,1); % bt(T,1) in the paper????
       Dt_star(t,1) = Dt(t,1) - Et(t,1)^2./St(t,1);
   else
       if St(t,1)==0
           Et(t,1)=Ft(t,1)\V_x(t,1) - Kt(t,1) * Rt(t,1);
           Rt(t-1,1)=mZ' * Et(t,1) + mT' * Rt(t,1);
           
           Dt_star(t,1) = Dt(t,1);
           Nt_star(t-1,1) = Nt(t-1,1);
           
           ut_star(t,1) = ut(t,1) + Et(t,1) * bt(t,1); % bt(T,1) in the paper????
           rt_star(t-1,1) = rt(t-1,1) + Rt(t-1,1) * bt(t,1);
       else
       
       Et(t,1)=Ft(t,1)\V_x(t,1) - Kt(t,1) * Rt(t,1);
       Rt(t-1,1)=mZ' * Et(t,1) + mT' * Rt(t,1);      
       
       %correction of disturbance smoother
       ut_star(t,1) = ut(t,1) + Et(t,1) * bt(t,1); % bt(T,1) in the paper????
       Dt_star(t,1) = Dt(t,1) - Et(t,1)^2./St(t,1);
       
       rt_star(t-1,1) = rt(t-1,1) + Rt(t-1,1) * bt(t,1);
       Nt_star(t-1,1) = Nt(t-1,1) - Rt(t-1,1)/St(t,1)* Rt(t-1,1)'; 
       end
   end
end    

 
%% Score function recalculation
score_lik(2,2)=1/2*sum(trace(ut_star*ut_star'-Dt_star) * Ht); %derivated by Ht
score_lik(1,2)=1/2*sum(trace(rt_star*rt_star'-Nt_star) * Qt); %derivated by Qt
end