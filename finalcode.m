clear all
load('nile.mat')
%% Kalman filter
% initial values
a = 0;
p = 10^7;
Ht = 15099;
Qt = 1469;
yt = data;
T = size(yt,1);

%Calculate all values of the kalman filter 
at(1,1) =  a;
Pt(1,1) = p;
vt(1,1) = yt(1,1) - at(1,1);
Ft(1,1) = Pt(1,1) + Ht;
Kt(1,1) = Pt(1,1) / Ft(1,1);

for t = 2:T
   at(t,1) = at(t-1,1) + Kt(t-1,1)*vt(t-1,1);  
   Pt(t,1) = Kt(t-1,1)*Ht + Qt;
   vt(t,1) = yt(t,1) - at(t,1);
   Ft(t,1) = Pt(t,1) + Ht;
   Kt(t,1) = Pt(t,1)/Ft(t,1);
end


confidence90Upper = at + sqrt(Pt)*norminv([0.9],0,1);
confidence90Lower = at - sqrt(Pt)*norminv([0.9],0,1);

%% plot book Fig. 2.1
figure('NumberTitle', 'off', 'Name', 'Figure 2.1');
t = (1872:1970)';
subplot(2,2,1)
q1 = plot(datenum(t,1,1),at(2:100));
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),yt(2:100));
set(q2,'Color','red','LineWidth',1.2);
confidence90Upperpper = plot(datenum(t,1,1),confidence90Upper(2:100),':');
confidence90Lower = plot(datenum(t,1,1),confidence90Lower(2:100),':');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 450 1400]);
legend('at','Data','90% cfi','Location','northeast')
title('(i)');

hold off

subplot(2,2,2)
t = (1871:1970)';
plot(datenum(t,1,1),Pt);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 5000 17500]);
title('(ii)');


subplot(2,2,3)
t = (1872:1970)';
z = zeros(99,1);
plot(datenum(t,1,1),vt(2:100));
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -450 400]);
title('(iii)');
hold off

subplot(2,2,4)
t = (1871:1970)';
plot(datenum(t,1,1),Ft);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 20000 32500]);
title('(iv)');

%% paragraph 2.4 state smoothing
Lt = ones(T,1) - Kt;
Ltrans = Lt';
r = 0;
N = 0;

rt(T,1)=r;
Nt(T,1) = N;

alphat(T,1) = at(T,1) + Pt(T,1)*rt(T,1);
Vt(T,1) = Pt(T,1) - Pt(T,1)*N*Pt(T,1);


for t = 1:T-1 
   rt(T-t,1) = inv(Ft(T-t+1,1))*vt(T-t+1,1) + Lt(T-t+1,1)*rt(T-t+1);
   Nt(T-t,1) = inv(Ft(T-t+1,1)) + Ltrans(1,T-t+1)*Nt(T-t+1,1)*Lt(T-t+1,1); 
    
   alphat(T-t,1) = at(T-t,1) + Pt(T-t,1)*rt(T-t,1);
   Vt(T-t,1) = Pt(T-t+1,1) - Pt(T-t+1,1)*Nt(T-t,1)*Pt(T-t+1,1);
   
   
end    
confidence90Uppers = alphat + sqrt(Vt)*norminv([0.9],0,1);
confidence90Lowers = alphat - sqrt(Vt)*norminv([0.9],0,1);


 %% Book Fig. 2.2
figure('NumberTitle', 'off', 'Name', 'Figure 2.2');
t = (1871:1970)';
subplot(2,2,1)
q1 = plot(datenum(t,1,1),alphat);
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),yt);
set(q2,'Color','red','LineWidth',1.2);
title('(i)');
plot(datenum(t,1,1),confidence90Uppers,':');
plot(datenum(t,1,1),confidence90Lowers,':');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 450 1400]);
legend('alpha_t','Data','90% cfi','Location','northeast')
title('(i)');
hold off

subplot(2,2,2)
t = (1871:1970)';
plot(datenum(t,1,1),Vt);
title('(ii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1971,1,1) 2300 4100]);
title('(ii)');


subplot(2,2,3)
z = zeros(101,1);
t = (1870:1970)';
plot(datenum(t,1,1),[r;rt]);
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1971,1,1) -0.036 0.024]);
title('(iii)');

hold off

subplot(2,2,4)
t = (1870:1970)';
plot(datenum(t,1,1),[N;Nt]);
title('(iv)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1971,1,1) 0.000025 0.00011]);
title('(iv)');


%% paragraph 2.5 disturbance smoothing
for t = 1:T
    ut(t,1) = inv(Ft(t,1))*vt(t,1) - Kt(t,1)*rt(t,1);
    Dt(t,1) = inv(Ft(t,1)) + Kt(t,1)^2*Nt(t,1);
    Finverse(t,1) = inv(Ft(t,1));
end

epsilont = Ht*ut;
varianceEpsilont = Ht*ones(size(Dt,1),1) - Ht^2*Dt;
etat = Qt*rt;
varianceEtat = Qt*ones(size(Nt,1),1) - Qt^2*Nt;

%% Figure 2.3
figure('NumberTitle', 'off', 'Name', 'Figure 2.3');
t = (1871:1970)';
subplot(2,2,1)
z = zeros(100,1);
plot(datenum(t,1,1),epsilont);
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -360 280]);
title('(i)');
hold off

subplot(2,2,2)
t = (1871:1970)';
plot(datenum(t,1,1),varianceEpsilont);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 2300 4100]);
title('(ii)');


subplot(2,2,3)
z = zeros(100,1);
t = (1871:1970)';
plot(datenum(t,1,1),etat);
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -50 35]);
title('(iii)');
hold off

subplot(2,2,4)
t = (1871:1970)';
plot(datenum(t,1,1),varianceEtat);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 1220 1475]);
title('(iv)');


%% paragraph 2.7 missing observations
for t=1:T
      if  (t>=0 && t<=20)
          dataM(t,1) = data(t,1);
      elseif (t>=41 && t<=60)    
          dataM(t,1) = data(t,1);
      elseif (t>=81)    
          dataM(t,1) = data(t,1);
      else    
      end    
     if (t>=21 && t<=40)
         dataN(t,1) = NaN;
     elseif (t>=61 && t<=80)
         dataN(t,1) = NaN;   
     else
         dataN(t,1) = data(t,1);
     end    
end 

a = 0;
p = 10^7;
for t = 1:size(data,1)
   if dataM(t,1) == 0
         muMt(t,1) = a;
         pMt(t,1) = p;
         p = p + Qt; 
   else    
   muMt(t,1) = a;
   pMt(t,1) = p;
   v = dataM(t,1) - a;
   k = p / (p + Ht);
   a = a + k*v;
   p = k*Ht + Qt;
   KMt(t,1) = k;
   vMt(t,1) = v;
   end
end

r = 0;
N = 0;
LMt = ones(size(KMt,1),1) - KMt;
LMtrans = LMt';
FMt = pMt + Ht;
for t = size(data,1):-1:1
   if dataM(t,1) == 0       
   else
        r = inv(FMt(t,1))*vMt(t,1) + LMt(t,1)*r;
        N = inv(FMt(t,1)) + LMtrans(1,t)*N*LMt(t,1);
   end
        alphaMt(t,1) = muMt(t,1) + pMt(t,1)*r;
        VMt(t,1) = pMt(t,1) - pMt(t,1)*N*pMt(t,1);   
end 

%% Figure 2.5
figure('NumberTitle', 'off', 'Name', 'Figure 2.5');
t = (1872:1970)';
subplot(2,2,1)
q1 = plot(datenum(t,1,1),muMt(2:100));
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),dataN(2:100));
set(q2,'Color','red','LineWidth',1.2);
title('(i)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 450 1400]);
hold off

subplot(2,2,2)
t = (1872:1970)';
plot(datenum(t,1,1),pMt(2:100));
title('(ii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 4000 36000]);

subplot(2,2,3)
t = (1871:1970)';
q1 = plot(datenum(t,1,1),alphaMt);
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),dataN);
set(q2,'Color','red','LineWidth',1.2);
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 450 1400]);
hold off

subplot(2,2,4)
t = (1871:1970)';
plot(datenum(t,1,1),VMt);
title('(iv)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 2000 10000]);
%% paragraph 2.8 forecasting
atF = at;
atF2 = [];
PtF = Pt;
FtF = [];
for t = 1:30
    atF = [atF;atF(end,1)];
    atF2 = [atF2;atF(end,1)];
    FtF = [FtF;PtF(end,1)+Ht];
    PtF = [PtF; PtF(end,1)+Qt];
end    
confUF = atF2 + sqrt(FtF)*norminv([0.75],0,1);
confDF = atF2 - sqrt(FtF)*norminv([0.75],0,1);

FtF = [PtF(1:100)+Ht;FtF];

%% Figure 2.6
figure('NumberTitle', 'off', 'Name', 'Figure 2.6');
t = (1872:2000)';
subplot(2,2,1)
plot(datenum(t,1,1),atF(2:end,1));
hold on
plot(datenum(t,1,1),[yt(2:end,1);NaN*ones(30,1)],'.');
title('(i)');
plot(datenum(t,1,1),[NaN*ones(99,1);confUF],':');
plot(datenum(t,1,1),[NaN*ones(99,1);confDF],':');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 450 1400]);
hold off

subplot(2,2,2)
t = (1872:2000)';
plot(datenum(t,1,1),PtF(2:end));
title('(ii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 4000 80000]);

subplot(2,2,3)
t = (1872:2000)';
plot(datenum(t,1,1),atF(2:end));
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 700 1200]);

subplot(2,2,4)
t = (1872:2000)';
plot(datenum(t,1,1),FtF(2:end));
title('(iv)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 19800 100000]);


%% paragraph 2.12 Diagnostic checking forecast errors
et = vt./sqrt(Ft);

%%  Figure 2.7
figure('NumberTitle', 'off', 'Name', 'Figure 2.7');
t = (1871:1970)';
z = zeros(100,1);
subplot(2,2,1)
plot(datenum(t,1,1),et);
hold on
plot(datenum(t,1,1),z,'k');
title('(i)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -2.8 2.8]);
hold off

subplot(2,2,2)
[f,x]=hist(et);
bar(x,f/trapz(x,f));
hold on
[ff,xi] = ksdensity(et);
plot(xi,ff);
title('(ii)');
axis([-4 4 0 0.55]);
hold off

subplot(2,2,3)
h=qqplot(et);
h(1).LineStyle = '-';
h(2).LineStyle = '-';
h(3).LineStyle = '-';
h(1).Color =[0,0,0];
title('(iii)');

subplot(2,2,4)
corr=autocorr(et);
bar(corr(2:end));
xlim([0 15]);
ylim([-1 1])
title('(iv)');

%% 2.12.2 Detection for ourliers and structural breaks
ssrut = Dt.^(-0.5).*ut;
ssrrt = Nt.^(-0.5).*rt;

%% Figure 2.8
figure('NumberTitle', 'off', 'Name', 'Figure 2.8');
subplot(2,2,1)
t = (1871:1970)';
plot(datenum(t,1,1),ssrut);
hold on
title('(i)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -3 2.5]);
hold off

subplot(2,2,2)
[f,x]=hist(ssrut,13);
bar(x,f/trapz(x,f));
hold on
[ff,xi] = ksdensity(ssrut);
plot(xi,ff);
title('(ii)');
axis([-4 3.5 0 0.55]);
hold off

subplot(2,2,3)
t = (1871:1970)';
z = zeros(100,1);
zu = z + 2;
zd = z - 2;
plot(datenum(t,1,1),ssrrt);
hold on
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -3 2.5]);
hold off

subplot(2,2,4)
[f,x]=hist(ssrrt,13);
bar(x,f/trapz(x,f));
hold on
[ff,xi] = ksdensity(ssrrt(1:98,1));
plot(xi,ff);
title('(iv)');
axis([-4.5 3 0 0.55]);
hold off

%% hraphs for AR(1) plus noise

%% Kalman filter
% initial values
a = 0;
phi=0.99998;%yt\at;
Ht = 15099;
Qt = 1469;
p =Qt/(1-phi^2);
yt = data;
T = size(yt,1);

%Calculate all values of the kalman filter 
at(1,1) =  a;
Pt(1,1) = p;
vt(1,1) = yt(1,1) - at(1,1);
Ft(1,1) = Pt(1,1) + Ht;
Kt(1,1) = phi.*Pt(1,1) / Ft(1,1);

for t = 2:T
   at(t,1) = phi*at(t-1,1) + Kt(t-1,1)*vt(t-1,1);  
   Pt(t,1) = phi^2*Kt(t-1,1)*Ht + Qt;
   vt(t,1) = yt(t,1) - at(t,1);
   Ft(t,1) = Pt(t,1) + Ht;
   Kt(t,1) = phi*Pt(t,1)/Ft(t,1);
end

confidence90Upper = at + sqrt(Pt)*norminv([0.9],0,1);
confidence90Lower = at - sqrt(Pt)*norminv([0.9],0,1);


%% plot book Fig. 2.1
figure('NumberTitle', 'off', 'Name', 'Figure 2.1');
t = (1872:1970)';
subplot(2,2,1)
q1 = plot(datenum(t,1,1),at(2:100));
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),yt(2:100));
set(q2,'Color','red','LineWidth',1.2);
confidence90Upperpper = plot(datenum(t,1,1),confidence90Upper(2:100),':');
confidence90Lower = plot(datenum(t,1,1),confidence90Lower(2:100),':');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 400 1400]);
legend('at','Data','90% cfi','Location','northeast')
title('(i)');

hold off

subplot(2,2,2)
t = (1871:1970)';
plot(datenum(t,1,1),Pt);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 4500 17500]);
title('(ii)');


subplot(2,2,3)
t = (1872:1970)';
z = zeros(99,1);
plot(datenum(t,1,1),vt(2:100));
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -450 450]);
title('(iii)');
hold off

subplot(2,2,4)
t = (1871:1970)';
plot(datenum(t,1,1),Ft);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 19000 32500]);
title('(iv)');

%% smoothing
Lt = phi*ones(T,1) - Kt;
Ltrans = Lt';
r = 0;
N = 0;

rt(T,1)=r;
Nt(T,1) = N;

alphat(T,1) = at(T,1) + Pt(T,1)*rt(T,1);
Vt(T,1) = Pt(T,1) - Pt(T,1)*N*Pt(T,1);


for t = 1:T-1 
   rt(T-t,1) = inv(Ft(T-t+1,1))*vt(T-t+1,1) + Lt(T-t+1,1)*rt(T-t+1);
   Nt(T-t,1) = inv(Ft(T-t+1,1)) + Ltrans(1,T-t+1)*Nt(T-t+1,1)*Lt(T-t+1,1); 
    
   alphat(T-t,1) = at(T-t,1) + Pt(T-t,1)*rt(T-t,1);
   Vt(T-t,1) = Pt(T-t+1,1) - Pt(T-t+1,1)*Nt(T-t,1)*Pt(T-t+1,1);
   
   
end    
confidence90Uppers = alphat + sqrt(Vt)*norminv([0.9],0,1);
confidence90Lowers = alphat - sqrt(Vt)*norminv([0.9],0,1);



 %% BooK Fig. 2.2
figure('NumberTitle', 'off', 'Name', 'Figure 2.2');
t = (1871:1970)';
subplot(2,2,1)
q1 = plot(datenum(t,1,1),alphat);
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),yt);
set(q2,'Color','red','LineWidth',1.2);
title('(i)');
plot(datenum(t,1,1),confidence90Uppers,':');
plot(datenum(t,1,1),confidence90Lowers,':');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 400 1400]);
legend('alpha_t','Data','90% cfi','Location','northeast')
title('(i)');
hold off

subplot(2,2,2)
t = (1871:1970)';
plot(datenum(t,1,1),Vt);
title('(ii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1971,1,1) 2100 4100]);
title('(ii)');

subplot(2,2,3)
z = zeros(101,1);
t = (1870:1970)';
plot(datenum(t,1,1),[r;rt]);
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1971,1,1) -0.036 0.04]);
title('(iii)');

hold off
subplot(2,2,4)
t = (1870:1970)';
plot(datenum(t,1,1),[N;Nt]);
title('(iv)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1971,1,1) 0.000025 0.00011]);
title('(iv)');


%% disturbance smoothing
for t = 1:T
    ut(t,1) = inv(Ft(t,1))*vt(t,1) - Kt(t,1)*rt(t,1);
    Dt(t,1) = inv(Ft(t,1)) + Kt(t,1)^2*Nt(t,1);
    Finverse(t,1) = inv(Ft(t,1));
end

epsilont = Ht*ut;
varianceEpsilont = Ht*ones(size(Dt,1),1) - Ht^2*Dt;
etat = Qt*rt;
varianceEtat = Qt*ones(size(Nt,1),1) - Qt^2*Nt;

%% Figure 2.3
figure('NumberTitle', 'off', 'Name', 'Figure 2.3');
t = (1871:1970)';
subplot(2,2,1)
z = zeros(100,1);
plot(datenum(t,1,1),epsilont);
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -360 280]);
title('(i)');
hold off

subplot(2,2,2)
t = (1871:1970)';
plot(datenum(t,1,1),varianceEpsilont);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 2200 4100]);
title('(ii)');

subplot(2,2,3)
z = zeros(100,1);
t = (1871:1970)';
plot(datenum(t,1,1),etat);
hold on
q1 = plot(datenum(t,1,1),z,':');
set(q1,'LineWidth',0.01);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -50 50]);
title('(iii)');
hold off

subplot(2,2,4)
t = (1871:1970)';
plot(datenum(t,1,1),varianceEtat);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 1220 1475]);
title('(iv)');


%% missing observations
for t=1:T
      if  (t>=0 && t<=20)
          dataM(t,1) = data(t,1);
      elseif (t>=41 && t<=60)    
          dataM(t,1) = data(t,1);
      elseif (t>=81)    
          dataM(t,1) = data(t,1);
      else    
      end    
     if (t>=21 && t<=40)
         dataN(t,1) = NaN;
     elseif (t>=61 && t<=80)
         dataN(t,1) = NaN;   
     else
         dataN(t,1) = data(t,1);
     end    
end 

a = 0;
p =Qt/(1-phi^2);
for t = 1:size(data,1)
   if dataM(t,1) == 0
         muMt(t,1) = phi*a;
         pMt(t,1) = p;
         p = phi^2*p + Qt; 
   else    
   muMt(t,1) = a;
   pMt(t,1) = p;
   v = dataM(t,1) - a;
   k = phi*p / (p + Ht);
   a = phi*a + k*v;
   p = phi^2*k*Ht + Qt;
   KMt(t,1) = k;
   vMt(t,1) = v;
   end
end

r = 0;
N = 0;
LMt = phi*ones(size(KMt,1),1) - KMt;
LMtrans = LMt';
FMt = pMt +Ht;
for t = size(data,1):-1:1
   if dataM(t,1) == 0       
   else
        r = inv(FMt(t,1))*vMt(t,1) + LMt(t,1)*r;
        N = inv(FMt(t,1)) + LMtrans(1,t)*N*LMt(t,1);
   end
        alfaMt(t,1) = muMt(t,1) + pMt(t,1)*r;
        VMt(t,1) = pMt(t,1) - pMt(t,1)*N*pMt(t,1);   
end 

% assignment week 2 - BooK Fig. 2.5
figure(5);
t = (1872:1970)';
subplot(2,2,1)
q1 = plot(datenum(t,1,1),muMt(2:100));
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),dataN(2:100),':');
set(q2,'Color','red','LineWidth',1.2);
title('(i)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 400 1400]);
hold off

subplot(2,2,2)
t = (1872:1970)';
plot(datenum(t,1,1),pMt(2:100));
title('(ii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 4000 36000]);

subplot(2,2,3)
t = (1871:1970)';
q1 = plot(datenum(t,1,1),alfaMt);
set(q1,'LineWidth',1.2);
hold on
q2 = plot(datenum(t,1,1),dataN,':');
set(q2,'Color','red','LineWidth',1.2);
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 450 1400]);
hold off

subplot(2,2,4)
t = (1871:1970)';
plot(datenum(t,1,1),VMt);
title('(iv)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 2000 10000]);

%% forecasting
atF = at;
atF2 = [];
PtF = Pt;
FtF = [];
for t = 1:30
    atF = [atF;atF(end,1)];
    atF2 = [atF2;atF(end,1)];
    FtF = [FtF;PtF(end,1)+Ht];
    PtF = [PtF; PtF(end,1)+Qt];
end    
confUF = atF2 + sqrt(FtF)*norminv([0.75],0,1);
confDF = atF2 - sqrt(FtF)*norminv([0.75],0,1);

FtF = [PtF(1:100)+Ht;FtF];

%% Figure 2.6
figure('NumberTitle', 'off', 'Name', 'Figure 2.6');
t = (1872:2000)';
subplot(2,2,1)
plot(datenum(t,1,1),atF(2:end,1));
hold on
plot(datenum(t,1,1),[yt(2:end,1);NaN*ones(30,1)],'.');
title('(i)');
plot(datenum(t,1,1),[NaN*ones(99,1);confUF],':');
plot(datenum(t,1,1),[NaN*ones(99,1);confDF],':');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 450 1400]);
hold off

subplot(2,2,2)
t = (1872:2000)';
plot(datenum(t,1,1),PtF(2:end));
title('(ii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 4000 80000]);

subplot(2,2,3)
t = (1872:2000)';
plot(datenum(t,1,1),atF(2:end));
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 600 1200]);

subplot(2,2,4)
t = (1872:2000)';
plot(datenum(t,1,1),FtF(2:end));
title('(iv)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(2000,1,1) 19800 100000]);


%% Diagnostic checking for forecast errors
et = vt./sqrt(Ft);

%%  Figure 2.7
figure('NumberTitle', 'off', 'Name', 'Figure 2.7');
t = (1871:1970)';
z = zeros(100,1);
subplot(2,2,1)
plot(datenum(t,1,1),et);
hold on
plot(datenum(t,1,1),z,'k');
title('(i)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -2.8 3.2]);
hold off

subplot(2,2,2)
[f,x]=hist(et);
bar(x,f/trapz(x,f));
hold on
[ff,xi] = ksdensity(et);
plot(xi,ff);
title('(ii)');
axis([-4 4 0 0.55]);
hold off

subplot(2,2,3)
h=qqplot(et);
h(1).LineStyle = '-';
h(2).LineStyle = '-';
h(3).LineStyle = '-';
h(1).Color =[0,0,0];
title('(iii)');
ylim([-2 2]);

subplot(2,2,4)
corr=autocorr(et);
bar(corr(2:end));
xlim([0 15]);
ylim([-1 1]);
title('(iv)');

%% Detection for ourliers and structural breaks
ssrut = Dt.^(-0.5).*ut;
ssrrt = Nt.^(-0.5).*rt;

%% Figure 2.8
figure('NumberTitle', 'off', 'Name', 'Figure 2.8');
subplot(2,2,1)
t = (1871:1970)';
plot(datenum(t,1,1),ssrut);
hold on
title('(i)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -3 2.5]);
hold off

subplot(2,2,2)
[f,x]=hist(ssrut,13);
bar(x,f/trapz(x,f));
hold on
[ff,xi] = ksdensity(ssrut);
plot(xi,ff);
title('(ii)');
axis([-4 3.5 0 0.55]);
hold off

subplot(2,2,3)
t = (1871:1970)';
z = zeros(100,1);
zu = z + 2;
zd = z - 2;
plot(datenum(t,1,1),ssrrt);
hold on
title('(iii)');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) -3 4]);
hold off

subplot(2,2,4)
[f,x]=hist(ssrrt,13);
bar(x,f/trapz(x,f));
hold on
[ff,xi] = ksdensity(ssrrt(1:98,1));
plot(xi,ff);
title('(iv)');
axis([-2 3 0 0.55]);
hold off
