

% simulation
% load example1_sys.mat

% n=1000;r=16;
tau=1;
x0=zeros(2*n,1);xr0=zeros(2*r,1);
tspan=0:0.001:1; %0:0.001:1; linspace(0,1,1000)
sol=dde23(@(t,x,x_tau) tilde_A*x+tilde_A_d*x_tau+tilde_B*exp(-1*t)*cos(15*t),tau,@(t)x0,tspan); % @(t)[ones(n,1)*sin(t); ones(n,1)*cos(t)]  % sinh(2*t)*sin(4*pi*t),sinh(2*t)*cos(3*pi*t+1),cos(10*pi*t)
% y=tilde_C*sol.y;
t=0:0.001:1;
x=deval(sol,t);
y=tilde_C*x;

solr1=dde23(@(t,xr1,xr1_tau) Ar1*xr1+Adr1*xr1_tau+Br1*exp(-1*t)*cos(15*t),tau,@(t)[Tl1 zeros(r,n);zeros(r,n) Tl1]*x0,tspan); % 不受影响的输入：exp(t),cos(t)
% yr=Cr*solr.y;
xr1=deval(solr1,t);
yr1=Cr1*xr1;

solr2=dde23(@(t,xr2,xr2_tau) Ar2*xr2+Adr2*xr2_tau+Br2*exp(-1*t)*cos(15*t),tau,@(t)[Tl2 zeros(r,n);zeros(r,n) Tl2]*x0,tspan); % 受影响的输入：exp(-t),sin(t)
% yr=Cr*solr.y;
xr2=deval(solr2,t);
yr2=Cr2*xr2;

solr3=dde23(@(t,xr3,xr3_tau) Ar3*xr3+Adr3*xr3_tau+Br3*exp(-1*t)*cos(15*t),tau,@(t)[Tl3 zeros(r,n);zeros(r,n) Tl3]*x0,tspan); % 其他输入：(heaviside(t-0.5)-heaviside(0.5-t))   
% yr=Cr*solr.y;
xr3=deval(solr3,t);
yr3=Cr3*xr3;

solr_lag=dde23(@(t,xr_lag,xr_lag_tau) Ar_lag*xr_lag+Adr_lag*xr_lag_tau+Br_lag*exp(-1*t)*cos(15*t),tau,@(t)[V_lag' zeros(r,n);zeros(r,n) V_lag']*x0,tspan); % sin(2*pi*t)*sinh(2*t),sinh(-2*t)*sin(1.5*pi*t),sinh(2*t)*sin(3*pi*t)
% yr_lag=Cr_lag*solr_lag.y;  %
% sinh(3*t)*sin(2*pi*t-0.5),sinh(-2*t)*sin(1.5*pi*t),sinh(-2*t)*cosh(-2*pi*t),sinh(-2*t)*sin(2*pi*t),sinh(-2*t+0.01)*sin(-2*pi*t),sinh(-3*t+0.05)*sin(-2*pi*t),cos(3*t)*cos(10*pi*t)
xr_lag=deval(solr_lag,t);
yr_lag=Cr_lag*xr_lag;

error1=abs(y-yr1)./abs(y);
error2=abs(y-yr2)./abs(y);
error3=abs(y-yr3)./abs(y);
error_lag=abs(y-yr_lag)./abs(y);

p=length(t);
k1=4;


figure(1)
plot(t(1:k1:p),yr_lag(1:k1:p),'rx',t(1+k1/4:k1:p),yr1(1+k1/4:k1:p),'b*',t(1+k1/2:k1:p),yr2(1+k1/2:k1:p),'m<',t(1+3*k1/4:k1:p),yr3(1+3*k1/4:k1:p),'gp',t,y(:),'k')
% plot(t(1+k1/4:k1:p),yr_lag(1+k1/4:k1:p),'rx',t(1:k1:p),yr1(1:k1:p),'b*',t(1+3*k1/4:k1:p),yr2(1+3*k1/4:k1:p),'m<',t(1+k1/2:k1:p),yr3(1+k1/2:k1:p),'gp',t,y(:),'k')
legend('Red S-1','Red S-2, m=70','Red S-2, m=80','Red S-2, m=90','Orig S')
% plot(t,y,'k',t,yr_lag,'r:',t,yr1,'m',t,yr2,'g-.',t,yr3,'b--','LineWidth',2)
% legend('Orig','Laguerre','50 samples','65 samples','80 samples')
% legend('Orig S','Red S-1','Red S-2, m=60','Red S-2, m=70','Red S-2, m=80')
xlabel('Time t (second)')
ylabel('y(t)')

% 
% figure(2)
% semilogy(t(1:k1:p),error_lag(1:k1:p),'r-x',t(1+k1/4:k1:p),error1(1+k1/4:k1:p),'b-*',t(1+k1/2:k1:p),error2(1+k1/2:k1:p),'m-<',t(1+3*k1/4:k1:p),error3(1+3*k1/4:k1:p),'g-p')
% legend('Red S-1','Red S-2, m=70','Red S-2, m=80','Red S-2, m=90')
% % semilogy(t,error_lag,'r:',t,error1,'m',t,error2,'g-.',t,error3,'b--','LineWidth',2)
% % legend('Red S-1','Red S-2, m=60','Red S-2, m=70','Red S-2, m=80')
% xlabel('Time t (second)')
% ylabel('Relative error')









