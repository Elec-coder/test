% clc;
% clear;
% 
% load('example3_174.mat');
n = 174; r=35;
tau = 1;
tilde_A=[zeros(n) eye(n); -inv(M)*K -inv(M)*D];
tilde_A_d=[zeros(n) zeros(n); -inv(M)*F -inv(M)*E];
tilde_B=[zeros(n,1); inv(M)*B];
tilde_C=[C zeros(1,n)];
x0=zeros(2*n,1);
x0r=zeros(2*r,1);
tspan=0:0.001:1;
t=0:0.001:1;
sol=dde23(@(t,x,x_tau) tilde_A*x + tilde_A_d*x_tau + tilde_B*exp(-1*t)*cos(15*t),tau,@(t)x0,tspan);
x=deval(sol,t);
y=tilde_C*x;

solr=dde23(@(t,xr,xr_tau) Ar_lag*xr + Adr_lag*xr_tau + Br_lag*exp(-1*t)*cos(15*t),tau,@(t)x0r,tspan);
xr=deval(solr,t);
yr=Cr_lag*xr;

plot(t,y,'k',t,yr,'b--')
legend('Orig','Proposed')
xlabel('Time t (second)')
ylabel('y(t)')