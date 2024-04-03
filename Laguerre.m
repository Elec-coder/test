clc;
clear;
load example3_174.mat;
n = 174; r=50;
tau = 1;
t=0:0.001:1;
alpha = 10; %
p=length(t);
k1=4;
x0=zeros(2*n,1);
xr0=zeros(2*r,1);
tspan=0:0.001:1;
tic
g0=sqrt(2*alpha)*exp(-alpha*tau)*laguerreL(0,2*alpha*tau);
g=zeros(r-3,1);
for j=1:r-3
    g(j)=sqrt(2*alpha)*exp(-alpha*tau)*laguerreL(j,2*alpha*tau);
end

delta_0=(M*alpha^2-alpha*D+K)+sqrt(alpha/2)*g0*E+sqrt(1/(2*alpha))*g0*F;
deltal_0=delta_0'; % two-sided
Delta=cell(r,1);Deltal=cell(r,1);
Delta{1}=(2*M*alpha^2-2*K)+sqrt(alpha/2)*(-g0+g(1))*E+sqrt(1/(2*alpha))*(-3*g0+g(1))*F;
Delta{2}=(M*alpha^2-alpha*D+K)+sqrt(alpha/2)*(-g0-g(1)+g(2))*E+sqrt(1/(2*alpha))*(3*g0-3*g(1)+g(2))*F;
Delta{3}=sqrt(alpha/2)*(g0-g(1)-g(2)+g(3))*E+sqrt(1/(2*alpha))*(-g0+3*g(1)-3*g(2)+g(3))*F;
Deltal{1}=Delta{1}';Deltal{2}=Delta{2}';Deltal{3}=Delta{3}';%
for j=4:r-3
    Delta{j}=sqrt(alpha/2)*(g(j-3)-g(j-2)-g(j-1)+g(j))*E+sqrt(1/(2*alpha))*(-g(j-3)+3*g(j-2)-3*g(j-1)+g(j))*F;
    Deltal{j}=Delta{j}';
end
Delta{r-2}=sqrt(alpha/2)*(g(r-5)-g(r-4)-g(r-3))*E+sqrt(1/(2*alpha))*(-g(r-5)+3*g(r-4)-3*g(r-3))*F;
Delta{r-1}=sqrt(alpha/2)*(g(r-4)-g(r-3))*E+sqrt(1/(2*alpha))*(-g(r-4)+3*g(r-3))*F;
Delta{r}=sqrt(alpha/2)*g(r-3)*E+sqrt(1/(2*alpha))*(-g(r-3))*F;
Deltal{r-2}=Delta{r-2}';Deltal{r-1}=Delta{r-1}';Deltal{r}=Delta{r}';%

l=inv(delta_0)*B; % inv(delta_0)*B,(sqrt(2*alpha))*inv(delta_0)*B
ll=inv(deltal_0)*C'; % inv(deltal_0)*C',(sqrt(2*alpha))*inv(deltal_0)*C'
AI=cell(r,1);
AIl=cell(r,1); %
for j=1:r
    AI{j}=inv(delta_0)*Delta{j};  %%--
    AIl{j}=inv(deltal_0)*Deltal{j}; %  %%--
end

q=zeros(n,r);q(:,1)=l/norm(l);
p=zeros(n,r-1);
T=zeros(r+1,r);
ql=zeros(n,r);ql(:,1)=ll/norm(ll); %
pl=zeros(n,r-1); %
Tll=zeros(r+1,r); %
for j=1:r
    v=AI{1}*q(:,j);
    vl=AIl{1}*ql(:,j); %
    for k=2:r
        v=v+AI{k}*p(:,r+1-k);
        vl=vl+AIl{k}*pl(:,r+1-k); %
    end
    for k=1:j
        T(k,j)=q(:,k)'*v;
        v=v-q(:,k)*T(k,j);
        Tll(k,j)=ql(:,k)'*vl;  %
        vl=vl-ql(:,k)*Tll(k,j);  %
    end
    T(j+1,j)=norm(v);
%     hat_T=[T;e(j,:)*norm(v)];
    q(:,j+1)=v/norm(v);
    Tll(j+1,j)=norm(vl);ql(:,j+1)=vl/norm(vl);  %
    for k=1:r-1
        p(:,r-k)=q(:,1:j+1)*[zeros(j,1) inv(T(2:j+1,1:j)); 0 zeros(1,j)]^k*[zeros(j,1);1];
        pl(:,r-k)=ql(:,1:j+1)*[zeros(j,1) inv(Tll(2:j+1,1:j)); 0 zeros(1,j)]^k*[zeros(j,1);1];  %
    end
end
V_lag=q(:,1:r);
W_lag=ql(:,1:r);  %

Mr_lag=W_lag'*M*V_lag;Dr_lag=W_lag'*D*V_lag;Kr_lag=W_lag'*K*V_lag;Er_lag=W_lag'*E*V_lag;Fr_lag=W_lag'*F*V_lag;Br_lag=W_lag'*B;Cr_lag=C*V_lag;
toc
save system_lag Mr_lag Dr_lag Kr_lag Er_lag Fr_lag Br_lag Cr_lag

Ar_lag=[zeros(r) eye(r); -inv(Mr_lag)*Kr_lag -inv(Mr_lag)*Dr_lag];
Adr_lag=[zeros(r) zeros(r); -inv(Mr_lag)*Fr_lag -inv(Mr_lag)*Er_lag];
Br_lag=[zeros(r,1);inv(Mr_lag)*Br_lag];
Cr_lag=[Cr_lag zeros(1,r)];


% solr_lag=dde23(@(t,xr_lag,xr_lag_tau) Ar_lag*xr_lag+Adr_lag*xr_lag_tau+Br_lag*exp(-1*t)*cos(15*t),tau,@(t)[V_lag' zeros(r,n);zeros(r,n) V_lag']*x0,tspan); % sin(2*pi*t)*sinh(2*t),sinh(-2*t)*sin(1.5*pi*t),sinh(2*t)*sin(3*pi*t)
% % yr_lag=Cr_lag*solr_lag.y;  %
% % sinh(3*t)*sin(2*pi*t-0.5),sinh(-2*t)*sin(1.5*pi*t),sinh(-2*t)*cosh(-2*pi*t),sinh(-2*t)*sin(2*pi*t),sinh(-2*t+0.01)*sin(-2*pi*t),sinh(-3*t+0.05)*sin(-2*pi*t),cos(3*t)*cos(10*pi*t)
% xr_lag=deval(solr_lag,t);
% yr_lag=Cr_lag*xr_lag;
% plot(t(1:k1:p),yr_lag(1:k1:p),'rx')
% % plot(t(1+k1/4:k1:p),yr_lag(1+k1/4:k1:p),'rx',t(1:k1:p),yr1(1:k1:p),'b*',t(1+3*k1/4:k1:p),yr2(1+3*k1/4:k1:p),'m<',t(1+k1/2:k1:p),yr3(1+k1/2:k1:p),'gp',t,y(:),'k')
% legend('Red S-1')
% % plot(t,y,'k',t,yr_lag,'r:',t,yr1,'m',t,yr2,'g-.',t,yr3,'b--','LineWidth',2)
% % legend('Orig','Laguerre','50 samples','65 samples','80 samples')
% % legend('Orig S','Red S-1','Red S-2, m=60','Red S-2, m=70','Red S-2, m=80')
% xlabel('Time t (second)')
% ylabel('y(t)')







t=-1:0.01:2;
w=10.^t;


 DelayRed_lag=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayRed(1).a=Adr_lag;
Gr_lag=delayss(Ar_lag,Br_lag,Cr_lag,0,DelayRed_lag);
[magdb_lag phase_lag]=bode(Gr_lag,w);
mag2_lag= squeeze(magdb_lag);
phase_lag = squeeze(phase_lag);
magdb_lag = 20*log10(mag2_lag);
semilogx(w,magdb_lag,'r:', 'LineWidth',1.5); 
% legend('Orig','50','65','80');
legend('Red S-1');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');












