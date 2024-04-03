
clear
clc

load example3_174.mat

n = 174; r=50;
tau = 1;
alpha = 10; %  Laguerre MOR (alpha=24.5,r=45),
% omega=1:m;

%---------------------------------------------------------------------------------- proposed MOR algorithm 50 samples
tic
t1 = linspace(-1,2,70);  % 觉得真正有用的区间只有 [10^(-3), 10]；因为取其他区间的样本点不增加 min(rank(Pa),rank(Qa))
omega=10.^t1;

m = length(omega);
Zc=zeros(n,m);
Zo=zeros(n,m);
parfor k=1:m
    R=inv((i*omega(k))^2*M+i*omega(k)*D+K+i*omega(k)*E*exp(-i*omega(k)*tau)+F*exp(-i*omega(k)*tau));
    zc=R*B;
    zo=R'*C';
    Zc(:,k)=zc;
    Zo(:,k)=zo;
end
save Frequency_sys1 Zc Zo
load Frequency_sys1.mat Zc Zo

Pa=Zc*Zc';Qa=Zo*Zo';
% r=min(rank(Pa),rank(Qa))

[U,S,V]=svd(Zc'*D'*Zo);
ss1=diag(S);
s1=zeros(r,1);
for j=1:r
    s1(j)=1/sqrt(ss1(j));
end
Tl1=diag(s1)*V(:,1:r)'*Zo';
Tr1=Zc*U(:,1:r)*diag(s1);

Mr1=Tl1*M*Tr1;Dr1=Tl1*D*Tr1;Kr1=Tl1*K*Tr1;Er1=Tl1*E*Tr1;Fr1=Tl1*F*Tr1;Br1=Tl1*B;Cr1=C*Tr1;
toc
%---------------------------------------------------------------------------------- proposed MOR algorithm 65 samples
tic
t2=linspace(-1,2,80);
omega=10.^t2;

m=length(omega);
Zc=zeros(n,m);
Zo=zeros(n,m);
parfor k=1:m
    R=inv((i*omega(k))^2*M+i*omega(k)*D+K+i*omega(k)*E*exp(-i*omega(k)*tau)+F*exp(-i*omega(k)*tau));
    zc=R*B;
    zo=R'*C';
    Zc(:,k)=zc;
    Zo(:,k)=zo;
end
save Frequency_sys2 Zc Zo
load Frequency_sys2.mat Zc Zo

Pa=Zc*Zc';Qa=Zo*Zo';
% r=min(rank(Pa),rank(Qa))

[U,S,V]=svd(Zc'*D'*Zo);
ss2=diag(S);
s1=zeros(r,1);
for j=1:r
    s1(j)=1/sqrt(ss2(j));
end
Tl2=diag(s1)*V(:,1:r)'*Zo';
Tr2=Zc*U(:,1:r)*diag(s1);

Mr2=Tl2*M*Tr2;Dr2=Tl2*D*Tr2;Kr2=Tl2*K*Tr2;Er2=Tl2*E*Tr2;Fr2=Tl2*F*Tr2;Br2=Tl2*B;Cr2=C*Tr2;
toc
%---------------------------------------------------------------------------------- proposed MOR algorithm 80 samples
tic
t3=linspace(-1,2,90);
omega=10.^t3;

m=length(omega);
Zc=zeros(n,m);
Zo=zeros(n,m);
parfor k=1:m
    R=inv((i*omega(k))^2*M+i*omega(k)*D+K+i*omega(k)*E*exp(-i*omega(k)*tau)+F*exp(-i*omega(k)*tau));
    zc=R*B;
    zo=R'*C';
    Zc(:,k)=zc;
    Zo(:,k)=zo;
end
save Frequency_sys3 Zc Zo
load Frequency_sys3.mat Zc Zo

Pa=Zc*Zc';Qa=Zo*Zo';
% r=min(rank(Pa),rank(Qa))

[U,S,V]=svd(Zc'*D'*Zo);
ss3=diag(S);
s1=zeros(r,1);
for j=1:r
    s1(j)=1/sqrt(ss3(j));
end
Tl3=diag(s1)*V(:,1:r)'*Zo';
Tr3=Zc*U(:,1:r)*diag(s1);

Mr3=Tl3*M*Tr3;Dr3=Tl3*D*Tr3;Kr3=Tl3*K*Tr3;Er3=Tl3*E*Tr3;Fr3=Tl3*F*Tr3;Br3=Tl3*B;Cr3=C*Tr3;
toc
%-------------------------------------------------------------------------------------------- Laguerre MOR
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
% load system_lag.mat Mr_lag Dr_lag Kr_lag Er_lag Fr_lag Br_lag Cr_lag

%------------------------------------------------------------------------------------------------ simulation

tilde_A=[zeros(n) eye(n); -inv(M)*K -inv(M)*D];
tilde_A_d=[zeros(n) zeros(n); -inv(M)*F -inv(M)*E];
tilde_B=[zeros(n,1); inv(M)*B];
tilde_C=[C zeros(1,n)];

Ar1=[zeros(r) eye(r); -inv(Mr1)*Kr1 -inv(Mr1)*Dr1];
Adr1=[zeros(r) zeros(r); -inv(Mr1)*Fr1 -inv(Mr1)*Er1];
Br1=[zeros(r,1);inv(Mr1)*Br1];
Cr1=[Cr1 zeros(1,r)];

Ar2=[zeros(r) eye(r); -inv(Mr2)*Kr2 -inv(Mr2)*Dr2];
Adr2=[zeros(r) zeros(r); -inv(Mr2)*Fr2 -inv(Mr2)*Er2];
Br2=[zeros(r,1);inv(Mr2)*Br2];
Cr2=[Cr2 zeros(1,r)];

Ar3=[zeros(r) eye(r); -inv(Mr3)*Kr3 -inv(Mr3)*Dr3];
Adr3=[zeros(r) zeros(r); -inv(Mr3)*Fr3 -inv(Mr3)*Er3];
Br3=[zeros(r,1);inv(Mr3)*Br3];
Cr3=[Cr3 zeros(1,r)];

Ar_lag=[zeros(r) eye(r); -inv(Mr_lag)*Kr_lag -inv(Mr_lag)*Dr_lag];
Adr_lag=[zeros(r) zeros(r); -inv(Mr_lag)*Fr_lag -inv(Mr_lag)*Er_lag];
Br_lag=[zeros(r,1);inv(Mr_lag)*Br_lag];
Cr_lag=[Cr_lag zeros(1,r)];

t=-1:0.01:2;
w=10.^t;

DelayOri=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayOri(1).a=tilde_A_d;
D=0;
G=delayss(tilde_A,tilde_B,tilde_C,0,DelayOri);
[magdb phase]=bode(G,w);
mag2= squeeze(magdb);
phase = squeeze(phase);
magdb = 20*log10(mag2);
save Frequency_orig magdb
load Frequency_orig.mat magdb

DelayRed1=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayRed1(1).a=Adr1;
Gr1=delayss(Ar1,Br1,Cr1,0,DelayRed1);
[magdbr1 phaser1]=bode(Gr1,w);
mag2r1= squeeze(magdbr1);
phaser1 = squeeze(phaser1);
magdbr1 = 20*log10(mag2r1);

DelayRed2=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayRed2(1).a=Adr2;
Gr2=delayss(Ar2,Br2,Cr2,0,DelayRed2);
[magdbr2 phaser2]=bode(Gr2,w);
mag2r2= squeeze(magdbr2);
phaser2 = squeeze(phaser2);
magdbr2 = 20*log10(mag2r2);

DelayRed3=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayRed3(1).a=Adr3;
Gr3=delayss(Ar3,Br3,Cr3,0,DelayRed3);
[magdbr3 phaser3]=bode(Gr3,w);
mag2r3= squeeze(magdbr3);
phaser3 = squeeze(phaser3);
magdbr3 = 20*log10(mag2r3);

DelayRed_lag=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayRed(1).a=Adr_lag;
Gr_lag=delayss(Ar_lag,Br_lag,Cr_lag,0,DelayRed_lag);
[magdb_lag phase_lag]=bode(Gr_lag,w);
mag2_lag= squeeze(magdb_lag);
phase_lag = squeeze(phase_lag);
magdb_lag = 20*log10(mag2_lag);

figure(1)
semilogx(w,magdb,'k',w,magdb_lag,'r:',w,magdbr1,'m',w,magdbr2,'g-.',w,magdbr3,'b--','LineWidth',1.5);  % ,'LineWidth',1.5
% legend('Orig','50','65','80');
legend('Orig S','Red S-1','Red S-2, m=60','Red S-2, m=75','Red S-2, m=90');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');

error1=abs(magdb-magdbr1)./abs(magdb);
error2=abs(magdb-magdbr2)./abs(magdb);
error3=abs(magdb-magdbr3)./abs(magdb);
error_lag=abs(magdb-magdb_lag)./abs(magdb);

figure(2)
loglog(w,error_lag,'r:',w,error1,'m',w,error2,'g-.',w,error3,'b--','LineWidth',1.5);
% legend('Proposed','Arnoldi','Laguerre')
legend('Red S-1','Red S-2, m=60','Red S-2, m=75','Red S-2, m=90');
xlabel('Frequency (rad/s)');
ylabel('Relative error');

% tic
% bode(G,'m',Gr,'b-.',Gr_kr,'g--',Gr_lag,'r:',{1e-4,1e0})
% legend('Orig','Proposed','Arnoldi','Laguerre');
% toc






