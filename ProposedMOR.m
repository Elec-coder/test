
clear
clc

load example3_174.mat

n=174;r=35;
tau=1;
% omega=1:m;

%---------------------------------------------------------------------------------- proposed MOR algorithm
tic
% t=linspace(-1,2,100);
t=-2:0.01:2;
omega=10.^t;
% omega=0:0.005:1;
m=length(omega);
Zc=zeros(n,m);
Zo=zeros(n,m);

% matlabpool open;
parfor k=1:m
    R=inv((i*omega(k))^2*M+i*omega(k)*D+K+i*omega(k)*E*exp(-i*omega(k)*tau)+F*exp(-i*omega(k)*tau));
    zc=R*B;
    zo=R'*C';
    Zc(:,k)=zc;
    Zo(:,k)=zo;
end
% matlabpool close;

Pa=Zc*Zc';Qa=Zo*Zo';
% r=min(rank(Pa),rank(Qa))

[U,S,V]=svd(Zc'*D'*Zo);
s=diag(S);
s1=zeros(r,1);
for j=1:r
    s1(j)=1/sqrt(s(j));
end
Tl=diag(s1)*V(:,1:r)'*Zo';
Tr=Zc*U(:,1:r)*diag(s1);

Mr=Tl*M*Tr;Dr=Tl*D*Tr;Kr=Tl*K*Tr;Er=Tl*E*Tr;Fr=Tl*F*Tr;Br=Tl*B;Cr=C*Tr;
toc

%------------------------------------------------------------------------------------------------ simulation

tilde_A=[zeros(n) eye(n); -inv(M)*K -inv(M)*D];
tilde_A_d=[zeros(n) zeros(n); -inv(M)*F -inv(M)*E];
tilde_B=[zeros(n,1); inv(M)*B];
tilde_C=[C zeros(1,n)];

Ar=[zeros(r) eye(r); -inv(Mr)*Kr -inv(Mr)*Dr];
Adr=[zeros(r) zeros(r); -inv(Mr)*Fr -inv(Mr)*Er];
Br=[zeros(r,1);inv(Mr)*Br];
Cr=[Cr zeros(1,r)];


t=-1:0.01:2;
w=10.^t;

% tic
% DelayOri=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
% DelayOri(1).a=tilde_A_d;
% D=0;
% G=delayss(tilde_A,tilde_B,tilde_C,D,DelayOri);
% [magdb phase]=bode(G,w);
% mag2= squeeze(magdb);
% phase = squeeze(phase);
% magdb = 20*log10(mag2);
% toc  %% 416.211765 seconds
% save Frequency_orig magdb
% load Frequency_orig.mat magdb
% 
% tic
% DelayRed=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
% DelayRed(1).a=Adr;
% Gr=delayss(Ar,Br,Cr,0,DelayRed);
% [magdbr phaser]=bode(Gr,w);
% mag2r= squeeze(magdbr);
% phaser = squeeze(phaser);
% magdbr = 20*log10(mag2r);
% toc
% 
% figure(1)
% semilogx(w,magdb,'m',w,magdbr,'b-.');  % ,'LineWidth',2
% legend('Orig','Proposed');
% xlabel('Frequency (rad/s)');
% ylabel('Magnitude (dB)');
% 
% error1=abs(magdb-magdbr)./abs(magdb);
% 
% 
% figure(2)
% loglog(w,error1,'b-');
% legend('Proposed')
% xlabel('Frequency (rad/s)');
% ylabel('Relative error');

%%------------------------------------------------------------------------------ Time Domain
x0=zeros(2*n,1);xr0=zeros(2*r,1);
tspan=0:0.001:1;
sol=dde23(@(t,x,x_tau) tilde_A*x+tilde_A_d*x_tau+tilde_B*sin(2.5*pi*t-1),tau,@(t)x0,tspan); % @(t)ones(2*n,1)*cos(t)
t=0:0.001:1;
x=deval(sol,t);
y=tilde_C*x;

solr=dde23(@(t,xr,xr_tau) Ar*xr+Adr*xr_tau+Br*sin(2.5*pi*t-1),tau,@(t)[Tl zeros(r,n);zeros(r,n) Tl]*x0,tspan); % @(t)ones(2*n,1)*cos(t)
xr=deval(solr,t);
yr=Cr*xr;

figure(1)
plot(t,y,'k',t,yr,'b--')
legend('Orig','Proposed')
xlabel('Time t (second)')
ylabel('y(t)')

error=abs(y-yr)./abs(y);
figure(2)
semilogy(t,error,'b--')
legend('Proposed')
xlabel('Time t (second)')
ylabel('Relative Error')





