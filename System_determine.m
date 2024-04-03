
% System determine


clc
clear 

load example1_sys.mat

n=2000;tau=1;


tilde_A=[zeros(n) eye(n); -inv(M)*K -inv(M)*D];
tilde_A_d=[zeros(n) zeros(n); -inv(M)*F -inv(M)*E];
tilde_B=[zeros(n,1); inv(M)*B];
tilde_C=[C zeros(1,n)];

% tilde_A=sparse(tilde_A);
% tilde_A_d=sparse(tilde_A_d);
% tilde_B=sparse(tilde_B);
% tilde_C=sparse(tilde_C);

t=-3:0.01:2;
w=10.^t;

DelayOri=struct('delay', {tau}, 'a',[],'b',[],'c',[], 'd',[]);
DelayOri(1).a=tilde_A_d;
% D=0;D_kr=0;
G=delayss(tilde_A,tilde_B,tilde_C,0,DelayOri);
% G=ss(tilde_A,tilde_B,tilde_C,0);  % second-order linear system
[magdb phase]=bode(G,w);
mag2= squeeze(magdb);
phase = squeeze(phase);
magdb = 20*log10(mag2);

figure(1)
semilogx(w,magdb);
legend('Orig');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');

















