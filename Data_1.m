
% Р§ 1 Ъ§Он

clear
clc

n=2000;q1=80;q2=q1;   % q1=q2=90,115


B=zeros(n,1);
B(1,1)=1;
C=B';
alpha=0.05;beta=0.05;   % alpha=0.05;beta=0.05;
M=zeros(n,n);
K=zeros(n,n);
for i=2:n
    M(i-1,i)=1;
    M(i,i-1)=1;
    M(i,i)=2/sqrt(1-alpha*beta);
    K(i-1,i)=-1;
    K(i,i-1)=-1;
    K(i,i)=2/sqrt(1-alpha*beta);
end
M(1,1)=(2+sqrt(1-alpha*beta))/sqrt(1-alpha*beta);
M(n,n)=M(1,1);
K(1,1)=(2-sqrt(1-alpha*beta))/sqrt(1-alpha*beta);
K(n,n)=K(1,1);
K=(alpha/beta)*K;
D=alpha*M+beta*K;


e1=0.1*[0;ones(q1-1,1)];  % -0.1,0.3,-0.1
e2=0.5*ones(q1,1);
e3=0.1*[ones(q1-1,1);0];
e=[e1 e2 e3];
d=[1;0;-1];
E1=spdiags(e,d,q1,q1);
E=[E1 zeros(q1,n-q1);zeros(n-q1,q1) zeros(n-q1)];
E=full(E);


e1=0.1*[0;ones(q2-1,1)];  % -0.1,0.5,0
e2=0.3*ones(q2,1);
e3=0.1*[ones(q2-1,1);0];
e=[e1 e2 e3];
d=[1;0;-1];
E1=spdiags(e,d,q2,q2);
F=[E1 zeros(q2,n-q2);zeros(n-q2,q2) zeros(n-q2)];
F=full(F);
F=2*F;


% F=[zeros(n-q) zeros(n-q,q);zeros(q,n-q) E1];
% F=full(F);
% f1=[0;-1*ones(n-1,1)];
% f2=2:n+1;%2*ones(n,1);
% f3=[-1*ones(n-1,1);0];
% f=[f1 f2' f3];
% d=[1;0;-1];
% F=spdiags(f,d,n,n);
% F=diag(1:n);
% F=full(F);
% E=5*F;

% load beam_348_2.mat
% M=full(M);D=full(D);K=full(K);B=B_2;C=C_1';


save example1_sys M D K E F B C




















