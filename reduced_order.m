
% ensure the reduced order


% s0=min(rank(Pa),rank(Qa));
% % s0=60;
% t=1:s0;
% sigma=1e-5;
% semilogy(t,s(t),'*',t,sigma*s(1)*ones(s0,1),'k')









s0=65;
t=1:s0;
sigma=1e-10;
semilogy(t,ss1(t),'b*',t,ss2(t),'m<',t,ss3(t),'gp')  % ,'MarkerSize',4
xlim([1,65])
ylim([1e-15,1e2])
legend('m=70','m=80','m=90')
xlabel('Number')
ylabel('Singular values')
hold on
% semilogy([37,37],[0,1e6],'r')
stem(r,1e2,'.')
hold off



