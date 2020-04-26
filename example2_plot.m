clear;
%%%%%%
%%% Requred packages 
%%% 1. CVX: http://cvxr.com/ 2. QETLAB: http://www.qetlab.com/Main_Page
%% set dimension
a=2;b=2;d=a*b;da=2;db=2;
dim=[a b];
%% plot
size=50;
for t=1:51
    r=(t-1)/size; 
    R(t)=r;
 
v1=[1 0 0 0 0 0 0 0 0]';
v1 = v1/norm(v1);
v2=[0 0 1 0 0 0 1 0 0]';
v2 = v2/norm(v2);
v3=[0 0 0 0 0 1 0 1 0]';
v3 = v3/norm(v3);
rho=r/2*v1*v1'+(1-r)/2*v2*v2' + 1/2*v3*v3';

%% compute the log negativity
EN(t)=log2(TraceNorm(PartialTranspose(rho,2)));

%% compute the kappa entanglement
cvx_begin sdp quiet
variable S(9, 9) hermitian ;
minimize ( trace(S)  ) ;
S >= 0 ;
- PartialTranspose(S,2) <= PartialTranspose(rho,2);
 PartialTranspose(S,2) >= PartialTranspose(rho,2) ;
cvx_end

EPPT(t) = log2(trace(S));
end
%  
 
bl = [0.05 0.32 0.78];
re = [0.8 0.28 0.08];
plot(R,EPPT,'LineWidth',2,'Color',bl)
grid on
hold on
plot(R,EN,'--','Color',re,'LineWidth',2)
%axis([0 1 0 1])
set(gca,'FontSize',12)
% leg = legend({'$R_{\max}(\mathcal{N}_{AD}^{(r )})$','$E_{PPT}(\mathcal{N}_{AD}^{(r)})$'});
leg = legend({'$E_{\rm{PPT}}(\omega_{p}) = E_{\kappa}(\omega_{p})$','$E_{N}({\omega_p})$'});
set(leg,'Interpreter','latex','FontSize',17,'Location','southwest');
xlabel('p from 0 to 1','FontSize',16)
ylabel('Ebit','FontSize',16)
