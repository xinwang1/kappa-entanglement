%%%%%%
%%% Requred packages 
%%% 1. CVX: http://cvxr.com/ 2. QETLAB: http://www.qetlab.com/Main_Page
%% input the biparite quantum state rhoAB
clear;
da = 4; %dimension of system A
db = 4; %dimension of system B¡®
rhoAB = RandomDensityMatrix(da*db,1,3);

%% compute the kappa entanglement

cvx_begin sdp quiet
variable SAB(da*db, da*db) hermitian ;
s = trace(SAB) 
    minimize s
    SAB >= 0 ;
    - PartialTranspose(SAB,2) <= PartialTranspose(rhoAB,2);
    PartialTranspose(SAB,2) >= PartialTranspose(rhoAB,2) ;
cvx_end

E_kappa_rho = log2(s)