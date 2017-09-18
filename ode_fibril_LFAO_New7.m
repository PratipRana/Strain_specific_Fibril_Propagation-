function dA_dt=ode_fibril_LFAO_New7(t,A ,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Defines the set of ODEs to be solved to simulate insulin fibrillation %
% A(1~n-1) is the vector of i-mer concentrations (i=1~n-1) %
% A(n) is fibril conc. and A(n+1) is natural hexamer conc. %
% t is time and dA_dt is the first order derivatives of A vector %
% n is the critical size of clusters and is assigned to be 6 %
% theta vector is the set of rate constants [knu1, kfb1, kfb_] %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dA_dt=zeros(size(A)); % First order derivatives of i-mer concentrations
% Flux of i-th fibrillation reaction
% The following is the list of rate constants
% Definitions of reaction fluxes Jfb
 % Reverse fibrillation rate constant
% kfb=theta(1);
kfb_=theta(2);
kel=theta(3); % First forward fibrillation rate constant
kel_=theta(4);
knu=theta(5);
knu_=theta(6);
kfb=ones(n,1)*theta(1); % First foward nucleation rate constants
for i=1:n-1
 kfb(i)=theta(1)/2*(1+i^(-1/3)); % Correct knu(i) by Stokes-Einstein Eq.
end
% Definitions of reaction fluxes Jfb
for i=1:n-1
 Jfb(i)=kfb(i)*A(n+1)*A(i)-kfb_*A(i+1); % The flux of i-mer nucleation rxn
 Jel(i)=kel*A(n)*A(i)-kel_*A(n); % The flux of i-mer elongation rxn
end

% There are n equations representing the conc. change of n species
% Derivative of monomer conc.
dA_dt(1)=-Jfb(1)-Jel(1);% Derivative of monomer conc.
% dA_dt(2)=-Jfb(2)+Jfb(1)+Jel(1);
for i=2:n-1 % from dimer to (n-1)-mer
 dA_dt(i)=-Jfb(i)+Jfb(i-1)-Jel(i); % Derivatives of oligomer concentrations
end

dA_dt(n)=Jfb(n-1);
dA_dt(n+1)=12*-sum(Jfb(1:n-1))-knu*A(n)*A(n+1)+knu_*A(n);
%dA_dt(n+1)=-12*(knu*A(n+1))+knu_*A(n+2)-kel*A(n)*A(n+1)+kel_*A(n);
%dA_dt(n+2)=knu*A(n+1)-knu_*A(n+2)-sum(Jfb(1:n-1));

% Derivative of fibril concentration
 % Derivative of insulin hexamer concentration 