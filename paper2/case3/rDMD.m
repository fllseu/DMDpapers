function [lambda,Phi,b,X2r]=rDMD(f,signal_number,system_rank,Q)


% 10/22/2019 


%% DMD operations
X1 = f(:,1:end-1);
X2 = f(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The system order
%r=ss%rank(X1);
ss = signal_number; 
r = system_rank; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,S,V]=svd(X1,'econ');
Ur=U(:,1:r); %To get the most effective if U
Vr=V(:,1:r); %To get the most effective of V
Sr=S(1:r,1:r); %To get the active eigenvalues
Atilda=Ur'*X2*Vr/Sr;

[W,D]=eig(Atilda); % to get the eigenvectors and eigenvalues
lambda=(diag(D)); % the eigenvalues are the diag of D
% omega=log(lambda)/dt; %to get the value of omega form (exp(omega*t))
%  Phi=Q*X2*Vr/Sr*W; %exact modes
  Phi = Q*Ur*W;
%P_freq=freq( freq>=0 ); % positive freq
x1=X1(:,2); % the intial condtion (x0)       <---  test it 
b=Phi\(Q*x1); %from the equation Phi*b=x1
%%%%%%%%%%%%%%%%%%%%%%


%% 01/15/2020
%%The IEEE paper steps Good for check

%     X2r=Ur*W*D*inv(W)*Sr*Vr';