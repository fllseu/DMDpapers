
function  [Q,B] = rsvd(f,r,q,p)

% 10/22/2019 

% I did this code, which is from 
% "Randomized Dynamic Mode Decomposition (2018)"
% Algorithm 3.1 Randomized QB decomposition.


m = size(f,2); % to get the number of columns
P = randn(m,r+p); %rank(f)=21 +p=10
Y = f*P;

for k=1:1
    [Q,~] = qr(Y,0);
    [Z,~] = qr(f'*Q,0);
    Y = f*Z;
end

[Q,~] = qr(Y,0);     % To get the obtimal basis
B=Q'*f;   %the small equivalent matrix


