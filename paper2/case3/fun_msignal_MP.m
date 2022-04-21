
function [y_hat,eig_s]= fun_msignal_MP(ya, M, dT)

 % ya is the inut signal. It should be a column vector for each signal. you
 % can give multiple signals. Each occupy one column.  
 % M is the assumed system order. 
 % h is the sampling time perod. 
 % output: eig_s gives identified eigenvalue of continuous system. Not in
 % discrete domian. 
 % Fig. 888 gives eigenvalue plot. 
 % Fig. 999 gives the match plot. 
 
% how many channels? set up subplot dimension. 
n_ch = size(ya, 2); 
N = size(ya,1)-1; t1 = 0:dT:N*dT;

if(n_ch<= 5)
    row_plot = n_ch;
    col_plot = 1; 
else
    if(mod(sqrt(n_ch),1)>0)
        row_plot= floor(sqrt(n_ch))+1; 
    else
        row_plot=sqrt(n_ch);
    end
    if (mod(n_ch/row_plot,1)>0)
        col_plot = floor(n_ch/row_plot) + 1; 
    else 
        col_plot = n_ch/row_plot; 
    end
    % e.g., 6 signals: 3*2
    % e.g.; 7 signals: 3*3
end

% figure(999);
% for i=1:n_ch
%     subplot(row_plot, col_plot, i); plot(t1, ya(:,i)); %legend('original signal');
%     hold on;
% end

D =[];
%N = size(ya,1)-1; t1 = 0:dT:N*dT;
L = floor(1/3*N);
%M=10; % order 

for k=1:size(ya,2) % visit each column
    % for each channel, build a Hankel matrix
    for i=1:L+2  
        H(i,:) = ya(i:i+N-L-1); 
    end
    D =[D, H];
end

[U,S,V] = svd(D); 
%figure(999)
%semilogy(diag(S));

U_prime= U(:,1:M);
U1 = U_prime(1:L+1, :);
U2 = U_prime(2:end, :); 

%Lambda = inv(U2'*U1)*(U2'*U2);
Lambda = inv(U1'*U1)*(U1'*U2);
z=eig(Lambda);
eig_s = log(z)/dT;

% figure(888);
% plot(real(eig_s), imag(eig_s)/2/pi,'+','Linewidth',2, 'Markersize',10);
% ylabel('Hz')
% title('eigenvalues'); grid on;



%% signal reconstruction 
for i1=1:N+1;
    for j1=1:M;
        Z(i1,j1)=z(j1)^(i1-1);
    end 
end
for i=1:size(ya,2) % for three signal, reconstruct
residue1 = pinv(Z)*ya(:,i);
y_hat(:,i)=Z*residue1;
end
% 
% figure(999); 
% for i=1:n_ch
%     subplot(row_plot, col_plot, i);      
%     %plot(t1, ya(:,i),'Linewidth',2); hold on;
%     plot(t1, real(y_hat(:,i)),'r','Linewidth',1);  
%     grid on
%     legend('Original', 'MP');
% end
