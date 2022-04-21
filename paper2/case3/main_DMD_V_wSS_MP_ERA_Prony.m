clear all, close all,clc


% 01/31/2020 : Case 3


%% Data Setting:-

load('case3row.mat') % fixed NaN measurments
T=num(:,1);
dt=T(2)-T(1);
% dt=T(3)-T(2);

start_time=51%130; % 70-70+10 gives 0.27
end_time=65%75%start_time+20;
ts=T(round(start_time/dt));
tf=T(round(end_time/dt));
ts_cell=round(start_time/dt);
tf_cell=round(end_time/dt);

tot=tf-ts; %total time

%% visualize the original signals
% figure(500)
%  plot(num(:,1),num(:,existed_signal)) 

%     ylim([59.98 60.06])

%%  % To get the needed measurments
y=[];
  mm=10+3; % freq=2 VM=3  ... check txt file 
    existed_signal=[mm:5:100];%176]%[mm:5:176];% To extract NaN measurments
%       existed_signal=existed_signal([3]);
 %  % the signal belong to the the following substation
 substation=[1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 8 8 8 8 9 9 10 10 10 11 11 12 12 12];   %check txt file 
%  subs=substation;
    
y=num(ts_cell:tf_cell,existed_signal); 
% y=num(round(ts/dt):round(tf/dt),existed_signal); % freq. starts w/ 2. Vm starts w/ 3. Vangle starts w/ 4
                                   % the first freq. signal is excluded bc.
                                   % it has zeros THAT's why we start with 7
                                   % 7

T=num(ts_cell:tf_cell,1);
% T=num(round(ts/dt):round(tf/dt),1);
ss=size(y,2); % Number of signals. For case 1 it's 31.
t = ts:dt:tf; 


%% Data pre-processing 

          y = interp1(T, y, t'); 


%       y=y';
% subtract mean-it helps.
%   [m,n]=size(y);   %  compute data size
%   mn=mean(y,1); %  compute mean for each row
%  y=y-repmat(mn,m,1);  %  repmat(mean,(repeat columns),(repeat rows))

%       y=detrend(y); % linearize the signal-it helps.

%       y=wdenoise(y); % linearize the signal-it helps.

%% Build The measurement matrix
 rows =round(0.20*length(y));% 2500 gives bad result  <------------------------------------------------------------
 m = length(y)-rows;

 %stack together
        y=y';
  f =[];
  tic
 for k = 1: rows
     f =[f; y(:,k:m+k)];
 end
 time.Hankel=toc;
 [U,S,V]=svd(f,'econ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The system order
% ss = 13; 
r = 20  %ss+200;%round(0.3*rank(f));%min(ss+20,round(0.7*rank(f)));  % the signals # is 31 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Standard DMD
tic
[lambda,Phi,b]=DMD(f,ss,r);
sDMD_time=toc;
time.sDMD_time=sDMD_time;
omega=log(lambda)/dt; %to get the value of omega form (exp(omega*t))


%% Randomized DMD

target_rank = r;
% the following values are suggested in the paper
q=2; % power iteratio. The paper suggests 2
p=10; %oversampling. The paper suggests 10 
tic
[Q,B] = rsvd(f,r,q,p);
rdecompostion_time=toc;
time.rdecompostion=rdecompostion_time;
tic
[lambdar,Phir,br]=rDMD(B,ss,target_rank,Q);
rDMD_time=toc;
time.rDMD=rDMD_time;
omegar=log(lambdar)/dt; %to get the value of omega form (exp(omega*t))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct the signals. 
 % % Continous
% The Signal reconstruction 
% 
   t2=[dt:dt:tf-ts+dt];
% %  time_dynamics_new=zeros(r,length(t));
% % for iter= 1:length(t2) 
% %  time_dynamics_new(:,iter)=(b.*exp(omega*t2(iter))); 
% % end
% %    X_dmd=Phi(1:r,:)*time_dynamics_new;
% 
%   

tic
  X_dmd = Phi*diag(b)*exp(omega*t2);
sDMD_reconsrection=toc;
time.sDMD_reconsrection=sDMD_reconsrection;
tic
  X_rdmd = Phir*diag(br)*exp(omegar*t2);
rDMD_reconsrection=toc;
time.rDMD_reconsrection=rDMD_reconsrection;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct the signals.
% % if we consider discrete-time :Works
% for iter= 1:length(t)
% % % X_dmd(:,iter)=Phi(1:r,:)*diag(lambda).^(iter-1)*b;
% % % X_dmd(:,iter)=Phi(1:r,:)*diag(lambda)^(iter-1)*b;
% % % X_dmd(:,iter)=Phi(1:r,:)*diag(b)*lambda.^(iter-1);
% LAMBDA(:,iter)=lambda.^(iter-1);
% end
% 
% X_dmd=Phi*diag(b)*LAMBDA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% ERA

 def=1;
 tic
 [X_ERA,A,B,C,D, eig_s, residue1, sysdisc,HSVs]=fun_mera(y,r,dt,def);
 ERA_time=toc;
time.ERA=ERA_time;
%% MP

tic
[X_MP,eig_s]=fun_msignal_MP(y', r, dt);
MP_time=toc;
time.MP=MP_time;
%% Prony

np=33;
tic
[eig_s, X_Prony] = fun_prony_rankreduction(y', dt, r, np);

Prony_time=toc;

time.Prony=Prony_time;


time
%% plots  reconstruction signals
figure(1)
 mycolor=[
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
          0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500];


for k=1:ss%2 %ss
    plot(t(1:length(y)),real(y(k,:)),'color',[mycolor(1,:)],'Linewidth',[1.2]) ,
% plot(t(1:length(y)-1),real(y(k,1:end-1))) ,
   hold on, 
plot(t,real(X_dmd(k,:)),'k-.','Linewidth',[1.2]) 
 hold on, 
    plot(t,real(X_rdmd(k,:)),'r-.','Linewidth',[1.2]) 
 hold on, 
    plot(t,real(X_MP(:,k)),'-.','color',[mycolor(4,:)],'Linewidth',[1.2]) 
 hold on, 
    plot(t,real(X_ERA(:,k)),'-.','color',[mycolor(5,:)],'Linewidth',[1.2]) 
 hold on, 
    plot(t,real(X_Prony(:,k)),'-.','color',[mycolor(3,:)],'Linewidth',[1.2]) 
if k==ss
   legend('Signals','DMD' ,'rDMD','MP','ERA','Prony')
end
end
%  legend('Measurements','DMD')
 axis([51.5 65.5 202.9 206.1])
title('Event 3: Algorithm Comparison')
xlabel('Time (s)'), ylabel('Voltage Magnitude (kV)')

%% Error


% area(years,expenses,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k')

Error_dmd=100*sum(abs(y-X_dmd(1:ss,:)))./sum(y);
Error_rdmd=100*sum(abs(y-X_rdmd(1:ss,:)))./sum(y);
Error_mp=100*sum(abs(y-X_MP'))./sum(y);
Error_era=100*sum(abs(y-X_ERA'))./sum(y);
Error_Prony=100*sum(abs(y-X_Prony'))./sum(y);

figure(5)


 area(t,Error_Prony)
   hold on, 
 area(t,Error_mp)
  hold on, 
 area(t,Error_era)
  hold on, 
 area(t,Error_rdmd)
  hold on, 
 area(t,Error_dmd)
 legend('Prony','MP','ERA','rDMD','DMD')
%  legend('Measurements','DMD')
  axis([51 66 0 0.03])
title('Event 3: The Difference of Actual and Reconstructed Signals')
xlabel('Time (s)'), ylabel('Error (%)')
alpha(.75)
%% plot signals seperetrly in a figure
% figure(100)
% 
% for k=1:ss%2 %ss
%     
% subplot(8,4,k), plot(t(1:length(y)),real(y(k,:)),'color',[mycolor(1,:)],'Linewidth',[1.2]) ,
% hold on, 
% subplot(8,4,k), plot(t,real(X_dmd(k,:)),'k-.','Linewidth',[1.2]) 
%  hold on, 
% subplot(8,4,k), plot(t,real(X_rdmd(k,:)),'r:','Linewidth',[1.2]) 
% end
% title('Event 1 Model Identification: With Shift-Stacking')
% xlabel('Time (s)'), ylabel('Signal Amplitude')


%% plot signals seperetrly in figures

% for k=1:2 %ss
%     figure(100+k)
%  plot(t(1:length(y)),real(y(k,:)),'color',[mycolor(1,:)],'Linewidth',[1.2]) ,
% hold on, 
% plot(t,real(X_dmd(k,:)),'k-.','Linewidth',[1.2]) 
%  hold on, 
% plot(t,real(X_rdmd(k,:)),'r:','Linewidth',[1.2]) 
% 
% % title('Event 1 Model Identification: With Shift-Stacking')
% 
% title(['Event 1, Frequency Signal:',num2str(k)])
% xlabel('Time (s)'), ylabel('Signal Amplitude')
% end

%%  Eigenvalues

%   figure(5)
%  
%   plot(real(omega), imag(omega)/2/pi,'kX','LineWidth',1); hold on;
%   plot(real(omegar), imag(omegar)/2/pi,'ro','LineWidth',1); hold on;
% %     plot(real(Dominant_mode), imag(Dominant_mode)/2/pi, 'bo','LineWidth',2,'MarkerSize',10); hold on;  
%   title('Eigenvalues');
%   xlabel('Damping'), ylabel('Frequency (Hz)');
%   legend('DMD','rDMD');
%   grid on;

%% Singular Values
% figure,
subplot(1,2,1), semilogy(diag(S),'-ok') , hold on, semilogy(diag(S(1:r,1:r)),'or')
xlabel('r'), ylabel('Singular Values');
title('Singular Values');
subplot(1,2,2), plot(cumsum(diag(S))/sum(diag(S)),'-ok'), hold on, plot(cumsum(diag(S(1:r,1:r)))/sum(diag(S)),'or')
xlabel('r'), ylabel('%');
title('Cumulative energy');

%% DMDvsFFT
% 
% bb=abs(b);
% bbr=abs(br);
% 
% 
% 
%  bb(find(bb>100))=0;
% % 
%  bbr(find(bbr>100))=0;
% 
% % bb=bb/max(bb);
% modeHz=imag(omega)/2/pi;
% rmodeHz=imag(omegar)/2/pi;
% 
% figure(6)
% stem(modeHz,2*bb,'k','filled','LineWidth',1)
% hold on
% stem(rmodeHz,2*bbr,'r-.','filled','LineWidth',1)
%  axis([-0.05 5 0 inf])
% % plot(freqs(1:maxFreq),xpower(1:maxFreq),'k','LineWidth',1.2)
% % grid on, hold on
% % scatter(abs(imag(DMDfreqs)),DMDpower,'r','LineWidth',1.2)
% % %    axis([-0.5 5 0 inf])
%   legend('DMD','rDMD')
%   title('DMD Spectrum');
%   xlabel('Frequency (Hz)'), ylabel('|b|');
%   
  
 %% Mode Shapes 
% %  
% %  for Nmode=1:1%:2*3 the dominant mode #
% %  [~,ind_bx]=max( bb);
% %  ind_b(Nmode)=ind_bx;
% %  end
% %    Phi_s = Phi(1:ss,ind_b);   
% 
% modeHz
% ind_b=3; % to get 1.58 Hz
% 
%  Phi_s = Phi(1:ss,ind_b);   
%  Mode = Phi_s; 
% 
%   figure(550)
% title('DMD Mode Shapes')
% h=compass(Mode(:,1));
% for m=1:ss
%   set(h(m),'color',[mycolor(substation(m),:)],'LineWidth',1.5)
% end
% xlabel(['Mode (Hz): ',num2str(imag(round(omega(ind_b)/(2*pi),2)))])
% 
% % textss={substation};
% %  legend('Sub:1:Ln:1','Sub:1:Ln:1-11','Sub:2:Ln:2','Sub:2:Ln:3','Sub:2:Ln:4','Sub:3:Ln:5','Sub:3:Ln:6','Sub:3:Ln:7','Sub:4:Ln:8','Sub:4:Ln:4-6','Sub:4:Ln:9','Sub:5:Ln:10','Sub:5:Ln:11','Sub:5:Ln:12','Sub:6:Ln:4-6','Sub:6:Ln:13','Sub:6:Ln:14','Sub:6:Gen:Gen1','Sub:7:Ln:15','Sub:7:Ln:16','Sub:7:Gen:Gen2','Sub:8:Ln:17','Sub:8:Ln:18','Sub:8:Ln:19','Sub:8:Gen:Gen1','Sub:9:Ln:20','Sub:9:Ln:21','Sub:11:Ln:25','Sub:11:Ln:1-11','Sub:12:Ln:26','Sub:12:Ln:27','Sub:12:Ln:28')
% %   legend('boxoff')
%   legend('sub: 1','sub:2','sub:2','sub:2','sub:3','sub:3','sub:3','sub:4','sub:4','sub:4','sub:5','sub:5','sub:5','sub:6','sub:6','sub:6','sub:6:Gen','sub:7','sub:7','sub:7:Gen2','sub:8','sub:8','sub:8','sub:8:Gen1','sub:9','sub:9','sub:10','sub:10','sub:10','sub:11','sub:11','sub:12','sub:12','sub:12')
%    legend('boxoff')