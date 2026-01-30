% --------- Plot results and generate pngs -------------

%% Pre-PD noise (Rates) 
R1_pre = load('results/r2/r2,16-ASK-0.2,S=4,C=0,Rb=1,n1=1,n2=0,n1m=0,n2m=0,P=cO,ps=RRC,a=0.01,Nsp=250,Npi=0,n=2048,L=4,Rs=300,Nr=1,Nb=16,I=250.mat'); 
R2_pre = load('results/r2/r2,16-ASK-1,S=4,C=0,Rb=1,n1=1,n2=0,n1m=0,n2m=0,P=cO,ps=RRC,a=0.01,Nsp=250,Npi=0,n=2048,L=4,Rs=300,Nr=1,Nb=16,I=250.mat'); 

figure(1);
clf(); 
plot(R1_pre.SNRdB+6, R1_pre.IqYX,'r.-'); %Shift by 6dB due to SNR definition 
hold on; 
grid on; 
grid minor; 
plot(R2_pre.SNRdB+6, R2_pre.IqYX,'b.-'); %Shift by 6dB due to SNR definition 
plot(R1_pre.SNRdB, 1/2*log2(1 + 10.^(R1_pre.SNRdB/10)),'k--'); %Coherent capacity
title('Pre-PD Noise')
xlabel('SNR [dB]')
ylabel('Rate [bpcu]')
xlim([-10,34])
ylim([0,4])

%Arrow
X = [24 24+5.6];
Y = [3.53, 3.53];
text(X(1)+0.5,Y(1)-0.15,'5.6 dB','FontSize',7)
plot(X,Y,'ks-','LineWidth',1,'MarkerFaceColor','auto');
legend('16-ASK-0.2','16-ASK-1','Coherent capacity (real)','','Location','northwest')

%Export 
set(gcf, 'Position', [100, 600, 600, 400]); % [left, bottom, width, height] 
saveas(gcf, 'png/prePD_noise-16-ASK-o_rates.png'); 

%% Pre-PD noise (Traces)
figure(3);
clf(); 
targetSNR = 24; %dB
idx_snr = find(abs(R1_pre.SNRdB - (targetSNR - 6)) < 0.5); 

R1_itermat = squeeze(R1_pre.RATEVAMP_iter_ten(1,idx_snr,:,1,:)).'; 
itvec = 1:size(R1_itermat,1); 
R1_itermat_ext = R1_itermat; 
for j = 1:size(R1_itermat,2)
    ii = find(R1_itermat(1:end,j) > 0 ,1,"last"); 
    R1_itermat_ext(ii-1:end,j) = R1_itermat_ext(ii-1,j); 
end 

PC = 5; %Percentile
up = prctile(R1_itermat_ext(:,:).',PC);
dw = prctile(R1_itermat_ext(:,:).',100-PC); 
fill([itvec fliplr(itvec)], ...
     [dw  fliplr(up)], ...
     [0.7 0.7 1], ...        % color
     'EdgeColor','none', ...
     'FaceAlpha',1);      % transparency
hold on; 

plot(R1_itermat_ext,'LineWidth',0.1,'Color', [1 0 0 0.2])
xlim([0,1.5*R1_pre.REQITER_vec(1,idx_snr)]); 
hold on; 
h = plot(squeeze(mean(R1_itermat_ext,2)).','k-','LineWidth',2)
grid on; 
grid minor; 
title(strcat("Traces of first SIC stage for 16-ASK-0.2, ",num2str(round(R1_pre.SNRdB(idx_snr)+6,1))," dB"));
xlabel("Iterations")
ylabel("Rate [bpcu]")
legend(h, 'mean', 'Location','northwest')
ylim([0,4])

%Export 
set(gcf, 'Position', [700, 600, 600, 400]); % [left, bottom, width, height] 
saveas(gcf, 'png/prePD_noise-16-ASK-o=0.2_traces.png'); 


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------

%% Post-PD noise (Rates) 
R1_post = load('results/r2/r2,16-ASK-0.2,S=4,C=0,Rb=0,n1=0,n2=1,n1m=0,n2m=1,P=cO,ps=RRC,a=0.01,Nsp=250,Npi=0,n=2048,L=4,Rs=300,Nr=1,Nb=16,I=250.mat'); 
R2_post = load('results/r2/r2,16-ASK-1,S=4,C=0,Rb=0,n1=0,n2=1,n1m=0,n2m=1,P=cO,ps=RRC,a=0.01,Nsp=250,Npi=0,n=2048,L=4,Rs=300,Nr=1,Nb=16,I=250.mat'); 

figure(2);
clf(); 
plot(R1_post.SNRdB, R1_post.IqYX,'r.-'); 
hold on; 
grid on; 
grid minor; 

plot(R2_post.SNRdB, R2_post.IqYX,'b.-'); 
title('Post-PD Noise')
xlabel('SNR [dB]')
ylabel('Rate [bpcu]')
xlim([-7,16]); 
ylim([0,4])

%Arrow
X = [10 10+2.7];
Y = [3.56, 3.56];
text(X(1)+0.0,Y(1)-0.15,'2.7 dB','FontSize',7)
plot(X,Y,'ks-','LineWidth',1,'MarkerFaceColor','auto');
legend('16-ASK-0.2','16-ASK-1','Location','northwest')


%Export 
set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height] 
saveas(gcf, 'png/postPD_noise-16-ASK-o_rates.png')

%% Pre-PD noise (Traces) 
figure(4);
clf(); 

targetSNR = 10; %dB
idx_snr = find(abs(R1_post.SNRdB - (targetSNR)) < 0.5); 

R1_itermat = squeeze(R1_post.RATEVAMP_iter_ten(1,idx_snr,:,1,:)).'; 
itvec = 1:size(R1_itermat,1); 
R1_itermat_ext = R1_itermat; 
for j = 1:size(R1_itermat,2)
    ii = find(R1_itermat(1:end,j) > 0 ,1,"last"); 
    R1_itermat_ext(ii-1:end,j) = R1_itermat_ext(ii-1,j); 
end 

PC = 5; %Percentile
up = prctile(R1_itermat_ext(:,:).',PC);
dw = prctile(R1_itermat_ext(:,:).',100-PC); 
fill([itvec fliplr(itvec)], ...
     [dw  fliplr(up)], ...
     [0.7 0.7 1], ...        % color
     'EdgeColor','none', ...
     'FaceAlpha',1);      % transparency
hold on; 

plot(R1_itermat_ext,'LineWidth',0.1,'Color', [1 0 0 0.2])
xlim([0,1.5*R1_post.REQITER_vec(1,idx_snr)]); 
hold on; 
h = plot(squeeze(mean(R1_itermat_ext,2)).','k-','LineWidth',2)
grid on; 
grid minor; 
title(strcat("Traces of first SIC stage for 16-ASK-0.2, ",num2str(targetSNR)," dB"));
xlabel("Iterations")
ylabel("Rate [bpcu]")
legend(h, 'mean', 'Location','northwest')
ylim([0,4])

%Export 
set(gcf, 'Position', [700, 100, 600, 400]); % [left, bottom, width, height] 
saveas(gcf, 'png/postPD_noise-16-ASK-o=0.2_traces.png');