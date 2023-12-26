clc
clear all
addpath("E:\")

set(groot,'defaultLineLineWidth',1)

filename_ftir = 'Biomass_Acid_FTIR.xlsx';
D_ftir= xlsread(filename_ftir);
T_ftir=D_ftir(1,2:end);
time_ftir=D_ftir(2,2:end);
X_ftir=D_ftir(3:end,:);
lam_ftir=X_ftir(:,1);
X11_ftir=X_ftir(:,2:end);
D_ftir=2-log10(X11_ftir);
%D_ftir=D_ftir-min(min(D_ftir));


filename_hnmr='Biomass_Acid_HNMR.xlsx';
D_hnmr= xlsread(filename_hnmr);  %data matrix
T_hnmr=D_hnmr(1,2:end); %temperature
time_hnmr=D_hnmr(2,2:end); %time
X_hnmr=D_hnmr(3:end,:);%intensity
conc_hnmr=X_hnmr(:,1);%chemical shift
X11_hnmr=X_hnmr(:,2:end); %intensity
D_hnmr=X11_hnmr;




Dback_ftir=msbackadj(lam_ftir,D_ftir,'STEP',100,'WINDOW',800);
Data_ftir=mssgolay(lam_ftir,Dback_ftir,'DEGREE',2,'Span',3);
Noise_ftir=Dback_ftir-X11_ftir;
Dback_hnmr=msbackadj(conc_hnmr,D_hnmr,'STEP',100,'WINDOW',800);
Data_hnmr=mssgolay(conc_hnmr,Dback_hnmr,'DEGREE',2,'Span',3);
Noise_hnmr=Dback_hnmr-Data_hnmr;

% scaling absorbances
% Data_ftir=(Data_ftir-min(Data_ftir(:)))./(max(Data_ftir(:))-min(Data_ftir(:)));

%% Initialization
global Z_ftir;  global P_ftir; global Z_hnmr; global P_hnmr; global R_f;
global R_h;

Z_ftir=NaN(3,9,1765); 

Z_ftir(1,1:3,:)=Data_ftir(:,1:3)';
Z_ftir(2,4:6,:)=Data_ftir(:,4:6)';
Z_ftir(3,7:9,:)=Data_ftir(:,7:9)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;
P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);    
Z_ftir=tensor(Z_ftir);



Z_hnmr =NaN(3,9,2048);
Z_hnmr(1,1:3,:)=Data_hnmr(:,1:3)';
Z_hnmr(2,4:6,:)=Data_hnmr(:,4:6)';
Z_hnmr(3,7:9,:)=Data_hnmr(:,7:9)';
Z_h=Z_hnmr;
Z_hnmr(isnan(Z_hnmr))=0;
P_hnmr=(Z_hnmr~=0);
P_hnmr=tensor(P_hnmr);
Z_hnmr=tensor(Z_hnmr);

%% Find rank of model
% % Using Concordia to get rank of ftir
%  LOF_F=[];CF_F=[];
% for i=1:8
%     [Factors_ftir,it_f,lof_f,cf_f]=parafac(Z_f,i);
%     %M_f = nmodel(Factors_ftir);
%     LOF_F=[LOF_F;lof_f];
%     cf_f = corcond(Z_f,Factors_ftir);
%     CF_F=[CF_F;cf_f];
% end
% 
% figure()
% subplot(1,2,1)
% plot(1:length(LOF_F),LOF_F,'-BX')
% axis tight
% xlabel('Number of components','fontweight','bold','FontSize',20)
% ylabel('Lack of fit (LOF)','fontweight','bold','FontSize',20)
% set(gca,'FontSize',20,'fontweight','bold')
% subplot(1,2,2)
% plot(1:length(CF_F),CF_F,'-BX')
% axis tight
% xlabel('Number of components','fontweight','bold','FontSize',20)
% ylabel('Core consistency','fontweight','bold','FontSize',20)
% set(gca,'FontSize',20,'fontweight','bold')
% [A_f,B_f,C_f]=fac2let(Factors_ftir);
% M_f = nmodel(Factors_ftir);
% 
% % %Rank has been found to be 4. Check for drop in core consistency
% 
% %
% %Use corcondia to det rank of HNMR
% LOF_H=[];CF_H=[];
% for i=1:8
%     [Factors_hnmr,it_h,lof_h,cf_h]=parafac(Z_h,i);
%     %M_f = nmodel(Factors_ftir);
%     LOF_H=[LOF_H;lof_h];
%     %cf_h = corcond(Z_h,Factors_hnmr);
%     CF_H=[CF_H;cf_h];
% end
% 
% figure()
% subplot(1,2,1)
% plot(1:length(LOF_H),LOF_H,'-BX')
% axis tight
% xlabel('Number of components','fontweight','bold','FontSize',20)
% ylabel('Lack of fit (LOF)','fontweight','bold','FontSize',20)
% set(gca,'FontSize',20,'fontweight','bold')
% subplot(1,2,2)
% plot(1:length(CF_H),CF_H,'-BX')
% axis tight
% xlabel('Number of components','fontweight','bold','FontSize',20)
% ylabel('Core consistency','fontweight','bold','FontSize',20)
% set(gca,'FontSize',20,'fontweight','bold')

 R_h=3;

 R_f=3;






%% Fitting A PARAFAC model. cp_wopt is used to add weighing matrix P

% % Initial guess
% [Ainit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,1)*tenmat(Z_ftir,1)'),Rf,[]);
% [Binit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,2)*tenmat(Z_ftir,2)'),Rf,[]);
% Cinit_ftir=double(tenmat(Z_ftir,3)*khatrirao(Binit_ftir,Ainit_ftir))*pinv((Binit_ftir'*Binit_ftir).*(Ainit_ftir'*Ainit_ftir));
% Minit_ftir={Ainit_ftir;Binit_ftir;Cinit_ftir};
% Minit_ftir = create_guess('Data',Z_ftir, 'Num_Factors', Rf, ...
%      'Factor_Generator', 'nvecs');
% %Call optimizer
SSE_Bioftir=[];M_Bio = struct('tensrs',[]);

for i=1:50
[M_Bioftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir,R_f,'lower',0);%,'init', Minit_ftir);
SSE_Bioftir=[SSE_Bioftir;sum(sum(sum((double((tensor(M_Bioftir))-Z_ftir)).^2)))];
M_Bio(i).tensrs=M_Bioftir;
end

csvwrite('SSE_Bioftir.csv',SSE_Bioftir);

[sse_minf,sse_locf]=min(SSE_Bioftir)

M_Bioftir=M_Bio(sse_locf).tensrs;


csvwrite('M_BioftirU1.csv',M_Bioftir.U{1});
csvwrite('M_BioftirU2.csv',M_Bioftir.U{2});
csvwrite('M_BioftirU3.csv',M_Bioftir.U{3});

%convert to a ttensor
R1 = length(M_Bioftir.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_Bioftir.lambda, repmat(R1,1,ndims(M_Bioftir))); %<-- Create a diagonal core.
Y_ftir = ttensor(core1, M_Bioftir.U) ;%<-- Assemble the ttensor.

%% Fit PARAFAC for Z_hnmr
SSE_Biohnmr=[];M_PEBHNMR = struct('tensrs',[]);
% load('M_DAOftirU1.csv');load('M_DAOftirU2.csv');load('M_DAOftirU3.csv');
% Spec_HNMR=double(tenmat(Z_hnmr,3)*khatrirao(M_DAOhnmrU2,M_DAOhnmrU1))*pinv((M_DAOhnmrU2'*M_DAOhnmrU2).*(M_DAOhnmrU1'*M_DAOhnmrU1));
% Minit_hnmr={M_DAOftirU1;M_DAOftirU2;Spec_HNMR};
for i=1:50
[M_Biohnmr,~,output_hnmr] = cp_wopt(Z_hnmr, P_hnmr, R_h,'lower',0);%,'init', Minit_hnmr);
SSE_Biohnmr=[SSE_Biohnmr;sum(sum(sum((double(tensor(M_Biohnmr)-Z_hnmr)).^2)))];
M_PEBHNMR(i).tensrs=M_Biohnmr;
end

csvwrite('SSE_Biohnmr.csv',SSE_Biohnmr);

[sse_minh,sse_loch]=min(SSE_Biohnmr)

M_Biohnmr=M_PEBHNMR(sse_loch).tensrs;


csvwrite('M_BiohnmrU1.csv',M_Biohnmr.U{1});
csvwrite('M_BiohnmrU2.csv',M_Biohnmr.U{2});
csvwrite('M_BiohnmrU3.csv',M_Biohnmr.U{3});

%convert to a ttensor
R2 = length(M_Biohnmr.lambda);  %<-- Number of factors in X.
core2 = tendiag(M_Biohnmr.lambda, repmat(R2,1,ndims(M_Biohnmr))); %<-- Create a diagonal core.
Y_hnmr = ttensor(core2, M_Biohnmr.U) ;%<-- Assemble the ttensor.
toc

%% Using LBFGSB to decompose the data
load('M_BioftirU1.csv');load('M_BioftirU2.csv');load('M_BioftirU3.csv');
load('M_BiohnmrU1.csv');load('M_BiohnmrU2.csv');load('M_BiohnmrU3.csv');
clc
global F; global H; global FH;
global alpha;global beta; global gamma; global lambda;

alpha=0;
beta=1;
gamma=0.001;
lambda=10;

ftir_train=[Data_ftir(:,2)  Data_ftir(:,4) Data_ftir(:,6)];
%ftir_train=[Data_ftir(:,3)  Data_ftir(:,12)  Data_ftir(:,16) Data_ftir(:,20) Data_ftir(:,24)];
F=cov(ftir_train);

hnmr_train=[Data_hnmr(:,2)  Data_hnmr(:,4)  Data_hnmr(:,6)];
%hnmr_train=[Data_hnmr(:,3)  Data_hnmr(:,10) Data_hnmr(:,14) Data_hnmr(:,17) Data_hnmr(:,20) ];
H=cov(hnmr_train);

fdev=ftir_train-(ones(size(ftir_train,1),size(ftir_train,1))*ftir_train)./size(ftir_train,1);
hdev=hnmr_train-(ones(size(hnmr_train,1),size(hnmr_train,1))*hnmr_train)./size(hnmr_train,1);
FH=(fdev*hdev')./size(hnmr_train',1); 


A=(M_BioftirU1+M_BiohnmrU1)./2;
B=(M_BioftirU2+M_BiohnmrU2)./2;
H1=M_BioftirU3;
H2=M_BiohnmrU3;

x0=[A(:);B(:);H1(:);H2(:)];
l=zeros(length(x0),1);
u=inf.*ones(length(x0),1);
opts.x0=x0;
fn=@(x)func_bio(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr);
gr=@(x)grad_bio(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr);
fabc= @(x)fminunc_wrapper(x,fn,gr); 

[x,fout,info] = lbfgsb(@(x)fabc(x),l,u,opts);
p=3*R_f;q=p+9*R_f;r=q+1765*R_f;s=r+2048*R_f;
Aoptb_Bio=reshape(x(1:p),[3,R_f]);
Boptb_Bio=reshape(x(p+1:q),[9,R_f]);
H1optb_Bio=reshape(x(q+1:r),[1765,R_f]);
H2optb_Bio=reshape(x(r+1:s),[2048,R_f]);

csvwrite('Aoptb_Bio.csv',Aoptb_Bio);
csvwrite('Boptb_Bio.csv',Boptb_Bio);
csvwrite('H1optb_Bio.csv',H1optb_Bio);
csvwrite('H2optb_Bio.csv',H2optb_Bio);


%%
% Plot PARAFAC Decomposition

load('M_BioftirU1.csv');load('M_BioftirU2.csv');load('M_BioftirU3.csv');
load('M_BiohnmrU1.csv');load('M_BiohnmrU2.csv');load('M_BiohnmrU3.csv');


figure()
subplot(2,2,1)
plot(1:length(M_BioftirU1),M_BioftirU1(:,1),'-k')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Temperature condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 1','FontSize',14,'fontweight','bold');


subplot(2,2,2)
plot(1:length(M_BioftirU1),M_BioftirU1(:,2),'-r')

set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Temperature condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 2','FontSize',14,'fontweight','bold');


subplot(2,2,3)
plot(1:length(M_BioftirU1),M_BioftirU1(:,3),'-b')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Temperature condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 3','FontSize',14,'fontweight','bold');

% Residence times

figure()
subplot(2,2,1)
plot(1:length(M_BioftirU2),M_BioftirU2(:,1),'-k')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Residence Time condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 1','FontSize',14,'fontweight','bold');


subplot(2,2,2)
plot(1:length(M_BioftirU2),M_BioftirU2(:,2),'-r')

set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Residence Time condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 2','FontSize',14,'fontweight','bold');


subplot(2,2,3)
plot(1:length(M_BioftirU2),M_BioftirU2(:,3),'-b')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Residence Time condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 3','FontSize',14,'fontweight','bold');

% FTIR
figure()
subplot(2,2,1)
plot(lam_ftir,M_BioftirU3(:,1),'-k')
set(gca,'FontSize',14)
set(gca, 'XDir','reverse')
axis tight
title('FTIR Spectrum of PC_1')
xlabel('Wavenumber (cm^{-1})','FontSize',14,'fontweight','bold')
ylabel('Absorbance','FontSize',14,'fontweight','bold')
subplot(2,2,2)
plot(lam_ftir,M_BioftirU3(:,2),'-r')
set(gca,'FontSize',14)
set(gca, 'XDir','reverse')
axis tight
title('FTIR Spectrum of PC_2')
xlabel('Wavenumber (cm^{-1})','FontSize',14,'fontweight','bold')
ylabel('Absorbance','FontSize',14,'fontweight','bold')
subplot(2,2,3)
plot(lam_ftir,M_BioftirU3(:,3),'-b')
set(gca,'FontSize',14)
set(gca, 'XDir','reverse')
axis tight
title('FTIR Spectrum of PC_3')
xlabel('Wavenumber (cm^{-1})','FontSize',14,'fontweight','bold')
ylabel('Absorbance','FontSize',14,'fontweight','bold')

% HNMR
figure()
subplot(2,2,1)
plot(1:length(M_BiohnmrU1),M_BiohnmrU1(:,1),'-k')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Temperature condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 1','FontSize',14,'fontweight','bold');


subplot(2,2,2)
plot(1:length(M_BiohnmrU1),M_BiohnmrU1(:,2),'-r')

set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Temperature condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 2','FontSize',14,'fontweight','bold');


subplot(2,2,3)
plot(1:length(M_BiohnmrU1),M_BiohnmrU1(:,3),'-b')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Temperature condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 3','FontSize',14,'fontweight','bold');


% Residence times
figure()
subplot(2,2,1)
plot(1:length(M_BiohnmrU2),M_BiohnmrU2(:,1),'-k')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Residence Time condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 1','FontSize',14,'fontweight','bold');


subplot(2,2,2)
plot(1:length(M_BiohnmrU2),M_BiohnmrU2(:,2),'-r')

set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Residence Time condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 2','FontSize',14,'fontweight','bold');


subplot(2,2,3)
plot(1:length(M_BiohnmrU2),M_BiohnmrU2(:,3),'-b')
set(gca,'FontSize',14,'fontweight','bold')
axis tight
xlabel('Residence Time condition','FontSize',14,'fontweight','bold');
ylabel('Pseudo-component concentration 3','FontSize',14,'fontweight','bold');

%Spectra
figure()
subplot(2,2,1)
plot(conc_hnmr,M_BiohnmrU3(:,1),'-k')
set(gca,'FontSize',14)
set(gca, 'XDir','reverse')
axis tight
title('HNMT Spectrum of PC_1')
xlabel('Chemical shift (ppm)','FontSize',14,'fontweight','bold')
ylabel('Absorbance','FontSize',14,'fontweight','bold')
subplot(2,2,2)
plot(conc_hnmr,M_BiohnmrU3(:,2),'-r')
set(gca,'FontSize',14)
set(gca, 'XDir','reverse')
axis tight
title('HNMR Spectrum of PC_2')
xlabel('Chemical shift (ppm)','FontSize',14,'fontweight','bold')
ylabel('Absorbance','FontSize',14,'fontweight','bold')
subplot(2,2,3)
plot(conc_hnmr,M_BiohnmrU3(:,3),'-b')
set(gca,'FontSize',14)
set(gca, 'XDir','reverse')
axis tight
title('HNMR Spectrum of PC_3')
xlabel('Chemicl shift (ppm)','FontSize',14,'fontweight','bold')
ylabel('Absorbance','FontSize',14,'fontweight','bold')

%% LBFGSB
Aoptb_Bio=load ('Aoptb_Bio_SCW.csv');
Boptb_Bio=load ('Boptb_Bio_SCW.csv');
H1optb_Bio=load ('H1optb_Bio_SCW.csv');
H2optb_Bio=load ('H2optb_Bio_SCW.csv');
figure()
subplot(4,1,1)
plot(1:length(Aoptb_Bio),Aoptb_Bio)
axis tight
ylabel('Conc');
xlabel('Temp');
legend('PC1','PC2','PC3','PC4','PC5','PC6','PC7')
subplot(4,1,2)
plot(1:length(Boptb_Bio),Boptb_Bio)
axis tight
ylabel('Conc');
legend('PC1','PC2','PC3','PC4','PC5','PC6','PC7')
subplot(4,1,3)
plot(lam_ftir,H1optb_Bio)
axis tight
ylabel('Absorbance');
xlabel('Wavenumbers');
legend('PC1','PC2','PC3','PC4','PC5','PC6','PC7')
subplot(4,1,4)
plot(conc_hnmr,H2optb_Bio)
axis tight
ylabel('Absorbance');
xlabel('Chem shift');
legend('PC1','PC2','PC3','PC4','PC5','PC6','PC7')


figure()
subplot(2,2,1)
plot(lam_ftir,H1optb_Bio(:,1),'-r','LineWidth',2)
set(gca, 'XDir','reverse')
axis tight
ylabel('Absorbance PC1','FontSize',21);
xlabel('Wavenumber(cm^{-1})','FontSize',21);


subplot(2,2,2)
plot(lam_ftir,H1optb_Bio(:,2),'-b','LineWidth',2)
set(gca, 'XDir','reverse')
axis tight
ylabel('Absorbance PC2','FontSize',21);
xlabel('Wavenumber(cm^{-1})','FontSize',21);


subplot(2,2,3)
plot(lam_ftir,H1optb_Bio(:,3),'-k','LineWidth',2)
set(gca, 'XDir','reverse')
axis tight
ylabel('Absorbance PC3','FontSize',21);
xlabel('Wavenumber(cm^{-1})','FontSize',21);


figure()
subplot(2,2,1)
plot(conc_hnmr,H2optb_Bio(:,1),'-r','LineWidth',2)
set(gca, 'XDir','reverse')
axis tight
ylabel('Absorbance PC1','FontSize',21);
xlabel('Chem shift (ppm)','FontSize',21);


subplot(2,2,2)
plot(conc_hnmr,H2optb_Bio(:,2),'-b','LineWidth',2)
set(gca, 'XDir','reverse')
axis tight
ylabel('Absorbance PC2','FontSize',21);
xlabel('Chem shift (ppm)','FontSize',21);

subplot(2,2,3)
plot(conc_hnmr,H2optb_Bio(:,3),'-k','LineWidth',2)
set(gca, 'XDir','reverse')
axis tight
ylabel('Absorbance PC3','FontSize',21);
xlabel('Chem shift (ppm)','FontSize',21);