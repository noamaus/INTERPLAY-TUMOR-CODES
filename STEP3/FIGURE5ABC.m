%%FIGURE 5 A,B and C

load('DFIN11.mat') %% load gastrointestinal and endometrial data with MSI
load('SIG3.mat') %%  MSI-Aneuploidy signature



subplot(1,3,1)
%%% 1. Plot KM survival analysis, based on aneuploidy and MSI

x1 = [DFIN11.OS(strcmp(DFIN11.MSI,'MSS')&DFIN11.AN>8),DFIN11.death(strcmp(DFIN11.MSI,'MSS')&DFIN11.AN>8)]; %% MSS aneuploid
x2 = [DFIN11.OS(strcmp(DFIN11.MSI,'MSI-H')),DFIN11.death(strcmp(DFIN11.MSI,'MSI-H'))]; %% MSI 
x3 = [DFIN11.OS(strcmp(DFIN11.MSI,'MSS')&DFIN11.AN<=8),DFIN11.death(strcmp(DFIN11.MSI,'MSS')&DFIN11.AN<=8)]; %% MSS dyploid

x1(isnan(x1))=0;
x2(isnan(x2))=0;
x3(isnan(x3))=0;

[ppp] = logrank(x1,x2,'Aneuploidy','MSI','r','b');
text(0,0,['rbP = ',num2str(ppp)])
[ppp] = logrank(x1,x3,'Aneuploidy','- -','r','g');
text(0,0.1,['rgP = ',num2str(ppp)])
hold on
[ppp] = logrank(x2,x3,'MSI','- -','b','g');
text(0,0.2,['bgP = ',num2str(ppp)])

text(4000,0.8,['Bn = ',num2str(length(x2))])
text(4000,0.9,['Rn = ',num2str(length(x1))])
text(4000,1,['Gn = ',num2str(length(x3))])

% 


s11=(sum(DFIN11.MUT(SIG3,:)));


subplot(1,3,2)
%%% 2. Plot KM survival analysis, based on MSI and aneuploid-MSI signature

x1 = [DFIN11.OS(strcmp(DFIN11.MSI,'MSS')&s11'>0),DFIN11.death(strcmp(DFIN11.MSI,'MSS')&s11'>0)];
x2 = [DFIN11.OS(strcmp(DFIN11.MSI,'MSI-H')&s11'>0),DFIN11.death(strcmp(DFIN11.MSI,'MSI-H')&s11'>0)];
x3 = [DFIN11.OS(strcmp(DFIN11.MSI,'MSS')&s11'<=0),DFIN11.death(strcmp(DFIN11.MSI,'MSS')&s11'<=0)];




x1(isnan(x1))=0;
x2(isnan(x2))=0;
x3(isnan(x3))=0;

[ppp] = logrank(x1,x2,'MSS&highSignature','MSI&highSignature','r','b');
text(0,0,['rbP = ',num2str(ppp)])

[ppp] = logrank(x1,x3,'MSS&highSignature','MSS&lowSignature','r','g');
text(0,0.1,['rgP = ',num2str(ppp)])

[ppp] = logrank(x2,x3,'MSI&highSignature','MSS&lowSignature','b','g');
text(0,0.2,['bgP = ',num2str(ppp)])


text(4000,0.8,['Bn = ',num2str(length(x2))])
text(4000,0.9,['Rn = ',num2str(length(x1))])
text(4000,1,['Gn = ',num2str(length(x3))])






subplot(1,3,3)

AN1=DFIN11.AN;

%%% 3. Plot KM survival analysis, based on aneuploidy and aneuploid-MSI signature
x1 = [DFIN11.OS(DFIN11.AN>8&s11'<=0),DFIN11.death(DFIN11.AN>8&s11'<=0)];%%AN
x2 = [DFIN11.OS(DFIN11.AN<=8&s11'>0),DFIN11.death(DFIN11.AN<=8&s11'>0)];%%SIG
x3 = [DFIN11.OS(DFIN11.AN<=8&s11'<=0),DFIN11.death(DFIN11.AN<=8&s11'<=0)];%%NON




x1(isnan(x1))=0;
x2(isnan(x2))=0;
x3(isnan(x3))=0;

[ppp] = logrank(x1,x2,'AN','SIG','r','b');
text(0,0,['rbP = ',num2str(ppp)])

[ppp] = logrank(x1,x3,'AN','NON','r','g');
text(0,0.1,['rgP = ',num2str(ppp)])

[ppp] = logrank(x2,x3,'SIG','NON','b','g');
text(0,0.2,['bgP = ',num2str(ppp)])

text(4000,0.8,['Bn = ',num2str(length(x2))])
text(4000,0.9,['Rn = ',num2str(length(x1))])
text(4000,1,['Gn = ',num2str(length(x3))])
% % 
