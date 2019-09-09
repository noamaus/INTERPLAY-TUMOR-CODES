%%%% Plots Figure 4C

load('DFIN11.mat') %% load gastrointestinal and endometrial data with MSI
load('SIG3.mat') %%  MSI-Aneuploidy signature
load('MGI_COAD.mat') %% MGI colon test data
load('MSK_UCEC.mat') %% MSK endometrial test data
load('DFCI_COAD.mat')   %% DFCI colon test data
load('PFZ_STAD.mat')    %% PFZ stomach test data


sm1=sum(DFIN11.MUT(SIG3,:));
MSI2=strcmp(DFIN11.MSI,'MSI-H');


[X,Y,T,AUCTR] = perfcurve(MSI2,sm1,1);

L11(1)=length(MSI2);
plot(X,Y)









gene1=DFIN11.gene(SIG3);

[a,b,c] = intersect(DFCI_COAD.gene,gene1);


s1 = sum(DFCI_COAD.MUT(b,:));
ll = strcmp(DFCI_COAD.MSI,'MSI-high')
[X,Y,T,AUC1] = perfcurve(ll,s1,1);
L11(2)=length(s1);

hold on
plot(X,Y)
xlabel('False positive rate')
ylabel('True positive rate')



[a,b,c] = intersect(MGI_COAD.gene,gene1);


s1 = sum(MGI_COAD.MUT(b,:));
ll = strcmp(MGI_COAD.MSI,'MSI')
L11(3)=length(ll);

[X,Y,T,AUC2] = perfcurve(ll,s1,1);
hold on
plot(X,Y)







[a,b,c] = intersect(gene1,PFZ_STAD.gene);

s1 = sum(PFZ_STAD.tab(c,:));

ll = double(strcmp(PFZ_STAD.MSI,'high-level microsatellite instabiliy'));

[X,Y,T,AUC3] = perfcurve(ll,s1,1);
L11(4)=length(ll);

plot(X,Y)






[a,b,c] = intersect(gene1,MSK_UCEC.gene);

s1 = sum(MSK_UCEC.tab(c,:));

ll = double(strcmp(MSK_UCEC.MSI,'MSI-H'));
L11(5)=length(ll);

[X,Y,T,AUC4] = perfcurve(ll,s1,1);
plot(X,Y)







xlabel('False positive rate')
ylabel('True positive rate')

legend(['TCGA (training), AUC = ',num2str(AUCTR)],['DFCI COADREAD, AUC = ',num2str(AUC1)],['MGI COADREAD, AUC = ',num2str(AUC2)],['Pfizer STAD, AUC = ',num2str(AUC3)],['MSK UCEC, AUC = ',num2str(AUC4)]);







