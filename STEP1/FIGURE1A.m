%%%PLOT FIGURE 1A

%%%Load datasets
load('DFIN.mat') %%PAN CANCER data
%%%DFIN.MUT - Mutation data
%%%DFIN.AN - aneuploidy levels
%%%DFIN.OS/death - overall survival/death
%%%DFIN.DFI/DFIIND - disease free interval/death
%%%DFIN.DSS/DFSSND - disease stable survival/death

load('UC1.mat') %% Cancer types order
load('PANDRIVE.mat') %% list of pan-cancer drivers


canc=DFIN.canc;
OS = DFIN.OS;
death = DFIN.death;
AN = DFIN.AN;



[a,b,c] = intersect(DFIN.gene,PANDRIVE);
DRC = sum(DFIN.MUT(b,:));
MDR = median(DRC); %% Median #drivers 


ps1=[];dAUC1=[];DDD=[];IND1=[];
for i = 1:length(UC1)
    cc1=strcmp(canc,UC1(i));
    OS1 = DFIN.OS(cc1==1);
    AN1=AN(cc1==1);
    
    death1 = DFIN.death(cc1==1);
    DRC1 = DRC(cc1==1);
    MDD(i) = median(DRC1);
    
    
    
    si1 = DRC1>=MDR;
    if nnz(DRC1>MDR)>0
        si1 = DRC1>MDR;
    end
    
    
    x1 = [OS1(si1==1),death1(si1==1)];
    x2 = [OS1(si1==0),death1(si1==0)];
    x1(isnan(x1))=0;
    x2(isnan(x2))=0;
    
    [ps1(i),dAUC1(i)] = logrank2(x1,x2,'High drivers','Low drivers'); %% KM survival based on driver load
    
    DDD=[DDD,DRC1];%%drive load
    IND1=[IND1,i*ones(size(DRC1))];%%cancer index
    
    [rh1(i),p1(i)] = corr(AN1,DRC1','type','Spearman'); %%Correlation of drivers and aneuploidy levles
    
    CN(IND1==i) = UC1(i);

end

%%%%Plot figures
subplot(2,1,1)
boxplot(DDD,CN,'Widths',0.8,'Symbol','go','Jitter',1)
ylabel('Driver mutations load')
xx=[dAUC1;rh1]
subplot(2,1,2)
bar(xx',1.1)
ylim([-0.85,0.85])
legend('Kaplan-Meier dAUC','Spearman rho correlation of drivers and aneuploidy')

tt={};
xp1=find(ps1<=0.1);
yp1=dAUC1(xp1);
yp1(yp1>=0) = yp1(yp1>=0)+0.025;
yp1(yp1<0) = yp1(yp1<0)-0.125;
xp1=xp1-0.35;
tt(xp1>0) = {'*'};
text(xp1,yp1,tt,'Color','blue','FontSize',25)

tt2={};
xp2=find(p1<=0.1);
yp2=rh1(xp2);
yp2(yp2>=0) = yp2(yp2>=0)+0.025;
yp2(yp2<0) = yp2(yp2<0)-0.125;
tt2(xp2>0) = {'*'};
text(xp2,yp2,tt2,'Color','red','FontSize',25)






