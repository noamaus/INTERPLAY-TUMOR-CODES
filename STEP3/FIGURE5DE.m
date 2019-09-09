%%% PLOT FIGURE 5D and E
load('DFIN11.mat') %% load gastrointestinal and endometrial data with MSI
load('SIG3.mat') %%  MSI-Aneuploidy signature
load('DRUG.mat') %% Drug response data
load('drugcons.mat') %% considered drugs (enough samples size)



MSI2=strcmp(DFIN11.MSI,'MSI-H');
an=DFIN11.AN;


X1=DFIN11.MUT(SIG3,:);
s11=(sum(DFIN11.MUT(SIG3,:)));

r = strcmp(DRUG.response,'Complete Response')|strcmp(DRUG.response,'Partial Response');

pp=[];
M1=[];M2=[];N2=[];AUC=[];c1=[];c2=[];PRR=[];PRR2=[];PRR3=[];PRR4=[];
clf
PRR=[];PRR2=[];PRR4=[];
%%% For each considered drug, find associations of the aneuploidy-MSI
%%% signatures with response, and of MSI and aneuploidy levels with
%%% response
for i = 1:length(drugcons)
    
    sm = DRUG.patient(strcmp(DRUG.drugn,drugcons(i)));
    r1=DRUG.response(strcmp(DRUG.drugn,drugcons(i)));
    r2=r(strcmp(DRUG.drugn,drugcons(i)));
    
    [a,b,c] = intersect(DFIN11.s,sm);
    s1=s11(b);
    a1=an(b);
    ms=MSI2(b);
    X2=X1(:,b);
    rsp = r1(c);
    r22 = r2(c);
    c1(i)=length(rsp);
    c2(i)=nnz(r22);
    
    [X,Y,T,AUC(i)] = perfcurve(double(r22),double(s1),1);
    
    
    [pp(i)] = ranksum(s1(r22==0),s1(r22>0),'tail','left');

    PRR(i,1) = nnz(r22&ms==1)/length(r22&ms==1); %%%MSI
    PRR(i,2) = nnz(r22&ms==0)/length(r22&ms==0); %%%MSS
    
    
    PRR2(i,1) = nnz(r22&a1>11)/length(r22&a1>11);
    PRR2(i,2) = nnz(r22&a1<=11)/length(r22&a1<=11);
    
    
    PRR3(i,1) = nnz(r22&s1'==0)/length(r22&s1'==0);
    PRR3(i,2) = nnz(r22&s1'>0)/length(r22&s1'>0);
    
    
    PRR4(i,1) = nnz(r22)/length(r22);
    
     PRR5(i,1) = nnz(r22&ms==0&s1'>0)/length(r22&ms==0&s1'>0); %%%MSI
     PRR5(i,2) = nnz(r22&ms==0&s1'==0)/length(r22&ms==0&s1'==0); %%%MSI
    
    
    PRR6(i,1) = nnz(r22&a1<=11&s1'>0)/length(r22&a1<=11&s1'>0);
    PRR6(i,2) = nnz(r22&a1<=11&s1'==0)/length(r22&a1<=11&s1'==1);
 
    subplot(1,6,i)
    boxplot(s1,r22)
    text([1,2],[10,10],{num2str(nnz(r22==0)),num2str(nnz(r22==1))})
    text(0.5,8,['p=',num2str(pp(i))])
    ylim([-0.5,12])
 
end

figure

load('mapG.mat')
subplot(1,3,1)
h1=heatmap({'MSI','MSS'},drugcons,round(PRR*1000)/10);
h1.Colormap=mapG
caxis([5,65])

subplot(1,3,2)
h1=heatmap({'Aneuploid','Diploid'},drugcons,round(PRR2*1000)/10);
h1.Colormap=mapG
caxis([5,65])



subplot(1,3,3)
h1=heatmap({'All'},drugcons,round(PRR4*1000)/10);
h1.Colormap=mapG
caxis([5,65])
