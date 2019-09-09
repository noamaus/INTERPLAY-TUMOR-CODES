%%% PLOT FIGURE 5F
load('DFIN11.mat') %% load gastrointestinal and endometrial data with MSI
load('SIG3.mat') %%  MSI-Aneuploidy signature
load('DRUG.mat') %% Drug response data
load('drugcons.mat') %% considered drugs (enough samples size)

MSI2=strcmp(DFIN11.MSI,'MSI-H');
an=DFIN11.AN;

X1=DFIN11.MUT(SIG3,:);
s11=(sum(DFIN11.MUT(SIG3,:)));


r = strcmp(DRUG.response,'Complete Response')|strcmp(DRUG.response,'Partial Response');

AUC2=[];PRR=[];PRR2=[];
M1=[];M2=[];N2=[];AUC=[];c1=[];c2=[];
clf

%%% For each considered drug, find associations of each feature in the aneuploidy-MSI
%%% signatures with response
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
    
    PRR(i,1) = nnz(r22&ms==0&s1>0)/length(r22&ms==0&s1>0);
    PRR(i,2) = nnz(r22&ms==1&s1==0)/length(r22&ms==1&s1==0);
    PRR(i,3) = nnz(r22&ms==0&s1==0)/length(r22&ms==0&s1==0);
    
    
    PRR2(i,1) = nnz(r22&s1>0&a1<=11)/length(r22&s1>0&a1<=11);
    PRR2(i,2) = nnz(r22&a1>11&s1==0)/length(r22&a1>11&s1==0);
    PRR2(i,3) = nnz(r22&s1==0&a1<=11)/length(r22&s1==0&a1<=11);

    for j = 1:length(SIG3)
        [XX,YY,T2,AUC2(i,j)] = perfcurve(double(r22),X2(j,:),1);   
    end
    

end


[svv,sii] = sort(median(AUC2),'descend');

ZZ=round(AUC2*100)/100;
h1=heatmap(DFIN11.gene(SIG3(sii)),drugcons,ZZ(:,sii));
load('map2.mat')
h1.Colormap=map2
caxis([0.3,0.7])

