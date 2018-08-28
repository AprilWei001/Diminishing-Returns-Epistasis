clear all
close all
clc
% step 1 in mapping
roundNum = 1; % use 1 if it is the first round
Round = strcat('round',num2str(roundNum)); 
[pheno,txt]=xlsread('BYxRM_Growth_log.xlsx');%input phenotype, change this name to i-th round of residual in (i+1)-th round of mapping
geno=importdata('BYxRM_geno_sorted2.txt');%input genotype
OutfileName = strcat('P_ttest_',Round,'.xlsx');%outputfile
star=1;
[row,column]=size(pheno);
loci=size(geno.data,1);
Result=zeros(loci,column);

for i=1:column;
    for j=1:loci;
        List0=[];
        List1=[];
        for k=1:row;
            if geno.data(j,k)==1;
                if isnan(pheno(k,i))==0;
                    List0=[List0,pheno(k,i)];
                end
            elseif isnan(pheno(k,i))==0;
                List1=[List1,pheno(k,i)];
            end
        end
        [h,p,ci,stats]=ttest2(List0,List1);
        Result(j,i)=p;
    end
end

sheet=1;
xlRange='A2';
xlswrite(OutfileName,geno.textdata(2:end,1),sheet,xlRange);
xlRange='A1';
xlswrite(OutfileName,txt(1,:),sheet,xlRange);
xlRange='B2';
xlswrite(OutfileName,Result,sheet,xlRange);

[p_value,txt]=xlsread(OutfileName);


[row,column]=size(p_value);
FDR=zeros(3,column);  %######   FDR=1%   and  FDR=5%
for i=1:column;
    all=sort(p_value(:,i));
    index=find(all>=0.5);
    num=2*length(index);
    bin=20;
    while(bin<row);
        sig=find(all<1/bin);
        percent=num/length(sig)/bin;
        if percent <=0.1;
            if FDR(1,i)==0;
                FDR(1,i)=1/bin;
            end
        end
        if percent <=0.05;
            if FDR(2,i)==0;
                FDR(2,i)=1/bin;
            end
        end
        if percent <=0.01;
            FDR(3,i)=1/bin;
            break
        end
        bin=bin+1;
    end
end
FDRFileName = strcat('FDR_',Round,'.xlsx');
sheet=1;
xlRange='A1';
xlswrite(FDRFileName,txt(1,:),sheet,xlRange);
xlRange='A2';
xlswrite(FDRFileName,{'FDR=10%';'FDR=5%';'FDR=1%'},sheet,xlRange);
xlRange='B2';
xlswrite(FDRFileName,FDR,sheet,xlRange);



[p_value,txt]=xlsread(OutfileName);

[FDR,txtFDR]=xlsread(FDRFileName);

[row,column]=size(p_value);
[row2,column]=size(FDR);

loci=txt(2:end,1);
chr='';
location=[];
for j=1:row;
    locus=char(loci(j));
    locus=strsplit(locus,'_');
    chr_next=char(locus(2));
    tf = strcmp(chr,chr_next);
    if tf==0;
        location=[location,j];
        chr=chr_next;
    end
end
location=[location,row+1];


outFileMarker = strcat('QTL_candidate_',Round,'.xlsx');
for r=1:row2;
    sheet=r;
    xlRange='A1';
    xlswrite(outFileMarker,txt(1,2:end)',sheet,xlRange);
    for i=1:column;
        cut=FDR(r,i);
        p_condition=p_value(:,i);
        Position=[];
        for j=1:16;
            p_chr=p_condition(location(j):location(j+1)-1);
            [M,I]=min(p_chr);
            if M<=cut;
                sig=find(p_chr==M);
                if length(sig)==1;
                    position=I+location(j)-1; % marker
                else
                    position=round(median(sig))+location(j)-1;
                end
                Position=[Position,position];
            end
        end
        xlRange=strcat('B',int2str(i));
        if length(Position)>0;
            xlswrite(outFileMarker,Position,sheet,xlRange);
        end
    end
    
end

[pheno,txt]=xlsread('BYxRM_Growth_log.xlsx');
geno=importdata('BYxRM_geno_sorted2.txt');
genotype=geno.data;

[row,column]=size(pheno);
total_loci=size(genotype,1);
if roundNum == 1
outFileMarkerFilt = strcat('Marker_Filt_',Round,'.xlsx');
for sheet=1:3;
    [Marker,txt2]=xlsread(outFileMarker,sheet);
    xlRange='A1';
    
    xlswrite(outFileMarkerFilt,txt2(:,1),sheet,xlRange);
    for i=1:column;
        marker=[];
        for j=1:size(Marker,2);
            if isnan(Marker(i,j))==0;
                marker=[marker,Marker(i,j)];
            else
                break
            end
        end
        Y=[];
        Y_index=[];
        for k=1:row;
            if pheno(k,i)>1.6;
                Y=[Y;pheno(k,i)];
                Y_index=[Y_index;k];
            end
        end
        X=[];
        for marker_No=1:length(marker);
            j=marker(marker_No);
            x=[];
            for k=1:length(Y_index);
                x=[x,genotype(j,Y_index(k))];
            end
            X=[X,x'];
        end
        if length(marker)>0;
            mdl=LinearModel.fit(X,Y);
            P=mdl.Coefficients.pValue(2:end);
        end
        marker_new=[];
        for marker_No=1:length(marker);
            if P(marker_No)<=0.05;
                marker_new=[marker_new,marker(marker_No)];
            else
                sheet
            end
        end
        xlRange=strcat('B',int2str(i));
        if length(marker_new)>0;
            xlswrite(outFileMarkerFilt,marker_new,sheet,xlRange);
        end
    end
end


else
    sheet = 1;
    outFileMarkerFiltPrev = strcat('Marker_Filt_round',num2str(roundNum -1),'.xlsx');;
    [Marker,txt2]=xlsread(outFileMarkerFiltPrev,sheet);
    [Marker2,txt3]=xlsread(outFileMarker,sheet);
	Marker2=[nan(8,size(Marker2,2));Marker2;nan(5,size(Marker2,2))];
    xlRange='A1';
    xlswrite(outFileMarkerFilt,txt2(:,1),1,xlRange);
    for i=1:column;
        marker=[];
        for j=1:size(Marker,2);
            if isnan(Marker(i,j))==0;
                marker=[marker,Marker(i,j)];
            else
                break
            end
        end
        for j=1:size(Marker2,2);
            if i<=size(Marker2,1);
                if isnan(Marker2(i,j))==0;
                    marker=[marker,Marker2(i,j)];
                else
                    break
                end
            end
        end
        Y=[];
        Y_index=[];
        for k=1:row;
            if isnan(pheno(k,i))==0;
                Y=[Y;pheno(k,i)];
                Y_index=[Y_index;k];
            end
        end
        X=[];
        for marker_No=1:length(marker);
            j=marker(marker_No);
            x=[];
            for k=1:length(Y_index);
                x=[x,genotype(j,Y_index(k))];
            end
            X=[X,x'];
        end
        if length(marker)>0;
            mdl=LinearModel.fit(X,Y);
            P=mdl.Coefficients.pValue(2:end);
        end
        marker_new=[];
        for marker_No=1:length(marker);
            if P(marker_No)<=0.05;
                marker_new=[marker_new,marker(marker_No)];
            else
                sheet
            end
        end
        xlRange=strcat('B',int2str(i));
        if length(marker_new)>0;
            xlswrite(outFileMarkerFilt,marker_new,1,xlRange);
        end
    end
end

[pheno,txt]=xlsread('BYxRM_Growth_log.xlsx');%no need to change
genotype=geno.data;
txt2=geno.textdata;
total_loci=size(genotype,1);
outResidualFile = strcat('Residuals_',Round,'.xlsx');
sheet=2;
[Marker,txt3]=xlsread(outFileMarkerFilt,sheet);
xlRange='A1';
xlswrite(outResidualFile,txt(:,1),1,xlRange);
xlRange='B1';
xlswrite(outResidualFile,txt(1,2:end),1,xlRange);

for i=1:column;
    marker=[];
    for j=1:size(Marker,2);
        if isnan(Marker(i,j))==0;
            marker=[marker,Marker(i,j)];
        else
            break
        end
    end
    if length(marker)>0;
        Y=[];
        Y_index=[];
        Geno_with_pheno=[];
        for k=1:row;
            if pheno(k,i)>1.6;
                Y=[Y;pheno(k,i)];
                Y_index=[Y_index;k];
                Geno_with_pheno=[Geno_with_pheno,genotype(:,k)];
            else
                pheno(k,i)=NaN(1);
            end
        end
        X=[];
        for marker_No=1:length(marker);
            j=marker(marker_No);
            x=Geno_with_pheno(j,:);
            X=[X,x'];
        end
        mdl=LinearModel.fit(X,Y);
        Y=mdl.Residuals.Raw;
        for k=1:length(Y_index);
            pheno(Y_index(k),i)=Y(k);
        end
    end
end

xlRange='B2';
xlswrite(outResidualFile,pheno,1,xlRange);