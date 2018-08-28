clear all
close all
clc

[Pheno txt1]=xlsread('BYxRM_Growth_log.xlsx');
A=importdata('BYxRM_geno_sorted2.txt');
Geno=A.data;
[row column]=size(Geno);
Result=zeros(47,4);
[QTL txt]=xlsread('QTLsForAverageGrowthRate.xlsx');
[FreqPVal txt]= xlsread('alleleFrequencyDeviationResult.xlsx');
PValCutOff = 0.05;
Fraction = 0.20;
row=1;
for i=1:47;
    pheno=Pheno(:,i);
    [pheno index]=sort(pheno);
    in=find(pheno<=1.6);
    pheno(1:length(in))=[];
    index(1:length(in))=[];
    ind=find(isnan(pheno)==1);
    index(ind)=[];
    pheno(ind)=[];
    geno=Geno(:,index);
    qtl=QTL(i,:);
    ind=find(isnan(qtl)==0);
    qtl=qtl(ind);
    Bootstrap=[mean(pheno)];
    for j=1:length(ind);
        numQTLOkay=0;
         if FreqPVal(qtl(j),i) >= PValCutOff;
             numQTLOkay = numQTLOkay + 1;
            genoOri = Geno(qtl(j),:);
            numBY = round(length(find(genoOri == -1))*Fraction);%number of BY genotype used at this SNP
            numRM = round(length(find(genoOri == 1))*Fraction);%number of RM genotype used at this SNP
            geno_j=geno(qtl(j),:);
            ind1=find(geno_j==1); % RM allele
            ind2=find(geno_j==-1); % BY allele
            diff_up=mean(pheno(ind1(end-numBY+1:end)))-mean(pheno(ind2(end-numRM+1:end)));
            diff_low=mean(pheno(ind1(1:numBY)))-mean(pheno(ind2(1:numRM)));
            diff=abs(diff_low)-abs(diff_up);
            Result(i,1)=mean(pheno);
            Diff_up=[];
            Diff_low=[];
            Diff_up_low=[];
            Count=0;
            for k=1:1000;
                diff_up=mean(datasample(pheno(ind1(end-numBY+1:end)),numBY))-mean(datasample(pheno(ind2(end-numRM+1:end)),numRM));
                diff_low=mean(datasample(pheno(ind1(1:numBY)),numBY))-mean(datasample(pheno(ind2(1:numRM)),numRM));
                Diff_up=[Diff_up;diff_up];
                Diff_low=[Diff_low;diff_low];
                if abs(diff_low)<abs(diff_up);
                    Count=Count+1;
                end
                %Diff_up_low=[Diff_up_low,abs(diff_low)-abs(diff_up)];
            end
            Bootstrap=[Bootstrap,Count/1000];
         end
    end
    if numQTLOkay > 0;
        xlswrite('QTLPvalueByBootStrap.xlsx',Bootstrap,1,strcat('A',num2str(row)));
        row=row+1;
    end
end

