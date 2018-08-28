%works under windows envioronment. 
%Mac environment, xlswrite may not work well, and a .csv file will be
%writen, therefore the input filename need to be changed

clear all
close all
clc
[Geno txt]=xlsread('Simulation_genotype_geometric.xlsx');%still use the simulation before
[Pheno txt]=xlsread('Simulation_phenotype_geometric.xlsx');%no need to change this file unless rerun the simulation part

[row1 col1]=size(Geno);
[row2 col2]=size(Pheno);
Equality=mean(Pheno,1)';
Rho=zeros(col1,2);% first column is for Q-s correlation, second column is for Q-s after controlling for fitness.
Med=0.9; % the fitness value chosen for Q-s correlation when controlling fitness.
%This value need to be relavent, because if it is too large or too small, one may not be able to get enough segregants for estimation of s'
for i=1:col1;
	geno=Geno(:,i);
	ind1=find(geno==1);
	ind2=find(geno==0);
	X=mean(Pheno(ind1,:),1)-mean(Pheno(ind2,:),1);
	[r p]=corr(abs(X'),Equality,'type','Spearman');
	Rho(i,1)=r;
end

Result=zeros(col2,col1);

Fraction = 0.4;

for i=1:col2;
	pheno=Pheno(:,i);
	[pheno index]=sort(pheno);%sorted from small to large
	for j=1:col1;
		geno=Geno(:,j);
		genoE=geno(index);
        ind_1=find(genoE==1);
        ind_2=find(genoE==0);
        pheno1=pheno(ind_1);
        pheno2=pheno(ind_2);
        num1 = round(length(ind_1)*Fraction);
        num2 = round(length(ind_2)*Fraction);
        ind_up1=find(pheno(ind_1)<Med);
        ind_up2=find(pheno(ind_2)<Med);
        lenUp1 = length(ind_up1)+round(num1/2);
        lenUp2 = length(ind_up2)+round(num2/2);
        diff1_2=0;
        diff2_1=0;
        if lenUp1 < length(ind_1);
            if lenUp1 > num1;
                if lenUp2 < length(ind_2);
                    if lenUp2 > num2;
                        mean1 = mean(pheno1(lenUp1-num1+1:lenUp1));
                        mean2 = mean(pheno2(lenUp2-num2+1:lenUp2));
                        percentileUp2 = round(lenUp1/length(ind_1)*length(ind_2));
                        percentileMean2 = mean(pheno2(percentileUp2-num2+1:percentileUp2));
                        percentileUp1 = round(lenUp2/length(ind_2)*length(ind_1));
                        percentileMean1 = mean(pheno1(percentileUp1-num1+1:percentileUp1));
                        diff1_2 = mean1 - percentileMean2;
                        diff2_1 = percentileMean1 - mean2;
                    end
                end
            end
        end
        diff=0;
        if diff1_2*diff2_1 ~= 0;
            diff= diff1_2 + diff2_1;
        end
        Result(i,j) = diff;
	end
end

for i=1:col1;
	result=Result(:,i);
	ind=find(result>0);
	[r p]=corr(result(ind),Equality(ind),'type','Spearman');
	Rho(i,2)=r;
end
length(find(Rho(:,2)<0))/100
xlswrite('Qs_from_simulation_geometric_controlFitness.xlsx',Rho,1,'A1');