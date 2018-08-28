clear all
close all
clc
all_P=[];
Sample=1000;
Mod=10;
numGenes = 10;
Gene_per_mod=ones(1,Mod)*numGenes;
No_environment=50;
E_effect=repmat([0.00:0.007:0.345]'*0.3+0.2,1,Mod)+normrnd(0,0.05,[No_environment,numGenes]);
Total_gene=sum(Gene_per_mod);
Gene_effect=repmat([11:20],1,numGenes)*0.008;
Genotype=randi([0,1],Sample,Total_gene);
Phenotype=zeros(Sample,No_environment);
for i=1:Sample;
    geno=Genotype(i,:);
    mod_eff=zeros(1,Mod);
    l=1;
    for j=1:Mod;
        no_mod=Gene_per_mod(j);
        for k=1:no_mod;
            mod_eff(j)=mod_eff(j)+geno(l)*Gene_effect(l);
            l=l+1;
        end
    end
    for j=1:No_environment;
        raw_pheno=E_effect(j,:)+mod_eff;
        ind=find(raw_pheno>1);
        raw_pheno(ind)=1;
        Phenotype(i,j)=sum(raw_pheno)/Mod+normrnd(0,0.008);
    end
end
xlswrite('Simulation_phenotype_additive.xlsx',Phenotype,1,'A1');
xlswrite('Simulation_genotype_additive.xlsx',Genotype,1,'A1');

Genotype=Genotype';
all_P=[all_P,p_mean_var];
Result_broad_narrow=zeros(No_environment,4);% result in the order of Q value, the g value when using absolute sH and sL
% the g value when considering only sH and sL are of the same sign, and the
% value when sH and sL are of the same sign, but sL < sH
for i=1:No_environment;
    pheno=Phenotype(:,i);
    [pheno index]=sort(pheno);
    geno=Genotype(:,index);
    for j=1:Total_gene;
        geno_j=geno(j,:);
        ind1=find(geno_j==1);
        ind2=find(geno_j==0);
        diff_up=mean(pheno(ind1(end-99:end)))-mean(pheno(ind2(end-99:end)));
        diff_low=mean(pheno(ind1(1:100)))-mean(pheno(ind2(1:100)));
        if diff_low<0;
            diff_low=-diff_low;
            diff_up=-diff_up;
        end
        if diff_low-diff_up>0;
            Result_broad_narrow(i,2)=Result_broad_narrow(i,2)+1;
        end
        if diff_up>0;
            if diff_low-diff_up>0;
                Result_broad_narrow(i,3)=Result_broad_narrow(i,3)+1;
            end
        end
        if diff_up<0;
            diff_up=-diff_up;
            diff_low=-diff_low;
        end
        if diff_low>0;
            if diff_up-diff_low>0;
                Result_broad_narrow(i,4)=Result_broad_narrow(i,4)+1;
            end
        end
        Result_broad_narrow(i,1)=mean(pheno);
    end
end
xlswrite('gValue_additiveModel.xlsx',Result_broad_narrow,1,'B1');

