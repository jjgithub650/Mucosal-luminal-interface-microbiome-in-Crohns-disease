library(phyloseq)
library(vioplot)
library(gridExtra)
library(FSA)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ape)
library(DESeq2)
library(qvalue)
library(caret)
library(pROC)
library(dplyr)
library(RColorBrewer)
library(scales)

#save.image(file='MLI_analysis.RData')
load('MLI_analysis.RData')

### Import 16S data and create analysis subsets

data<-read.csv(choose.files(),header=T,row.names=1)  # ASV count table 
meta<-read.csv(choose.files(),header=T,row.names=1)  # mapping file
tdata=t(data)
length(which(rownames(t(data))==rownames(meta)))

s1<-factor(meta$AdequateSequenceDepth)
split.m <- split.data.frame(meta, s1)
dim(split.m$'Y')
dim(split.m$'N')
newmeta=split.m$'Y'
split.s <- split.data.frame(tdata, s1)
dim(split.s$'Y')
dim(split.s$'N')
newdata=t(as.matrix(split.s$'Y'))
length(which(colnames(newdata)==rownames(newmeta)))

# Add alpha diversity metrics to metadata file
OTU=otu_table(newdata, taxa_are_rows = TRUE)
taxmat=as.matrix(read.csv(choose.files(),header=T,row.names=1))
TAX = tax_table(taxmat)
biom = phyloseq(OTU,TAX)
map=sample_data(newmeta)
data=merge_phyloseq(biom,map)
print(data)
sample_sums(data)
min(sample_sums(data))
rarified = rarefy_even_depth(data,sample.size=50862,replace=FALSE)
sample_sums(rarified)
alphadiversity=estimate_richness(rarified)
length(which(rownames(newmeta)==rownames(alphadiversity)))
newmeta$Chao1<-alphadiversity$Chao1
newmeta$Shannon<-alphadiversity$Shannon

s2<-factor(newmeta$Site)
split.m2 <- split.data.frame(newmeta, s2)
dim(split.m2$'Sigmoid')
dim(split.m2$'Cecum')
newmeta_sigmoid=split.m2$'Sigmoid'
newmeta_cecum=split.m2$'Cecum'
tdata2=t(newdata)
split.s2 <- split.data.frame(tdata2, s2)
dim(split.s2$'Sigmoid')
dim(split.s2$'Cecum')
newdata_sigmoid=t(as.matrix(split.s2$'Sigmoid'))
newdata_cecum=t(as.matrix(split.s2$'Cecum'))
length(which(colnames(newdata_sigmoid)==rownames(newmeta_sigmoid)))
length(which(colnames(newdata_cecum)==rownames(newmeta_cecum)))

s3<-factor(newmeta_sigmoid$Disease)
split.m3 <- split.data.frame(newmeta_sigmoid, s3)
dim(split.m3$'Non-IBD')
dim(split.m3$'CD')
newmeta_sigmoid_nonIBD=split.m3$'Non-IBD'
newmeta_sigmoid_CD=split.m3$'CD'
tdata3=t(newdata_sigmoid)
split.s3 <- split.data.frame(tdata3, s3)
dim(split.s3$'Non-IBD')
dim(split.s3$'CD')
newdata_sigmoid_nonIBD=t(as.matrix(split.s3$'Non-IBD'))
newdata_sigmoid_CD=t(as.matrix(split.s3$'CD'))
length(which(colnames(newdata_sigmoid_nonIBD)==rownames(newmeta_sigmoid_nonIBD)))
length(which(colnames(newdata_sigmoid_CD)==rownames(newmeta_sigmoid_CD)))

s4<-factor(newmeta_cecum$Disease)
split.m4 <- split.data.frame(newmeta_cecum, s4)
dim(split.m4$'Non-IBD')
dim(split.m4$'CD')
newmeta_cecum_nonIBD=split.m4$'Non-IBD'
newmeta_cecum_CD=split.m4$'CD'
tdata4=t(newdata_cecum)
split.s4 <- split.data.frame(tdata4, s4)
dim(split.s4$'Non-IBD')
dim(split.s4$'CD')
newdata_cecum_nonIBD=t(as.matrix(split.s4$'Non-IBD'))
newdata_cecum_CD=t(as.matrix(split.s4$'CD'))
length(which(colnames(newdata_cecum_nonIBD)==rownames(newmeta_cecum_nonIBD)))
length(which(colnames(newdata_cecum_CD)==rownames(newmeta_cecum_CD)))

s7<-factor(newmeta_sigmoid_nonIBD$BMI_known)
split.m7 <- split.data.frame(newmeta_sigmoid_nonIBD, s7)
dim(split.m7$'Y')
dim(split.m7$'N')
newmeta_sigmoid_nonIBD_BMI=split.m7$'Y'
tdata7=t(newdata_sigmoid_nonIBD)
split.s7 <- split.data.frame(tdata7, s7)
dim(split.s7$'Y')
dim(split.s7$'N')
newdata_sigmoid_nonIBD_BMI=t(as.matrix(split.s7$'Y'))
length(which(colnames(newdata_sigmoid_nonIBD_BMI)==rownames(newmeta_sigmoid_nonIBD_BMI)))

s8<-factor(newmeta_cecum_nonIBD$BMI_known)
split.m8 <- split.data.frame(newmeta_cecum_nonIBD, s8)
dim(split.m8$'Y')
dim(split.m8$'N')
newmeta_cecum_nonIBD_BMI=split.m8$'Y'
tdata8=t(newdata_cecum_nonIBD)
split.s8 <- split.data.frame(tdata8, s8)
dim(split.s8$'Y')
dim(split.s8$'N')
newdata_cecum_nonIBD_BMI=t(as.matrix(split.s8$'Y'))
length(which(colnames(newdata_cecum_nonIBD_BMI)==rownames(newmeta_cecum_nonIBD_BMI)))

s9<-factor(newmeta_cecum_nonIBD$Genotype)
split.m9 <- split.data.frame(newmeta_cecum_nonIBD, s9)
dim(split.m9$'Y')
dim(split.m9$'N')
newmeta_cecum_nonIBD_Genotype=split.m9$'Y'
tdata9=t(newdata_cecum_nonIBD)
split.s9 <- split.data.frame(tdata9, s9)
dim(split.s9$'Y')
dim(split.s9$'N')
newdata_cecum_nonIBD_Genotype=t(as.matrix(split.s9$'Y'))
length(which(colnames(newdata_cecum_nonIBD_Genotype)==rownames(newmeta_cecum_nonIBD_Genotype)))

s10<-factor(newmeta_sigmoid_nonIBD$Genotype)
split.m10 <- split.data.frame(newmeta_sigmoid_nonIBD, s10)
dim(split.m10$'Y')
dim(split.m10$'N')
newmeta_sigmoid_nonIBD_Genotype=split.m10$'Y'
tdata10=t(newdata_sigmoid_nonIBD)
split.s10 <- split.data.frame(tdata10, s10)
dim(split.s10$'Y')
dim(split.s10$'N')
newdata_sigmoid_nonIBD_Genotype=t(as.matrix(split.s10$'Y'))
length(which(colnames(newdata_sigmoid_nonIBD_Genotype)==rownames(newmeta_sigmoid_nonIBD_Genotype)))

s11<-factor(newmeta_sigmoid$B_known)
split.m11 <- split.data.frame(newmeta_sigmoid, s11)
dim(split.m11$'Y')
dim(split.m11$'N')
newmeta_sigmoid_CD_B=split.m11$'Y'
tdata11=t(newdata_sigmoid)
split.s11 <- split.data.frame(tdata11, s11)
dim(split.s11$'Y')
dim(split.s11$'N')
newdata_sigmoid_CD_B=t(as.matrix(split.s11$'Y'))
length(which(colnames(newdata_sigmoid_CD_B)==rownames(newmeta_sigmoid_CD_B)))

s12<-factor(newmeta_cecum$B_known)
split.m12 <- split.data.frame(newmeta_cecum, s12)
dim(split.m12$'Y')
dim(split.m12$'N')
newmeta_cecum_CD_B=split.m12$'Y'
tdata12=t(newdata_cecum)
split.s12 <- split.data.frame(tdata12, s12)
dim(split.s12$'Y')
dim(split.s12$'N')
newdata_cecum_CD_B=t(as.matrix(split.s12$'Y'))
length(which(colnames(newdata_cecum_CD_B)==rownames(newmeta_cecum_CD_B)))

s13<-factor(newmeta_cecum$Outcome_B)
split.m13 <- split.data.frame(newmeta_cecum, s13)
dim(split.m13$'Y')
dim(split.m13$'N')
newmeta_cecum_CD_progression=split.m13$'Y'
tdata13=t(newdata_cecum)
split.s13 <- split.data.frame(tdata13, s13)
dim(split.s13$'Y')
dim(split.s13$'N')
newdata_cecum_CD_progression=t(as.matrix(split.s13$'Y'))
length(which(colnames(newdata_cecum_CD_progression)==rownames(newmeta_cecum_CD_progression)))

s14<-factor(newmeta_sigmoid$Outcome_B)
split.m14 <- split.data.frame(newmeta_sigmoid, s14)
dim(split.m14$'Y')
dim(split.m14$'N')
newmeta_sigmoid_CD_progression=split.m14$'Y'
tdata14=t(newdata_sigmoid)
split.s14 <- split.data.frame(tdata14, s14)
dim(split.s14$'Y')
dim(split.s14$'N')
newdata_sigmoid_CD_progression=t(as.matrix(split.s14$'Y'))
length(which(colnames(newdata_sigmoid_CD_progression)==rownames(newmeta_sigmoid_CD_progression)))

s15<-factor(newmeta_cecum$Outcome)
split.m15 <- split.data.frame(newmeta_cecum, s15)
dim(split.m15$'Y')
dim(split.m15$'N')
newmeta_cecum_CD_progression_only=split.m15$'Y'
tdata15=t(newdata_cecum)
split.s15 <- split.data.frame(tdata15, s15)
dim(split.s15$'Y')
dim(split.s15$'N')
newdata_cecum_CD_progression_only=t(as.matrix(split.s15$'Y'))
length(which(colnames(newdata_cecum_CD_progression_only)==rownames(newmeta_cecum_CD_progression_only)))

s16<-factor(newmeta_sigmoid$Outcome)
split.m16 <- split.data.frame(newmeta_sigmoid, s16)
dim(split.m16$'Y')
dim(split.m16$'N')
newmeta_sigmoid_CD_progression_only=split.m16$'Y'
tdata16=t(newdata_sigmoid)
split.s16 <- split.data.frame(tdata16, s16)
dim(split.s16$'Y')
dim(split.s16$'N')
newdata_sigmoid_CD_progression_only=t(as.matrix(split.s16$'Y'))
length(which(colnames(newdata_sigmoid_CD_progression_only)==rownames(newmeta_sigmoid_CD_progression_only)))

s29<-factor(newmeta_cecum_CD$Genotype)
split.m29 <- split.data.frame(newmeta_cecum_CD, s29)
dim(split.m29$'Y')
dim(split.m29$'N')
newmeta_cecum_CD_Genotype=split.m29$'Y'
tdata29=t(newdata_cecum_CD)
split.s29 <- split.data.frame(tdata29, s29)
dim(split.s29$'Y')
dim(split.s29$'N')
newdata_cecum_CD_Genotype=t(as.matrix(split.s29$'Y'))
length(which(colnames(newdata_cecum_CD_Genotype)==rownames(newmeta_cecum_CD_Genotype)))

s30<-factor(newmeta_sigmoid_CD$Genotype)
split.m30 <- split.data.frame(newmeta_sigmoid_CD, s30)
dim(split.m30$'Y')
dim(split.m30$'N')
newmeta_sigmoid_CD_Genotype=split.m30$'Y'
tdata30=t(newdata_sigmoid_CD)
split.s30 <- split.data.frame(tdata30, s30)
dim(split.s30$'Y')
dim(split.s30$'N')
newdata_sigmoid_CD_Genotype=t(as.matrix(split.s30$'Y'))
length(which(colnames(newdata_sigmoid_CD_Genotype)==rownames(newmeta_sigmoid_CD_Genotype)))

s31<-factor(newmeta_cecum_CD_Genotype$B_known)
split.m31 <- split.data.frame(newmeta_cecum_CD_Genotype, s31)
dim(split.m31$'Y')
dim(split.m31$'N')
newmeta_cecum_CD_Genotype_B=split.m31$'Y'
tdata31=t(newdata_cecum_CD_Genotype)
split.s31 <- split.data.frame(tdata31, s31)
dim(split.s31$'Y')
dim(split.s31$'N')
newdata_cecum_CD_Genotype_B=t(as.matrix(split.s31$'Y'))
length(which(colnames(newdata_cecum_CD_Genotype_B)==rownames(newmeta_cecum_CD_Genotype_B)))

s32<-factor(newmeta_sigmoid_CD_Genotype$B_known)
split.m32 <- split.data.frame(newmeta_sigmoid_CD_Genotype, s32)
dim(split.m32$'Y')
dim(split.m32$'N')
newmeta_sigmoid_CD_Genotype_B=split.m32$'Y'
tdata32=t(newdata_sigmoid_CD_Genotype)
split.s32 <- split.data.frame(tdata32, s32)
dim(split.s32$'Y')
dim(split.s32$'N')
newdata_sigmoid_CD_Genotype_B=t(as.matrix(split.s32$'Y'))
length(which(colnames(newdata_sigmoid_CD_Genotype_B)==rownames(newmeta_sigmoid_CD_Genotype_B)))


### Alpha diversity analysis

kruskal.test(Chao1 ~ Disease, data=newmeta_sigmoid)
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + Disease, data = newmeta_sigmoid)  
anova(fit.data)
coef(fit.data)
kruskal.test(Chao1 ~ Disease, data=newmeta_cecum)
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + Disease, data = newmeta_cecum)  
anova(fit.data)
coef(fit.data)

kruskal.test(Shannon ~ Disease, data=newmeta_sigmoid)
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + Disease, data = newmeta_sigmoid)  
anova(fit.data)
coef(fit.data)
kruskal.test(Shannon ~ Disease, data=newmeta_cecum)
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + Disease, data = newmeta_cecum)  
anova(fit.data)
coef(fit.data)

par(mfrow = c(1,2))
newmeta_sigmoid$Disease_ordered<-factor(newmeta_sigmoid$Disease,levels=c("Non-IBD","CD"))
p1<-ggplot(data=newmeta_sigmoid,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_sigmoid,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

par(mfrow = c(1,2))
newmeta_cecum$Disease_ordered<-factor(newmeta_cecum$Disease,levels=c("Non-IBD","CD"))
p1<-ggplot(data=newmeta_cecum,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

kruskal.test(Chao1 ~ B_combine, data=newmeta_sigmoid_CD_B)
dunnTest(Chao1 ~ B_combine, data=newmeta_sigmoid_CD_B,method="bh") 
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_sigmoid_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")

kruskal.test(Chao1 ~ B_combine, data=newmeta_cecum_CD_B)
dunnTest(Chao1 ~ B_combine, data=newmeta_cecum_CD_B,method="bh") 
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_cecum_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")

kruskal.test(Shannon ~ B_combine, data=newmeta_sigmoid_CD_B)
dunnTest(Shannon ~ B_combine, data=newmeta_sigmoid_CD_B,method="bh") 
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_sigmoid_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")

kruskal.test(Shannon ~ B_combine, data=newmeta_cecum_CD_B)
dunnTest(Shannon ~ B_combine, data=newmeta_cecum_CD_B,method="bh")
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_cecum_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")

par(mfrow = c(1,2))
newmeta_sigmoid_CD_B$Disease_ordered<-factor(newmeta_sigmoid_CD_B$B_combine,levels=c("B1","B2_B3"))
p1<-ggplot(data=newmeta_sigmoid_CD_B,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_sigmoid_CD_B,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

par(mfrow = c(1,2))
newmeta_cecum_CD_B$Disease_ordered<-factor(newmeta_cecum_CD_B$B_combine,levels=c("B1","B2_B3"))
p1<-ggplot(data=newmeta_cecum_CD_B,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum_CD_B,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

kruskal.test(Chao1 ~ Obesity, data=newmeta_sigmoid_nonIBD_BMI)
dunnTest(Chao1 ~ Obesity, data=newmeta_sigmoid_nonIBD_BMI,method="bh") 
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity, data = newmeta_sigmoid_nonIBD_BMI)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Obesity")

kruskal.test(Chao1 ~ Obesity, data=newmeta_cecum_nonIBD_BMI)
dunnTest(Chao1 ~ Obesity, data=newmeta_cecum_nonIBD_BMI,method="bh") 
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity, data = newmeta_cecum_nonIBD_BMI)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Obesity")

kruskal.test(Shannon ~ Obesity, data=newmeta_sigmoid_nonIBD_BMI)
dunnTest(Shannon ~ Obesity, data=newmeta_sigmoid_nonIBD_BMI,method="bh") 
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity, data = newmeta_sigmoid_nonIBD_BMI)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Obesity")

kruskal.test(Shannon ~ Obesity, data=newmeta_cecum_nonIBD_BMI)
dunnTest(Shannon ~ Obesity, data=newmeta_cecum_CD_nonIBD_BMI,method="bh")
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity, data = newmeta_cecum_nonIBD_BMI)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Obesity")

par(mfrow = c(1,2))
newmeta_sigmoid_nonIBD_BMI$Disease_ordered<-factor(newmeta_sigmoid_nonIBD_BMI$Obesity,levels=c("Normal_weight","Overweight","Obese"))
p1<-ggplot(data=newmeta_sigmoid_nonIBD_BMI,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_sigmoid_nonIBD_BMI,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

par(mfrow = c(1,2))
newmeta_cecum_nonIBD_BMI$Disease_ordered<-factor(newmeta_cecum_nonIBD_BMI$Obesity,levels=c("Normal_weight","Overweight","Obese"))
p1<-ggplot(data=newmeta_cecum_nonIBD_BMI,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum_nonIBD_BMI,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

kruskal.test(Chao1 ~ Quartile_CD_all, data=newmeta_sigmoid_nonIBD_Genotype)
dunnTest(Chao1 ~ Quartile_CD_all, data=newmeta_sigmoid_nonIBD_Genotype,method="bh") 
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_sigmoid_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
termplot(fit.data, terms="GRS_CD_all")
qplot(GRS_CD_all, Chao1, data = newmeta_sigmoid_nonIBD_Genotype, geom=c("point", "smooth"))+theme_pubr()
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_sigmoid_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Chao1 ~ Quartile_CD_all, data=newmeta_cecum_nonIBD_Genotype)
dunnTest(Chao1 ~ Quartile_CD_all, data=newmeta_cecum_nonIBD_Genotype,method="bh") 
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_cecum_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_cecum_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Shannon ~ Quartile_CD_all, data=newmeta_sigmoid_nonIBD_Genotype)
dunnTest(Shannon ~ Quartile_CD_all, data=newmeta_sigmoid_nonIBD_Genotype,method="bh") 
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_sigmoid_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_sigmoid_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Shannon ~ Quartile_CD_all, data=newmeta_cecum_nonIBD_Genotype)
dunnTest(Shannon ~ Quartile_CD_all, data=newmeta_cecum_CD_nonIBD_Genotype,method="bh")
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_cecum_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_cecum_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Chao1 ~ Quartile_CD_all, data=newmeta_sigmoid_CD_Genotype)
dunnTest(Chao1 ~ Quartile_CD_all, data=newmeta_sigmoid_CD_Genotype,method="bh") 
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
qplot(GRS_CD_all, Chao1, data = newmeta_sigmoid_CD_Genotype, geom=c("point", "smooth"))+theme_pubr()
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Chao1 ~ Quartile_CD_all, data=newmeta_cecum_CD_Genotype)
dunnTest(Chao1 ~ Quartile_CD_all, data=newmeta_cecum_CD_Genotype,method="bh") 
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)
qplot(GRS_CD_all, Chao1, data = newmeta_cecum_CD_Genotype, geom=c("point", "smooth"))+theme_pubr()
fit.data <- aov(Chao1 ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Shannon ~ Quartile_CD_all, data=newmeta_sigmoid_CD_Genotype)
dunnTest(Shannon ~ Quartile_CD_all, data=newmeta_sigmoid_CD_Genotype,method="bh") 
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

kruskal.test(Shannon ~ Quartile_CD_all, data=newmeta_cecum_CD_Genotype)
dunnTest(Shannon ~ Quartile_CD_all, data=newmeta_cecum_CD_Genotype,method="bh")
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(Shannon ~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Quartile_CD_all")

m_CD = lm(Chao1 ~ GRS_CD_all, data=newmeta_sigmoid_CD_Genotype)
m_nonIBD = lm(Chao1 ~ GRS_CD_all, data=newmeta_sigmoid_nonIBD_Genotype)
dd_m_CD = data.frame(x=newmeta_sigmoid_CD_Genotype$GRS_CD_all, y=predict(m_CD, newmeta_sigmoid_CD_Genotype), Disease=newmeta_sigmoid_CD_Genotype$Disease)
dd_m_nonIBD = data.frame(x=newmeta_sigmoid_nonIBD_Genotype$GRS_CD_all, y=predict(m_nonIBD, newmeta_sigmoid_nonIBD_Genotype), Disease=newmeta_sigmoid_nonIBD_Genotype$Disease)
ggplot(newmeta_sigmoid_genotype) + geom_point(aes(GRS_CD_all, Chao1, colour=Disease)) + geom_line(data=dd_m_CD, aes(x, y, colour=Disease)) + geom_line(data=dd_m_nonIBD, aes(x, y, colour=Disease))+theme_pubr()

m_CD = lm(Chao1 ~ GRS_CD_all, data=newmeta_cecum_CD_Genotype)
m_nonIBD = lm(Chao1 ~ GRS_CD_all, data=newmeta_cecum_nonIBD_Genotype)
dd_m_CD = data.frame(x=newmeta_cecum_CD_Genotype$GRS_CD_all, y=predict(m_CD, newmeta_cecum_CD_Genotype), Disease=newmeta_cecum_CD_Genotype$Disease)
dd_m_nonIBD = data.frame(x=newmeta_cecum_nonIBD_Genotype$GRS_CD_all, y=predict(m_nonIBD, newmeta_cecum_nonIBD_Genotype), Disease=newmeta_cecum_nonIBD_Genotype$Disease)
ggplot(newmeta_cecum_genotype) + geom_point(aes(GRS_CD_all, Chao1, colour=Disease)) + geom_line(data=dd_m_CD, aes(x, y, colour=Disease)) + geom_line(data=dd_m_nonIBD, aes(x, y, colour=Disease))+theme_pubr()

m_CD = lm(Shannon ~ GRS_CD_all, data=newmeta_sigmoid_CD_Genotype)
m_nonIBD = lm(Shannon ~ GRS_CD_all, data=newmeta_sigmoid_nonIBD_Genotype)
dd_m_CD = data.frame(x=newmeta_sigmoid_CD_Genotype$GRS_CD_all, y=predict(m_CD, newmeta_sigmoid_CD_Genotype), Disease=newmeta_sigmoid_CD_Genotype$Disease)
dd_m_nonIBD = data.frame(x=newmeta_sigmoid_nonIBD_Genotype$GRS_CD_all, y=predict(m_nonIBD, newmeta_sigmoid_nonIBD_Genotype), Disease=newmeta_sigmoid_nonIBD_Genotype$Disease)
ggplot(newmeta_sigmoid_genotype) + geom_point(aes(GRS_CD_all, Shannon, colour=Disease)) + geom_line(data=dd_m_CD, aes(x, y, colour=Disease)) + geom_line(data=dd_m_nonIBD, aes(x, y, colour=Disease))+theme_pubr()

m_CD = lm(Shannon ~ GRS_CD_all, data=newmeta_cecum_CD_Genotype)
m_nonIBD = lm(Shannon ~ GRS_CD_all, data=newmeta_cecum_nonIBD_Genotype)
dd_m_CD = data.frame(x=newmeta_cecum_CD_Genotype$GRS_CD_all, y=predict(m_CD, newmeta_cecum_CD_Genotype), Disease=newmeta_cecum_CD_Genotype$Disease)
dd_m_nonIBD = data.frame(x=newmeta_cecum_nonIBD_Genotype$GRS_CD_all, y=predict(m_nonIBD, newmeta_cecum_nonIBD_Genotype), Disease=newmeta_cecum_nonIBD_Genotype$Disease)
ggplot(newmeta_cecum_genotype) + geom_point(aes(GRS_CD_all, Shannon, colour=Disease)) + geom_line(data=dd_m_CD, aes(x, y, colour=Disease)) + geom_line(data=dd_m_nonIBD, aes(x, y, colour=Disease))+theme_pubr()

kruskal.test(Chao1 ~ Progression, data=newmeta_sigmoid_CD_progression_only)
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + Progression, data = newmeta_sigmoid_CD_progression)  
anova(fit.data)
coef(fit.data)
kruskal.test(Chao1 ~ Progression, data=newmeta_cecum_CD_progression_only)
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + Progression, data = newmeta_cecum_CD_progression)  
anova(fit.data)
coef(fit.data)

kruskal.test(Shannon ~ Progression, data=newmeta_sigmoid_CD_progression_only)
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + Progression, data = newmeta_sigmoid_CD_progression)  
anova(fit.data)
coef(fit.data)
kruskal.test(Shannon ~ Progression, data=newmeta_cecum_CD_progression_only)
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + Progression, data = newmeta_cecum_CD_progression)  
anova(fit.data)
coef(fit.data)

kruskal.test(Chao1 ~ B, data=newmeta_sigmoid_CD_B)
dunnTest(Chao1 ~ B, data=newmeta_sigmoid_CD_B,method="bh") 
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + B, data = newmeta_sigmoid_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B")

kruskal.test(Chao1 ~ B, data=newmeta_cecum_CD_B)
dunnTest(Chao1 ~ B, data=newmeta_cecum_CD_B,method="bh") 
fit.data <- aov(Chao1 ~ Batch + Age_at_Collection + Gender + Obesity + B, data = newmeta_cecum_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B")

kruskal.test(Shannon ~ B, data=newmeta_sigmoid_CD_B)
dunnTest(Shannon ~ B, data=newmeta_sigmoid_CD_B,method="bh") 
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + B, data = newmeta_sigmoid_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B")

kruskal.test(Shannon ~ B, data=newmeta_cecum_CD_B)
dunnTest(Shannon ~ B, data=newmeta_cecum_CD_B,method="bh")
fit.data <- aov(Shannon ~ Batch + Age_at_Collection + Gender + Obesity + B, data = newmeta_cecum_CD_B)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B")

par(mfrow = c(1,2))
newmeta_sigmoid_CD_B$Disease_ordered<-factor(newmeta_sigmoid_CD_B$B,levels=c("B1","B2","B3"))
p1<-ggplot(data=newmeta_sigmoid_CD_B,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_sigmoid_CD_B,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

par(mfrow = c(1,2))
newmeta_cecum_CD_B$Disease_ordered<-factor(newmeta_cecum_CD_B$B,levels=c("B1","B2","B3"))
p1<-ggplot(data=newmeta_cecum_CD_B,aes(x=Disease_ordered,y=Chao1,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum_CD_B,aes(x=Disease_ordered,y=Shannon,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)


### Beta diversity analyis

filtereddata_sigmoid<-newdata_sigmoid[ rowSums(newdata_sigmoid > 0) >= 20, ] 
dim(filtereddata_sigmoid)
# write.table(filtereddata_sigmoid, "Sigmoid_s20_ASV_count_table.txt", sep="\t", col.names=NA)
# data<-read.table(choose.files(),check.names=FALSE)   ## to import distance matrix after running DEICODE
# bray_dist_sigmoid <- as.dist(as(data, "matrix"))   
inputdata<-t(filtereddata_sigmoid)  # below commands to run Bray-Curtis
bray_dist_sigmoid<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid ~ Batch + Age_at_Collection, data=newmeta_sigmoid, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid ~ Batch + Gender, data=newmeta_sigmoid, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid ~ Batch + Obesity, data=newmeta_sigmoid, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid ~ Batch + Disease, data=newmeta_sigmoid, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid~ Batch + Age_at_Collection + Gender + Obesity + Disease, data=newmeta_sigmoid, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("CD"="red","Non-IBD"="blue")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_sigmoid$Disease)],pch=16,xlab = "PC1 (16%)", ylab = "PC2 (11%)", axes = TRUE, main = "PCoA")
plot=as.data.frame(pcoa_output$vectors)
ggplot(plot,aes(x=Axis.1,y=Axis.2, colour=col.list[paste(newmeta_sigmoid$Disease)])) + geom_point() + stat_ellipse() + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity()+theme_pubr()

filtereddata_cecum<-newdata_cecum[ rowSums(newdata_cecum > 0) >= 20, ] 
dim(filtereddata_cecum)
inputdata<-t(filtereddata_cecum)  
bray_dist_cecum<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum~ Batch + Age_at_Collection, data=newmeta_cecum, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum ~ Batch + Gender, data=newmeta_cecum, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum ~ Batch + Obesity, data=newmeta_cecum, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum ~ Batch + Disease, data=newmeta_cecum, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum ~ Batch + Age_at_Collection + Gender + Obesity + Disease, data=newmeta_cecum, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("CD"="red","Non-IBD"="blue")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_cecum$Disease)],pch=16,xlab = "PC1 (15%)", ylab = "PC2 (10%)", axes = TRUE, main = "PCoA")
plot=as.data.frame(pcoa_output$vectors)
ggplot(plot,aes(x=Axis.1,y=Axis.2, colour=col.list[paste(newmeta_cecum$Disease)])) + geom_point() + stat_ellipse() + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity()+theme_pubr()

filtereddata_sigmoid_CD_B<-newdata_sigmoid_CD_B[ rowSums(newdata_sigmoid_CD_B > 0) >= 8, ] 
dim(filtereddata_sigmoid_CD_B)
inputdata<-t(filtereddata_sigmoid_CD_B)  # below commands to run Bray-Curtis
bray_dist_sigmoid_CD_B<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_CD_B ~ Batch + B_combine, data=newmeta_sigmoid_CD_B, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_B~ Batch + Age_at_Collection + Gender + Obesity + B_combine, data=newmeta_sigmoid_CD_B, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid_CD_B)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("B2_B3"="red","B1"="blue")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_sigmoid_CD_B$B_combine)],pch=16,xlab = "PC1 (19%)", ylab = "PC2 (15%)", axes = TRUE, main = "PCoA")

filtereddata_cecum_CD_B<-newdata_cecum_CD_B[ rowSums(newdata_cecum_CD_B > 0) >= 8, ] 
dim(filtereddata_cecum_CD_B)
inputdata<-t(filtereddata_cecum_CD_B)  # below commands to run Bray-Curtis
bray_dist_cecum_CD_B<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_CD_B ~ Batch + B_combine, data=newmeta_cecum_CD_B, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_B ~ Batch + Age_at_Collection + Gender + Obesity + B_combine, data=newmeta_cecum_CD_B, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum_CD_B)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("B2_B3"="red","B1"="blue")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_cecum_CD_B$B_combine)],pch=16,xlab = "PC1 (17%)", ylab = "PC2 (12%)", axes = TRUE, main = "PCoA")

filtereddata_sigmoid_nonIBD_BMI<-newdata_sigmoid_nonIBD_BMI[ rowSums(newdata_sigmoid_nonIBD_BMI > 0) >= 10, ] 
dim(filtereddata_sigmoid_nonIBD_BMI)
inputdata<-t(filtereddata_sigmoid_nonIBD_BMI)  # below commands to run Bray-Curtis
bray_dist_sigmoid_nonIBD_BMI<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_nonIBD_BMI ~ Age_at_Collection, data=newmeta_sigmoid_nonIBD_BMI, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_nonIBD_BMI ~ Gender, data=newmeta_sigmoid_nonIBD_BMI, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_nonIBD_BMI ~ Obesity, data=newmeta_sigmoid_nonIBD_BMI, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_nonIBD_BMI~ Age_at_Collection + Gender + Obesity, data=newmeta_sigmoid_nonIBD_BMI, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid_nonIBD_BMI)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("Obese"="red","Normal_weight"="blue","Overweight"="grey")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_sigmoid_nonIBD_BMI$Obesity)],pch=16,xlab = "PC1 (15%)", ylab = "PC2 (11%)", axes = TRUE, main = "PCoA")
plot=as.data.frame(pcoa_output$vectors)
ggplot(plot,aes(x=Axis.1,y=Axis.2, colour=col.list[paste(newmeta_sigmoid_nonIBD_BMI$Obesity)])) + geom_point() + stat_ellipse() + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity()+theme_pubr()

filtereddata_cecum_nonIBD_BMI<-newdata_cecum_nonIBD_BMI[ rowSums(newdata_cecum_nonIBD_BMI > 0) >= 10, ] 
dim(filtereddata_cecum_nonIBD_BMI)
inputdata<-t(filtereddata_cecum_nonIBD_BMI)  # below commands to run Bray-Curtis
bray_dist_cecum_nonIBD_BMI<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_nonIBD_BMI ~ Age_at_Collection, data=newmeta_cecum_nonIBD_BMI, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_nonIBD_BMI ~ Gender, data=newmeta_cecum_nonIBD_BMI, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_nonIBD_BMI ~ Obesity, data=newmeta_cecum_nonIBD_BMI, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_nonIBD_BMI~ Age_at_Collection + Gender + Obesity, data=newmeta_cecum_nonIBD_BMI, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum_nonIBD_BMI)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("Obese"="red","Normal_weight"="blue","Overweight"="grey")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_cecum_nonIBD_BMI$Obesity)],pch=16,xlab = "PC1 (14%)", ylab = "PC2 (10%)", axes = TRUE, main = "PCoA")
plot=as.data.frame(pcoa_output$vectors)
ggplot(plot,aes(x=Axis.1,y=Axis.2, colour=col.list[paste(newmeta_cecum_nonIBD_BMI$Obesity)])) + geom_point() + stat_ellipse() + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity()+theme_pubr()

filtereddata_sigmoid_nonIBD_Genotype<-newdata_sigmoid_nonIBD_Genotype[ rowSums(newdata_sigmoid_nonIBD_Genotype> 0) >= 10, ] 
dim(filtereddata_sigmoid_nonIBD_Genotype)
inputdata<-t(filtereddata_sigmoid_nonIBD_Genotype)  # below commands to run Bray-Curtis
bray_dist_sigmoid_nonIBD_Genotype<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_nonIBD_Genotype ~ GRS_CD_all, data=newmeta_sigmoid_nonIBD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_nonIBD_Genotype ~ Quartile_CD_all, data=newmeta_sigmoid_nonIBD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_nonIBD_Genotype~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data=newmeta_sigmoid_nonIBD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_nonIBD_Genotype~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data=newmeta_sigmoid_nonIBD_Genotype, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid_nonIBD_Genotype)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
newmeta_sigmoid_nonIBD_Genotype$Color=3+newmeta_sigmoid_nonIBD_Genotype$GRS_CD_all
newmeta_sigmoid_nonIBD_Genotype$PC1=pcoa_output$vectors[,1]
newmeta_sigmoid_nonIBD_Genotype$PC2=pcoa_output$vectors[,2]
ggplot(newmeta_sigmoid_nonIBD_Genotype,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Color))+ scale_colour_gradient2()+theme_pubr()

filtereddata_cecum_nonIBD_Genotype<-newdata_cecum_nonIBD_Genotype[ rowSums(newdata_cecum_nonIBD_Genotype> 0) >= 10, ] 
dim(filtereddata_cecum_nonIBD_Genotype)
inputdata<-t(filtereddata_cecum_nonIBD_Genotype)  # below commands to run Bray-Curtis
bray_dist_cecum_nonIBD_Genotype<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_nonIBD_Genotype ~ GRS_CD_all, data=newmeta_sigmoid_nonIBD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_nonIBD_Genotype ~ Quartile_CD_all, data=newmeta_cecum_nonIBD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_nonIBD_Genotype~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data=newmeta_cecum_nonIBD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_nonIBD_Genotype~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data=newmeta_cecum_nonIBD_Genotype, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum_nonIBD_Genotype)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
newmeta_cecum_nonIBD_Genotype$Color=3+newmeta_cecum_nonIBD_Genotype$GRS_CD_all
newmeta_cecum_nonIBD_Genotype$PC1=pcoa_output$vectors[,1]
newmeta_cecum_nonIBD_Genotype$PC2=pcoa_output$vectors[,2]
ggplot(newmeta_cecum_nonIBD_Genotype,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Color))+ scale_colour_gradient2()+theme_pubr()

filtereddata_sigmoid_CD_Genotype<-newdata_sigmoid_CD_Genotype[ rowSums(newdata_sigmoid_CD_Genotype> 0) >= 8, ] 
dim(filtereddata_sigmoid_CD_Genotype)
inputdata<-t(filtereddata_sigmoid_CD_Genotype)  # below commands to run Bray-Curtis
bray_dist_sigmoid_CD_Genotype<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_CD_Genotype ~ GRS_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_CD_Genotype ~ Quartile_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_Genotype~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_Genotype~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_Genotype~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_Genotype ~ Batch + Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid_CD_Genotype)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
newmeta_sigmoid_CD_Genotype$Color=3+newmeta_sigmoid_CD_Genotype$GRS_CD_all
newmeta_sigmoid_CD_Genotype$PC1=pcoa_output$vectors[,1]
newmeta_sigmoid_CD_Genotype$PC2=pcoa_output$vectors[,2]
ggplot(newmeta_sigmoid_CD_Genotype,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Color))+ scale_colour_gradient2()+theme_pubr()

filtereddata_cecum_CD_Genotype<-newdata_cecum_CD_Genotype[ rowSums(newdata_cecum_CD_Genotype> 0) >= 8, ] 
dim(filtereddata_cecum_CD_Genotype)
inputdata<-t(filtereddata_cecum_CD_Genotype)  # below commands to run Bray-Curtis
bray_dist_cecum_CD_Genotype<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_CD_Genotype ~ GRS_CD_all, data=newmeta_sigmoid_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_CD_Genotype ~ Quartile_CD_all, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_Genotype~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_Genotype~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_CD_Genotype ~ Batch, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_Genotype~ Batch + Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_Genotype~ Age_at_Collection + Gender + Obesity + Quartile_CD_all, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_Genotype~ Age_at_Collection + Gender + Obesity + B_combine + Quartile_CD_all, data=newmeta_cecum_CD_Genotype, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum_CD_Genotype)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
newmeta_cecum_CD_Genotype$Color=3+newmeta_cecum_CD_Genotype$GRS_CD_all
newmeta_cecum_CD_Genotype$PC1=pcoa_output$vectors[,1]
newmeta_cecum_CD_Genotype$PC2=pcoa_output$vectors[,2]
ggplot(newmeta_cecum_CD_Genotype,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Color))+ scale_colour_gradient2()+theme_pubr()

filtereddata_sigmoid_CD_progression<-newdata_sigmoid_CD_progression_only[ rowSums(newdata_sigmoid_CD_progression_only > 0) >= 8, ] 
dim(filtereddata_sigmoid_CD_progression)
inputdata<-t(filtereddata_sigmoid_CD_progression)  # below commands to run Bray-Curtis
bray_dist_sigmoid_CD_progression<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_CD_progression ~ Batch, data=newmeta_sigmoid_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_CD_progression ~ Age_at_Collection, data=newmeta_sigmoid_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_CD_progression ~ Gender, data=newmeta_sigmoid_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_CD_progression ~ Obesity, data=newmeta_sigmoid_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_CD_progression ~ Progression, data=newmeta_sigmoid_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_progression ~ Batch + Age_at_Collection + Gender + Obesity + Progression, data=newmeta_sigmoid_CD_progression_only, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid_CD_progression)

filtereddata_cecum_CD_progression<-newdata_cecum_CD_progression_only[ rowSums(newdata_cecum_CD_progression_only > 0) >= 8, ] 
dim(filtereddata_cecum_CD_progression)
inputdata<-t(filtereddata_cecum_CD_progression)  # below commands to run Bray-Curtis
bray_dist_cecum_CD_progression<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_CD_progression ~ Batch, data=newmeta_cecum_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_CD_progression ~ Age_at_Collection, data=newmeta_cecum_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_CD_progression ~ Gender, data=newmeta_cecum_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_CD_progression ~ Obesity, data=newmeta_cecum_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_CD_progression ~ Progression, data=newmeta_cecum_CD_progression_only, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_progression ~ Batch + Age_at_Collection + Gender + Obesity + Progression, data=newmeta_cecum_CD_progression_only, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum_CD_progression)

filtereddata_sigmoid_CD_B<-newdata_sigmoid_CD_B[ rowSums(newdata_sigmoid_CD_B > 0) >= 8, ] 
dim(filtereddata_sigmoid_CD_B)
inputdata<-t(filtereddata_sigmoid_CD_B)  # below commands to run Bray-Curtis
bray_dist_sigmoid_CD_B<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_CD_B ~ Batch + B, data=newmeta_sigmoid_CD_B, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_CD_B~ Batch + Age_at_Collection + Gender + Obesity + B, data=newmeta_sigmoid_CD_B, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_sigmoid_CD_B)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("B3"="red","B1"="blue","B2"="grey")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_sigmoid_CD_B$B)],pch=16,xlab = "PC1 (19%)", ylab = "PC2 (15%)", axes = TRUE, main = "PCoA")

filtereddata_cecum_CD_B<-newdata_cecum_CD_B[ rowSums(newdata_cecum_CD_B > 0) >= 8, ] 
dim(filtereddata_cecum_CD_B)
inputdata<-t(filtereddata_cecum_CD_B)  # below commands to run Bray-Curtis
bray_dist_cecum_CD_B<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_CD_B ~ Batch + B, data=newmeta_cecum_CD_B, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_CD_B ~ Batch + Age_at_Collection + Gender + Obesity + B, data=newmeta_cecum_CD_B, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_dist_cecum_CD_B)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained
col.list=c("B3"="red","B1"="blue","B2"="grey")
plot(pcoa_output$vectors[,1],pcoa_output$vectors[,2],col=col.list[paste(newmeta_cecum_CD_B$B)],pch=16,xlab = "PC1 (17%)", ylab = "PC2 (12%)", axes = TRUE, main = "PCoA")


### Phylum level taxa summary plots by disease status

# First need to create phyloseq objects as in DESeq2 section
data_sigmoid_phylum <- tax_glom(data_sigmoid, "Phylum")
ps0 <- transform_sample_counts(data_sigmoid_phylum, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Disease")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")+ theme_pubr()
data_cecum_phylum <- tax_glom(data_cecum, "Phylum")
ps0 <- transform_sample_counts(data_cecum_phylum, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Disease")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")+ theme_pubr()
# Use DESeq2 for significance testing
diagdds=phyloseq_to_deseq2(data_sigmoid_phylum, ~ Batch + Gender + Obesity + Disease)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=FALSE, contrast=c("Disease","CD","Non-IBD"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid)[rownames(CDvsNorm), ], "matrix"))
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid CD vs Non-IBD adjusted batch gender obesity - PHYLUM.csv") 
diagdds=phyloseq_to_deseq2(data_cecum_phylum, ~ Batch + Gender + Obesity + Disease)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=FALSE, contrast=c("Disease","CD","Non-IBD"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum)[rownames(CDvsNorm), ], "matrix"))
write.csv(CDvsNormMatrix, "DESEq2 - Cecum CD vs Non-IBD adjusted batch gender obesity - PHYLUM.csv") 

# Summary for Obesity analysis
data_sigmoid_nonIBD_BMI_phylum <- tax_glom(data_sigmoid_nonIBD_BMI, "Phylum")
ps0 <- transform_sample_counts(data_sigmoid_nonIBD_BMI_phylum, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Obesity")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")+ theme_pubr()
data_cecum_nonIBD_BMI_phylum <- tax_glom(data_cecum_nonIBD_BMI, "Phylum")
ps0 <- transform_sample_counts(data_cecum_nonIBD_BMI_phylum, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Obesity")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")+ theme_pubr()
# Use DESeq2 for significance testing
diagdds=phyloseq_to_deseq2(data_sigmoid_nonIBD_BMI_phylum, ~ Gender + Obesity)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=FALSE, contrast=c("Obesity","Obese","Normal_weight"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_nonIBD_BMI)[rownames(CDvsNorm), ], "matrix"))
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid Obese vs. Normal weight Non-IBD adjusted gender - PHYLUM.csv") 
diagdds=phyloseq_to_deseq2(data_cecum_nonIBD_BMI_phylum, ~ Gender + Obesity)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=FALSE, contrast=c("Obesity","Obese","Normal_weight"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_nonIBD_BMI)[rownames(CDvsNorm), ], "matrix"))
write.csv(CDvsNormMatrix, "DESEq2 - Cecum Obese vs. Normal weight Non-IBD adjusted gender - PHYLUM.csv") 



### DESeq2 differential abundance analysis

alpha = 0.05 # set Q value filtering threshold for plots
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}

filtereddata_sigmoid2<-newdata_sigmoid[ rowSums(newdata_sigmoid > 0) >= 50, ] 
OTU_sigmoid=otu_table(filtereddata_sigmoid2, taxa_are_rows = TRUE)
biom_sigmoid = phyloseq(OTU_sigmoid,TAX)
map_sigmoid=sample_data(newmeta_sigmoid)
data_sigmoid=merge_phyloseq(biom_sigmoid,map_sigmoid)
print(data_sigmoid)
diagdds=phyloseq_to_deseq2(data_sigmoid, ~ Batch + Gender + Obesity + Disease)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Disease","CD","Non-IBD"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid CD vs Non-IBD adjusted batch gender obesity - ASVs s50.csv") 
DESEQ2_data_sigmoid = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid = DESEQ2_data_sigmoid[DESEQ2_data_sigmoid$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid = DESEQ2_data_sigmoid[!is.na(DESEQ2_data_sigmoid$Order), ] 
y = tapply(DESEQ2_data_sigmoid$log2FoldChange, DESEQ2_data_sigmoid$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid$Genus= factor(as.character(DESEQ2_data_sigmoid$Genus), levels = names(y))
dev.new(width=8, height=10)
ggplot(DESEQ2_data_sigmoid, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum2<-newdata_cecum[ rowSums(newdata_cecum> 0) >= 50, ] 
OTU_cecum=otu_table(filtereddata_cecum2, taxa_are_rows = TRUE)
biom_cecum = phyloseq(OTU_cecum,TAX)
map_cecum=sample_data(newmeta_cecum)
data_cecum=merge_phyloseq(biom_cecum,map_cecum)
print(data_cecum)
diagdds=phyloseq_to_deseq2(data_cecum, ~ Batch + Gender + Obesity + Disease)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Disease","CD","Non-IBD"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum CD vs Non-IBD adjusted batch gender obesity - ASVs s50.csv") 
DESEQ2_data_cecum = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum = DESEQ2_data_cecum[DESEQ2_data_cecum$Abundance >= 0.00001, ] 
DESEQ2_data_cecum = DESEQ2_data_cecum[!is.na(DESEQ2_data_cecum$Order), ] 
y = tapply(DESEQ2_data_cecum$log2FoldChange, DESEQ2_data_cecum$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum$Genus= factor(as.character(DESEQ2_data_cecum$Genus), levels = names(y))
dev.new(width=8, height=10)
ggplot(DESEQ2_data_cecum, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_sigmoid_CD_B2<-newdata_sigmoid_CD_B[ rowSums(newdata_sigmoid_CD_B > 0) >= 20, ] 
OTU_sigmoid_CD_B=otu_table(filtereddata_sigmoid_CD_B2, taxa_are_rows = TRUE)
biom_sigmoid_CD_B = phyloseq(OTU_sigmoid_CD_B,TAX)
map_sigmoid_CD_B=sample_data(newmeta_sigmoid_CD_B)
data_sigmoid_CD_B=merge_phyloseq(biom_sigmoid_CD_B,map_sigmoid_CD_B)
print(data_sigmoid_CD_B)
diagdds=phyloseq_to_deseq2(data_sigmoid_CD_B, ~ Gender + Obesity + B_combine)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("B_combine","B2_B3","B1"))  # rerun with different terms, e.g. "B3" vs. "B1" 
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_CD_B)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid B2+B3 vs. B1 adjusted batch gender obesity - ASVs s20.csv") 
DESEQ2_data_sigmoid_CD_B = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid_CD_B = DESEQ2_data_sigmoid_CD_B[DESEQ2_data_sigmoid_CD_B$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid_CD_B = DESEQ2_data_sigmoid_CD_B[!is.na(DESEQ2_data_sigmoid_CD_B$Order), ]
y = tapply(DESEQ2_data_sigmoid_CD_B$log2FoldChange, DESEQ2_data_sigmoid_CD_B$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid_CD_B$Genus= factor(as.character(DESEQ2_data_sigmoid_CD_B$Genus), levels = names(y))
dev.new(width=8, height=6)
ggplot(DESEQ2_data_sigmoid_CD_B, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum_CD_B2<-newdata_cecum_CD_B[ rowSums(newdata_cecum_CD_B > 0) >= 20, ] 
OTU_cecum_CD_B=otu_table(filtereddata_cecum_CD_B2, taxa_are_rows = TRUE)
biom_cecum_CD_B = phyloseq(OTU_cecum_CD_B,TAX)
map_cecum_CD_B=sample_data(newmeta_cecum_CD_B)
data_cecum_CD_B=merge_phyloseq(biom_cecum_CD_B,map_cecum_CD_B)
print(data_cecum_CD_B)
diagdds=phyloseq_to_deseq2(data_cecum_CD_B, ~ Gender + Obesity + B_combine)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("B_combine","B2_B3","B1"))  # rerun with different terms, e.g. "B3" vs. "B1" 
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_CD_B)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum B2+B3 vs. B1 adjusted batch gender obesity - ASVs s20.csv") 
DESEQ2_data_cecum_CD_B = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum_CD_B = DESEQ2_data_cecum_CD_B[DESEQ2_data_cecum_CD_B$Abundance >= 0.00001, ] 
DESEQ2_data_cecum_CD_B = DESEQ2_data_cecum_CD_B[!is.na(DESEQ2_data_cecum_CD_B$Order), ] 
y = tapply(DESEQ2_data_cecum_CD_B$log2FoldChange, DESEQ2_data_cecum_CD_B$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum_CD_B$Genus= factor(as.character(DESEQ2_data_cecum_CD_B$Genus), levels = names(y))
dev.new(width=8, height=5)
ggplot(DESEQ2_data_cecum_CD_B, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_sigmoid_CD_progression2<-newdata_sigmoid_CD_progression_only[ rowSums(newdata_sigmoid_CD_progression_only > 0) >= 18, ] 
OTU_sigmoid_CD_progression=otu_table(filtereddata_sigmoid_CD_progression2, taxa_are_rows = TRUE)
biom_sigmoid_CD_progression = phyloseq(OTU_sigmoid_CD_progression,TAX)
map_sigmoid_CD_progression=sample_data(newmeta_sigmoid_CD_progression_only)
data_sigmoid_CD_progression=merge_phyloseq(biom_sigmoid_CD_progression,map_sigmoid_CD_progression)
print(data_sigmoid_CD_progression)
diagdds=phyloseq_to_deseq2(data_sigmoid_CD_progression, ~ Gender + Obesity + Progression)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Progression","Progressor","Non-progressor"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_CD_progression)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid CD Progressor vs Non-progressor adjusted gender obesity - ASVs s18.csv") 
DESEQ2_data_sigmoid_CD_progression = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid_CD_progression = DESEQ2_data_sigmoid_CD_progression[DESEQ2_data_sigmoid_CD_progression$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid_CD_progression = DESEQ2_data_sigmoid_CD_progression[!is.na(DESEQ2_data_sigmoid_CD_progression$Order), ] 
y = tapply(DESEQ2_data_sigmoid_CD_progression$log2FoldChange, DESEQ2_data_sigmoid_CD_progression$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid_CD_progression$Genus= factor(as.character(DESEQ2_data_sigmoid_CD_progression$Genus), levels = names(y))
dev.new(width=8, height=5)
ggplot(DESEQ2_data_sigmoid_CD_progression, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum_CD_progression2<-newdata_cecum_CD_progression_only[ rowSums(newdata_cecum_CD_progression_only > 0) >= 18, ] 
OTU_cecum_CD_progression=otu_table(filtereddata_cecum_CD_progression2, taxa_are_rows = TRUE)
biom_cecum_CD_progression = phyloseq(OTU_cecum_CD_progression,TAX)
map_cecum_CD_progression=sample_data(newmeta_cecum_CD_progression_only)
data_cecum_CD_progression=merge_phyloseq(biom_cecum_CD_progression,map_cecum_CD_progression)
print(data_cecum_CD_progression)
diagdds=phyloseq_to_deseq2(data_cecum_CD_progression, ~ Gender + Obesity + Progression)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Progression","Progressor","Non-progressor"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_CD_progression)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,'(f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,'(o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum CD Progressor vs Non-progressor adjusted gender obesity - ASVs s18.csv") 
DESEQ2_data_cecum_CD_progression = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum_CD_progression = DESEQ2_data_cecum_CD_progression[DESEQ2_data_cecum_CD_progression$Abundance >= 0.00001, ] 
DESEQ2_data_cecum_CD_progression = DESEQ2_data_cecum_CD_progression[!is.na(DESEQ2_data_cecum_CD_progression$Order), ] 
y = tapply(DESEQ2_data_cecum_CD_progression$log2FoldChange, DESEQ2_data_cecum_CD_progression$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum_CD_progression$Genus= factor(as.character(DESEQ2_data_cecum_CD_progression$Genus), levels = names(y))
dev.new(width=8, height=4)
ggplot(DESEQ2_data_cecum_CD_progression, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_sigmoid_nonIBD_BMI2<-newdata_sigmoid_nonIBD_BMI[ rowSums(newdata_sigmoid_nonIBD_BMI > 0) >= 24, ] 
OTU_sigmoid_nonIBD_BMI=otu_table(filtereddata_sigmoid_nonIBD_BMI2, taxa_are_rows = TRUE)
biom_sigmoid_nonIBD_BMI = phyloseq(OTU_sigmoid_nonIBD_BMI,TAX)
map_sigmoid_nonIBD_BMI=sample_data(newmeta_sigmoid_nonIBD_BMI)
data_sigmoid_nonIBD_BMI=merge_phyloseq(biom_sigmoid_nonIBD_BMI,map_sigmoid_nonIBD_BMI)
print(data_sigmoid_nonIBD_BMI)
diagdds=phyloseq_to_deseq2(data_sigmoid_nonIBD_BMI, ~ Gender + Obesity)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Obesity","Obese","Normal_weight"))    # Overweight
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_nonIBD_BMI)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid Non-IBD Obese vs. Overweight adjusted gender - ASVs s24.csv") 
DESEQ2_data_sigmoid_nonIBD_BMI = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid_nonIBD_BMI = DESEQ2_data_sigmoid_nonIBD_BMI[DESEQ2_data_sigmoid_nonIBD_BMI$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid_nonIBD_BMI = DESEQ2_data_sigmoid_nonIBD_BMI[!is.na(DESEQ2_data_sigmoid_nonIBD_BMI$Order), ] 
y = tapply(DESEQ2_data_sigmoid_nonIBD_BMI$log2FoldChange, DESEQ2_data_sigmoid_nonIBD_BMI$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid_nonIBD_BMI$Genus= factor(as.character(DESEQ2_data_sigmoid_nonIBD_BMI$Genus), levels = names(y))
dev.new(width=8, height=3.5)
ggplot(DESEQ2_data_sigmoid_nonIBD_BMI, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum_nonIBD_BMI2<-newdata_cecum_nonIBD_BMI[ rowSums(newdata_cecum_nonIBD_BMI > 0) >= 24, ] 
OTU_cecum_nonIBD_BMI=otu_table(filtereddata_cecum_nonIBD_BMI2, taxa_are_rows = TRUE)
biom_cecum_nonIBD_BMI = phyloseq(OTU_cecum_nonIBD_BMI,TAX)
map_cecum_nonIBD_BMI=sample_data(newmeta_cecum_nonIBD_BMI)
data_cecum_nonIBD_BMI=merge_phyloseq(biom_cecum_nonIBD_BMI,map_cecum_nonIBD_BMI)
print(data_cecum_nonIBD_BMI)
diagdds=phyloseq_to_deseq2(data_cecum_nonIBD_BMI, ~ Gender + Obesity)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Obesity","Obese","Normal_weight"))    #  Overweight
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_nonIBD_BMI)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum Non-IBD Overweight vs. Normal adjusted gender - ASVs s24.csv") 
DESEQ2_data_cecum_nonIBD_BMI = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum_nonIBD_BMI = DESEQ2_data_cecum_nonIBD_BMI[DESEQ2_data_cecum_nonIBD_BMI$Abundance >= 0.00001, ] 
DESEQ2_data_cecum_nonIBD_BMI = DESEQ2_data_cecum_nonIBD_BMI[!is.na(DESEQ2_data_cecum_nonIBD_BMI$Order), ] 
y = tapply(DESEQ2_data_cecum_nonIBD_BMI$log2FoldChange, DESEQ2_data_cecum_nonIBD_BMI$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum_nonIBD_BMI$Genus= factor(as.character(DESEQ2_data_cecum_nonIBD_BMI$Genus), levels = names(y))
dev.new(width=8, height=5)
ggplot(DESEQ2_data_cecum_nonIBD_BMI, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_sigmoid_nonIBD_Genotype2<-newdata_sigmoid_nonIBD_Genotype[ rowSums(newdata_sigmoid_nonIBD_Genotype > 0) >= 24, ] 
OTU_sigmoid_nonIBD_Genotype=otu_table(filtereddata_sigmoid_nonIBD_Genotype2, taxa_are_rows = TRUE)
biom_sigmoid_nonIBD_Genotype = phyloseq(OTU_sigmoid_nonIBD_Genotype,TAX)
map_sigmoid_nonIBD_Genotype=sample_data(newmeta_sigmoid_nonIBD_Genotype)
data_sigmoid_nonIBD_Genotype=merge_phyloseq(biom_sigmoid_nonIBD_Genotype,map_sigmoid_nonIBD_Genotype)
print(data_sigmoid_nonIBD_Genotype)
diagdds=phyloseq_to_deseq2(data_sigmoid_nonIBD_Genotype, ~ Gender + Obesity + GRS_CD_all)   # Note: DESeq2 models did not converge with age included
# diagdds=phyloseq_to_deseq2(data_sigmoid_nonIBD_Genotype, ~ Gender + Obesity + Combine_CD_all)  #   Quartile_CD_all
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, name="GRS_CD_all") 
# CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Combine_CD_all","Four","OneandTwo"))  #  "Quartile_CD_all","Four","One"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_nonIBD_Genotype)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid Non-IBD GRS CD_all adjusted gender obesity - ASVs s24.csv") 
# write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid Non-IBD GRS CD_all Quartile 4 vs. 1_and_2_combined adjusted gender obesity - ASVs s24.csv") 
alpha=0.1
DESEQ2_data_sigmoid_nonIBD_Genotype = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid_nonIBD_Genotype = DESEQ2_data_sigmoid_nonIBD_Genotype[DESEQ2_data_sigmoid_nonIBD_Genotype$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid_nonIBD_Genotype = DESEQ2_data_sigmoid_nonIBD_Genotype[!is.na(DESEQ2_data_sigmoid_nonIBD_Genotype$Order), ] 
y = tapply(DESEQ2_data_sigmoid_nonIBD_Genotype$log2FoldChange, DESEQ2_data_sigmoid_nonIBD_Genotype$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid_nonIBD_Genotype$Genus= factor(as.character(DESEQ2_data_sigmoid_nonIBD_Genotype$Genus), levels = names(y))
dev.new(width=8, height=4)
ggplot(DESEQ2_data_sigmoid_nonIBD_Genotype, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum_nonIBD_Genotype2<-newdata_cecum_nonIBD_Genotype[ rowSums(newdata_cecum_nonIBD_Genotype > 0) >= 24, ] 
OTU_cecum_nonIBD_Genotype=otu_table(filtereddata_cecum_nonIBD_Genotype2, taxa_are_rows = TRUE)
biom_cecum_nonIBD_Genotype = phyloseq(OTU_cecum_nonIBD_Genotype,TAX)
map_cecum_nonIBD_Genotype=sample_data(newmeta_cecum_nonIBD_Genotype)
data_cecum_nonIBD_Genotype=merge_phyloseq(biom_cecum_nonIBD_Genotype,map_cecum_nonIBD_Genotype)
print(data_cecum_nonIBD_Genotype)
diagdds=phyloseq_to_deseq2(data_cecum_nonIBD_Genotype, ~ Gender + Obesity + GRS_CD_all)   # Note: DESeq2 models did not converge with age included
# diagdds=phyloseq_to_deseq2(data_cecum_nonIBD_Genotype, ~ Gender + Obesity + Quartile_CD_all)  #   Combine_CD_all
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, name="GRS_CD_all") 
# CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Quartile_CD_all","Four","One"))  # "Combine_CD_all","Four","OneandTwo" 
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_nonIBD_Genotype)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum Non-IBD GRS CD_all adjusted gender obesity - ASVs s24.csv") 
# write.csv(CDvsNormMatrix, "DESEq2 - Cecum Non-IBD GRS CD_all Quartile 4 vs. 1 adjusted gender obesity - ASVs s24.csv") 
DESEQ2_data_cecum_nonIBD_Genotype = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum_nonIBD_Genotype = DESEQ2_data_cecum_nonIBD_Genotype[DESEQ2_data_cecum_nonIBD_Genotype$Abundance >= 0.00001, ] 
DESEQ2_data_cecum_nonIBD_Genotype = DESEQ2_data_cecum_nonIBD_Genotype[!is.na(DESEQ2_data_cecum_nonIBD_Genotype$Order), ] 
y = tapply(DESEQ2_data_cecum_nonIBD_Genotype$log2FoldChange, DESEQ2_data_cecum_nonIBD_Genotype$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum_nonIBD_Genotype$Genus= factor(as.character(DESEQ2_data_cecum_nonIBD_Genotype$Genus), levels = names(y))
dev.new(width=8, height=4)
ggplot(DESEQ2_data_cecum_nonIBD_Genotype, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_sigmoid_CD_B2<-newdata_sigmoid_CD_B[ rowSums(newdata_sigmoid_CD_B > 0) >= 20, ] 
OTU_sigmoid_CD_B=otu_table(filtereddata_sigmoid_CD_B2, taxa_are_rows = TRUE)
biom_sigmoid_CD_B = phyloseq(OTU_sigmoid_CD_B,TAX)
map_sigmoid_CD_B=sample_data(newmeta_sigmoid_CD_B)
data_sigmoid_CD_B=merge_phyloseq(biom_sigmoid_CD_B,map_sigmoid_CD_B)
print(data_sigmoid_CD_B)
diagdds=phyloseq_to_deseq2(data_sigmoid_CD_B, ~ Gender + Obesity + B)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("B","B2","B1"))  # rerun with different terms, e.g. "B3" vs. "B1" 
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_CD_B)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid B2 vs. B1 adjusted batch gender obesity - ASVs s20.csv") 
DESEQ2_data_sigmoid_CD_B = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid_CD_B = DESEQ2_data_sigmoid_CD_B[DESEQ2_data_sigmoid_CD_B$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid_CD_B = DESEQ2_data_sigmoid_CD_B[!is.na(DESEQ2_data_sigmoid_CD_B$Order), ] 
y = tapply(DESEQ2_data_sigmoid_CD_B$log2FoldChange, DESEQ2_data_sigmoid_CD_B$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid_CD_B$Genus= factor(as.character(DESEQ2_data_sigmoid_CD_B$Genus), levels = names(y))
dev.new(width=8, height=5)
ggplot(DESEQ2_data_sigmoid_CD_B, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum_CD_B2<-newdata_cecum_CD_B[ rowSums(newdata_cecum_CD_B > 0) >= 20, ] 
OTU_cecum_CD_B=otu_table(filtereddata_cecum_CD_B2, taxa_are_rows = TRUE)
biom_cecum_CD_B = phyloseq(OTU_cecum_CD_B,TAX)
map_cecum_CD_B=sample_data(newmeta_cecum_CD_B)
data_cecum_CD_B=merge_phyloseq(biom_cecum_CD_B,map_cecum_CD_B)
print(data_cecum_CD_B)
diagdds=phyloseq_to_deseq2(data_cecum_CD_B, ~ Gender + Obesity + B)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("B","B2","B1"))  # rerun with different terms, e.g. "B3" vs. "B1" 
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_CD_B)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum B3 vs. B2 adjusted batch gender obesity - ASVs s20.csv") 
DESEQ2_data_cecum_CD_B = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum_CD_B = DESEQ2_data_cecum_CD_B[DESEQ2_data_cecum_CD_B$Abundance >= 0.00001, ] 
DESEQ2_data_cecum_CD_B = DESEQ2_data_cecum_CD_B[!is.na(DESEQ2_data_cecum_CD_B$Order), ] 
y = tapply(DESEQ2_data_cecum_CD_B$log2FoldChange, DESEQ2_data_cecum_CD_B$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum_CD_B$Genus= factor(as.character(DESEQ2_data_cecum_CD_B$Genus), levels = names(y))
dev.new(width=8, height=2.5)
ggplot(DESEQ2_data_cecum_CD_B, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_sigmoid_CD_Genotype2<-newdata_sigmoid_CD_Genotype_B[ rowSums(newdata_sigmoid_CD_Genotype_B > 0) >= 18, ] 
OTU_sigmoid_CD_Genotype=otu_table(filtereddata_sigmoid_CD_Genotype2, taxa_are_rows = TRUE)
biom_sigmoid_CD_Genotype = phyloseq(OTU_sigmoid_CD_Genotype,TAX)
map_sigmoid_CD_Genotype=sample_data(newmeta_sigmoid_CD_Genotype_B)
data_sigmoid_CD_Genotype=merge_phyloseq(biom_sigmoid_CD_Genotype,map_sigmoid_CD_Genotype)
print(data_sigmoid_CD_Genotype)
diagdds=phyloseq_to_deseq2(data_sigmoid_CD_Genotype, ~ Gender + Obesity + B_combine + GRS_CD_all)   # Note: DESeq2 models did not converge with age included
# diagdds=phyloseq_to_deseq2(data_sigmoid_CD_Genotype, ~ Gender + Obesity + Combine_CD_all)  #   Quartile_CD_all
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, name="GRS_CD_all") 
# CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Combine_CD_all","Four","OneandTwo"))  #  "Quartile_CD_all","Four","One"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_sigmoid_CD_Genotype)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid CD GRS CD_all adjusted gender obesity B_combine - ASVs s24.csv") 
# write.csv(CDvsNormMatrix, "DESEq2 - Sigmoid CD GRS CD_all Quartile 4 vs. 1_and_2_combined adjusted gender obesity B_combine - ASVs s24.csv") 
DESEQ2_data_sigmoid_CD_Genotype = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_sigmoid_CD_Genotype = DESEQ2_data_sigmoid_CD_Genotype[DESEQ2_data_sigmoid_CD_Genotype$Abundance >= 0.00001, ] 
DESEQ2_data_sigmoid_CD_Genotype = DESEQ2_data_sigmoid_CD_Genotype[!is.na(DESEQ2_data_sigmoid_CD_Genotype$Order), ] 
y = tapply(DESEQ2_data_sigmoid_CD_Genotype$log2FoldChange, DESEQ2_data_sigmoid_CD_Genotype$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_sigmoid_CD_Genotype$Genus= factor(as.character(DESEQ2_data_sigmoid_CD_Genotype$Genus), levels = names(y))
dev.new(width=8, height=6)
ggplot(DESEQ2_data_sigmoid_CD_Genotype, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_cecum_CD_Genotype2<-newdata_cecum_CD_Genotype_B[ rowSums(newdata_cecum_CD_Genotype_B > 0) >= 18, ] 
OTU_cecum_CD_Genotype=otu_table(filtereddata_cecum_CD_Genotype2, taxa_are_rows = TRUE)
biom_cecum_CD_Genotype = phyloseq(OTU_cecum_CD_Genotype,TAX)
map_cecum_CD_Genotype=sample_data(newmeta_cecum_CD_Genotype_B)
data_cecum_CD_Genotype=merge_phyloseq(biom_cecum_CD_Genotype,map_cecum_CD_Genotype)
print(data_cecum_CD_Genotype)
diagdds=phyloseq_to_deseq2(data_cecum_CD_Genotype, ~ Gender + Obesity + B_combine + GRS_CD_all)   # Note: DESeq2 models did not converge with age included
# diagdds=phyloseq_to_deseq2(data_cecum_CD_Genotype, ~ Gender + Obesity + Quartile_CD_all)  #   Combine_CD_all
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, name="GRS_CD_all") 
# CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Quartile_CD_all","Four","One"))  # "Combine_CD_all","Four","OneandTwo" 
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_cecum_CD_Genotype)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Cecum CD GRS CD_all adjusted gender obesity B_combine - ASVs s24.csv") 
# write.csv(CDvsNormMatrix, "DESEq2 - Cecum CD GRS CD_all Quartile 4 vs. 1 adjusted gender obesity B_combine - ASVs s24.csv") 
DESEQ2_data_cecum_CD_Genotype = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_cecum_CD_Genotype = DESEQ2_data_cecum_CD_Genotype[DESEQ2_data_cecum_CD_Genotype$Abundance >= 0.00001, ] 
DESEQ2_data_cecum_CD_Genotype = DESEQ2_data_cecum_CD_Genotype[!is.na(DESEQ2_data_cecum_CD_Genotype$Order), ] 
y = tapply(DESEQ2_data_cecum_CD_Genotype$log2FoldChange, DESEQ2_data_cecum_CD_Genotype$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_cecum_CD_Genotype$Genus= factor(as.character(DESEQ2_data_cecum_CD_Genotype$Genus), levels = names(y))
dev.new(width=8, height=4)
ggplot(DESEQ2_data_cecum_CD_Genotype, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")



### Random forests classifiers

## Sigmoid CD vs Non-IBD

sigmoid_index = createDataPartition(newmeta_sigmoid$Disease, p = 0.6, list = FALSE)
newmeta_sigmoid_train = newmeta_sigmoid[sigmoid_index, ]
newmeta_sigmoid_test = newmeta_sigmoid[-sigmoid_index, ]
dim(newmeta_sigmoid_train)
dim(newmeta_sigmoid_test)
write.csv(newmeta_sigmoid_train, "RF_sigmoid_train_meta.csv")
write.csv(newmeta_sigmoid_test, "RF_sigmoid_test_meta.csv")
relativedata=sweep(newdata_sigmoid,2,colSums(newdata_sigmoid),"/")
list=rownames(DESEQ2_data_sigmoid)
inputdata=relativedata[list,]
tnewdata_sigmoid_train = t(inputdata)[sigmoid_index, ]
tnewdata_sigmoid_test = t(inputdata)[-sigmoid_index, ]
dim(tnewdata_sigmoid_train)
dim(tnewdata_sigmoid_test)
length(which(rownames(tnewdata_sigmoid_train)==rownames(newmeta_sigmoid_train)))
newmeta_sigmoid_train$Disease=factor(newmeta_sigmoid_train$Disease,levels=c("Non-IBD","CD"))
Disease=as.character(newmeta_sigmoid_train$Disease)
Disease[Disease == "Non-IBD"]="NonIBD"
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_train,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_sigmoid_train,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$CD >= 2,])
inputdata=relativedata[newlist,]
tnewdata_sigmoid_train2 = t(inputdata)[sigmoid_index, ]
tnewdata_sigmoid_test2 = t(inputdata)[-sigmoid_index, ]
dim(tnewdata_sigmoid_train2)
dim(tnewdata_sigmoid_test2)
length(which(rownames(tnewdata_sigmoid_train2)==rownames(newmeta_sigmoid_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_train2,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_sigmoid=train(tnewdata_sigmoid_train2,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_sigmoid
selectedIndices <- rfFit_sigmoid$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_sigmoid$pred$obs[selectedIndices],rfFit_sigmoid$pred$CD[selectedIndices])  #Indicate category to predict
rfImportance_sigmoid=varImp(rfFit_sigmoid,scale=FALSE)
write.csv(rfImportance_sigmoid$importance, "RF_1001trees_Importance_Scores_Sigmoid_CD_vs_non-IBD.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_sigmoid, newdata = tnewdata_sigmoid_test2,type="raw")
probability_predictions=predict(rfFit_sigmoid, newdata = tnewdata_sigmoid_test2,type="prob")
newmeta_sigmoid_test$Disease[newmeta_sigmoid_test$Disease == "Non-IBD"]="NonIBD"
newmeta_sigmoid_test$Disease<-factor(newmeta_sigmoid_test$Disease,levels=c("NonIBD","CD"))
confusionMatrix(category_predictions,newmeta_sigmoid_test[,'Disease'],positive="CD")
plot.roc(newmeta_sigmoid_test$Disease,probability_predictions[[2]]) 
ROC_sigmoid <- roc(newmeta_sigmoid_test$Disease,probability_predictions[[2]])  
auc(ROC_sigmoid)
ci.auc(ROC_sigmoid, method="bootstrap", conf.level=0.95)
CI_plot1 <- plot.roc(newmeta_sigmoid_test$Disease,probability_predictions[[2]])  
plot(CI_plot1)
ROC1_CI <- ci.sp(CI_plot1, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)

## Cecum CD vs Non-IBD

cecum_index = createDataPartition(newmeta_cecum$Disease, p = 0.6, list = FALSE)
newmeta_cecum_train = newmeta_cecum[cecum_index, ]
newmeta_cecum_test = newmeta_cecum[-cecum_index, ]
dim(newmeta_cecum_train)
dim(newmeta_cecum_test)
write.csv(newmeta_cecum_train, "RF_cecum_train_meta.csv")
write.csv(newmeta_cecum_test, "RF_cecum_test_meta.csv")
relativedata=sweep(newdata_cecum,2,colSums(newdata_cecum),"/")
list=rownames(DESEQ2_data_cecum)
inputdata=relativedata[list,]
tnewdata_cecum_train = t(inputdata)[cecum_index, ]
tnewdata_cecum_test = t(inputdata)[-cecum_index, ]
dim(tnewdata_cecum_train)
dim(tnewdata_cecum_test)
length(which(rownames(tnewdata_cecum_train)==rownames(newmeta_cecum_train)))
newmeta_cecum_train$Disease=factor(newmeta_cecum_train$Disease,levels=c("Non-IBD","CD"))
Disease=as.character(newmeta_cecum_train$Disease)
Disease[Disease == "Non-IBD"]="NonIBD"
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_cecum_train,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_cecum_train,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$CD >= 2,])
inputdata=relativedata[newlist,]
tnewdata_cecum_train2 = t(inputdata)[cecum_index, ]
tnewdata_cecum_test2 = t(inputdata)[-cecum_index, ]
dim(tnewdata_cecum_train2)
dim(tnewdata_cecum_test2)
length(which(rownames(tnewdata_cecum_train2)==rownames(newmeta_cecum_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_cecum_train2,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_cecum=train(tnewdata_cecum_train2,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_cecum
selectedIndices <- rfFit_cecum$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_cecum$pred$obs[selectedIndices],rfFit_cecum$pred$CD[selectedIndices])  #Indicate category to predict
rfImportance_cecum=varImp(rfFit_cecum,scale=FALSE)
write.csv(rfImportance_cecum$importance, "RF_1001trees_Importance_Scores_Cecum_CD_vs_non-IBD.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_cecum, newdata = tnewdata_cecum_test2,type="raw")
probability_predictions=predict(rfFit_cecum, newdata = tnewdata_cecum_test2,type="prob")
newmeta_cecum_test$Disease[newmeta_cecum_test$Disease == "Non-IBD"]="NonIBD"
newmeta_cecum_test$Disease<-factor(newmeta_cecum_test$Disease,levels=c("NonIBD","CD"))
confusionMatrix(category_predictions,newmeta_cecum_test[,'Disease'],positive="CD")
plot.roc(newmeta_cecum_test$Disease,probability_predictions[[2]]) 
ROC_cecum <- roc(newmeta_cecum_test$Disease,probability_predictions[[2]])  
auc(ROC_cecum)
ci.auc(ROC_cecum, method="bootstrap", conf.level=0.95)
plot(CI_plot1)
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)
CI_plot2 <- plot.roc(newmeta_cecum_test$Disease,probability_predictions[[2]], add=TRUE)  
plot(CI_plot2, add=TRUE)
ROC2_CI <- ci.sp(CI_plot2, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC2_CI, type="shape",col="red",no.roc=FALSE)
roc.test(ROC_sigmoid,ROC_cecum,method="bootstrap")  #  p-value = 0.7453

## Sigmoid B2+B3 vs. B1

sigmoid_CD_B_index = createDataPartition(newmeta_sigmoid_CD_B$B_combine, p = 0.6, list = FALSE)
newmeta_sigmoid_CD_B_train = newmeta_sigmoid_CD_B[sigmoid_CD_B_index, ]
newmeta_sigmoid_CD_B_test = newmeta_sigmoid_CD_B[-sigmoid_CD_B_index, ]
dim(newmeta_sigmoid_CD_B_train)
dim(newmeta_sigmoid_CD_B_test)
write.csv(newmeta_sigmoid_CD_B_train, "RF_sigmoid_CD_B_combine_train_meta.csv")
write.csv(newmeta_sigmoid_CD_B_test, "RF_sigmoid_CD_B_combine_test_meta.csv")
relativedata=sweep(newdata_sigmoid_CD_B,2,colSums(newdata_sigmoid_CD_B),"/")
list=rownames(DESEQ2_data_sigmoid_CD_B)
inputdata=relativedata[list,]
tnewdata_sigmoid_CD_B_train = t(inputdata)[sigmoid_CD_B_index, ]
tnewdata_sigmoid_CD_B_test = t(inputdata)[-sigmoid_CD_B_index, ]
dim(tnewdata_sigmoid_CD_B_train)
dim(tnewdata_sigmoid_CD_B_test)
length(which(rownames(tnewdata_sigmoid_CD_B_train)==rownames(newmeta_sigmoid_CD_B_train)))
newmeta_sigmoid_CD_B_train$B_combine=factor(newmeta_sigmoid_CD_B_train$B_combine,levels=c("B1","B2_B3"))
B_combine=newmeta_sigmoid_CD_B_train$B_combine
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_CD_B_train,B_combine,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_sigmoid_CD_B_train,B_combine,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$B2_B3 >= 0,])
inputdata=relativedata[newlist,]
tnewdata_sigmoid_CD_B_train2 = t(inputdata)[sigmoid_CD_B_index, ]
tnewdata_sigmoid_CD_B_test2 = t(inputdata)[-sigmoid_CD_B_index, ]
dim(tnewdata_sigmoid_CD_B_train2)
dim(tnewdata_sigmoid_CD_B_test2)
length(which(rownames(tnewdata_sigmoid_CD_B_train2)==rownames(newmeta_sigmoid_CD_B_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_CD_B_train2,B_combine,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_sigmoid_CD_B=train(tnewdata_sigmoid_CD_B_train2,B_combine,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_sigmoid_CD_B
selectedIndices <- rfFit_sigmoid_CD_B$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_sigmoid_CD_B$pred$obs[selectedIndices],rfFit_sigmoid_CD_B$pred$B2_B3[selectedIndices])  #Indicate category to predict
rfImportance_sigmoid_CD_B=varImp(rfFit_sigmoid_CD_B,scale=FALSE)
write.csv(rfImportance_sigmoid_CD_B$importance, "RF_1001trees_Importance_Scores_Sigmoid_CD_B2+B3_vs_B1.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_sigmoid_CD_B, newdata = tnewdata_sigmoid_CD_B_test2,type="raw")
probability_predictions=predict(rfFit_sigmoid_CD_B, newdata = tnewdata_sigmoid_CD_B_test2,type="prob")
newmeta_sigmoid_CD_B_test$B_combine=factor(newmeta_sigmoid_CD_B_test$B_combine,levels=c("B1","B2_B3"))
confusionMatrix(category_predictions,newmeta_sigmoid_CD_B_test[,'B_combine'],positive="B2_B3")
plot.roc(newmeta_sigmoid_CD_B_test$B_combine,probability_predictions[[2]]) 
ROC_sigmoid_CD_B <- roc(newmeta_sigmoid_CD_B_test$B_combine,probability_predictions[[2]])  
auc(ROC_sigmoid_CD_B)
ci.auc(ROC_sigmoid_CD_B, method="bootstrap", conf.level=0.95)
CI_plot1 <- plot.roc(newmeta_sigmoid_CD_B_test$B_combine,probability_predictions[[2]])  
plot(CI_plot1)
ROC1_CI <- ci.sp(CI_plot1, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)

## Cecum B2+B3 vs. B1

cecum_CD_B_index = createDataPartition(newmeta_cecum_CD_B$B_combine, p = 0.6, list = FALSE)
newmeta_cecum_CD_B_train = newmeta_cecum_CD_B[cecum_CD_B_index, ]
newmeta_cecum_CD_B_test = newmeta_cecum_CD_B[-cecum_CD_B_index, ]
dim(newmeta_cecum_CD_B_train)
dim(newmeta_cecum_CD_B_test)
write.csv(newmeta_cecum_CD_B_train, "RF_cecum_CD_B_combine_train_meta.csv")
write.csv(newmeta_cecum_CD_B_test, "RF_cecum_CD_B_combine_test_meta.csv")
relativedata=sweep(newdata_cecum_CD_B,2,colSums(newdata_cecum_CD_B),"/")
list=rownames(DESEQ2_data_cecum_CD_B)
inputdata=relativedata[list,]
tnewdata_cecum_CD_B_train = t(inputdata)[cecum_CD_B_index, ]
tnewdata_cecum_CD_B_test = t(inputdata)[-cecum_CD_B_index, ]
dim(tnewdata_cecum_CD_B_train)
dim(tnewdata_cecum_CD_B_test)
length(which(rownames(tnewdata_cecum_CD_B_train)==rownames(newmeta_cecum_CD_B_train)))
newmeta_cecum_CD_B_train$B_combine=factor(newmeta_cecum_CD_B_train$B_combine,levels=c("B1","B2_B3"))
B_combine=newmeta_cecum_CD_B_train$B_combine
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_cecum_CD_B_train,B_combine,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_cecum_CD_B_train,B_combine,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$B2_B3 >= 0,])
inputdata=relativedata[newlist,]
tnewdata_cecum_CD_B_train2 = t(inputdata)[cecum_CD_B_index, ]
tnewdata_cecum_CD_B_test2 = t(inputdata)[-cecum_CD_B_index, ]
dim(tnewdata_cecum_CD_B_train2)
dim(tnewdata_cecum_CD_B_test2)
length(which(rownames(tnewdata_cecum_CD_B_train2)==rownames(newmeta_cecum_CD_B_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_cecum_CD_B_train2,B_combine,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_cecum_CD_B=train(tnewdata_cecum_CD_B_train2,B_combine,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_cecum_CD_B
selectedIndices <- rfFit_cecum_CD_B$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_cecum_CD_B$pred$obs[selectedIndices],rfFit_cecum_CD_B$pred$B2_B3[selectedIndices])  #Indicate category to predict
rfImportance_cecum_CD_B=varImp(rfFit_cecum_CD_B,scale=FALSE)
write.csv(rfImportance_cecum_CD_B$importance, "RF_1001trees_Importance_Scores_Cecum_CD_B2+B3_vs_B1.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_cecum_CD_B, newdata = tnewdata_cecum_CD_B_test2,type="raw")
probability_predictions=predict(rfFit_cecum_CD_B, newdata = tnewdata_cecum_CD_B_test2,type="prob")
newmeta_cecum_CD_B_test$B_combine=factor(newmeta_cecum_CD_B_test$B_combine,levels=c("B1","B2_B3"))
confusionMatrix(category_predictions,newmeta_cecum_CD_B_test[,'B_combine'],positive="B2_B3")
plot.roc(newmeta_cecum_CD_B_test$B_combine,probability_predictions[[2]]) 
ROC_cecum_CD_B <- roc(newmeta_cecum_CD_B_test$B_combine,probability_predictions[[2]])  
auc(ROC_cecum_CD_B)
ci.auc(ROC_cecum_CD_B, method="bootstrap", conf.level=0.95)
plot(CI_plot1)
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)
CI_plot2 <- plot.roc(newmeta_cecum_CD_B_test$B_combine,probability_predictions[[2]], add=TRUE)
plot(CI_plot2, add=TRUE)
ROC2_CI <- ci.sp(CI_plot2, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC2_CI, type="shape",col="red",no.roc=FALSE)
roc.test(ROC_sigmoid_CD_B,ROC_cecum_CD_B,method="bootstrap")  # p-value = 0.9331

## Sigmoid CD Progressor vs. Non-progressor

sigmoid_CD_progression_index = createDataPartition(newmeta_sigmoid_CD_progression_only$Progression, p = 0.6, list = FALSE)
newmeta_sigmoid_CD_progression_train = newmeta_sigmoid_CD_progression_only[sigmoid_CD_progression_index, ]
newmeta_sigmoid_CD_progression_test = newmeta_sigmoid_CD_progression_only[-sigmoid_CD_progression_index, ]
dim(newmeta_sigmoid_CD_progression_train)
dim(newmeta_sigmoid_CD_progression_test)
write.csv(newmeta_sigmoid_CD_progression_train, "RF_sigmoid_CD_progression_train_meta.csv")
write.csv(newmeta_sigmoid_CD_progression_test, "RF_sigmoid_CD_progression_test_meta.csv")
relativedata=sweep(newdata_sigmoid_CD_progression_only,2,colSums(newdata_sigmoid_CD_progression_only),"/")
list=rownames(DESEQ2_data_sigmoid_CD_progression)
inputdata=relativedata[list,]
tnewdata_sigmoid_CD_progression_train = t(inputdata)[sigmoid_CD_progression_index, ]
tnewdata_sigmoid_CD_progression_test = t(inputdata)[-sigmoid_CD_progression_index, ]
dim(tnewdata_sigmoid_CD_progression_train)
dim(tnewdata_sigmoid_CD_progression_test)
length(which(rownames(tnewdata_sigmoid_CD_progression_train)==rownames(newmeta_sigmoid_CD_progression_train)))
newmeta_sigmoid_CD_progression_train$Progression=factor(newmeta_sigmoid_CD_progression_train$Progression,levels=c("Non-progressor","Progressor"))
Progression=as.character(newmeta_sigmoid_CD_progression_train$Progression)
Progression[Progression== "Non-progressor"]="Nonprogressor"
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,3,4,5,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_sigmoid_CD_progression_train,Progression,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$Progressor >= 0,])
inputdata=relativedata[newlist,]
tnewdata_sigmoid_CD_progression_train2 = t(inputdata)[sigmoid_CD_progression_index, ]
tnewdata_sigmoid_CD_progression_test2 = t(inputdata)[-sigmoid_CD_progression_index, ]
dim(tnewdata_sigmoid_CD_progression_train2)
dim(tnewdata_sigmoid_CD_progression_test2)
length(which(rownames(tnewdata_sigmoid_CD_progression_train2)==rownames(newmeta_sigmoid_CD_progression_train)))
mtry=c(2,3,4,5,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_sigmoid_CD_progression=train(tnewdata_sigmoid_CD_progression_train2,Progression,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_sigmoid_CD_progression
selectedIndices <- rfFit_sigmoid_CD_progression$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_sigmoid_CD_progression$pred$obs[selectedIndices],rfFit_sigmoid_CD_progression$pred$Progressor[selectedIndices])  #Indicate category to predict
rfImportance_sigmoid_CD_progression=varImp(rfFit_sigmoid_CD_progression,scale=FALSE)
write.csv(rfImportance_sigmoid_CD_progression$importance, "RF_1001trees_Importance_Scores_Sigmoid_CD_Progressor_vs_non-progressor.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_sigmoid_CD_progression, newdata = tnewdata_sigmoid_CD_progression_test2,type="raw")
probability_predictions=predict(rfFit_sigmoid_CD_progression, newdata = tnewdata_sigmoid_CD_progression_test2,type="prob")
newmeta_sigmoid_CD_progression_test$Progression[newmeta_sigmoid_CD_progression_test$Progression== "Non-progressor"]="Nonprogressor"
newmeta_sigmoid_CD_progression_test$Progression<-factor(newmeta_sigmoid_CD_progression_test$Progression,levels=c("Nonprogressor","Progressor"))
confusionMatrix(category_predictions,newmeta_sigmoid_CD_progression_test[,'Progression'],positive="Progressor")
plot.roc(newmeta_sigmoid_CD_progression_test$Progression,probability_predictions[[2]]) 
ROC_sigmoid_CD_progression <- roc(newmeta_sigmoid_CD_progression_test$Progression,probability_predictions[[2]])  
auc(ROC_sigmoid_CD_progression)
ci.auc(ROC_sigmoid_CD_progression, method="bootstrap", conf.level=0.95)
CI_plot1 <- plot.roc(newmeta_sigmoid_CD_progression_test$Progression,probability_predictions[[2]])  
plot(CI_plot1)
ROC1_CI <- ci.sp(CI_plot1, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)
# Bar graph of taxa contributing to the RF classifier
newlist=rownames(rfImportance_sigmoid_CD_progression$importance)
sigmoid_CD_progression_RF_taxonomy=DESEQ2_data_sigmoid_CD_progression[newlist,]
length(which(rownames(rfImportance_sigmoid_CD_progression$importance)==rownames(sigmoid_CD_progression_RF_taxonomy)))
bar_plot_input=cbind(rownames(rfImportance_sigmoid_CD_progression$importance),rfImportance_sigmoid_CD_progression$importance[,1],sigmoid_CD_progression_RF_taxonomy)
bar_plot_input$Direction <- ifelse(bar_plot_input$log2FoldChange >=0, "Enriched", "Depleted")
colnames(bar_plot_input)[1]="ASV"
colnames(bar_plot_input)[2]="Importance"
dev.new(width=7, height=4)
ggplot(bar_plot_input,aes(x=reorder(Genus,Importance),y=Importance,fill=Phylum))+geom_bar(stat='identity')+coord_flip()+theme_pubr()+theme(axis.text.y = element_text(colour = ifelse(bar_plot_input$Direction == "Enriched", "darkgreen", "darkred")),axis.title.y=element_blank())

## Cecum CD Progressor vs. Non-progressor

cecum_CD_progression_index = createDataPartition(newmeta_cecum_CD_progression$Progression, p = 0.6, list = FALSE)
newmeta_cecum_CD_progression_train = newmeta_cecum_CD_progression[cecum_CD_progression_index, ]
newmeta_cecum_CD_progression_test = newmeta_cecum_CD_progression[-cecum_CD_progression_index, ]
dim(newmeta_cecum_CD_progression_train)
dim(newmeta_cecum_CD_progression_test)
write.csv(newmeta_cecum_CD_progression_train, "RF_cecum_CD_progression_train_meta.csv")
write.csv(newmeta_cecum_CD_progression_test, "RF_cecum_CD_progression_test_meta.csv")
relativedata=sweep(newdata_cecum_CD_progression,2,colSums(newdata_cecum_CD_progression),"/")
list=rownames(DESEQ2_data_cecum_CD_progression)
inputdata=relativedata[list,]
tnewdata_cecum_CD_progression_train = t(inputdata)[cecum_CD_progression_index, ]
tnewdata_cecum_CD_progression_test = t(inputdata)[-cecum_CD_progression_index, ]
dim(tnewdata_cecum_CD_progression_train)
dim(tnewdata_cecum_CD_progression_test)
length(which(rownames(tnewdata_cecum_CD_progression_train)==rownames(newmeta_cecum_CD_progression_train)))
newmeta_cecum_CD_progression_train$Progression=factor(newmeta_cecum_CD_progression_train$Progression,levels=c("Non-progressor","Progressor"))
Progression=as.character(newmeta_cecum_CD_progression_train$Progression)
Progression[Progression== "Non-progressor"]="Nonprogressor"
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_cecum_CD_progression_train,Progression,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$Progressor >= 0,])
inputdata=relativedata[newlist,]
tnewdata_cecum_CD_progression_train2 = t(inputdata)[cecum_CD_progression_index, ]
tnewdata_cecum_CD_progression_test2 = t(inputdata)[-cecum_CD_progression_index, ]
dim(tnewdata_cecum_CD_progression_train2)
dim(tnewdata_cecum_CD_progression_test2)
length(which(rownames(tnewdata_cecum_CD_progression_train2)==rownames(newmeta_cecum_CD_progression_train)))
mtry=c(2,3,4,6,8,10)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_cecum_CD_progression=train(tnewdata_cecum_CD_progression_train2,Progression,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_cecum_CD_progression
selectedIndices <- rfFit_cecum_CD_progression$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_cecum_CD_progression$pred$obs[selectedIndices],rfFit_cecum_CD_progression$pred$Progressor[selectedIndices])  #Indicate category to predict
rfImportance_cecum_CD_progression=varImp(rfFit_cecum_CD_progression,scale=FALSE)
write.csv(rfImportance_cecum_CD_progression$importance, "RF_1001trees_Importance_Scores_Cecum_CD_Progressor_vs_non-progressor.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_cecum_CD_progression, newdata = tnewdata_cecum_CD_progression_test2,type="raw")
probability_predictions=predict(rfFit_cecum_CD_progression, newdata = tnewdata_cecum_CD_progression_test2,type="prob")
newmeta_cecum_CD_progression_test$Progression[newmeta_cecum_CD_progression_test$Progression== "Non-progressor"]="Nonprogressor"
newmeta_cecum_CD_progression_test$Progression<-factor(newmeta_cecum_CD_progression_test$Progression,levels=c("Nonprogressor","Progressor"))
confusionMatrix(category_predictions,newmeta_cecum_CD_progression_test[,'Progression'],positive="Progressor")
plot.roc(newmeta_cecum_CD_progression_test$Progression,probability_predictions[[2]]) 
ROC_cecum_CD_progression <- roc(newmeta_cecum_CD_progression_test$Progression,probability_predictions[[2]])  
auc(ROC_cecum_CD_progression)
ci.auc(ROC_cecum_CD_progression, method="bootstrap", conf.level=0.95)
plot(CI_plot1)
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)
CI_plot2 <- plot.roc(newmeta_cecum_CD_progression_test$Progression,probability_predictions[[2]], add=TRUE)  
plot(CI_plot2, add=TRUE)
ROC2_CI <- ci.sp(CI_plot2, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC2_CI, type="shape",col="red",no.roc=FALSE)
roc.test(ROC_sigmoid_CD_progression,ROC_cecum_CD_progression,method="bootstrap")  # p-value =  0.7836
# Bar graph of taxa contributing to the RF classifier
newlist=rownames(rfImportance_cecum_CD_progression$importance)
cecum_CD_progression_RF_taxonomy=DESEQ2_data_cecum_CD_progression[newlist,]
length(which(rownames(rfImportance_cecum_CD_progression$importance)==rownames(cecum_CD_progression_RF_taxonomy)))
bar_plot_input=cbind(rownames(rfImportance_cecum_CD_progression$importance),rfImportance_cecum_CD_progression$importance[,1],cecum_CD_progression_RF_taxonomy)
bar_plot_input$Direction <- ifelse(bar_plot_input$log2FoldChange >=0, "Enriched", "Depleted")
colnames(bar_plot_input)[1]="ASV"
colnames(bar_plot_input)[2]="Importance"
dev.new(width=7, height=3)
ggplot(bar_plot_input,aes(x=reorder(Genus,Importance),y=Importance,fill=Phylum))+geom_bar(stat='identity')+coord_flip()+theme_pubr()+theme(axis.text.y = element_text(colour = ifelse(bar_plot_input$Direction == "Enriched", "darkgreen", "darkred")),axis.title.y=element_blank())

## Sigmoid Non-IBD Obese vs. Normal weight

s27<-factor(newmeta_sigmoid_nonIBD_BMI$Obesity)
split.m27 <- split.data.frame(newmeta_sigmoid_nonIBD_BMI, s27)
dim(split.m27$'Obese')
dim(split.m27$'Normal_weight')
newmeta_sigmoid_nonIBD_obese_normal=rbind(split.m27$'Obese',split.m27$'Normal_weight')
dim(newmeta_sigmoid_nonIBD_obese_normal)
tdata27=t(newdata_sigmoid_nonIBD_BMI)
split.s27 <- split.data.frame(tdata27, s27)
dim(split.s27$'Obese')
dim(split.s27$'Normal_weight')
newdata_sigmoid_nonIBD_obese_normal=t(as.matrix(rbind(split.s27$'Obese',split.s27$'Normal_weight')))
dim(newdata_sigmoid_nonIBD_obese_normal)
length(which(colnames(newdata_sigmoid_nonIBD_obese_normal)==rownames(newmeta_sigmoid_nonIBD_obese_normal)))

# Do pairwise Adonis
filtereddata_sigmoid_nonIBD_obese_normal<-newdata_sigmoid_nonIBD_obese_normal[ rowSums(newdata_sigmoid_nonIBD_obese_normal > 0) >= 6, ] 
dim(filtereddata_sigmoid_nonIBD_obese_normal)
inputdata<-t(filtereddata_sigmoid_nonIBD_obese_normal)  # below commands to run Bray-Curtis
bray_dist_sigmoid_nonIBD_obese_normal<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_sigmoid_nonIBD_obese_normal ~ Age_at_Collection, data=newmeta_sigmoid_nonIBD_obese_normal, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_nonIBD_obese_normal ~ Gender, data=newmeta_sigmoid_nonIBD_obese_normal, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_sigmoid_nonIBD_obese_normal ~ Obesity, data=newmeta_sigmoid_nonIBD_obese_normal, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_sigmoid_nonIBD_obese_normal ~ Age_at_Collection + Gender + Obesity, data=newmeta_sigmoid_nonIBD_obese_normal, permutations=10000)
data.adonis

sigmoid_obese_normal_index = createDataPartition(newmeta_sigmoid_nonIBD_obese_normal$Obesity, p = 0.6, list = FALSE)
newmeta_sigmoid_obese_normal_train = newmeta_sigmoid_nonIBD_obese_normal[sigmoid_obese_normal_index, ]
newmeta_sigmoid_obese_normal_test = newmeta_sigmoid_nonIBD_obese_normal[-sigmoid_obese_normal_index, ]
dim(newmeta_sigmoid_obese_normal_train)
dim(newmeta_sigmoid_obese_normal_test)
write.csv(newmeta_sigmoid_obese_normal_train, "RF_sigmoid_obese_normal_train_meta.csv")
write.csv(newmeta_sigmoid_obese_normal_test, "RF_sigmoid_obese_normal_test_meta.csv")
relativedata=sweep(newdata_sigmoid_nonIBD_obese_normal,2,colSums(newdata_sigmoid_nonIBD_obese_normal),"/")
list=rownames(DESEQ2_data_sigmoid_nonIBD_BMI)
inputdata=relativedata[list,]
tnewdata_sigmoid_obese_normal_train = t(inputdata)[sigmoid_obese_normal_index, ]
tnewdata_sigmoid_obese_normal_test = t(inputdata)[-sigmoid_obese_normal_index, ]
dim(tnewdata_sigmoid_obese_normal_train)
dim(tnewdata_sigmoid_obese_normal_test)
length(which(rownames(tnewdata_sigmoid_obese_normal_train)==rownames(newmeta_sigmoid_obese_normal_train)))
newmeta_sigmoid_obese_normal_train$Obesity=factor(newmeta_sigmoid_obese_normal_train$Obesity,levels=c("Normal_weight","Obese"))
Disease=as.character(newmeta_sigmoid_obese_normal_train$Obesity)
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_obese_normal_train,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_sigmoid_obese_normal_train,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$Obese >= 1,])
inputdata=relativedata[newlist,]
tnewdata_sigmoid_obese_normal_train2 = t(inputdata)[sigmoid_obese_normal_index, ]
tnewdata_sigmoid_obese_normal_test2 = t(inputdata)[-sigmoid_obese_normal_index, ]
dim(tnewdata_sigmoid_obese_normal_train2)
dim(tnewdata_sigmoid_obese_normal_test2)
length(which(rownames(tnewdata_sigmoid_obese_normal_train2)==rownames(newmeta_sigmoid_obese_normal_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_obese_normal_train2,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_sigmoid_obese_normal=train(tnewdata_sigmoid_obese_normal_train2,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_sigmoid_obese_normal
selectedIndices <- rfFit_sigmoid_obese_normal$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_sigmoid_obese_normal$pred$obs[selectedIndices],rfFit_sigmoid_obese_normal$pred$Obese[selectedIndices])  #Indicate category to predict
rfImportance_sigmoid_obese_normal=varImp(rfFit_sigmoid_obese_normal,scale=FALSE)
write.csv(rfImportance_sigmoid_obese_normal$importance, "RF_1001trees_Importance_Scores_Sigmoid_non-IBD_Obese_vs_Normalweight.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_sigmoid_obese_normal, newdata = tnewdata_sigmoid_obese_normal_test2,type="raw")
probability_predictions=predict(rfFit_sigmoid_obese_normal, newdata = tnewdata_sigmoid_obese_normal_test2,type="prob")
newmeta_sigmoid_obese_normal_test$Obesity<-factor(newmeta_sigmoid_obese_normal_test$Obesity,levels=c("Normal_weight","Obese"))
confusionMatrix(category_predictions,newmeta_sigmoid_obese_normal_test[,'Obesity'],positive="Obese")
plot.roc(newmeta_sigmoid_obese_normal_test$Obesity,probability_predictions[[2]]) 
ROC_sigmoid_obese_normal <- roc(newmeta_sigmoid_obese_normal_test$Obesity,probability_predictions[[2]])  
auc(ROC_sigmoid_obese_normal)
CI_plot1 <- plot.roc(newmeta_sigmoid_obese_normal_test$Obesity,probability_predictions[[2]])  
plot(CI_plot1)
ROC1_CI <- ci.sp(CI_plot1, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)

# Cecum pairwise Adonis
s28<-factor(newmeta_cecum_nonIBD_BMI$Obesity)
split.m28 <- split.data.frame(newmeta_cecum_nonIBD_BMI, s28)
dim(split.m28$'Obese')
dim(split.m28$'Normal_weight')
newmeta_cecum_nonIBD_obese_normal=rbind(split.m28$'Obese',split.m28$'Normal_weight')
dim(newmeta_cecum_nonIBD_obese_normal)
tdata28=t(newdata_cecum_nonIBD_BMI)
split.s28 <- split.data.frame(tdata28, s28)
dim(split.s28$'Obese')
dim(split.s28$'Normal_weight')
newdata_cecum_nonIBD_obese_normal=t(as.matrix(rbind(split.s28$'Obese',split.s28$'Normal_weight')))
dim(newdata_cecum_nonIBD_obese_normal)
length(which(colnames(newdata_cecum_nonIBD_obese_normal)==rownames(newmeta_cecum_nonIBD_obese_normal)))

# Do pairwise Adonis
filtereddata_cecum_nonIBD_obese_normal<-newdata_cecum_nonIBD_obese_normal[ rowSums(newdata_cecum_nonIBD_obese_normal > 0) >= 6, ] 
dim(filtereddata_cecum_nonIBD_obese_normal)
inputdata<-t(filtereddata_cecum_nonIBD_obese_normal)  # below commands to run Bray-Curtis
bray_dist_cecum_nonIBD_obese_normal<-vegdist(inputdata, method="bray")
data.adonis=adonis(bray_dist_cecum_nonIBD_obese_normal ~ Age_at_Collection, data=newmeta_cecum_nonIBD_obese_normal, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_nonIBD_obese_normal ~ Gender, data=newmeta_cecum_nonIBD_obese_normal, permutations=10000)
data.adonis
data.adonis=adonis(bray_dist_cecum_nonIBD_obese_normal ~ Obesity, data=newmeta_cecum_nonIBD_obese_normal, permutations=10000)
data.adonis
data.adonis=adonis2(bray_dist_cecum_nonIBD_obese_normal ~ Age_at_Collection + Gender + Obesity, data=newmeta_cecum_nonIBD_obese_normal, permutations=10000)
data.adonis


### CD Dysbiosis Index

## Sigmoid

# RF_ASVs_sigmoid=rownames(t(tnewdata_sigmoid_train2))   # Use these three lines if selecting differential ASVs from DESEQ2 that were included in the final RF classifier
# RF_ASV_DESEQ2_sigmoid=DESEQ2_data_sigmoid[RF_ASV_names_sigmoid,]
# RF_ASVs_sigmoid_CD_enriched=RF_ASV_DESEQ2_sigmoid[RF_ASV_DESEQ2_sigmoid$log2FoldChange > 0,]
ASVs_sigmoid_CD_enriched=rownames(DESEQ2_data_sigmoid[DESEQ2_data_sigmoid$log2FoldChange > 0,])
ASVs_sigmoid_CD_depleted=rownames(DESEQ2_data_sigmoid[DESEQ2_data_sigmoid$log2FoldChange < 0,])
enriched=newdata_sigmoid[ASVs_sigmoid_CD_enriched, ]
enriched_sum=colSums(enriched) 
depleted=newdata_sigmoid[ASVs_sigmoid_CD_depleted, ]
depleted_sum=colSums(depleted) 
depleted_sum[depleted_sum==0]=1
CD_dysbiosis_index_sigmoid=log10(enriched_sum/depleted_sum)
newmeta_sigmoid$CD_dysbiosis_index=CD_dysbiosis_index_sigmoid

s17<-factor(newmeta_sigmoid$BMI_known)
split.m17 <- split.data.frame(newmeta_sigmoid, s17)
dim(split.m17$'Y')
dim(split.m17$'N')
newmeta_sigmoid_BMI=split.m17$'Y'
newmeta_sigmoid_BMI$Obesity_CD=factor(newmeta_sigmoid_BMI$Obesity_CD,levels=c("Normal_weight","Overweight","Obese","CD"))
kruskal.test(CD_dysbiosis_index ~ Obesity_CD, data=newmeta_sigmoid_BMI)
dunnTest(CD_dysbiosis_index ~ Obesity_CD, data=newmeta_sigmoid_BMI,method="bh") 
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity_CD, data = newmeta_sigmoid_BMI)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Obesity_CD")

s18<-factor(newmeta_sigmoid$Genotype)
split.m18 <- split.data.frame(newmeta_sigmoid, s18)
dim(split.m18$'Y')
dim(split.m18$'N')
newmeta_sigmoid_genotype=split.m18$'Y'
newmeta_sigmoid_genotype$GRS_CD_all=as.numeric(newmeta_sigmoid_genotype$GRS_CD_all)
newmeta_sigmoid_genotype$Disease=factor(newmeta_sigmoid_genotype$Disease,levels=c("Non-IBD","CD"))
qplot(GRS_CD_all, CD_dysbiosis_index, data = newmeta_sigmoid_genotype, color =Disease, geom=c("point", "smooth"))+theme_pubr()
qplot(GRS_CD_all, Chao1, data = newmeta_sigmoid_genotype, color =Disease, geom=c("point", "smooth"))+theme_pubr()
qplot(GRS_CD_all, Shannon, data = newmeta_sigmoid_genotype, color =Disease, geom=c("point", "smooth"))+theme_pubr()

## Cecum

ASVs_cecum_CD_enriched=rownames(DESEQ2_data_cecum[DESEQ2_data_cecum$log2FoldChange > 0,])
ASVs_cecum_CD_depleted=rownames(DESEQ2_data_cecum[DESEQ2_data_cecum$log2FoldChange < 0,])
enriched=newdata_cecum[ASVs_cecum_CD_enriched, ]
enriched_sum=colSums(enriched) 
depleted=newdata_cecum[ASVs_cecum_CD_depleted, ]
depleted_sum=colSums(depleted)
depleted_sum[depleted_sum==0]=1
CD_dysbiosis_index_cecum=log10(enriched_sum/depleted_sum)
newmeta_cecum$CD_dysbiosis_index=CD_dysbiosis_index_cecum

s19<-factor(newmeta_cecum$BMI_known)
split.m19 <- split.data.frame(newmeta_cecum, s19)
dim(split.m19$'Y')
dim(split.m19$'N')
newmeta_cecum_BMI=split.m19$'Y'
newmeta_cecum_BMI$Obesity_CD=factor(newmeta_cecum_BMI$Obesity_CD,levels=c("Normal_weight","Overweight","Obese","CD"))
kruskal.test(CD_dysbiosis_index ~ Obesity_CD, data=newmeta_cecum_BMI)
dunnTest(CD_dysbiosis_index ~ Obesity_CD, data=newmeta_cecum_BMI,method="bh") 
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity_CD, data = newmeta_cecum_BMI)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Obesity_CD")

par(mfrow = c(1,2))
p1<-ggplot(data=newmeta_sigmoid_BMI,aes(x=Obesity_CD,y=CD_dysbiosis_index,fill=Obesity_CD))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum_BMI,aes(x=Obesity_CD,y=CD_dysbiosis_index,fill=Obesity_CD))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

s20<-factor(newmeta_cecum$Genotype)
split.m20 <- split.data.frame(newmeta_cecum, s20)
dim(split.m20$'Y')
dim(split.m20$'N')
newmeta_cecum_genotype=split.m20$'Y'
newmeta_cecum_genotype$GRS_CD_all=as.numeric(newmeta_cecum_genotype$GRS_CD_all)
newmeta_cecum_genotype$Disease=factor(newmeta_cecum_genotype$Disease,levels=c("Non-IBD","CD"))
qplot(GRS_CD_all, CD_dysbiosis_index, data = newmeta_cecum_genotype, color =Disease, geom=c("point", "smooth"))+theme_pubr()

## Assess CD dysbiosis index by B category
#rerun prior split of newmeta_cecum and newmeta_sigmoid into CD with known B by B_known, then combine with non-IBD
newmeta_sigmoid_nonIBD_Bknown=rbind(newmeta_sigmoid_nonIBD,newmeta_sigmoid_CD_B)
newmeta_cecum_nonIBD_Bknown=rbind(newmeta_cecum_nonIBD,newmeta_cecum_CD_B)
kruskal.test(CD_dysbiosis_index ~ B_combine, data=newmeta_sigmoid_nonIBD_Bknown)
dunnTest(CD_dysbiosis_index ~ B_combine, data=newmeta_sigmoid_nonIBD_Bknown,method="bh") 
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_sigmoid_nonIBD_Bknown)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")
kruskal.test(CD_dysbiosis_index ~ B_combine, data=newmeta_cecum_nonIBD_Bknown)
dunnTest(CD_dysbiosis_index ~ B_combine, data=newmeta_cecum_nonIBD_Bknown,method="bh") 
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_cecum_nonIBD_Bknown)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")
par(mfrow = c(1,2))
p1<-ggplot(data=newmeta_sigmoid_nonIBD_Bknown,aes(x=B_combine,y=CD_dysbiosis_index,fill=B_combine))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum_nonIBD_Bknown,aes(x=B_combine,y=CD_dysbiosis_index,fill=B_combine))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

## Assess CD dysbiosis index by CD disease progression
#rerun prior split of newmeta_cecum and newmeta_sigmoid into CD by Outcome, then combine with non-IBD
newmeta_sigmoid_nonIBD_CD_progression=rbind(newmeta_sigmoid_nonIBD,newmeta_sigmoid_CD_progression_only)
newmeta_cecum_nonIBD_CD_progression=rbind(newmeta_cecum_nonIBD,newmeta_cecum_CD_progression_only)
kruskal.test(CD_dysbiosis_index ~ Progression, data=newmeta_sigmoid_nonIBD_CD_progression)
dunnTest(CD_dysbiosis_index ~ Progression, data=newmeta_sigmoid_nonIBD_CD_progression,method="bh") 
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + B_combine + Progression, data = newmeta_sigmoid_nonIBD_CD_progression)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Progression")
kruskal.test(CD_dysbiosis_index ~ Progression, data=newmeta_cecum_nonIBD_CD_progression)
dunnTest(CD_dysbiosis_index ~ Progression, data=newmeta_cecum_nonIBD_CD_progression,method="bh") 
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + B_combine + Progression, data = newmeta_cecum_nonIBD_CD_progression)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "Progression")
par(mfrow = c(1,2))
p1<-ggplot(data=newmeta_sigmoid_nonIBD_Bknown,aes(x=Progression,y=CD_dysbiosis_index,fill=Progression))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum_nonIBD_Bknown,aes(x=Progression,y=CD_dysbiosis_index,fill=Progression))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

## Assess CD dysbiosis index relationship with CD GRS
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_sigmoid_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + GRS_CD_all, data = newmeta_cecum_nonIBD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data = newmeta_sigmoid_CD_Genotype)  
anova(fit.data)
coef(fit.data)
fit.data <- aov(CD_dysbiosis_index ~ Age_at_Collection + Gender + Obesity + B_combine + GRS_CD_all, data = newmeta_cecum_CD_Genotype)  
anova(fit.data)
coef(fit.data)

m_CD = lm(CD_dysbiosis_index ~ GRS_CD_all, data=newmeta_sigmoid_CD_Genotype)
m_nonIBD = lm(CD_dysbiosis_index ~ GRS_CD_all, data=newmeta_sigmoid_nonIBD_Genotype)
dd_m_CD = data.frame(x=newmeta_sigmoid_CD_Genotype$GRS_CD_all, y=predict(m_CD, newmeta_sigmoid_CD_Genotype), Disease=newmeta_sigmoid_CD_Genotype$Disease)
dd_m_nonIBD = data.frame(x=newmeta_sigmoid_nonIBD_Genotype$GRS_CD_all, y=predict(m_nonIBD, newmeta_sigmoid_nonIBD_Genotype), Disease=newmeta_sigmoid_nonIBD_Genotype$Disease)
ggplot(newmeta_sigmoid_genotype) + geom_point(aes(GRS_CD_all, CD_dysbiosis_index, colour=Disease)) + geom_line(data=dd_m_CD, aes(x, y, colour=Disease)) + geom_line(data=dd_m_nonIBD, aes(x, y, colour=Disease))+theme_pubr()

m_CD = lm(CD_dysbiosis_index ~ GRS_CD_all, data=newmeta_cecum_CD_Genotype)
m_nonIBD = lm(CD_dysbiosis_index ~ GRS_CD_all, data=newmeta_cecum_nonIBD_Genotype)
dd_m_CD = data.frame(x=newmeta_cecum_CD_Genotype$GRS_CD_all, y=predict(m_CD, newmeta_cecum_CD_Genotype), Disease=newmeta_cecum_CD_Genotype$Disease)
dd_m_nonIBD = data.frame(x=newmeta_cecum_nonIBD_Genotype$GRS_CD_all, y=predict(m_nonIBD, newmeta_cecum_nonIBD_Genotype), Disease=newmeta_cecum_nonIBD_Genotype$Disease)
ggplot(newmeta_cecum_genotype) + geom_point(aes(GRS_CD_all, CD_dysbiosis_index, colour=Disease)) + geom_line(data=dd_m_CD, aes(x, y, colour=Disease)) + geom_line(data=dd_m_nonIBD, aes(x, y, colour=Disease))+theme_pubr()


## Plot of GRS_CD_all by B category

ggplot(data=newmeta_sigmoid_genotype,aes(x=B_combine,y=GRS_CD_all,fill=B_combine))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
fit.data <- aov(GRS_CD_all ~ Age_at_Collection + Gender + Obesity + B_combine, data = newmeta_sigmoid_nonIBD_Bknown)  
anova(fit.data)
coef(fit.data)
TukeyHSD(fit.data, "B_combine")
kruskal.test(GRS_CD_all ~ B_combine, data=newmeta_cecum_nonIBD_Bknown)
dunnTest(GRS_CD_all ~ B_combine, data=newmeta_cecum_nonIBD_Bknown,method="bh") 


### Cecum vs. Sigmoid analysis

s22=factor(newmeta$Both_sites)
split.m22 <- split.data.frame(newmeta, s22)
dim(split.m22$'Y')
dim(split.m22$'N')
newmeta_bothsites=split.m22$'Y'
tdata22=t(newdata)
split.s22 <- split.data.frame(tdata22, s22)
dim(split.s22$'Y')
dim(split.s22$'N')
newdata_bothsites=t(as.matrix(split.s22$'Y'))
length(which(colnames(newdata_bothsites)==rownames(newmeta_bothsites)))

## Alpha diversity

s23=factor(newmeta_bothsites$Site)
split.m23 <- split.data.frame(newmeta_bothsites, s23)
dim(split.m23$'Sigmoid')
dim(split.m23$'Cecum')
newmeta_bothsites_sigmoid=split.m23$'Sigmoid'
newmeta_bothsites_cecum=split.m23$'Cecum'
tdata23=t(newdata_bothsites)
split.s23 <- split.data.frame(tdata23, s23)
dim(split.s23$'Sigmoid')
dim(split.s23$'Cecum')
newdata_bothsites_sigmoid=t(as.matrix(split.s23$'Sigmoid'))
newdata_bothsites_cecum=t(as.matrix(split.s23$'Cecum'))
length(which(colnames(newdata_bothsites_sigmoid)==rownames(newmeta_bothsites_sigmoid)))
length(which(colnames(newdata_bothsites_cecum)==rownames(newmeta_bothsites_cecum)))
length(which(newmeta_bothsites_sigmoid$Subject==newmeta_bothsites_cecum$Subject))
newmeta_bothsites_cecum$Chao1_difference=newmeta_bothsites_cecum$Chao1 - newmeta_bothsites_sigmoid$Chao1
newmeta_bothsites_cecum$Shannon_difference=newmeta_bothsites_cecum$Shannon - newmeta_bothsites_sigmoid$Shannon

s24=factor(newmeta_bothsites_sigmoid$Disease)
split.m24 <- split.data.frame(newmeta_bothsites_sigmoid, s24)
dim(split.m24$'CD')
dim(split.m24$'Non-IBD')
newmeta_bothsites_sigmoid_CD=split.m24$'CD'
newmeta_bothsites_sigmoid_nonIBD=split.m24$'Non-IBD'
tdata24=t(newdata_bothsites_sigmoid)
split.s24 <- split.data.frame(tdata24, s24)
dim(split.s24$'CD')
dim(split.s24$'Non-IBD')
newdata_bothsites_sigmoid_CD=t(as.matrix(split.s24$'CD'))
newdata_bothsites_sigmoid_nonIBD=t(as.matrix(split.s24$'Non-IBD'))
s25=factor(newmeta_bothsites_cecum$Disease)
split.m25 <- split.data.frame(newmeta_bothsites_cecum, s25)
dim(split.m25$'CD')
dim(split.m25$'Non-IBD')
newmeta_bothsites_cecum_CD=split.m25$'CD'
newmeta_bothsites_cecum_nonIBD=split.m25$'Non-IBD'
tdata25=t(newdata_bothsites_cecum)
split.s25 <- split.data.frame(tdata25, s25)
dim(split.s25$'CD')
dim(split.s25$'Non-IBD')
newdata_bothsites_cecum_CD=t(as.matrix(split.s25$'CD'))
newdata_bothsites_cecum_nonIBD=t(as.matrix(split.s25$'Non-IBD'))

wilcox.test(newmeta_bothsites_cecum_CD$Chao1, newmeta_bothsites_sigmoid_CD$Chao1, paired=TRUE)
wilcox.test(newmeta_bothsites_cecum_nonIBD$Chao1, newmeta_bothsites_sigmoid_nonIBD$Chao1, paired=TRUE)
wilcox.test(newmeta_bothsites_cecum_CD$Shannon, newmeta_bothsites_sigmoid_CD$Shannon, paired=TRUE)
wilcox.test(newmeta_bothsites_cecum_nonIBD$Shannon, newmeta_bothsites_sigmoid_nonIBD$Shannon, paired=TRUE)

wilcox.test(newmeta_bothsites_cecum_nonIBD$Chao1_difference, newmeta_bothsites_cecum_CD$Chao1_difference, paired=FALSE)
wilcox.test(newmeta_bothsites_cecum_nonIBD$Shannon_difference, newmeta_bothsites_cecum_CD$Shannon_difference, paired=FALSE)

par(mfrow = c(1,2))
newmeta_bothsites_cecum$Disease_ordered<-factor(newmeta_bothsites_cecum$Disease,levels=c("Non-IBD","CD"))
p1<-ggplot(data=newmeta_bothsites_cecum,aes(x=Disease_ordered,y=Chao1_difference,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_bothsites_cecum,aes(x=Disease_ordered,y=Shannon_difference,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

s36<-factor(newmeta_bothsites_sigmoid$BMI_known)
split.m36 <- split.data.frame(newmeta_bothsites_sigmoid, s36)
dim(split.m36$'Y')
dim(split.m36$'N')
newmeta_bothsites_sigmoid_BMI=split.m36$'Y'
newmeta_bothsites_sigmoid_BMI$Obesity_CD=factor(newmeta_bothsites_sigmoid_BMI$Obesity_CD,levels=c("Normal_weight","Overweight","Obese","CD"))
s37<-factor(newmeta_bothsites_cecum$BMI_known)
split.m37 <- split.data.frame(newmeta_bothsites_cecum, s37)
dim(split.m37$'Y')
dim(split.m37$'N')
newmeta_bothsites_cecum_BMI=split.m37$'Y'
newmeta_bothsites_cecum_BMI$Obesity_CD=factor(newmeta_bothsites_cecum_BMI$Obesity_CD,levels=c("Normal_weight","Overweight","Obese","CD"))
kruskal.test(Chao1_difference ~ Obesity_CD, data=newmeta_bothsites_cecum_BMI)
dunnTest(Chao1_difference ~ Obesity_CD, data=newmeta_bothsites_cecum_BMI,method="bh") 
kruskal.test(Shannon_difference ~ Obesity_CD, data=newmeta_bothsites_cecum_BMI)
dunnTest(Shannon_difference ~ Obesity_CD, data=newmeta_bothsites_cecum_BMI,method="bh") 
par(mfrow = c(1,2))
p1<-ggplot(data=newmeta_bothsites_cecum_BMI,aes(x=Obesity_CD,y=Chao1_difference,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_bothsites_cecum_BMI,aes(x=Obesity_CD,y=Shannon_difference,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

## Beta diversity

filtereddata_bothsites=newdata_bothsites[ rowSums(newdata_bothsites > 0) >= 39, ] 
dim(filtereddata_bothsites)
inputdata=t(filtereddata_bothsites)
bray_bothsites=vegdist(inputdata, method="bray")
bray_select=as.matrix(bray_bothsites)[which(newmeta_bothsites$Site == 'Cecum'),which(newmeta_bothsites$Site == 'Sigmoid')]
newmeta_bothsites_cecum$pairwise_bray=diag(bray_select)
# rerun s24 and s25 factorization from alpha diversity section

s26=factor(newmeta_bothsites$Disease)
split.m26 <- split.data.frame(newmeta_bothsites, s26)
dim(split.m26$'CD')
dim(split.m26$'Non-IBD')
newmeta_bothsites_CD=split.m26$'CD'
newmeta_bothsites_nonIBD=split.m26$'Non-IBD'
tdata26=t(newdata_bothsites)
split.s26 <- split.data.frame(tdata26, s26)
dim(split.s26$'CD')
dim(split.s26$'Non-IBD')
newdata_bothsites_CD=t(as.matrix(split.s26$'CD'))
newdata_bothsites_nonIBD=t(as.matrix(split.s26$'Non-IBD'))

wilcox.test(newmeta_bothsites_cecum_nonIBD$pairwise_bray, newmeta_bothsites_cecum_CD$pairwise_bray, paired=FALSE)
par(mfrow = c(1,2))
newmeta_bothsites_cecum$Disease_ordered<-factor(newmeta_bothsites_cecum$Disease,levels=c("Non-IBD","CD"))
p1<-ggplot(data=newmeta_bothsites_cecum,aes(x=Disease_ordered,y=pairwise_bray,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_bothsites_cecum,aes(x=Disease_ordered,y=pairwise_bray,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

kruskal.test(pairwise_bray ~ Obesity_CD, data=newmeta_bothsites_cecum_BMI)
dunnTest(pairwise_bray ~ Obesity_CD, data=newmeta_bothsites_cecum_BMI,method="bh") 
par(mfrow = c(1,2))
p1<-ggplot(data=newmeta_bothsites_cecum_BMI,aes(x=Obesity_CD,y=pairwise_bray,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_bothsites_cecum_BMI,aes(x=Obesity_CD,y=pairwise_bray,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)

data.adonis=adonis(bray_bothsites ~ Batch + Disease + Subject + Site, data=newmeta_bothsites, permutations=10000)
data.adonis
pcoa_output<-pcoa(bray_bothsites)
variance_explained=(pcoa_output$values[,1])/pcoa_output$trace   # use this to identify percent of variance explained
variance_explained

pcoa_coordinates=pcoa_output$vectors
pcoa_output_CD=pcoa_coordinates[which(newmeta_bothsites$Disease == 'CD'),]
PC1=pcoa_output_CD[,1]
PC2=pcoa_output_CD[,2]
plotdata=cbind(PC1,PC2,newmeta_bothsites_CD)
col.list=c("Cecum"="blue","Sigmoid"="red")
pch.list=c("Cecum"=1,"Sigmoid"=16) 
p <- ggplot(plotdata,aes(x=PC1,y=PC2,group=Subject))
p + theme_bw() + geom_point(aes(size='0.1',colour=col.list[paste(plotdata$Site)],shape=pch.list[paste(plotdata$Site)])) + geom_path(size=0.5) + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity() + labs(x = "PC1 (15%)", y = "PC2 (11%)")
p + theme_bw() + geom_point(aes(size='0.1',colour=col.list[paste(plotdata$Site)],shape=16)) + geom_path(size=0.5) + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity() + labs(x = "PC1 (15%)", y = "PC2 (11%)")

pcoa_output_nonIBD=pcoa_coordinates[which(newmeta_bothsites$Disease == 'Non-IBD'),]
PC1=pcoa_output_nonIBD[,1]
PC2=pcoa_output_nonIBD[,2]
plotdata=cbind(PC1,PC2,newmeta_bothsites_nonIBD)
col.list=c("Cecum"="blue","Sigmoid"="red")
pch.list=c("Cecum"=1,"Sigmoid"=16) 
p <- ggplot(plotdata,aes(x=PC1,y=PC2,group=Subject))
p + theme_bw() + geom_point(aes(size='0.1',colour=col.list[paste(plotdata$Site)],shape=pch.list[paste(plotdata$Site)])) + geom_path(size=0.5) + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity() + labs(x = "PC1 (15%)", y = "PC2 (11%)")
p + theme_bw() + geom_point(aes(size='0.1',colour=col.list[paste(plotdata$Site)],shape=16)) + geom_path(size=0.5) + theme(aspect.ratio=1/1.12) + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+scale_color_identity()+scale_shape_identity() + labs(x = "PC1 (15%)", y = "PC2 (11%)")

## Differential taxa

filtereddata_bothsites_CD<-newdata_bothsites_CD[ rowSums(newdata_bothsites_CD > 0) >= 43, ] 
OTU_bothsites_CD=otu_table(filtereddata_bothsites_CD, taxa_are_rows = TRUE)
biom_bothsites_CD = phyloseq(OTU_bothsites_CD,TAX)
map_bothsites_CD=sample_data(newmeta_bothsites_CD)
data_bothsites_CD=merge_phyloseq(biom_bothsites_CD,map_bothsites_CD)
print(data_bothsites_CD)
diagdds=phyloseq_to_deseq2(data_bothsites_CD, ~ Subject + Site)
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Site","Cecum","Sigmoid"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_bothsites_CD)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - CD Cecum vs. Sigmoid adjusted batch gender obesity - ASVs s43.csv") 
DESEQ2_data_bothsites_CD = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_bothsites_CD = DESEQ2_data_bothsites_CD[DESEQ2_data_bothsites_CD$Abundance >= 0.00001, ] 
DESEQ2_data_bothsites_CD = DESEQ2_data_bothsites_CD[!is.na(DESEQ2_data_bothsites_CD$Order), ] 
y = tapply(DESEQ2_data_bothsites_CD$log2FoldChange, DESEQ2_data_bothsites_CD$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_bothsites_CD$Genus= factor(as.character(DESEQ2_data_bothsites_CD$Genus), levels = names(y))
dev.new(width=8, height=4)
ggplot(DESEQ2_data_bothsites_CD, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")

filtereddata_bothsites_nonIBD<-newdata_bothsites_nonIBD[ rowSums(newdata_bothsites_nonIBD > 0) >= 54, ] 
OTU_bothsites_nonIBD=otu_table(filtereddata_bothsites_nonIBD, taxa_are_rows = TRUE)
biom_bothsites_nonIBD = phyloseq(OTU_bothsites_nonIBD,TAX)
map_bothsites_nonIBD=sample_data(newmeta_bothsites_nonIBD)
data_bothsites_nonIBD=merge_phyloseq(biom_bothsites_nonIBD,map_bothsites_nonIBD)
print(data_bothsites_nonIBD)
diagdds=phyloseq_to_deseq2(data_bothsites_nonIBD, ~ Subject + Site)  
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
CDvsNorm = results(diagdds, independentFiltering=TRUE, contrast=c("Site","Cecum","Sigmoid"))
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix = cbind(as(CDvsNorm, "data.frame"), as(tax_table(data_bothsites_nonIBD)[rownames(CDvsNorm), ], "matrix"))   # need to update data file name
CDvsNormMatrix$Abundance=(CDvsNormMatrix$baseMean)/sum(CDvsNormMatrix$baseMean)
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Genus), paste0(CDvsNormMatrix$Family,' (f)'), paste0(CDvsNormMatrix$Genus))
CDvsNormMatrix$Genus = ifelse(is.na(CDvsNormMatrix$Family), paste0(CDvsNormMatrix$Order,' (o)'), paste0(CDvsNormMatrix$Genus))
qobj=qvalue(CDvsNormMatrix[,"pvalue"],pi0.method="bootstrap",robust=FALSE)
CDvsNormMatrix=cbind(CDvsNormMatrix,qobj$qvalues)
names(CDvsNormMatrix)[names(CDvsNormMatrix) == 'qobj$qvalues'] <- 'Qvalue'
write.csv(CDvsNormMatrix, "DESEq2 - Non-IBD Cecum vs. Sigmoid adjusted subject - ASVs s54.csv") 
DESEQ2_data_bothsites_nonIBD = CDvsNormMatrix[(CDvsNormMatrix$Qvalue < alpha), ]  
DESEQ2_data_bothsites_nonIBD = DESEQ2_data_bothsites_nonIBD[DESEQ2_data_bothsites_nonIBD$Abundance >= 0.00001, ] 
DESEQ2_data_bothsites_nonIBD = DESEQ2_data_bothsites_nonIBD[!is.na(DESEQ2_data_bothsites_nonIBD$Order), ] 
y = tapply(DESEQ2_data_bothsites_nonIBD$log2FoldChange, DESEQ2_data_bothsites_nonIBD$Genus, function(y) max(y))
y = sort(y, FALSE)   #switch to TRUE to reverse direction
DESEQ2_data_bothsites_nonIBD$Genus= factor(as.character(DESEQ2_data_bothsites_nonIBD$Genus), levels = names(y))
dev.new(width=8, height=4)
ggplot(DESEQ2_data_bothsites_nonIBD, aes(x = log2FoldChange, y = Genus, color = Phylum)) + geom_point(aes(size = sqrt(Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.00001),sqrt(0.3)),breaks=c(sqrt(0.0001),sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")




#### Predicted metagenomics

# Performed in QIIME2-2019.10
qiime tools import --input-path Deblur_merged_final.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path table.qza
qiime fragment-insertion sepp --i-representative-sequences rep-seqs.qza --p-threads 4 --i-reference-database picrust2_default_sepp_ref.qza --output-dir output_directory
qiime picrust2 custom-tree-pipeline --i-table table.qza --i-tree output_directory/tree.qza --output-dir q2-picrust2_output --p-threads 4 --p-hsp-method mp --p-max-nsti 2 --verbose
qiime tools export --input-path q2-picrust2_output/ko_metagenome.qza --output-path q2-picrust2_output/ko_metagenome
biom convert -i q2-picrust2_output/ko_metagenome/feature-table.biom -o q2-picrust2_output/ko_metagenome/feature-table.txt --to-tsv

KO_data<-read.csv(choose.files(),header=T,row.names=1)  
tKO_data=t(KO_data)
length(which(rownames(t(KO_data))==rownames(meta)))

split.s27 <- split.data.frame(tKO_data, s1)
dim(split.s27$'Y')
dim(split.s27$'N')
KO_newdata=t(as.matrix(split.s27$'Y'))
length(which(colnames(KO_newdata)==rownames(newmeta)))

tdata28=t(KO_newdata)
split.s28 <- split.data.frame(tdata28, s2)
dim(split.s28$'Sigmoid')
dim(split.s28$'Cecum')
KO_newdata_sigmoid=t(as.matrix(split.s28$'Sigmoid'))
KO_newdata_cecum=t(as.matrix(split.s28$'Cecum'))
length(which(colnames(KO_newdata_sigmoid)==rownames(newmeta_sigmoid)))
length(which(colnames(KO_newdata_cecum)==rownames(newmeta_cecum)))

KEGG_bile<-read.csv(choose.files(),header=T,row.names=1)  
KEGG_bile_filter=rownames(KEGG_bile)
dim(KO_newdata_sigmoid)
relativeKOsigmoiddata<-as.data.frame(sweep(KO_newdata_sigmoid,2,colSums(KO_newdata_sigmoid),"/"))
KEGG_bile_sigmoid<-t(relativeKOsigmoiddata[KEGG_bile_filter,])
dim(KEGG_bile_sigmoid)
write.csv(KEGG_bile_sigmoid, "Sigmoid - bile acid KO relative abundances.csv") 
dim(KO_newdata_cecum)
relativeKOcecumdata<-as.data.frame(sweep(KO_newdata_cecum,2,colSums(KO_newdata_cecum),"/"))
KEGG_bile_cecum<-t(relativeKOcecumdata[KEGG_bile_filter,])
dim(KEGG_bile_cecum)
write.csv(KEGG_bile_cecum, "Cecum - bile acid KO relative abundances.csv") 

length(which(rownames(KEGG_bile_sigmoid)==rownames(newmeta_sigmoid)))
newmeta_sigmoid$Bai=rowSums(KEGG_bile_sigmoid)
write.csv(newmeta_sigmoid, "Sigmoid - Bai bile acid gene relative abundances.csv") 
kruskal.test(Bai ~ Disease, data=newmeta_sigmoid)
length(which(rownames(KEGG_bile_cecum)==rownames(newmeta_cecum)))
newmeta_cecum$Bai=rowSums(KEGG_bile_cecum)
write.csv(newmeta_cecum, "Cecum - Bai bile acid gene relative abundances.csv") 
kruskal.test(Bai ~ Disease, data=newmeta_cecum)
p1<-ggplot(data=newmeta_sigmoid,aes(x=Disease,y=Bai,fill=Disease))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=newmeta_cecum,aes(x=Disease,y=Bai,fill=Disease))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
par(mfrow = c(1,2))
grid.arrange(p1, p2, nrow = 1)



#### Metabolomics

library("MetaboAnalystR")
source("D:/Dropbox/R/R scripts/KNN.obs.sel.r")
library(sva)
library("limma")
library("qvalue")
library(ape)
library(vegan)
library(caret)
library(pROC)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(DESeq2)
library(WGCNA)
library(reshape)

#save.image(file='Metabolomics_MLI_analysis.RData')
load('Metabolomics_MLI_analysis.RData')

meta<-read.csv(choose.files(),header=T,row.names=1)  #mapping file
rawdata_POS<-read.csv(choose.files(),header=T,row.names=1) 
dim(rawdata_POS)
filteredrawdata_POS<-rawdata_POS[ rowSums(is.na(rawdata_POS)) <= 347, ]  # filter metabolites detected in less than 10% of samples
dim(filteredrawdata_POS)
dat=t(filteredrawdata_POS) # algorithm requires metabolites to be in columns (referred to as "variables") and samples as rows (referred to as "observations")
results <- impute.knn.obs.sel(dat, K=10)
imputed_data_POS=t(results)
write.csv(imputed_data_POS, "Positive_ESI_metabolomics_data__filtered_10percent_KNN_obs_sel_K10_imputation.csv")
s1<-factor(meta$POS)
split.m <- split.data.frame(meta, s1)
dim(split.m$'Y')
dim(split.m$'N')
meta_POS=split.m$'Y'
length(which(rownames(t(imputed_data_POS))==rownames(meta_POS)))
rawdata_NEG<-read.csv(choose.files(),header=T,row.names=1) 
dim(rawdata_NEG)
filteredrawdata_NEG<-rawdata_NEG[ rowSums(is.na(rawdata_NEG)) <= 349, ]  # filter metabolites detected in less than 10% of samples
dim(filteredrawdata_NEG)
dat=t(filteredrawdata_NEG) # algorithm requires metabolites to be in columns (referred to as "variables") and samples as rows (referred to as "observations")
results <- impute.knn.obs.sel(dat, K=10)
imputed_data_NEG=t(results)
write.csv(imputed_data_NEG, "Negative_ESI_metabolomics_data_filtered_10percent_KNN_obs_sel_K10_imputation.csv")
s2<-factor(meta$NEG)
split.m2 <- split.data.frame(meta, s2)
dim(split.m2$'Y')
dim(split.m2$'N')
meta_NEG=split.m2$'Y'
length(which(rownames(t(imputed_data_NEG))==rownames(meta_NEG)))

## Batch correction by ESI mode

batch_input_POS=cbind(meta_POS$Disease,t(imputed_data_POS))
write.csv(batch_input_POS, "temp.csv")

metabo_batch_POS <- InitDataObjects("pktable", "stat", FALSE)
metabo_batch_POS<-Read.TextData(metabo_batch_POS, "temp.csv", "rowu", "disc")  # "disc"=discrete variable, can also make continuous; rowu= samples in rows and unpaired
metabo_batch_POS<-SanityCheckData(metabo_batch_POS)
metabo_batch_POS<-ReplaceMin(metabo_batch_POS)
metabo_batch_POS<-FilterVariable(metabo_batch_POS, "iqr", "F", 25)
metabo_batch_POS<-PreparePrenormData(metabo_batch_POS)
metabo_batch_POS<-Normalization(metabo_batch_POS, "QuantileNorm", "LogNorm", "NULL", ref=NULL, ratio=FALSE, ratioNum=20)  #  normalize by quantiles, log transform data, then mean center
normalized_set <- metabo_batch_POS[["dataSet"]][["norm"]]
ordered_normalized_set <- normalized_set[order(row.names(normalized_set)), ]
new_normalized_set <- cbind(meta_POS$Disease, meta_POS$Batch_POS, meta_POS$Injection_Order_POS, ordered_normalized_set)
write.csv(new_normalized_set,file = "Positive_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized.csv") ## make sure to add name in upper left cell
metabo_batch_POS<- PCA.Anal(metabo_batch_POS)
metabo_batch_POS<- PlotPCAPairSummary(metabo_batch_POS, "pca_pair_0_", "png", 72, width=NA, 2)  # change file name as needed
metabo_batch_POS<- PlotPCAScree(metabo_batch_POS, "pca_scree_0_", "png", 72, width=NA, 2)
metabo_batch_POS<- PlotPCA2DScore(metabo_batch_POS, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
metabo_batch_POS<- PlotPCALoading(metabo_batch_POS, "pca_loading_0_", "png", 72, width=NA, 1,2)
metabo_batch_POS<- PlotPCABiplot(metabo_batch_POS, "pca_biplot_0_", "png", 72, width=NA, 1,2)

metabo_norm_batch_POS <- InitDataObjects("pktable", "utils", FALSE)
metabo_norm_batch_POS <- Read.BatchDataTB(metabo_norm_batch_POS, "Positive_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized.csv", "row")
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="Combat")  # Can't do auto due to error
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="EigenMS")  
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="QC RLSC")  
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="ANCOVA") 
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="RUV random")
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="RUV 2")
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="RUV")
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="CCMN")
metabo_norm_batch_POS$dataSet$interbatch_dis
metabo_norm_batch_POS <- PerformBatchCorrection(metabo_norm_batch_POS, Method="WaveICA")
metabo_norm_batch_POS$dataSet$interbatch_dis
data_corrected <- read.csv("MetaboAnalyst_batch_data.csv")
data_corrected_new <- data_corrected[,-3]
write.csv(data_corrected_new,file = "Positive_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized_Combat_batch_corrected.csv",row.names =  F)
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "Positive_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized_Combat_batch_corrected.csv", "row", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20) # No need to normalize again
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)

rownames(data_corrected_new)<-rownames(t(imputed_data_POS))
s3<-factor(meta_POS$Site)
split.m3 <- split.data.frame(meta_POS, s3)
dim(split.m3$'Cecum')
dim(split.m3$'Sigmoid')
meta_cecum_POS=split.m3$'Cecum'
meta_sigmoid_POS=split.m3$'Sigmoid'
split.s3 <- split.data.frame(data_corrected_new, s3)
dim(split.s3$'Cecum')
dim(split.s3$'Sigmoid')
data_corrected_sigmoid_POS=as.matrix(split.s3$'Sigmoid')
data_corrected_sigmoid_POS<- data_corrected_sigmoid_POS[,-c(1,2)]
data_corrected_cecum_POS=as.matrix(split.s3$'Cecum')
data_corrected_cecum_POS<- data_corrected_cecum_POS[,-c(1,2)]
length(which(rownames(data_corrected_cecum_POS)==rownames(meta_cecum_POS)))
length(which(rownames(data_corrected_sigmoid_POS)==rownames(meta_sigmoid_POS)))

batch_input_NEG=cbind(meta_NEG$Disease,t(imputed_data_NEG))
write.csv(batch_input_NEG, "temp.csv")

metabo_batch_NEG <- InitDataObjects("pktable", "stat", FALSE)
metabo_batch_NEG<-Read.TextData(metabo_batch_NEG, "temp.csv", "rowu", "disc")  # "disc"=discrete variable, can also make continuous; rowu= samples in rows and unpaired
metabo_batch_NEG<-SanityCheckData(metabo_batch_NEG)
metabo_batch_NEG<-ReplaceMin(metabo_batch_NEG)
metabo_batch_NEG<-FilterVariable(metabo_batch_NEG, "iqr", "F", 25)
metabo_batch_NEG<-PreparePrenormData(metabo_batch_NEG)
metabo_batch_NEG<-Normalization(metabo_batch_NEG, "QuantileNorm", "LogNorm", "NULL", ref=NULL, ratio=FALSE, ratioNum=20)  #  normalize by quantiles, log transform data, then mean center
normalized_set <- metabo_batch_NEG[["dataSet"]][["norm"]]
ordered_normalized_set <- normalized_set[order(row.names(normalized_set)), ]
new_normalized_set <- cbind(meta_NEG$Disease, meta_NEG$Batch_NEG, meta_NEG$Injection_Order_NEG, ordered_normalized_set)
write.csv(new_normalized_set,file = "Negative_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized.csv") ## make sure to add name in upper left cell
metabo_batch_NEG<- PCA.Anal(metabo_batch_NEG)
metabo_batch_NEG<- PlotPCAPairSummary(metabo_batch_NEG, "pca_pair_0_", "png", 72, width=NA, 2)  # change file name as needed
metabo_batch_NEG<- PlotPCAScree(metabo_batch_NEG, "pca_scree_0_", "png", 72, width=NA, 2)
metabo_batch_NEG<- PlotPCA2DScore(metabo_batch_NEG, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
metabo_batch_NEG<- PlotPCALoading(metabo_batch_NEG, "pca_loading_0_", "png", 72, width=NA, 1,2)
metabo_batch_NEG<- PlotPCABiplot(metabo_batch_NEG, "pca_biplot_0_", "png", 72, width=NA, 1,2)

metabo_norm_batch_NEG <- InitDataObjects("pktable", "utils", FALSE)
metabo_norm_batch_NEG <- Read.BatchDataTB(metabo_norm_batch_NEG, "Negative_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized.csv", "row")
metabo_norm_batch_NEG <- PerformBatchCorrection(metabo_norm_batch_NEG, Method="Combat")  # Can't do auto due to error
metabo_norm_batch_NEG$dataSet$interbatch_dis
data_corrected <- read.csv("MetaboAnalyst_batch_data.csv")
data_corrected_new <- data_corrected[,-3]
write.csv(data_corrected_new,file = "Negative_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized_Combat_batch_corrected.csv",row.names =  F)
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "Negative_ESI_metabolomics_data_filtered_10percent_KNN_quantile_log_normalized_Combat_batch_corrected.csv", "row", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20) # No need to normalize again
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)

s4<-factor(meta_NEG$Site)
split.m4 <- split.data.frame(meta_NEG, s4)
dim(split.m4$'Cecum')
dim(split.m4$'Sigmoid')
meta_sigmoid_NEG=split.m4$'Sigmoid'
meta_cecum_NEG=split.m4$'Cecum'
split.s4 <- split.data.frame(data_corrected_new, s4)
dim(split.s4$'Sigmoid')
dim(split.s4$'Cecum')
data_corrected_sigmoid_NEG=as.matrix(split.s4$'Sigmoid')
data_corrected_sigmoid_NEG<- data_corrected_sigmoid[,-c(1,2)]
data_corrected_cecum_NEG=as.matrix(split.s4$'Cecum')
data_corrected_cecum_NEG<- data_corrected_cecum_NEG[,-c(1,2)]
length(which(rownames(data_corrected_sigmoid_NEG)==rownames(meta_sigmoid_NEG)))
length(which(rownames(data_corrected_cecum_NEG)==rownames(meta_cecum_NEG)))


## PCoA, Differential abundance, and pathway analysis - Cecum

s5<-factor(meta_cecum_NEG$POS_NEG)
split.m5 <- split.data.frame(meta_cecum_NEG, s5)
dim(split.m5$'Y')
dim(split.m5$'N')
meta_cecum_NEG_Y=split.m5$'Y'
split.s5 <- split.data.frame(data_corrected_cecum_NEG, s5)
dim(split.s5$'Y')
dim(split.s5$'N')
data_corrected_cecum_NEG_Y=split.s5$'Y'
s6<-factor(meta_cecum_POS$POS_NEG)
split.m6 <- split.data.frame(meta_cecum_POS, s6)
dim(split.m6$'Y')
dim(split.m6$'N')
meta_cecum_POS_Y=split.m6$'Y'
split.s6 <- split.data.frame(data_corrected_cecum_POS, s6)
dim(split.s6$'Y')
dim(split.s6$'N')
data_corrected_cecum_POS_Y=split.s6$'Y'
length(which(rownames(data_corrected_cecum_NEG_Y)==rownames(data_corrected_cecum_POS_Y)))
data_corrected_cecum_NEGandPOS = cbind(data_corrected_cecum_NEG_Y, data_corrected_cecum_POS_Y)
dim(data_corrected_cecum_NEGandPOS)
write.csv(data_corrected_cecum_NEGandPOS, "Cecum_batch_corrected_NEGandPOS.csv")

pca_input=read.csv("Cecum_batch_corrected_NEGandPOS.csv",header=T,row.names=1) 
pca_output<-prcomp(pca_input)
PoV <- pca_output$sdev^2/sum(pca_output$sdev^2)
PoV
col.list=c("CD"="red","Non-IBD"="blue")
plot(pca_output$rotation[,1],pca_output$rotation[,2],col=col.list[paste(meta_sigmoid_POS_Y$Disease)],pch=16,xlab = "PC1 (19%)", ylab = "PC2 (7%)", axes = TRUE, main = "PCoA")

euclidean_dist_cecum_NEGandPOS <-dist(data_corrected_cecum_NEGandPOS, method="euclidean")
data.adonis=adonis2(euclidean_dist_cecum_NEGandPOS~ Age_at_Collection + Gender + Obesity + Disease, data=meta_cecum_POS_Y, permutations=10000)
data.adonis
data.adonis=adonis2(euclidean_dist_cecum_NEGandPOS~ Gender + Obesity + Disease, data=meta_cecum_POS_Y, permutations=10000)
data.adonis

Disease=factor(meta_cecum_POS_Y$Disease,levels=c("Non-IBD","CD"))
Obesity=meta_cecum_POS_Y$Obesity
Gender=meta_cecum_POS_Y$Gender
input=t(data_corrected_cecum_NEGandPOS)
write.csv(input, "temp.csv")
inputdata<-read.csv("temp.csv",header=T,row.names=1) 

design <- model.matrix(~ 0 + Gender + Obesity + Disease)
fit <- lmFit(inputdata, design)
fit2 <- eBayes(fit,robust="TRUE")
fit2

topTable(fit2, coef=6, adjust="BH")  #make sure to look at fit2 to select the proper coef#
a<-fit2$p.value
b<-fit2$coefficients
c<-fit2$Amean
d<-a[,6]                 
e<-b[,6]                 
f<-fit2$t
g<-f[,6]
q1=qvalue(d,pi0.method="bootstrap",robust=FALSE)
Qvalues=q1$qvalues
# Qvalues=p.adjust(d, method="BH")
Export=cbind(c,e,g,d,Qvalues) 
colnames(Export)=c("Mean","Log2FC","T-statistic","P-value","Q-value")
write.csv(Export, "Limma_CD_vs_nonIBD_NEGandPOS_Cecum_adjusted_SexObesity.csv")

# create txt file with mz in first column, p-values from limma in second column, retention time in 3rd column, and t scores in 4th column

metabo_norm_batch_mummichog_NEGandPOS<-InitDataObjects("mass_all", "mummichog", FALSE)
SetPeakFormat("mprt")
metabo_norm_batch_mummichog_NEGandPOS<-UpdateInstrumentParameters(metabo_norm_batch_mummichog_NEGandPOS, 10, "mixed");
metabo_norm_batch_mummichog_NEGandPOS<-Read.PeakListData(metabo_norm_batch_mummichog_NEGandPOS, "Mummichog_input_cecum_NEGandPOS.txt");
metabo_norm_batch_mummichog_NEGandPOS<-SanityCheckMummichogData(metabo_norm_batch_mummichog_NEGandPOS)
metabo_norm_batch_mummichog_NEGandPOS<-SetPeakEnrichMethod(metabo_norm_batch_mummichog_NEGandPOS, "mum", version="v2")   # "integ"
metabo_norm_batch_mummichog_NEGandPOS<-SetMummichogPval(metabo_norm_batch_mummichog_NEGandPOS, 0.04)  # adjust p-value below 0.05 such that approximately 10% of peaks are significant
metabo_norm_batch_mummichog_NEGandPOS<-PerformPSEA(metabo_norm_batch_mummichog_NEGandPOS, "hsa_mfn", "current", permNum =100)  # hsa_kegg, eco_kegg 
metabo_norm_batch_mummichog_NEGandPOS<-PlotIntegPaths(metabo_norm_batch_mummichog_NEGandPOS, "peaks_to_paths_", "png", 300)
metabo_norm_batch_mummichog_NEGandPOS$integ.resmat
# Use this function to view a table of the significant / non-significant compound hits and m/z matching details in a selected pathway
metabo_norm_batch_mummichog_NEGandPOS <- GetMummichogPathSetDetails(metabo_norm_batch_mummichog_NEGandPOS, "Bile acid biosynthesis")  # creates new file with hits
# Use this function to view a table of matching details (i.e. adducts, t-scores, m.z) for a compound
metabo_norm_batch_mummichog_NEGandPOS <- GetCompoundDetails(metabo_norm_batch_mummichog_NEGandPOS, "C02528")
# Note p-values from mummichog are permutation-based estimates of FDR using permutation of dataset


## PCoA, Differential abundance, and pathway analysis - Sigmoid

s7<-factor(meta_sigmoid_NEG$POS_NEG)
split.m7 <- split.data.frame(meta_sigmoid_NEG, s7)
dim(split.m7$'Y')
dim(split.m7$'N')
meta_sigmoid_NEG_Y=split.m7$'Y'
split.s7 <- split.data.frame(data_corrected_sigmoid_NEG, s7)
dim(split.s7$'Y')
dim(split.s7$'N')
data_corrected_sigmoid_NEG_Y=split.s7$'Y'
s8<-factor(meta_sigmoid_POS$POS_NEG)
split.m8 <- split.data.frame(meta_sigmoid_POS, s8)
dim(split.m8$'Y')
dim(split.m8$'N')
meta_sigmoid_POS_Y=split.m8$'Y'
split.s8 <- split.data.frame(data_corrected_sigmoid_POS, s8)
dim(split.s8$'Y')
dim(split.s8$'N')
data_corrected_sigmoid_POS_Y=split.s8$'Y'
length(which(rownames(data_corrected_sigmoid_NEG_Y)==rownames(data_corrected_sigmoid_POS_Y)))
data_corrected_sigmoid_NEGandPOS = cbind(data_corrected_sigmoid_NEG_Y, data_corrected_sigmoid_POS_Y)
dim(data_corrected_sigmoid_NEGandPOS)
write.csv(data_corrected_sigmoid_NEGandPOS, "Sigmoid_batch_corrected_NEGandPOS.csv")

pca_input=read.csv("Sigmoid_batch_corrected_NEGandPOS.csv",header=T,row.names=1) 
pca_output<-prcomp(pca_input)
PoV <- pca_output$sdev^2/sum(pca_output$sdev^2)
PoV
col.list=c("CD"="red","Non-IBD"="blue")
plot(pca_output$rotation[,1],pca_output$rotation[,2],col=col.list[paste(meta_sigmoid_POS_Y$Disease)],pch=16,xlab = "PC1 (19%)", ylab = "PC2 (7%)", axes = TRUE, main = "PCoA")

euclidean_dist_sigmoid_NEGandPOS <-dist(data_corrected_sigmoid_NEGandPOS, method="euclidean")
data.adonis=adonis2(euclidean_dist_sigmoid_NEGandPOS~ Age_at_Collection + Gender + Obesity + Disease, data=meta_sigmoid_POS_Y, permutations=10000)
data.adonis
data.adonis=adonis2(euclidean_dist_sigmoid_NEGandPOS~ Gender + Obesity + Disease, data=meta_sigmoid_POS_Y, permutations=10000)
data.adonis

Disease=factor(meta_sigmoid_POS_Y$Disease,levels=c("Non-IBD","CD"))
Obesity=meta_sigmoid_POS_Y$Obesity
Gender=meta_sigmoid_POS_Y$Gender
input=t(data_corrected_sigmoid_NEGandPOS)
write.csv(input, "temp.csv")
inputdata<-read.csv("temp.csv",header=T,row.names=1) 

design <- model.matrix(~ 0 + Gender + Obesity + Disease)
fit <- lmFit(inputdata, design)
fit2 <- eBayes(fit,robust="TRUE")
fit2

topTable(fit2, coef=6, adjust="BH")  #make sure to look at fit2 to select the proper coef#
a<-fit2$p.value
b<-fit2$coefficients
c<-fit2$Amean
d<-a[,6]                 
e<-b[,6]                  
f<-fit2$t
g<-f[,6]
q1=qvalue(d,pi0.method="bootstrap",robust=FALSE)
Qvalues=q1$qvalues
# Qvalues=p.adjust(d, method="BH")
Limma_sigmoid_NEGandPOS=cbind(c,e,g,d,Qvalues)
colnames(Limma_sigmoid_NEGandPOS)=c("Mean","Log2FC","T-statistic","P-value","Qvalue")
Limma_sigmoid_NEGandPOS=as.data.frame(Limma_sigmoid_NEGandPOS)
write.csv(Limma_sigmoid_NEGandPOS, "Limma_CD_vs_nonIBD_NEGandPOS_Sigmoid_adjusted_SexObesity.csv")

# create txt file with mz in first column, p-values from limma in second column, retention time in 3rd column, and t scores in 4th column

metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-InitDataObjects("mass_all", "mummichog", FALSE)
SetPeakFormat("mprt")
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-UpdateInstrumentParameters(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, 10, "mixed");
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-Read.PeakListData(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, "Mummichog_input_sigmoid_NEGandPOS.txt");
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-SanityCheckMummichogData(metabo_norm_batch_mummichog_sigmoid_NEGandPOS)
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-SetPeakEnrichMethod(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, "mum", version="v2")
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-SetMummichogPval(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, 0.02)
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-PerformPSEA(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, "hsa_mfn", "current", permNum =100)
metabo_norm_batch_mummichog_sigmoid_NEGandPOS<-PlotIntegPaths(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, "peaks_to_paths_", "png", 300)
metabo_norm_batch_mummichog_sigmoid_NEGandPOS$integ.resmat
# Use this function to view a table of the significant / non-significant compound hits and m/z matching details in a selected pathway
metabo_norm_batch_mummichog_sigmoid_NEGandPOS <- GetMummichogPathSetDetails(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, "Bile acid biosynthesis")  # creates new file with hits
# Use this function to view a table of matching details (i.e. adducts, t-scores, m.z) for a compound
metabo_norm_batch_mummichog_sigmoid_NEGandPOS <- GetCompoundDetails(metabo_norm_batch_mummichog_sigmoid_NEGandPOS, "C02528")


## Comparison of cecum-sigmoid bile acid distances between CD and non-IBD

selecttemp=as.data.frame(data_corrected_cecum_NEGandPOS)
drop <- c("N_514.3238336_2.402429147")
data_corrected_cecum_NEGandPOS_new=selecttemp[ , !(names(selecttemp) %in% drop)]
length(which(colnames(data_corrected_cecum_NEGandPOS_new)==colnames(data_corrected_sigmoid_NEGandPOS)))
data_corrected_NEGandPOS=rbind(data_corrected_sigmoid_NEGandPOS,data_corrected_cecum_NEGandPOS_new)
data_corrected_NEGandPOS <- data_corrected_NEGandPOS[ order(row.names(data_corrected_NEGandPOS)), ]
rownames(data_corrected_NEGandPOS)
length(which(colnames(meta_sigmoid_POS_Y)==colnames(meta_cecum_POS_Y)))
meta_combined_NEGandPOS=rbind(meta_sigmoid_POS_Y,meta_cecum_POS_Y)
meta_combined_NEGandPOS <- meta_combined_NEGandPOS[ order(row.names(meta_combined_NEGandPOS)), ]
length(which(rownames(meta_combined_NEGandPOS)==rownames(data_corrected_NEGandPOS)))

s11=factor(meta_combined_NEGandPOS$BothSites)
split.m11 <- split.data.frame(meta_combined_NEGandPOS, s11)
dim(split.m11$'Y')
dim(split.m11$'N')
meta_combined_NEGandPOS_bothsites=split.m11$'Y'
split.s11 <- split.data.frame(data_corrected_NEGandPOS, s11)
dim(split.s11$'Y')
dim(split.s11$'N')
data_corrected_NEGandPOS_bothsites=split.s11$'Y'

s12=factor(meta_combined_NEGandPOS_bothsites$Site)
split.m12 <- split.data.frame(meta_combined_NEGandPOS_bothsites, s12)
dim(split.m12$'Cecum')
dim(split.m12$'Sigmoid')
meta_combined_NEGandPOS_bothsites_cecum=split.m12$'Cecum'

write.csv(meta_combined_NEGandPOS, "temp.csv")
meta_combined_NEGandPOS<-read.csv(choose.files(),header=T,row.names=1)  

euclidean_dist_combined_NEGandPOS <-dist(data_corrected_NEGandPOS_bothsites, method="euclidean")
euclidean_select=as.matrix(euclidean_dist_combined_NEGandPOS)[which(meta_combined_NEGandPOS_bothsites$Site == 'Cecum'),which(meta_combined_NEGandPOS_bothsites$Site == 'Sigmoid')]
meta_combined_NEGandPOS_bothsites_cecum$pairwise_euclidean=diag(euclidean_select)

s13=factor(meta_combined_NEGandPOS_bothsites_cecum$Disease)
split.m13 <- split.data.frame(meta_combined_NEGandPOS_bothsites_cecum, s13)
dim(split.m13$'CD')
dim(split.m13$'Non-IBD')
meta_combined_NEGandPOS_bothsites_cecum_CD=split.m13$'CD'
meta_combined_NEGandPOS_bothsites_cecum_nonIBD=split.m13$'Non-IBD'

wilcox.test(as.numeric(meta_combined_NEGandPOS_bothsites_cecum_nonIBD$pairwise_euclidean), as.numeric(meta_combined_NEGandPOS_bothsites_cecum_CD$pairwise_euclidean), paired=FALSE)        # p=0.2245

par(mfrow = c(1,2))
meta_combined_NEGandPOS_bothsites_cecum$Disease_ordered<-factor(meta_combined_NEGandPOS_bothsites_cecum$Disease,levels=c("Non-IBD","CD"))
p1<-ggplot(data=meta_combined_NEGandPOS_bothsites_cecum,aes(x=Disease_ordered,y=pairwise_euclidean,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=meta_combined_NEGandPOS_bothsites_cecum,aes(x=Disease_ordered,y=pairwise_euclidean,fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)


## Cholic acid / Deoxycholic acid log ratio

# P_391.2869_8.12 = cholic acid
# P_393.2998_8.13 = deoxycholic acid
# removed batch 2013_05_09
a=as.numeric(data_corrected_sigmoid_NEGandPOS[,'P_391.2869_8.12'])
b=as.numeric(data_corrected_sigmoid_NEGandPOS[,'P_393.2998_8.13'])
logratio_sigmoid= a - b
length(which(rownames(data_corrected_sigmoid_NEGandPOS)==rownames(meta_sigmoid_POS_Y)))
logratio_bydisease_sigmoid=cbind(meta_sigmoid_POS_Y$Disease,logratio_sigmoid)
rownames(logratio_bydisease_sigmoid)<-rownames(meta_sigmoid_POS_Y)
colnames(logratio_bydisease_sigmoid)<-c("Disease","Logratio")
logratio_bydisease_sigmoid=as.data.frame(logratio_bydisease_sigmoid)

s9=factor(logratio_bydisease_sigmoid$Disease)
split.s9 <- split.data.frame(logratio_bydisease_sigmoid, s9)
dim(split.s9$'CD')
dim(split.s9$'Non-IBD')
logratio_bydisease_sigmoid_CD=split.s9$'CD'
logratio_bydisease_sigmoid_nonIBD=split.s9$'Non-IBD'

a=as.numeric(data_corrected_cecum_NEGandPOS[,'P_391.2869_8.12'])
b=as.numeric(data_corrected_cecum_NEGandPOS[,'P_393.2998_8.13'])
logratio_cecum= a - b
length(which(rownames(data_corrected_cecum_NEGandPOS)==rownames(meta_cecum_POS_Y)))
logratio_bydisease_cecum=cbind(meta_cecum_POS_Y$Disease,logratio_cecum)
rownames(logratio_bydisease_cecum)<-rownames(meta_cecum_POS_Y)
colnames(logratio_bydisease_cecum)<-c("Disease","Logratio")
logratio_bydisease_cecum=as.data.frame(logratio_bydisease_cecum)

s10=factor(logratio_bydisease_cecum$Disease)
split.s10 <- split.data.frame(logratio_bydisease_cecum, s10)
dim(split.s10$'CD')
dim(split.s10$'Non-IBD')
logratio_bydisease_cecum_CD=split.s10$'CD'
logratio_bydisease_cecum_nonIBD=split.s10$'Non-IBD'

wilcox.test(as.numeric(logratio_bydisease_sigmoid_CD$Logratio), as.numeric(logratio_bydisease_sigmoid_nonIBD$Logratio), paired=FALSE)
wilcox.test(as.numeric(logratio_bydisease_cecum_CD$Logratio), as.numeric(logratio_bydisease_cecum_nonIBD$Logratio), paired=FALSE)

par(mfrow = c(1,2))
logratio_bydisease_sigmoid$Disease_ordered<-factor(logratio_bydisease_sigmoid$Disease,levels=c("Non-IBD","CD"))
logratio_bydisease_cecum$Disease_ordered<-factor(logratio_bydisease_cecum$Disease,levels=c("Non-IBD","CD"))
p1<-ggplot(data=logratio_bydisease_sigmoid,aes(x=Disease_ordered,y=as.numeric(Logratio),fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
p2<-ggplot(data=logratio_bydisease_cecum,aes(x=Disease_ordered,y=as.numeric(Logratio),fill=Disease_ordered))+scale_fill_viridis_d(option="D")+geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(shape=21,size=2,position=position_jitter(width=0.25),color="black",alpha=1)+theme_pubr()
grid.arrange(p1, p2, nrow = 1)


## Random forests classifier CD vs. Non-IBD - Sigmoid

sigmoid_index = createDataPartition(meta_sigmoid_POS_Y$Disease, p = 0.6, list = FALSE)
newmeta_sigmoid_train = meta_sigmoid_POS_Y[sigmoid_index, ]
newmeta_sigmoid_test = meta_sigmoid_POS_Y[-sigmoid_index, ]
dim(newmeta_sigmoid_train)
dim(newmeta_sigmoid_test)
write.csv(newmeta_sigmoid_train, "RF_sigmoid_train_meta.csv")
write.csv(newmeta_sigmoid_test, "RF_sigmoid_test_meta.csv")

FDR025_sigmoid = Limma_sigmoid_NEGandPOS[(Limma_sigmoid_NEGandPOS$Qvalue < 0.25), ]  
list=rownames(FDR025_sigmoid)
inputdata=data_corrected_sigmoid_NEGandPOS[,list]
tnewdata_sigmoid_train = inputdata[sigmoid_index, ]
tnewdata_sigmoid_test = inputdata[-sigmoid_index, ]
dim(tnewdata_sigmoid_train)
dim(tnewdata_sigmoid_test)
length(which(rownames(tnewdata_sigmoid_train)==rownames(newmeta_sigmoid_train)))
newmeta_sigmoid_train$Disease=factor(newmeta_sigmoid_train$Disease,levels=c("Non-IBD","CD"))
Disease=as.character(newmeta_sigmoid_train$Disease)
Disease[Disease == "Non-IBD"]="NonIBD"
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_train,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_sigmoid_train,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$CD >= 2,])
inputdata=data_corrected_sigmoid_NEGandPOS[,newlist]
tnewdata_sigmoid_train2 = inputdata[sigmoid_index, ]
tnewdata_sigmoid_test2 = inputdata[-sigmoid_index, ]
dim(tnewdata_sigmoid_train2)
dim(tnewdata_sigmoid_test2)
length(which(rownames(tnewdata_sigmoid_train2)==rownames(newmeta_sigmoid_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_sigmoid_train2,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_sigmoid=train(tnewdata_sigmoid_train2,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_sigmoid
selectedIndices <- rfFit_sigmoid$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_sigmoid$pred$obs[selectedIndices],rfFit_sigmoid$pred$CD[selectedIndices])  #Indicate category to predict
rfImportance_sigmoid=varImp(rfFit_sigmoid,scale=FALSE)
write.csv(rfImportance_sigmoid$importance, "RF_1001trees_Importance_Scores_Sigmoid_CD_vs_non-IBD.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_sigmoid, newdata = tnewdata_sigmoid_test2,type="raw")
probability_predictions=predict(rfFit_sigmoid, newdata = tnewdata_sigmoid_test2,type="prob")
newmeta_sigmoid_test$Disease[newmeta_sigmoid_test$Disease == "Non-IBD"]="NonIBD"
newmeta_sigmoid_test$Disease<-factor(newmeta_sigmoid_test$Disease,levels=c("NonIBD","CD"))
confusionMatrix(category_predictions,newmeta_sigmoid_test[,'Disease'],positive="CD")
plot.roc(newmeta_sigmoid_test$Disease,probability_predictions[[2]]) 
ROC_sigmoid <- roc(newmeta_sigmoid_test$Disease,probability_predictions[[2]])  
auc(ROC_sigmoid)
ci.auc(ROC_sigmoid, method="bootstrap", conf.level=0.95)	
CI_plot1 <- plot.roc(newmeta_sigmoid_test$Disease,probability_predictions[[2]])  
plot(CI_plot1)
ROC1_CI <- ci.sp(CI_plot1, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)


## Random forests classifier CD vs. Non-IBD - Cecum

cecum_index = createDataPartition(meta_cecum_POS_Y$Disease, p = 0.6, list = FALSE)
newmeta_cecum_train = meta_cecum_POS_Y[cecum_index, ]
newmeta_cecum_test = meta_cecum_POS_Y[-cecum_index, ]
dim(newmeta_cecum_train)
dim(newmeta_cecum_test)
write.csv(newmeta_cecum_train, "RF_cecum_train_meta.csv")
write.csv(newmeta_cecum_test, "RF_cecum_test_meta.csv")

FDR025_cecum = Limma_cecum_NEGandPOS[(Limma_cecum_NEGandPOS$Qvalue < 0.25), ]  
list=rownames(FDR025_cecum)
inputdata=data_corrected_cecum_NEGandPOS[,list]
tnewdata_cecum_train = inputdata[cecum_index, ]
tnewdata_cecum_test = inputdata[-cecum_index, ]
dim(tnewdata_cecum_train)
dim(tnewdata_cecum_test)
length(which(rownames(tnewdata_cecum_train)==rownames(newmeta_cecum_train)))
newmeta_cecum_train$Disease=factor(newmeta_cecum_train$Disease,levels=c("Non-IBD","CD"))
Disease=as.character(newmeta_cecum_train$Disease)
Disease[Disease == "Non-IBD"]="NonIBD"
fitControl<-trainControl(method="cv",number=5,classProbs=T,summaryFunction=twoClassSummary,savePredictions=T)
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_cecum_train,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit=train(tnewdata_cecum_train,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfImportance=varImp(rfFit,scale=FALSE)
newlist=rownames(rfImportance$importance[rfImportance$importance$CD >= 2,])
inputdata=data_corrected_cecum_NEGandPOS[,newlist]
tnewdata_cecum_train2 = inputdata[cecum_index, ]
tnewdata_cecum_test2 = inputdata[-cecum_index, ]
dim(tnewdata_cecum_train2)
dim(tnewdata_cecum_test2)
length(which(rownames(tnewdata_cecum_train2)==rownames(newmeta_cecum_train)))
mtry=c(2,4,6,8,10,12,14,16,18,20)
mtryGrid=as.data.frame(mtry)
RFtune=train(tnewdata_cecum_train2,Disease,method="rf",trControl=fitControl,ntree=251,tuneGrid=mtryGrid)  
RFtune
mtry=c(2,4,6)   #refine around best mtry from initial set
mtrySelect=as.data.frame(mtry)
rfFit_cecum=train(tnewdata_cecum_train2,Disease,method="rf",trControl=fitControl,ntree=1001,tuneGrid=mtrySelect,importance=TRUE)
rfFit_cecum
selectedIndices <- rfFit_cecum$pred$mtry == 2  # select mtry parameter with best AUC
plot.roc(rfFit_cecum$pred$obs[selectedIndices],rfFit_cecum$pred$CD[selectedIndices])  #Indicate category to predict
rfImportance_cecum=varImp(rfFit_cecum,scale=FALSE)
write.csv(rfImportance_cecum$importance, "RF_1001trees_Importance_Scores_Cecum_CD_vs_non-IBD.csv")
# Evaluate test characteristics using test dataset
category_predictions=predict(rfFit_cecum, newdata = tnewdata_cecum_test2,type="raw")
probability_predictions=predict(rfFit_cecum, newdata = tnewdata_cecum_test2,type="prob")
newmeta_cecum_test$Disease[newmeta_cecum_test$Disease == "Non-IBD"]="NonIBD"
newmeta_cecum_test$Disease<-factor(newmeta_cecum_test$Disease,levels=c("NonIBD","CD"))
confusionMatrix(category_predictions,newmeta_cecum_test[,'Disease'],positive="CD")
plot.roc(newmeta_cecum_test$Disease,probability_predictions[[2]]) 
ROC_cecum <- roc(newmeta_cecum_test$Disease,probability_predictions[[2]])  
auc(ROC_cecum)
ci.auc(ROC_cecum, method="bootstrap", conf.level=0.95)vg txbb	
plot(CI_plot1)
plot(ROC1_CI, type="shape",col="blue",no.roc=FALSE)
CI_plot2 <- plot.roc(newmeta_cecum_test$Disease,probability_predictions[[2]], add=TRUE)
plot(CI_plot2, add=TRUE)
ROC2_CI <- ci.sp(CI_plot2, sensitivities=seq(0, 1, .01), boot.n=1000)  
plot(ROC2_CI, type="shape",col="red",no.roc=FALSE)
roc.test(ROC_sigmoid,ROC_cecum,method="bootstrap")  # p-value = 0.7455


## Correlations between fecal microbes and metabolites from multivariate models

# CECUM
include_mmvec_cecum=intersect(rownames(data_corrected_cecum_NEGandPOS),colnames(filtereddata_cecum2))
export=t(data_corrected_cecum_NEGandPOS)
cecum_metabolites_for_interomic=export[,include_mmvec_cecum]
dim(cecum_metabolites_for_interomic)
cecum_mapping_for_interomic=meta_cecum_POS_Y[include_mmvec_cecum,]
dim(cecum_mapping_for_interomic)
length(which(rownames(cecum_mapping_for_interomic)==colnames(cecum_metabolites_for_interomic)))
cecum_microbes_for_interomic=filtereddata_cecum2[,include_mmvec_cecum]
dim(cecum_microbes_for_interomic)
length(which(colnames(cecum_metabolites_for_interomic)==colnames(cecum_microbes_for_interomic)))
Disease=factor(cecum_mapping_for_interomic$Disease,levels=c("Non-IBD","CD"))
Obesity=cecum_mapping_for_interomic$Obesity
Gender=cecum_mapping_for_interomic$Gender
write.csv(cecum_metabolites_for_interomic, "temp.csv")
inputdata<-read.csv("temp.csv",header=T,row.names=1) 
design <- model.matrix(~ 0 + Gender + Obesity + Disease)
fit <- lmFit(inputdata, design)
fit2 <- eBayes(fit,robust="TRUE")
fit2
metres_corrected_cecum_NEGandPOS=residuals.MArrayLM(fit2,inputdata)

OTU_cecum_for_interomic=otu_table(cecum_microbes_for_interomic, taxa_are_rows = TRUE)
biom_cecum_for_interomic = phyloseq(OTU_cecum_for_interomic,TAX)
map_cecum_for_interomic=sample_data(cecum_mapping_for_interomic)
data_cecum_for_interomic=merge_phyloseq(biom_cecum_for_interomic,map_cecum_for_interomic)
print(data_cecum_for_interomic)
diagdds=phyloseq_to_deseq2(data_cecum_for_interomic, ~ Batch + Gender + Obesity + Disease)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
dds = DESeq(diagdds, fitType="local")
fitted.common.scale = t( t( assays(dds)[["mu"]] ) / sizeFactors(dds) )
micres_cecum=counts(dds, normalized=TRUE) - fitted.common.scale 
dim(micres_cecum)
micres_cecum_differential_taxa=micres_cecum[rownames(DESEQ2_data_cecum),]
dim(micres_cecum_differential_taxa)
length(which(rownames(micres_cecum_differential_taxa)==rownames(DESEQ2_data_cecum)))
write.csv(DESEQ2_data_cecum$Genus, "temp.csv")   # enter numbers to prevent duplicates, e.g. s1, s2, etc.
cecum_microbes_unique_names<-read.csv("List_of_cecum_microbes_for_bile_acid_heatmap.csv",header=F) 
rownames(micres_cecum_differential_taxa)=cecum_microbes_unique_names$V1

write.csv(DESEQ2_data_cecum$log2FoldChange, "temp.csv")   # change to Enriched or Depleted
cecum_microbe_CD_differential_annotation<-read.csv("temp.csv",header=F) 
rownames(cecum_microbe_CD_differential_annotation)=rownames(micres_cecum_differential_taxa)
colnames(cecum_microbe_CD_differential_annotation)[1]='Differential'

bile_acids=c('P_391.2869_8.12','P_393.2998_8.13')
dim(metres_corrected_cecum_NEGandPOS)
metres=metres_corrected_cecum_NEGandPOS[bile_acids,]
dim(metres)
# Note: need samples to be in rows that match; use transpose if samples are in columns
met<-t(metres)
mic<-t(micres_cecum_differential_taxa)
mic_met_shared_samples=intersect(rownames(mic),rownames(met))
final_mic=mic[mic_met_shared_samples,]
dim(final_mic)
final_met=met[mic_met_shared_samples,]
dim(final_met)
length(which(rownames(final_mic)==rownames(final_met)))

enableWGCNAThreads()
moduleTraitCor = cor(final_mic, final_met, use = "p", method="spearman"); 
nSamples = nrow(final_mic)  
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue[is.na(moduleTraitPvalue)]<-1
moduleTraitQvalueset = qvalue(moduleTraitPvalue, pi0.method="bootstrap", robust=FALSE)
moduleTraitQvalue = moduleTraitQvalueset$qvalues
#write.csv(moduleTraitCor, "Metabolites and metabolic pathways correlation.csv")
#write.csv(moduleTraitPvaluena, "Metabolites and metabolic pathways P-values.csv")
#write.csv(moduleTraitQvalue, "Metabolites and metabolic pathways Q-values.csv")
dim(moduleTraitPvalue)
dim(moduleTraitQvalue)
# BHFDR=p.adjust(moduleTraitPvalue,method="BH")
correlations.upperTriangle<-moduleTraitQvalue #take a copy of the original cor-mat
correlations_melted<-melt(correlations.upperTriangle, value.name ="Qvalue") 
corr.upperTriangle<-moduleTraitCor #take a copy of the original cor-mat
corr_melted<-melt(corr.upperTriangle, value.name ="Correlation")
combine<-cbind(correlations_melted,corr_melted[,3])
colnames(combine)<-c("OTUs", "Metabolites", "Qvalue", "Correlation")
combinesig_cecum = combine[(combine$Qvalue < 0.05), ]   
write.csv(combinesig_cecum, "Bile acid and Microbes correlations of residuals q less than 0.05 - Cecum.csv")
correlations_cecum=moduleTraitCor

pheatmap(correlations_cecum, show_rownames = F, show_colnames = T, breaks = seq(-0.3, 0.3, length = 101), annotation_row=cecum_microbe_CD_differential_annotation, annotation_colors = list(Differential=c(`Enriched` = "red", `Depleted` = "blue")))
pheatmap(correlations_cecum, show_rownames = T, show_colnames = T, breaks = seq(-0.3, 0.3, length = 101), annotation_row=cecum_microbe_CD_differential_annotation, annotation_colors = list(Differential=c(`Enriched` = "red", `Depleted` = "blue")))

# SIGMOID
include_mmvec_sigmoid=intersect(rownames(data_corrected_sigmoid_NEGandPOS),colnames(filtereddata_sigmoid2))
export=t(data_corrected_sigmoid_NEGandPOS)
sigmoid_metabolites_for_interomic=export[,include_mmvec_sigmoid]
dim(sigmoid_metabolites_for_interomic)
sigmoid_mapping_for_interomic=meta_sigmoid_POS_Y[include_mmvec_sigmoid,]
dim(sigmoid_mapping_for_interomic)
length(which(rownames(sigmoid_mapping_for_interomic)==colnames(sigmoid_metabolites_for_interomic)))
sigmoid_microbes_for_interomic=filtereddata_sigmoid2[,include_mmvec_sigmoid]
dim(sigmoid_microbes_for_interomic)
length(which(colnames(sigmoid_metabolites_for_interomic)==colnames(sigmoid_microbes_for_interomic)))
Disease=factor(sigmoid_mapping_for_interomic$Disease,levels=c("Non-IBD","CD"))
Obesity=sigmoid_mapping_for_interomic$Obesity
Gender=sigmoid_mapping_for_interomic$Gender
write.csv(sigmoid_metabolites_for_interomic, "temp.csv")
inputdata<-read.csv("temp.csv",header=T,row.names=1) 
design <- model.matrix(~ 0 + Gender + Obesity + Disease)
fit <- lmFit(inputdata, design)
fit2 <- eBayes(fit,robust="TRUE")
fit2
metres_corrected_sigmoid_NEGandPOS=residuals.MArrayLM(fit2,inputdata)

OTU_sigmoid_for_interomic=otu_table(sigmoid_microbes_for_interomic, taxa_are_rows = TRUE)
biom_sigmoid_for_interomic = phyloseq(OTU_sigmoid_for_interomic,TAX)
map_sigmoid_for_interomic=sample_data(sigmoid_mapping_for_interomic)
data_sigmoid_for_interomic=merge_phyloseq(biom_sigmoid_for_interomic,map_sigmoid_for_interomic)
print(data_sigmoid_for_interomic)
diagdds=phyloseq_to_deseq2(data_sigmoid_for_interomic, ~ Batch + Gender + Obesity + Disease)   # Note: DESeq2 models did not converge with age included
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
dds = DESeq(diagdds, fitType="local")
fitted.common.scale = t( t( assays(dds)[["mu"]] ) / sizeFactors(dds) )
micres_sigmoid=counts(dds, normalized=TRUE) - fitted.common.scale 
dim(micres_sigmoid)
micres_sigmoid_differential_taxa=micres_sigmoid[rownames(DESEQ2_data_sigmoid),]
dim(micres_sigmoid_differential_taxa)
length(which(rownames(micres_sigmoid_differential_taxa)==rownames(DESEQ2_data_sigmoid)))
write.csv(DESEQ2_data_sigmoid$Genus, "temp.csv")   # enter numbers to prevent duplicates, e.g. s1, s2, etc.
sigmoid_microbes_unique_names<-read.csv("List_of_sigmoid_microbes_for_bile_acid_heatmap.csv",header=F) 
rownames(micres_sigmoid_differential_taxa)=sigmoid_microbes_unique_names$V1

write.csv(DESEQ2_data_sigmoid$log2FoldChange, "temp.csv")   # change to Enriched or Depleted
sigmoid_microbe_CD_differential_annotation<-read.csv("temp.csv",header=F) 
rownames(sigmoid_microbe_CD_differential_annotation)=rownames(micres_sigmoid_differential_taxa)
colnames(sigmoid_microbe_CD_differential_annotation)[1]='Differential'

bile_acids=c('P_391.2869_8.12','P_393.2998_8.13')
dim(metres_corrected_sigmoid_NEGandPOS)
metres=metres_corrected_sigmoid_NEGandPOS[bile_acids,]
dim(metres)
# Note: need samples to be in rows that match; use transpose if samples are in columns
met<-t(metres)
mic<-t(micres_sigmoid_differential_taxa)
mic_met_shared_samples=intersect(rownames(mic),rownames(met))
final_mic=mic[mic_met_shared_samples,]
dim(final_mic)
final_met=met[mic_met_shared_samples,]
dim(final_met)
length(which(rownames(final_mic)==rownames(final_met)))

enableWGCNAThreads()
moduleTraitCor = cor(final_mic, final_met, use = "p", method="spearman"); 
nSamples = nrow(final_mic)  
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue[is.na(moduleTraitPvalue)]<-1
moduleTraitQvalueset = qvalue(moduleTraitPvalue, pi0.method="bootstrap", robust=FALSE)
moduleTraitQvalue = moduleTraitQvalueset$qvalues
#write.csv(moduleTraitCor, "Metabolites and metabolic pathways correlation.csv")
#write.csv(moduleTraitPvaluena, "Metabolites and metabolic pathways P-values.csv")
#write.csv(moduleTraitQvalue, "Metabolites and metabolic pathways Q-values.csv")
dim(moduleTraitPvalue)
dim(moduleTraitQvalue)
# BHFDR=p.adjust(moduleTraitPvalue,method="BH")
correlations.upperTriangle<-moduleTraitQvalue #take a copy of the original cor-mat
correlations_melted<-melt(correlations.upperTriangle, value.name ="Qvalue") 
corr.upperTriangle<-moduleTraitCor #take a copy of the original cor-mat
corr_melted<-melt(corr.upperTriangle, value.name ="Correlation")
combine<-cbind(correlations_melted,corr_melted[,3])
colnames(combine)<-c("OTUs", "Metabolites", "Qvalue", "Correlation")
combinesig_sigmoid = combine[(combine$Qvalue < 0.05), ]   
write.csv(combinesig_sigmoid, "Bile acid and Microbes correlations of residuals q less than 0.05 - Sigmoid.csv")
correlations_sigmoid=moduleTraitCor

pheatmap(correlations_sigmoid, show_rownames = F, show_colnames = T, breaks = seq(-0.3, 0.3, length = 101), annotation_row=sigmoid_microbe_CD_differential_annotation, annotation_colors = list(Differential=c(`Enriched` = "red", `Depleted` = "blue")))
pheatmap(correlations_sigmoid, show_rownames = T, show_colnames = T, breaks = seq(-0.3, 0.3, length = 101), annotation_row=sigmoid_microbe_CD_differential_annotation, annotation_colors = list(Differential=c(`Enriched` = "red", `Depleted` = "blue")))


