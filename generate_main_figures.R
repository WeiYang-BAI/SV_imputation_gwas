setwd('/home/bwy/R_figure')

library(circlize)
library(ComplexHeatmap)
library(data.table)
library(ggblanket)
library(ggmosaic)
library(ggplot2)
library(ggpubr) 
library(ggthemr)
library(gridExtra)
library(ImageGP)
library(patchwork)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))
library(cowplot)
library(see)


man_theme <- ggthemr('fresh')

theme_set(
  man_theme$theme +
    theme(legend.title = element_text(size = 8))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.text.x = element_text(size = 8))+
    theme(axis.title.y = element_text(size = 8))+
    theme(axis.title.x = element_text(size = 8))+
    theme(axis.line = element_line(colour = "grey70"))+
    theme(axis.ticks = element_line(color = 'grey70', size = 0.3))+
    theme(axis.ticks.length = unit(0.8, "mm"))+
    theme(panel.grid = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.major.y = element_blank())+
    theme(strip.text.x = element_text(size = 8, colour = "Black"))+
    theme(legend.text = element_text(size = 8))+
    theme(legend.key.size = unit(0.8, "lines"))+
    theme(legend.title = element_blank())
)



# SV N each haplotype
data=fread('J1.JOINT_make_stage2-3_var/plot/VAR_count_264_assemblies_stage2_SV_del_ins.tsv')

data=subset(data,data$SOURCE != "1KCP")
data$SAMPLE=factor(data$SAMPLE, c('NA18906_h1','HG02717_h1','NA19240_h1','HG02723_h2','HG03486_h2','HG01891_h2','NA19240_h2','HG03579_h2','HG03453_h1','HG02622_h1','HG03516_h1','NA18906_h2','HG03098_h1','HG03098_h2','HG02622_h2','HG02055_h1','HG03579_h1','HG03453_h2','HG03486_h1','HG02723_h1','HG03065_h2','HG02886_h1','HG02257_h1','HG02486_h1','HG02717_h2','HG03516_h2','HG02630_h1','HG02559_h1','HG03540_h1','HG02630_h2','HG02572_h2','HG02559_h2','NA19238_h2','HG01891_h1','HG02886_h2','HG02109_h1','HG02011_h1','NA19239_h2','HG02055_h2','NA20129_h1','HG02145_h2','HG02572_h1','HG03540_h2','NA19239_h1','HG02145_h1','HG02587_h2','HG02818_h2','HG03125_h1','HG02109_h2','HG03371_h1','HG03065_h1','HG03125_h2','NA19983_h1','HG03371_h2','NA21309_h1','NA19238_h1','HG02257_h2','HG02818_h1','HG02587_h1','HG02486_h2','HG02011_h2','NA21309_h2','NA19983_h2','NA20129_h2','HG01109_h2','HG01243_h1','HG01243_h2','HG01928_h2','HG02148_h1','HG01952_h1','HG01928_h1','HG01952_h2','HG01978_h1','HG01109_h1','HG01361_h2','HG00741_h2','HG01123_h2','HG01978_h2','HG01358_h2','HG01106_h1','HG02148_h2','HG00735_h1','HG01106_h2','HG00735_h2','HG01123_h1','HG01361_h1','HG00733_h2','HG00733_h1','HG01358_h1','HG01175_h2','HG00741_h1','HG01071_h1','HG00732_h2','HG01258_h2','HG01175_h1','HG00732_h1','HG01258_h1','HG01071_h2','NA19650_h1','HG00731_h1','HG00731_h2','HG01114_h1','HG01114_h2','NA19650_h2','HG00438_h2','HG005_h1','HG005_h2','HG00621_h2','HG00621_h1','HG00438_h1','HG00673_h2','HG00673_h1','HG02080_h1','HG02080_h2','HG00513_h2','HG00514_h1','HG00864_h1','NA18534_h1','HG00513_h1','NA18939_h1','NA18534_h2','HG00512_h1','NA18939_h2','HG00514_h2','HG00512_h2','HG00864_h2','HG01596_h1','HG01596_h2','HG002_h1','HG002_h2','NA20509_h2','HG00171_h1','HG00096_h1','NA12329_h2','NA12329_h1','HG01505_h1','NA12878_h2','HG00171_h2','NA20509_h1','HG00096_h2','NA12878_h1','HG01505_h2','HG03683_h1','HG03683_h2','HG03492_h1','HG02492_h2','HG02492_h1','HG03492_h2','HG03732_h2','HG03732_h1','HG03009_h2','HG03009_h1','NA20847_h1','NA20847_h2'))

data_long<-melt(
  data,                       #待转换的数据集名称
  id.vars=c("SAMPLE","GENDER","MODE","POP","SBP","OKG","SOURCE","VTYPE","MAPPER","CALLER","TOTAL"),  #要保留的主字段
  variable.name="svTYPE",         #转换后的分类字段名称（维度）
  value.name="count"#转换后的度量值名称
)

a=ggplot(data_long, aes(x=SAMPLE, y=count/1000, fill=svTYPE)) + # show pop
  geom_bar(stat="identity", position = "stack", width=1, size=0.2) +
  labs(x = "Super populations",y = "Number of SVs (K)") + scale_fill_manual(values = c("#FFCF8F","#7C9FB0"))+ 
  scale_y_continuous(expand = c(0,1),limits = c(0, 19))+
  theme(legend.position = c(0.8,0.9))+
  theme(legend.direction = "horizontal")+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())
a


# SV N in population level
data=fread('09.ANNOTATE/plot/population_level_stat.tsv')

data$vTYPE<-factor(data$vTYPE, c('Insertion', 'Deletion'))

b=ggplot(data, aes(vTYPE, count/1000, fill=length)) + 
  geom_bar(stat="identity", position = "stack", width=0.7) +
  labs(x = "VAR",y = "Number of variants (K)") + 
  scale_fill_manual(values = c("#48cae4","#415a77"))+
  facet_wrap(.~gTYPE)+
  theme(legend.position = c(0.8,0.9))+
  theme(axis.ticks.x = element_blank())
b



# SV length distribution
data = read.table('09.ANNOTATE/plot/HGSVC2_HPRCY1.chrAuto.stage2.SV.s5.cr50.AnnotSV.id-chr-pos1-pos2-len-type.tsv', sep = ' ', header = T)

# length
minor_ticks <- c(seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000), seq(200000, 1000000, by = 100000))
c1=ggplot(data, aes(x = length, color = type, fill=type)) +
  geom_histogram(binwidth=0.002) +
  scale_x_log10(
    breaks = c(50, 100, 1000, 10000, 100000, 1000000, minor_ticks),  # Specify major and minor tick positions
    labels = c(50, 100, 1000, 10000, "100000", 1000000, rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  ) +
  scale_color_manual(values = c("#FFCF8F","#7C9FB0"))+ 
  scale_fill_manual(values = c("#FFCF8F","#7C9FB0"))+ 
  labs(x='SV length (bp)',y='Number of SVs (K)', title = '')+
  theme(legend.position = c(0.8,0.2))
c1


# length boxplot
data = read.table('09.ANNOTATE/plot/HGSVC2_HPRCY1.chrAuto.stage2.SV.s5.cr50.AnnotSV.id-chr-pos1-pos2-len-type.tsv', sep = ' ', header = T)

minor_ticks <- c(seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000), seq(200000, 1000000, by = 100000))
c2=ggplot(data, aes(x = type, y = length, fill=type)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.1)+
  geom_boxplot(fill = "white", color = "black",outlier.size = 0.05, outlier.alpha = 0.05, width=0.08, size=0.2) +
  scale_fill_manual(values = c("#FFCF8F","#7C9FB0"))+
  labs(x = "", y = "Length (bp)")+
  scale_y_log10(
    breaks = c(50, 100, 1000, 10000, 100000, 1000000, minor_ticks),  # Specify major and minor tick positions
    labels = c(50, 100, 1000, 10000, "100000", 1000000, rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    #expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  )+
  theme(legend.position = 'none')
c2

c=ggdraw(c1)+draw_plot(c2,x=0.1,y=0.1,scale = 0.6, width = 1,height = 1);c




## basic anno SV
data=fread('09.ANNOTATE/plot/basic_anno_SV.txt')

data$type<-factor(data$type,c('Intergenic', 'Intra-intron', 'Intra-exon', 'Splicing site', '5UTR', '3UTR', 'Whole transcript'))
data$proportion <- data$num / sum(data$num)

#去掉背景网格
fig1 <- ggplot(data,aes(x=type,y=num/1000,fill =type )) + geom_bar(stat='identity',position=position_dodge(), width = 0.75) +
  labs(x=NULL,y=NULL)+    #使用labs自定义x轴和y轴标签名字
  coord_cartesian(ylim = c(0,0.8)) +  #使用ylim设置下面一半的范围
  geom_text(aes(label = paste0(round(proportion * 100, digits=1), "%")), 
            vjust = -0.5, color = "black", size = 3) +
  scale_fill_manual(values=brewer.pal(7,'Set2')) +
  theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = 'none')
fig1 #下面半部分

fig2 <- ggplot(data,aes(x=type,y=num/1000,fill =type )) + geom_bar(stat='identity',position=position_dodge(), width = 0.75) +
  labs(x=NULL,y="Number of SVs (K)") +   #上面半部分标签不需要添加
  coord_cartesian(ylim = c(2.850,3.350)) +  #使用ylim设置上面一半的范围
  geom_text(aes(label = paste0(round(proportion * 100, digits=1), "%")), 
            vjust = -0.5, color = "black", size = 3) +
  #scale_y_continuous(breaks = c(17,20,1)) + 
  scale_fill_manual(values=brewer.pal(7,'Set2')) +
  theme(legend.position = 'none')+
  theme(axis.line.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())
fig2

fig3 <- ggplot(data,aes(x=type,y=num/1000,fill =type )) + geom_bar(stat='identity',position=position_dodge(), width = 0.75) +
  labs(x=NULL,y=NULL) +   #上面半部分标签不需要添加
  coord_cartesian(ylim = c(120,140)) +  #使用ylim设置上面一半的范围
  geom_text(aes(label = paste0(round(proportion * 100, digits=1), "%")), 
            vjust = -0.5, color = "black", size = 3) +
  #scale_y_continuous(breaks = c(17,20,1)) + 
  scale_fill_manual(values=brewer.pal(7,'Set2')) +
  theme(legend.position = 'none')+
  theme(axis.line.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())
fig3


#利用ggarrange（）可以组合不同的图，align为对齐参数， align = "v"为垂直对齐， align = "h"为水平对齐
d=ggarrange(fig3,fig2,fig1,heights=c(1/5, 1/5, 2/5),ncol = 1, nrow = 3, align = "v")+theme(legend.position = "none")

d


#### VNTR allele_N
data=fread('09.ANNOTATE/plot/VNTR_ale_N.tsv')

e1<-ggplot(data, aes(x = ALE_N)) +
  geom_bar(color='black',size=0.01,fill='#778da9',width = 0.7)+
  labs(x = "Alternative allele number of VNTR", y = "Number of VNTR loci (K)")
e1



### motif length
data=fread('09.ANNOTATE/plot/HGSVC2_HPRCY1.HG38.M24sl5F.VNTR.chrAuto.summary')

minor_ticks <- c(seq(20, 100, by = 10), seq(200, 1000, by = 100))
e2=ggplot(data, aes(x = "", y = REF_SoP)) +
  geom_boxplot(fill = "#778da9", color = "black",outlier.size = 0.05, outlier.alpha = 0.3,) +
  labs(x = "Repeat motif", y = "Length (bp)")+
  scale_y_log10(
    breaks = c(10, 100, 1000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    #expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  )+
  theme(legend.position = 'none')
e2

e=ggdraw(e1)+draw_plot(e2,x=0.1,y=0.1,scale = 0.6, width = 1,height = 1);e



## DEL INS trDEL trINS length
svlen = read.table('09.ANNOTATE/plot/zSV_trSV_len.tsv', sep='\t', header = T)
svlen$type<-factor(svlen$type,c('non-VNTR-INS', 'VNTR-INS','non-VNTR-DEL', 'VNTR-DEL'))
my_comparisons <- list( c("non-VNTR-DEL", "VNTR-DEL"), c("non-VNTR-INS", "VNTR-INS"))

minor_ticks <- c(seq(200, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
f <- ggplot(svlen, aes(x=type, y=length, fill=type)) +
  geom_violin(position=position_dodge(1),width=0.8, size=0.1, adjust=1.5)+
  geom_boxplot(width=0.08, cex=0.2, outlier.size = 0.1, outlier.alpha = 0.1, fill = "white")+
  scale_fill_manual(values=c("#76B6EA","#A2D5F2","#F3C677","#F9E8B2")) +
  scale_y_log10(
    breaks = c(100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    labels = c(100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  labs(x = "VAR", y = "Length (bp)")+
  theme(axis.ticks.x = element_blank())
f



nVD=subset(svlen,svlen$type=='non-VNTR-DEL')
VD=subset(svlen,svlen$type=='VNTR-DEL')
nVI=subset(svlen,svlen$type=='non-VNTR-INS')
VI=subset(svlen,svlen$type=='VNTR-INS')

test_result <- wilcox.test(nVD$length, VD$length)
test_result$p.value

test_result <- wilcox.test(nVI$length, VI$length)
test_result$p.value


layout <- "
AAAAAAAAAABB##CCC
DDDD#EEEEE#FFFFF#
"


p=a + b + c + d + e + f  + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))


ggsave('Figure_R_1.pdf', p, width = 15, height = 8)



data_long=fread('J3.JOINT_impu_performance/plot/HG002/HG002_M23s5F_M4_md_default.txt')
data_long$array<-factor(data_long$array,c('UKBA', 'MEGA', 'OMNI-2.5M', 'Pre-imputed', 'STA','Short-read WGS'))
data_long$metric<-factor(data_long$metric,c('Recall', 'Precision', 'F1 score', 'GCR'))
data_long$bench<-factor(data_long$bench,c('Tier1', 'CMRG'))

#colors1=c("#d9ed92","#99d98c","#52b69a","#168aad","#184e77")
colors1=c("#cee5f2","#accbe1","#7c98b3","#637081","#536b78")
colors2=c("#c9ada7","#735d78","#415a77")

p1=ggplot(data_long, aes(metric, value,group=array,fill=array)) + 
  geom_bar(stat="identity", position = position_dodge(0.8), width=0.7) +
  scale_fill_manual(values = colors1)+
  facet_grid(.~bench, scales = "free")+
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(x = "",y = "", title = "", subtitle = "", caption = "", tag = "")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5, vjust = 1))
p1

## temp
tmp=subset(data_long, data_long$bench=='Tier1')
fwrite(tmp,'tmp.txt')
tmp=fread('tmp.txt')
tmp$array<-factor(tmp$array,c('UKB Axion array (0.8M SNPs)', 'Multi-Ethnic Global array (1.8M SNPs)', 'Infinium Omni 2.5M array (2.4M SNPs)', 'Pre-imputed (3.2M SNPs)', 'Short-read WGS (3.4M SNPs)'))
tmp$metric<-factor(tmp$metric,c('Recall', 'Precision', 'F1 score', 'GT concordance'))
ggplot(tmp, aes(metric, value,group=array,fill=array)) + 
  geom_bar(stat="identity", position = position_dodge(0.8), width=0.7) +
  scale_fill_manual(values=brewer.pal(5,'Set2'))+
  facet_grid(.~bench, scales = "free")+
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(x = "",y = "", title = "", subtitle = "", caption = "", tag = "")+
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1))
## 




data=fread('J3.JOINT_impu_performance/plot/loo/zGROUP_77_s5_2t_impu_vs_call_loo_region5.txt')

data=subset(data,data$array!='STA')
data=subset(data,data$array!='STArsq08')

data$array<-factor(data$array,c('UKBA', 'MEGA', 'OMNI-2.5M', 'Pre-imputed','Short-read WGS'))
data$pop<-factor(data$pop,c('AFR', 'AMR', 'EUR', 'SAS', 'EAS'))
data$group<-factor(data$group,c('HGSVC2_HPRCY1'))
data$bench<-factor(data$bench,c('Easy', 'Whole genome', 'Tandem repeat', 'Low mappability', 'Difficult'))

data=subset(data,data$model=='M23s5F')
data=subset(data,data$threshold=='default')
data=subset(data,data$adj=='N')
data=subset(data,data$group=='HGSVC2_HPRCY1') # change group

colors=c("#444a5b","#cc5856","#78a4a1","#dfaf6a")

data=data.frame(data$sample,data$array,data$bench,data$model,data$group,data$precision,data$recall,data$f1,data$gt_concordance)
names(data) <- c("sample", "array", "bench", "model", "group", "Precision", "Recall", "F1 score", "GCR")

data_long<-melt(
  data,                       #待转换的数据集名称
  id.vars=c('sample','array','bench','model','group'),  #要保留的主字段
  variable.name="metric",         #转换后的分类字段名称（维度）
  value.name="value"#转换后的度量值名称
)

data_long=subset(data_long,data_long$bench!='Difficult')
data_long$metric<-factor(data_long$metric,c('Recall', 'Precision', 'F1 score', 'GCR'))


p2=ggplot(data_long, aes(array, value, color=bench)) + 
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.15, 
               notch = F, width=0.5, size=0.4, fill="white")+
  scale_y_continuous(limits = c(0.4, 1), breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_color_manual(values = colors, breaks = c("Whole genome", "Easy", "Tandem repeat", "Low mappability"))+
  facet_grid( .~metric,scale = "free")+
  labs(x = "", y = "", color='Region',title = "", subtitle = "", caption = "", tag = "")
  #theme(legend.position = c(0.90,0.3))
p2



data=fread('J3.JOINT_impu_performance/plot/estiR2/zGROUP_77_pop5_s5_2t_estiR2.txt')
data=subset(data, data$MafBin!='< 1%')

data$MafBin<-factor(data$MafBin,c("< 1%","1-5%","5-10%","10-20%","20-50%"))
data$array<-factor(data$array,c('UKBA', 'MEGA', 'OMNI-2.5M', 'Pre-imputed', 'STArsq08', 'Short-read WGS'))
data$pop<-factor(data$pop,c('AFR', 'AMR', 'EUR', 'SAS', 'EAS'))

data=subset(data,data$model=='without INDEL')
data=subset(data,data$array!='STArsq08')
data=subset(data,data$group=='HGSVC2_HPRCY1')
colors=c("#bc4749","#f77f00","#457b9d","#a98467","#adb5bd")

p3=ggplot(data, aes(x=MafBin, y=R2_mean,group=pop,color=pop,fill=pop)) + # show pop
  geom_line(size=0.35) +
  geom_point(group="pop", size = 0.5)+
  scale_y_continuous(limits = c(0.2, 0.95), breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid( .~array)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  labs(x = "MAF",y = "Mean imputation Rsq",title = "", subtitle = "", caption = "", tag = "")
  #theme(legend.position = c(0.95,0.3))
p3




data=fread('J8.Real_Rsq_VNTR/plot/GP1_ARR5_M1_CATE6_DY')

data$array<-factor(data$array,c('UKBA', 'MEGA', 'OMNI-2.5M', 'Pre-imputed', 'Short-read WGS'))
data$cate<-factor(data$cate,c('all imputed loci', 'rmImpuSGT', 'rmImpuDBT', 'omit singleton repeats', 'rmUniq_rmImpuSGT', 'rmUniq_rmImpuDBT'))

data=subset(data,data$DY=='0')
data=subset(data,data$cate!='rmUniq_rmImpuDBT')
data=subset(data,data$cate!='rmUniq_rmImpuSGT')
data=subset(data,data$cate!='rmImpuDBT')
data=subset(data,data$cate!='rmImpuSGT')

data=subset(data,data$group=='HGSVC2_HPRCY1')


p4=ggplot(data, aes(x=array, y=RsqDOS,color=cate, fill=cate)) + # show pop
  geom_violin(position=position_dodge(1),width=0.88, size=0.1)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.15, 
               notch = F, width=0.13, size=0.4, fill="white")+
  #geom_point(size=1) +
  labs(x = "",y = "rsq (Imputed vs. assembly-based)") + 
  scale_y_continuous(limits = c(0, 1), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  #scale_color_viridis_d(direction = -1) + 
  #facet_wrap( .~array, scales = "free")+
  #facet_grid( .~array)+
  scale_fill_manual(values=c("#f2cc8f","#6d6a75"))+
  scale_color_manual(values=c("#f2cc8f","#6d6a75"))+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5, vjust = 1))+
  theme(legend.spacing.y = unit(5, 'cm'))
  #theme(legend.position = "top")
p4

layout <- "
A
B
C
D
"
p=p1 + p2 + p3 + p4 + plot_layout(design = layout) #+ plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

ggsave('Figure_R_2.pdf', p, width = 9, height = 16)





#################################### s2 volin plot ###########################################
data=fread('04.GREML/res_plot/s2/GREML_meta_bench_cont.tsv.SS10K.rename')
data=subset(data,data$mTYPE=="marginal")
data$hTYPE<-factor(data$hTYPE,c('SGVs', 'SVs'))

SGV_margin=subset(data,data$hTYPE=="SGVs")
SV_margin=subset(data,data$hTYPE=="SVs")

p1 <- ggplot(data,aes(hTYPE,hsq))+ 
  geom_violin(aes(fill=hTYPE),width=0.5,cex=0.1)+ 
  geom_text(data = labels, aes(label = label1, x = group, y = 0.2, vjust = -4, hjust= -0.2), size = 3) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.2, vjust = -2, hjust= -0.2), size = 3) +
  labs(x = "",y = "Variance explained", subtitle = "Marginal analysis") + 
  scale_fill_manual(values = c("#959595","#bc6c25"))+
  geom_boxplot(width=0.08,cex=0.2, fill="white",outlier.size = 1, outlier.alpha = 0.15)+ 
  theme(legend.position = 'none')
p1


data=fread('04.GREML/res_plot/s2/GREML_meta_bench_cont.tsv.SS10K.rename')
data=subset(data,data$mTYPE=="joint")
data$hTYPE<-factor(data$hTYPE,c('SGVs', 'SVs', 'Sum'))

SGV_joint=subset(data,data$hTYPE=="SGVs")
SV_joint=subset(data,data$hTYPE=="SVs")
SUM_joint=subset(data,data$hTYPE=="Sum")


p2 <- ggplot(data,aes(hTYPE,hsq))+ 
  geom_violin(aes(fill=hTYPE),width=1.2,cex=0.1)+ 
  geom_text(data = labels, aes(label = label1, x = group, y = 0.2, vjust = -4, hjust= -0.2), size = 3) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.2, vjust = -2, hjust= -0.2), size = 3) +
  #geom_text(data = labels, aes(label = label3, x = group, y = 0.2, vjust = -0, hjust= -0.2), size = 3) +
  labs(x = "",y = "Variance explained", subtitle = "Joint analysis") + 
  scale_fill_manual(values = c("#ffba08","#ff5a5f","#0077b6"))+
  geom_boxplot(width=0.08,cex=0.2, fill="white",outlier.size = 1, outlier.alpha = 0.15)+ 
  theme(legend.position = 'none')
p2


layout <- "
AABBB
"
p=p1 + p2 + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

p

ggsave('Figure_R_3.pdf', p, width = 9, height = 4)



# pleiotropy

data=fread('05.GWAS/summ/bench_s5_2t/zBench.EUR.gwas.cojo.condi.lda.fm.anno.tsv.oriSig.pleoN.base')


colors1=c("#FFCF8F","#7C9FB0")
colors2=c("#F4A460","#3B5998")

a1=ggplot(data,aes(x =traits_N, group=type, fill=type, color=type))+
  geom_histogram(aes(y=..count..), alpha=1, binwidth = 1, color="white")+
  labs(x = "Number of traits assoiated with SVs",y = "Number of SVs")+
  scale_fill_manual(values = colors1)+
  scale_color_manual(values = colors2)+
  scale_x_continuous(limits = c(0, 20))+
  theme(legend.position = 'none')
#geom_density(adjust = 1.8, alpha=0) # 密度曲线


a2=ggplot(data,aes(x =traits_N, group=type, fill=type, color=type))+
  geom_histogram(aes(y=..count..), alpha=1, binwidth = 0.35)+
  labs(x = "",y = "")+
  scale_fill_manual(values = colors1)+
  scale_color_manual(values = colors1)+
  scale_x_continuous(limits = c(20, 70), breaks = c(20, 45, 70))+
  theme(legend.position = c(1,0.5))
#geom_density(adjust = 1.8, alpha=0) # 密度曲线
a2

a=ggdraw(a1)+draw_plot(a2,x=0.1,y=0.1,scale = 0.6, width = 1,height = 1);a





# annotSV
colors=c("#FFCF8F","#7C9FB0")
colors_sig=c("#bc4749","#457b9d")
data=fread('05.GWAS/summ/bench_s5_2t/zAnnotSV_for_plot/zBench.Merged.EUR.maf0.01.hwe-6.chrAuto.fastGWA.type3.SV.uniq.annotSV.length',sep = '\t')
data$significant<-factor(data$significant,c('sig', 'nonSig'))

# length
minor_ticks <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
b=ggplot(data, aes(x=SV_type, y=SV_length, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.05, #outlier.shape = NA,
               notch = F, width=0.1, size=0.2, fill="white") +
  labs(x = "",y = "SV length (bp)",title = "", subtitle = "", caption = "", tag = "") +
  #theme(legend.position = 'none') +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, 10000, 100000, rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  scale_fill_manual(values = colors_sig, labels = c("Significant SVs", "Other SVs"))+
  scale_color_manual(values = c("black","black"), labels = c("Significant SVs", "Other SVs"))+ stat_compare_means(aes(group = significant,label = paste0("p =", ..p.format..)))
b

# nearest tss
data=fread('05.GWAS/summ/bench_s5_2t/zAnnotSV_for_plot/zBench.Merged.EUR.maf0.01.hwe-6.chrAuto.fastGWA.type3.SV.uniq.annotSV.tss_dist',sep = '\t')
data$significant<-factor(data$significant,c('sig', 'nonSig'))
minor_ticks <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
c=ggplot(data, aes(x=SV_type, y=Dist_nearest_TSS, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.05,
               notch = F, width=0.1, size=0.2, fill="white") +
  labs(x = "",y = "Distance to nearest TSS (bp)",title = "", subtitle = "", caption = "", tag = "") +
  theme(legend.position = 'none') +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  scale_fill_manual(values = colors_sig)+
  scale_color_manual(values = c("black","black"))+ stat_compare_means(aes(group = significant,label = paste0("p", ..p.format..)))
c

nVD=subset(data,data$significant=='sig' & data$SV_type=='DEL')
VD=subset(data,data$significant=='nonSig' & data$SV_type=='DEL')
nVI=subset(data,data$significant=='sig' & data$SV_type=='INS')
VI=subset(data,data$significant=='nonSig' & data$SV_type=='INS')

test_result <- wilcox.test(nVD$Dist_nearest_TSS, VD$Dist_nearest_TSS)
test_result$p.value

test_result <- wilcox.test(nVI$Dist_nearest_TSS, VI$Dist_nearest_TSS)
test_result$p.value





# nearest ss
data=fread('05.GWAS/summ/bench_s5_2t/zAnnotSV_for_plot/zBench.Merged.EUR.maf0.01.hwe-6.chrAuto.fastGWA.type3.SV.uniq.annotSV.ss_dist',sep = '\t')
data$significant<-factor(data$significant,c('sig', 'nonSig'))
minor_ticks <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
d=ggplot(data, aes(x=SV_type, y=Dist_nearest_SS, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.05,
               notch = F, width=0.1, size=0.2, fill="white") +
  labs(x = "",y = "Distance to nearest splicing site (bp)",title = "", subtitle = "", caption = "", tag = "") +
  theme(legend.position = 'none') +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  scale_fill_manual(values = colors_sig)+
  scale_color_manual(values = c("black","black"))+ stat_compare_means(aes(group = significant,label = paste0("p", ..p.format..)))
d


nVD=subset(data,data$significant=='sig' & data$SV_type=='DEL')
VD=subset(data,data$significant=='nonSig' & data$SV_type=='DEL')
nVI=subset(data,data$significant=='sig' & data$SV_type=='INS')
VI=subset(data,data$significant=='nonSig' & data$SV_type=='INS')

test_result <- wilcox.test(nVD$Dist_nearest_SS, VD$Dist_nearest_SS)
test_result$p.value

test_result <- wilcox.test(nVI$Dist_nearest_SS, VI$Dist_nearest_SS)
test_result$p.value




# constrain
data=fread('05.GWAS/summ/bench_s5_2t/zAnnotSV_for_plot/zBench.Merged.EUR.maf0.01.hwe-6.chrAuto.fastGWA.type3.SV.uniq.annotSV.constrain.sv',sep = '\t')

data_long<-melt(
  data,                       #待转换的数据集名称
  id.vars=c('significant','SV_chrom','SV_start','SV_end','SV_length','SV_type','ID','CytoBand','Gene_name','loc'),  #要保留的主字段
  variable.name="metric",         #转换后的分类字段名称（维度）
  value.name="value"#转换后的度量值名称
)


# LOEUF_bin
data_long_sub=subset(data_long,data_long$metric=='LOEUF_bin' & data_long$loc=='exonic')
data_long_sub$significant<-factor(data_long_sub$significant,c('sig', 'nonSig'))

e=ggplot(data_long_sub, aes(x=SV_type, y=value, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1.5, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.15,
               notch = F, width=0.1, size=0.2, fill="white") +
  labs(x = "",y = "LOEUF bin",title = "", subtitle = "", caption = "", tag = "") +
  scale_y_continuous(limits = c(0, 9.5), breaks = c(0, 2.5, 5.0, 7.5))+
  theme(legend.position = 'none') +
  scale_fill_manual(values = colors_sig)+
  scale_color_manual(values = c("black","black")) + stat_compare_means(aes(group = significant,label = paste0("p = ", ..p.format..)), 
                                                                       method = "wilcox.test",na.rm = T, label.y = 9.5)
e

# DDD_HI_percent
data_long_sub=subset(data_long,data_long$metric=='DDD_HI_percent' & data_long$loc=='exonic')
data_long_sub$significant<-factor(data_long_sub$significant,c('sig', 'nonSig'))

f=ggplot(data_long_sub, aes(x=SV_type, y=value, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1.5, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.15, #show.legend = FALSE,
               notch = F, width=0.1, size=0.2, fill="white") +
  #geom_jitter(size=0.2, position=position_dodge(1))+
  labs(x = "",y = "DDD haploinsufficiency percent",title = "", subtitle = "", caption = "", tag = "") +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25, 50, 75, 100))+
  theme(legend.position = 'none') +
  scale_fill_manual(values = colors_sig, labels = c("Significant SVs", "Other SVs"))+
  scale_color_manual(values = c("black","black"), labels = c("Significant SVs", "Other SVs")) + stat_compare_means(aes(group = significant,label = paste0("p = ", ..p.format..)),
                                                                                                                   method = "wilcox.test",na.rm = T, label.y = 105)

f



annotation_row = data.frame(
  Class = factor(rep(c('AD','ADIP','BLD','BONE','BRN','CA','EDCR','EDTH','EPTH','ESC','ESDV','EYE','GI','HRT','HSC','IPSC','KDNY','LCL','LNG','LVR','MSC','MUSC','MYO','NSPH','OTHR','PANC','PCTA','PNS','REPR','SMTH','SPLN','STRM','THYM','URNY'), 
                     c(1,2,9,5,26,74,8,8,25,4,14,2,16,13,9,15,5,14,2,4,4,10,4,2,5,3,7,2,3,5,2,18,2,3))))
name <- readLines("14.ENRICHMENT/plot/cs_cell_names.txt")
rownames(annotation_row) = name
####

data=fread('14.ENRICHMENT/plot/chromatin_state_ov.SV.logistic.noRsqDss.indep.lr.mini350')
data=subset(data,data$Tissue=='BLD'|data$Tissue=='BONE'|data$Tissue=='REPR'|data$Tissue=='MUSC'|data$Tissue=='BRN'|data$Tissue=='EYE'|data$Tissue=='KDNY'|data$Tissue=='LVR'|data$Tissue=='LNG'|data$Tissue=='PANC'|data$Tissue=='REPR'|data$Tissue=='SPLN')

data$dire=NA
data$dire[data$OR<1]=0.5
data$dire[data$OR>=1]=2


data_OR=data.frame(data$cell,data$state,data$OR)

OR_matrix <- acast(data_OR, data.cell ~ data.state, value.var = "data.OR")

data_P=data.frame(data$cell,data$state,data$P)
P_matrix <- acast(data_P, data.cell ~ data.state, value.var = "data.P")


g1=pheatmap(OR_matrix,
            annotation_row = annotation_row,
            #scale = "row",
            border_color ="white",
            number_color="white",
            fontsize_number=14,
            fontsize_row=9,
            fontsize_col=9,
            cellwidth=15,
            cellheight=10,
            cluster_rows=T,
            treeheight_row=14,
            cluster_cols=F,
            color=c(colorRampPalette(colors = c("#1d3557","#f1faee"))(105),colorRampPalette(colors = c("#f1faee","#e63946"))(300)),
            #display_numbers = p2,
            #show_rownames=T，
            display_numbers = display
)


g11 <- cowplot::as_gtable(g1$gtable)

#ggsave('chromatin_state.pdf', g, width = 10, height = 15)


data=fread('14.ENRICHMENT/plot/histone_modi_ov.SV.logistic.noRsqDss.indep.lr')
data=subset(data,data$histone=='H3K4me1'|data$histone=='H3K4me3'|data$histone=='H3K36me3'|data$histone=='H3K27me3'|data$histone=='H3K9me3'|data$histone=='H3K27ac'|data$histone=='H3K9ac')
data=subset(data,data$cell=='astrocyte'|data$cell=='CD14-positive_monocyte'|data$cell=='DND-41'|data$cell=='endodermal_cell'|data$cell=='endothelial_cell_of_umbilical_vein'|data$cell=='fibroblast_of_dermis'|data$cell=='GM12878'|data$cell=='H1'|data$cell=='H9'|data$cell=='HCT116'|data$cell=='HeLa-S3'|data$cell=='hepatocyte'|data$cell=='HepG2'|data$cell=='HUES64'|data$cell=='IMR-90'|data$cell=='K562'|data$cell=='Karpas-422'|data$cell=='keratinocyte'|data$cell=='mammary_epithelial_cell'|data$cell=='MCF-7'|data$cell=='mesenchymal_stem_cell'|data$cell=='MM.1S'|data$cell=='myotube'|data$cell=='NCI-H929'|data$cell=='neural_progenitor_cell'|data$cell=='OCI-LY1'|data$cell=='OCI-LY3'|data$cell=='PC-3'|data$cell=='PC-9'|data$cell=='SK-N-SH'|data$cell=='skeletal_muscle_myoblast'|data$cell=='smooth_muscle_cell'|data$cell=='trophoblast_cell')


data_OR=data.frame(data$cell,data$histone,data$OR)
OR_matrix <- acast(data_OR, data.cell ~ data.histone, value.var = "data.OR")
missing_count <- apply(OR_matrix, 1, function(x) sum(is.na(x)))
OR_matrix <- OR_matrix[missing_count <= 3,]


data_P=data.frame(data$cell,data$histone,data$P)
P_matrix <- acast(data_P, data.cell ~ data.histone, value.var = "data.P")
missing_count <- apply(P_matrix, 1, function(x) sum(is.na(x)))
P_matrix <- P_matrix[missing_count <= 3,]


g2=pheatmap(OR_matrix,
            scale = "none",
            border_color ="white",
            number_color="white",
            fontsize_number=14,
            fontsize_row=9,
            fontsize_col=9,
            cellwidth=15,
            cellheight=10,
            cluster_rows=T,
            treeheight_row=14,
            cluster_cols=T,
            treeheight_col=14,
            color=c(colorRampPalette(colors = c("#1d3557","#f1faee"))(105),colorRampPalette(colors = c("#f1faee","#e63946"))(240)),
            #color = mycol,
            #display_numbers = p2,
            #show_rownames=T,
            display_numbers = display
)

g22 <- cowplot::as_gtable(g2$gtable)
#ggsave('histone_modi.pdf', g, width = 5, height = 10)


data=fread('14.ENRICHMENT/plot/TAD_boundary_ov.SV.logistic.noRsqDss.indep.lr')
data=data[order(data$OR),]
data$expriment=factor(data$expriment,levels = data$expriment)

g3 <- ggplot(data, aes(x = OR, y = expriment,color=Sig)) + 
  scale_color_manual(values = c("#457b9d","#bc4749"))+
  geom_vline(aes(xintercept = 1), linewidth = .45, linetype = "dashed", color="black") + 
  geom_errorbarh(aes(xmax = `CI_1`, xmin = `CI_2`), size = .6, 
                 height = 0) +
  geom_point(size = 2.5,shape=19) +
  ylab("") + xlab("Odds ratio") +
  #facet_wrap( ~ resolution, scales = "free_x")+
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.8,0.1))

#ggsave('enrichment_TAD_boundary.pdf', g, width = 5, height = 7)



#############################################################  RBP

data=fread('14.ENRICHMENT/plot/RBP_eCLIP_allRBP_ov.allRBP.SV.logistic.lr')
data=data[order(data$OR),]

data$Sig=NA
data=data[order(data$OR),]
data$cell=factor(data$cell,levels = data$cell)

g4 <- ggplot(data, aes(x = OR, y = cell,color=Sig)) + 
  scale_color_manual(values = c("#bc4749","#457b9d"))+
  geom_vline(aes(xintercept = 1), linewidth = .45, linetype = "dashed", color="black") + 
  geom_errorbarh(aes(xmax = `CI_1`, xmin = `CI_2`), size = .6, 
                 height = 0) +
  geom_point(size = 2.5,shape=19) +
  ylab("") + xlab("Odds ratio") +
  #facet_wrap( ~ resolution, scales = "free_x")+
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.8,0.1)) +
  theme(legend.position = 'none')

#ggsave('enrichment_RBP.pdf', g, width = 5, height = 2)
g4



layout <- "
AAAA##########
AAAA##########
AAAA##########
AAAA##########
AAAA##########
AAAA##########
AAAA##########
BBCCFFF#######
BBCCFFF#######
BBCCFFF#######
DDEEFFF#######
DDEEFFFGGG####
DDEEFFFGGG####
"


p=a + c + d + e + f + g3 + g4 + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))


ggsave('Figure_R_4.pdf', p, width = 15, height = 13)




# pleiotropic  D8_HGSVC2_HPRCY1_381799
data=fread('08.eQTL/zCherry_GWAS_exp/D8_HGSVC2_HPRCY1_381799.copy.sid.exp.HPR.Brain_Three_tissues')

data_summary <- data %>%
  group_by(tissue, copy) %>%
  summarise(lower = quantile(exp, 0.25),
            middle = median(exp),
            upper = quantile(exp, 0.75))

p1=ggplot(data, aes(x = copy, y = exp, group=copy)) +
  geom_violin(aes(fill=tissue),width=0.5,cex=0)+
  #geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.15, width=0.1, size=0.5, fill="white",color="white")+
  geom_errorbar(data = data_summary, aes(ymin = lower, ymax = lower, y = middle), size = 0.4, width=0.5, color="white") +
  geom_errorbar(data = data_summary, aes(ymin = upper, ymax = upper, y = middle), size = 0.4, width=0.5, color="white") +
  geom_point(data = data_summary, aes(y=middle), size = 2,shape=19, color="#faf0ca") +
  scale_x_continuous(limits = c(-0.5,2.5), breaks = c(0,1,2)) +
  labs(x = "No. of risk allele", y = "Normalized expression of HPR")+
  facet_grid( .~tissue)+
  #scale_fill_manual(values = brewer.pal(3,"Set1"))+
  scale_fill_manual(values = c('#7FCDBB','#1D91C0','#225EA8'))+
  theme(legend.position = 'none')
p1




# pleiotropic  D8_HGSVC2_HPRCY1_108290
data=fread('08.eQTL/zCherry_GWAS_exp/D8_HGSVC2_HPRCY1_108290.copy.sid.exp.multi-gets')
#data=subset(data,data$tissue=='Liver'|data$tissue=='Small_Intestine_Terminal_Ileum'|data$tissue=='Spleen')
data=subset(data,data$gene=='UGT2B17')
data$tissue=factor(data$tissue,c('Liver','Spleen','Testis','Lung','Whole_Blood','Small_Intestine_Terminal_Ileum','Colon_Transverse','Colon_Sigmoid','Cells_EBV-transformed_lymphocytes'))

data_summary <- data %>%
  group_by(tissue, copy) %>%
  summarise(lower = quantile(exp, 0.25),
            middle = median(exp),
            upper = quantile(exp, 0.75))

p2=ggplot(data, aes(x = copy, y = exp, group=copy)) +
  geom_violin(aes(fill=tissue),width=0.8,cex=0)+
  #geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.15, width=0.1, size=0.5, fill="white",color="white")+
  geom_errorbar(data = data_summary, aes(ymin = lower, ymax = lower, y = middle), size = 0.4, width=0.8, color="white") +
  geom_errorbar(data = data_summary, aes(ymin = upper, ymax = upper, y = middle), size = 0.4, width=0.8, color="white") +
  geom_point(data = data_summary, aes(y=middle), size = 2,shape=19, color="#faf0ca") +
  scale_x_continuous(limits = c(-0.5,2.5), breaks = c(0,1,2)) +
  labs(x = "No. of risk allele", y = "Normalized expression of UGT2B17")+
  facet_grid( .~tissue)+
  #scale_fill_manual(values = c(colorRampPalette(colors = c("#708d81","#bf0603"))(9)))+
  scale_fill_manual(values = c(brewer.pal(8,"Dark2"), "lightblue"))+
  #scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#333533"))+
  theme(legend.position = 'none')
p2


layout <- "
AAABBBBBBBBB
"
p=p1 + p2 + plot_layout(design = layout)

p

ggsave('Figure_R_5.pdf', p, width = 12, height = 3)







data = read.table('06.VNTR/fine_plot/candidate/three_sig_pltE100/formatted.HGSVC2_HPRCY1_HG38_chr22_24605024_24606149_7.tsv', sep='\t', header = T)
data=subset(data,data$trait=='Gamma Glutamyltransferase')
(p1 <- ggplot(data, aes(x=genoCN, y=phenoAVG, color=as.factor(trait)))+
    geom_segment(aes(x = genoCN, y = frq/0.2 - 0.6, xend = genoCN, yend = -0.6), size = 4.2, colour =  "#A9B4C0")+
    geom_line(linetype=1,cex=0.7)+
    geom_point(aes(x=genoCN, y=phenoAVG),size=1.5)+
    geom_errorbar(aes(x=genoCN,y = phenoAVG, ymin = ci_up, ymax = ci_down), width = 0.3,linetype=1,cex=0.5)+
    scale_y_continuous(sec.axis = sec_axis(~ (. + 0.6) * 0.2, name = 'Alele frequency'))+
    #scale_x_continuous(breaks = c(0,2,4,6,8))+
    labs(x='GGT1 VNTR length',y='Normalized phenotypes', title = '', col='Trait')+
    scale_color_manual(values = "#b71c1c")+
    man_theme$theme + #主题基本大小
    theme(axis.line = element_line(colour = "grey70"),
          axis.text = element_text(color = 'black'),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.ticks = element_line(color = 'grey70', size = 0.3),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"),
          legend.position = c(0.28,0.98),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          #panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(0.5,0.5,0.5,0.5,'cm'),
          plot.title = element_text(size=13,hjust = 0.5),
          strip.text.x = element_text(size = 8, colour = "Black"))
)

# brewer.pal(8,"Accent")
# "#7FC97F" "#BEAED4" "#FDC086" "#FFFF99" "#386CB0" "#F0027F" "#BF5B17" "#666666"
# brewer.pal(8,"Dark2")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"


data = read.table('06.VNTR/fine_plot/candidate/three_sig/formatted.HGSVC2_HPRCY1_HG38_chr1_247870261_247871078_75.tsv', sep='\t', header = T)
data=subset(data,data$trait!='Red Blood Cell Distribution Width' & data$trait!='Platelet Cout' & data$trait!='Mean Platelet Volume' & data$trait!='Platelet Distribution Width')
(p2 <- ggplot(data, aes(x=genoCN, y=phenoAVG, color=as.factor(trait)))+
    geom_segment(aes(x = genoCN, y = frq/1.1 - 0.07, xend = genoCN, yend = -0.07), size = 4.2, colour =  "#A9B4C0")+
    geom_line(linetype=1,cex=0.7)+
    geom_point(aes(x=genoCN, y=phenoAVG),size=1.5)+
    geom_errorbar(aes(x=genoCN,y = phenoAVG, ymin = ci_up, ymax = ci_down), width = 0.3,linetype=1,cex=0.5)+
    scale_y_continuous(sec.axis = sec_axis(~ (. + 0.07) * 1.1, name = 'Alele frequency'))+
    #scale_x_continuous(breaks = c(0,2,4,6,8))+
    labs(x='TRIM58 VNTR length',y='Normalized phenotypes', title = '', col='Trait')+
    scale_color_manual(values = c("#1B9E77","#BF5B17"))+
    man_theme$theme + #主题基本大小
    theme(axis.line = element_line(colour = "grey70"),
          axis.text = element_text(color = 'black'),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.ticks = element_line(color = 'grey70', size = 0.3),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"),
          legend.position = c(0.35,0.98),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          #panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(0.5,0.5,0.5,0.5,'cm'),
          plot.title = element_text(size=13,hjust = 0.5),
          strip.text.x = element_text(size = 8, colour = "Black"))
)



data = read.table('15.MOTIF/summ/plot/correlation/GGT1_RBMS/Small_Intestine_Terminal_Ileum.GGT1_RBMS', sep='\t', header = T)
(c1=ggplot(data, aes(x=GGT1, y=RBMS1)) +
    geom_abline(slope = -1, colour = "Black", linetype="dashed")+
    geom_point(size=1.5, colour = "grey20",shape = 21,stroke = 0.2, fill="#ddbea8") +
    labs(x = "Normalized expression of GGT1",y = "Normalized expression of RBMS1", title = "Small intestine terminal ileum")+
    annotate("text", x = -1, y = -2.2, label = "r = -0.82, P = 1.9e-35", color = "black", size=2)
)

(c2=ggplot(data, aes(x=GGT1, y=RBMS2)) +
    geom_abline(slope = -1, colour = "Black", linetype="dashed")+
    geom_point(size=1.5, colour = "grey20",shape = 21,stroke = 0.2, fill="#ddbea8") +
    labs(x = "Normalized expression of GGT1",y = "Normalized expression of RBMS2",title = "Small intestine terminal ileum")+
    annotate("text", x = -1, y = -2.2, label = "r = -0.55, P = 2.1e-12", color = "black", size=2)
)

(c3=ggplot(data, aes(x=GGT1, y=RBMS3)) +
    geom_abline(slope = -1, colour = "Black", linetype="dashed")+
    geom_point(size=1.5, colour = "grey20",shape = 21,stroke = 0.2, fill="#ddbea8") +
    labs(x = "Normalized expression of GGT1",y = "Normalized expression of RBMS3",title = "Small intestine terminal ileum")+
    annotate("text", x = -1, y = -2.2, label = "r = -0.72, P = 6.0e-24", color = "black", size=2)
)

layout <- "
A
B
C
"


c=c1+c2+c3+ plot_layout(design = layout)

c


data=fread('08.eQTL/zCherry_GWAS_exp/HGSVC2_HPRCY1_HG38_chr22_24605024_24606149_7.copy.sid.exp.GGT1.TS6')
data$copyINT=round(data$copy, digits = 0)
#data=subset(data,data$copyINT>-17)
data=subset(data,data$tissue=='Thyroid')

data_summary <- data %>%
  group_by(copyINT) %>%
  summarise(lower = quantile(exp, 0.25),
            middle = median(exp),
            upper = quantile(exp, 0.75))

(p3=ggplot(data, aes(x = copyINT, y = exp, group = copyINT)) +
    #geom_point(size=1, colour="grey60")+
    geom_jitter(size=1.3, width=0.5, colour="#d8e2dc")+
    #geom_beeswarm(dodge.width=0.5, size=1, colour="grey")+
    geom_errorbar(data = data_summary, aes(ymin = lower, ymax = upper, y = middle), size = 0.6, width=0.3, color="#3B5998") +
    geom_point(data = data_summary, aes(y=middle), size = 1.5, shape=19, color="#3B5998") +
    #scale_x_continuous(limits = c(-16,3)) +
    labs(x = "GGT1 VNTR length", y = "Normalized expression of GGT1")+
    annotate("text", x = -21, y = 3, label = "Thyroid", color = "black", size=3)
  #facet_grid( tissue~.)+
  #scale_color_manual(values = brewer.pal(3,"Accent"))
)


# OR2T8

data=fread('08.eQTL/zCherry_GWAS_exp/HGSVC2_HPRCY1_HG38_chr1_247870261_247871078_75.copy.sid.exp.OR2T8.Thyroid')
data$copyINT=round(data$copy, digits = 0)

data_summary <- data %>%
  group_by(copyINT) %>%
  summarise(lower = quantile(exp, 0.25),
            middle = median(exp),
            upper = quantile(exp, 0.75))

(p4=ggplot(data, aes(x = copyINT, y = exp, group = copyINT)) +
    #geom_point(aes(x = copy, y = exp), size=1, colour="grey60")+
    geom_jitter(size=1.3, width=0.5, colour="#d8e2dc")+
    #geom_beeswarm(dodge.width=0.5, size=1, colour="grey")+
    geom_errorbar(data = data_summary, aes(ymin = lower, ymax = upper, y = middle), size = 0.6, width=0.3, color="#3B5998") +
    geom_point(data = data_summary, aes(y=middle), size = 1.5, shape=19, color="#3B5998") +
    labs(x = "TRIM58 VNTR length", y = "Normalized expression of OR2T8") +
    annotate("text", x = 12, y = 3, label = "Thyroid", color = "black", size=3)
  #facet_wrap( .~tissue)+
  #scale_color_manual(values = brewer.pal(3,"Accent"))+
)


layout <- "
AAAAAACCCDDDDDD
AAAAAACCCDDDDDD
AAAAAACCCDDDDDD
BBBBBBCCCEEEEEE
BBBBBBCCCEEEEEE
"


p=p1+p3+c+p2+p4 + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

ggsave('Figure_R_6.pdf', p, width = 12, height = 7)



