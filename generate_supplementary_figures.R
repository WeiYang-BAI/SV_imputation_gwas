setwd('/home/bwy/R_figure')

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(plyr)
library(ImageGP)
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
library(devtools)
library(tidyverse)
library(ggthemr)
library(ggmosaic)
library(cowplot)
library(see)
library(ggmosaic)
library(ggplot2)
library(ggpubr) 
library(ggthemr)
library(gridExtra)
library(ImageGP)
library(patchwork)
library(ggpointdensity)

man_theme <- ggthemr('fresh')

theme_set(
  man_theme$theme +
    theme(legend.title = element_text(size = 10, colour="black"))+
    theme(axis.text.y = element_text(size = 10, colour="black"))+
    theme(axis.text.x = element_text(size = 10, colour="black"))+
    theme(axis.title.y = element_text(size = 10, colour="black"))+
    theme(axis.title.x = element_text(size = 10, colour="black"))+
    theme(axis.line = element_line(colour = "black"))+
    theme(axis.ticks = element_line(color = 'black', size = 0.3))+
    #theme(axis.ticks.length = unit(0.8, "mm"))+
    theme(panel.grid = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.major.y = element_blank())+
    theme(strip.text.x = element_text(size = 10, colour = "Black"))+
    theme(legend.text = element_text(size = 10, colour = "Black"))+
    theme(legend.key.size = unit(0.8, "lines"))+
    theme(legend.title = element_blank())
)


# base calling, reci benchmarking, HG002 calling bench, mergeing threshold

data = read.table('J1.JOINT_make_stage2-3_var/plot/VAR_count_264_assemblies_md0908.tsv', sep = '\t', header = T)

data=subset(data,data$SOURCE != "1KCP")

data$VTYPE <- factor(data$VTYPE, levels = c('SV', 'INDEL', 'SNP'))
data$CONFIDENCE <- factor(data$CONFIDENCE, levels = c('low', 'high'))

#### base count
a <- ggplot(data,aes(CALLER,NUM/1000))+ 
  geom_violin(aes(fill=CALLER),width=0.6,cex=0.1)+ 
  facet_grid(VTYPE ~ CONFIDENCE, scales = "free", labeller = labeller(CONFIDENCE = c(low = "Raw set", high = "Filtered set")))+
  labs(x = "",y = "Number of variants (K)") + 
  geom_boxplot(width=0.04,cex=0.2, fill="white",outlier.size = 0.5, outlier.alpha = 0.1,)+
  scale_fill_manual(values = c("#264a5f","#d15034","#f1ae2d"))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1))
#scale_fill_brewer(palette = "Spectral")
#scale_fill_manual(values = brewer.pal(3,"Set2"))
a



# three tools reci benchmark
data<-read.table("J1.JOINT_make_stage2-3_var/plot/reco_bmk_132_samples.tsv", head=T, sep='\t')
data$CALLER <- factor(data$CALLER, levels = c('dipcall vs. SVIM-asm', 'SVIM-asm vs. dipcall', 
                                              'dipcall vs. PAV', 'PAV vs. dipcall', 
                                              'SVIM-asm vs. PAV', 'PAV vs. SVIM-asm'))
Data=subset(data,data$SOURCE != "1KCP")

b = ggplot(Data, aes(x=RECALL, y=PRECISION, shape=CONFIDENCE, color=CALLER))+
  geom_point(size=2,stroke = 0.4)+
  scale_color_manual(values=c('#56B4E9','#0072B2','#FFCC00', '#FF9900', '#66CC99', '#009E73'))+
  scale_shape_manual(values=c(0,2))+
  labs(x = "Recall", y = "Precision")
b



# called HG002 vs GIAB 
data=fread('J1.JOINT_make_stage2-3_var/plot/called_bench_HG002_v0.6.txt')
data$caller<-factor(data$caller,c('dipcall', 'PAV', 'SVIM-asm', 'Merged'))
data$metric<-factor(data$metric,c('Recall', 'Precision', 'F1', 'GT concordance'))

colors2=c("#264a5f","#d15034","#f1ae2d","grey")

c=ggplot(data, aes(metric, value, group=caller, fill=caller)) + 
  geom_bar(stat="identity", position = position_dodge(0.8), width=0.7) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_fill_manual(values = c("#264a5f","#d15034","#f1ae2d","grey"))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.title.x = element_blank())
#scale_fill_brewer(palette = "Spectral")
#scale_fill_manual(values = brewer.pal(4,"Set2"))

c



# merge thres

data=fread('J2.JOINT_merge_make_panel/plot/HGSVC2_HPRCY1.chrAuto.evalu_merge_paras_plot.txt')

#data=subset(data,data$feature=='simlarity')

d=ggplot(data, aes(x=threshold, y=last_diff, group=feature, color=feature)) +
  #geom_col(position = 'dodge')+
  geom_line(size=0.5,linetype = "solid") +
  geom_point(size=0.5) +
  labs(x = "Merging threshold (1 - SV similarity)",y = "step-wise redundancy") + 
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  scale_color_manual(values = c("#457b9d","#e07a5f"))
#scale_color_brewer(palette = "Spectral")
#scale_color_manual(values = brewer.pal(4,"Set2"))

d


layout <- "
AAABBB
CCCDDD
"


p=a + b + c + d + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))


ggsave('SF_base_call.pdf', p, width = 12, height = 10)


######

## VNTR base anno

data=fread('09.ANNOTATE/plot/basic_anno_VNTR.txt')

data$type<-factor(data$type,c('Intergenic', 'Intra-intron', 'Intra-exon', 'Splicing site', '5UTR', '3UTR', 'Whole transcript'))
data$proportion <- data$num / sum(data$num)

#去掉背景网格
fig1 <- ggplot(data,aes(x=type,y=num,fill =type )) + geom_bar(stat='identity',position=position_dodge(), width = 0.75) +
  labs(x=NULL,y=NULL)+    #使用labs自定义x轴和y轴标签名字
  coord_cartesian(ylim = c(0,250)) +  #使用ylim设置下面一半的范围
  geom_text(aes(label = paste0(round(proportion * 100, digits=2), "%")), 
            vjust = -0.5, color = "black", size = 3) +
  scale_fill_manual(values=brewer.pal(7,'Set2')) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5))+
  theme(axis.ticks.x = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())
fig1 #下面半部分


fig2 <- ggplot(data,aes(x=type,y=num,fill =type )) + geom_bar(stat='identity',position=position_dodge(), width = 0.75) +
  labs(x=NULL,y="Number of VNTR") +   #上面半部分标签不需要添加
  coord_cartesian(ylim = c(6500,8500)) +  #使用ylim设置上面一半的范围
  geom_text(aes(label = paste0(round(proportion * 100, digits=2), "%")), 
            vjust = -0.5, color = "black", size = 3) +
  #scale_y_continuous(breaks = c(6500,7500)) + 
  scale_fill_manual(values=brewer.pal(7,'Set2')) +
  theme_classic()+
  theme(axis.line.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())
fig2

#利用ggarrange（）可以组合不同的图，align为对齐参数， align = "v"为垂直对齐， align = "h"为水平对齐
p=ggarrange(fig2,fig1,heights=c(2/5, 3/5),ncol = 1, nrow = 2, common.legend = T,legend="right",align = "v")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
p
ggsave('SF_basic_anno_VNTR.pdf', p, width = 6, height = 5)


###

# Rsq with or without INDELs

data=fread('J3.JOINT_impu_performance/plot/estiR2/zGROUP_77_pop5_s5_2t_estiR2.txt')

data$MafBin<-factor(data$MafBin,c("< 1%","1-5%","5-10%","10-20%","20-50%"))
data$array<-factor(data$array,c('UKBA', 'MEGA', 'OMNI-2.5M', 'Imputed_A', 'STArsq08', 'Short-read WGS'))
data$pop<-factor(data$pop,c('AFR', 'AMR', 'EUR', 'SAS', 'EAS'))

data=subset(data,data$array=='Short-read WGS')
data=subset(data,data$model=='with INDELs' | data$model=='without INDELs')
data$model<-factor(data$model,c('without INDELs', 'with INDELs'))
#data=subset(data,data$pop=='EUR')
data=subset(data,data$array!='STArsq08')
data=subset(data,data$array!='MEGA')
data=subset(data,data$array!='OMNI-2.5M')
data=subset(data,data$MafBin!='< 1%')
data=subset(data,data$group=='HGSVC2_HPRCY1')

colors=c("#d00000","#ffba08","#0077b6","#bc6c25","#959595")
#colors=c("#e07a5f","#457b9d","#81b29a","#f2cc8f")
man_theme <- ggthemr('pale')
p1=ggplot(data, aes(x=MafBin, y=R2_mean,group=model,color=model)) + # show pop
  geom_line(size=0.65) +
  #geom_point(size=1) +
  labs(x = "MAF",y = "Mean imputation Rsq") + 
  scale_y_continuous(limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
  #scale_color_viridis_d(direction = -1) + 
  #facet_wrap( .~array, scales = "free")+
  facet_grid( .~pop)+
  man_theme$theme + scale_color_manual(values = colors)+
  theme(legend.title = element_text(size = 11))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 9, angle = 45, vjust = 0.5))+
  theme(axis.title.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 10))+
  theme(panel.grid.major = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color = 'black'))+
  theme(legend.title = element_blank())
#theme(legend.position = 'none')
p1
ggsave('SF_R2_M18s5F_M23s5F_WGS_meanRsq.pdf', p1, width = 8, height = 5)


### 
# Full MAF Rsq evaluations

data=fread('J3.JOINT_impu_performance/plot/estiR2/zGROUP_77_pop5_s5_2t_estiR2.txt')
#data=subset(data, data$MafBin!='< 1%')

data$MafBin<-factor(data$MafBin,c("< 1%","1-5%","5-10%","10-20%","20-50%"))
data$array<-factor(data$array,c('UKBA', 'MEGA', 'OMNI-2.5M', 'Pre-imputed', 'STArsq08', 'Short-read WGS'))
data$pop<-factor(data$pop,c('AFR', 'AMR', 'EUR', 'SAS', 'EAS'))

data=subset(data,data$model=='without INDELs')
data=subset(data,data$array!='STArsq08')
data=subset(data,data$group=='HGSVC2_HPRCY1')
colors=c("#bc4749","#f77f00","#457b9d","#a98467","#adb5bd")

p=ggplot(data, aes(x=MafBin, y=R2_mean,group=pop,color=pop,fill=pop)) + # show pop
  geom_line(size=0.35, linetype="dashed") +
  geom_point(group="pop", size = 1)+
  scale_y_continuous(limits = c(0.2, 0.95), breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid( .~array)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  labs(x = "MAF",y = "Mean imputation Rsq",title = "", subtitle = "", caption = "", tag = "")+
  man_theme$theme +
  theme(legend.title = element_text(size = 8))+
  #theme(text=element_text(family="Arial"))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(axis.line = element_line(colour = "grey70"))+
  theme(axis.ticks = element_line(color = 'grey70', size = 0.3))+
  #theme(axis.ticks.length = unit(0.8, "mm"))+
  #theme(panel.grid = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.major.y = element_blank())+
  theme(strip.text.x = element_text(size = 8, colour = "Black"))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.key.size = unit(0.8, "lines"))+
  theme(legend.title = element_blank())

ggsave('SF_full_MAF_meanRsq.pdf', p, width = 9, height = 5)



### 

# Variance explained for AMR AFR EUR

data=fread('04.GREML/res_plot/s5_2t/aframr/P3.comapre.hsq', col.names = c('pop','tt','phe','hTYPE','hsq','se','p'))
data$hsq=round(data$hsq, 6)
data=subset(data,data$pop=="AMR")

data$hTYPE<-factor(data$hTYPE,c('SGVs (marginal)', 'SVs (marginal)', 'SGVs (joint)', 'SVs (joint)', 'Sum (joint)'))

SGV_margin=subset(data,data$hTYPE=="SGVs (marginal)")
SV_margin=subset(data,data$hTYPE=="SVs (marginal)")
SGV_joint=subset(data,data$hTYPE=="SGVs (joint)")
SV_joint=subset(data,data$hTYPE=="SVs (joint)")
SUM_joint=subset(data,data$hTYPE=="Sum (joint)")


a=ggplot(data,aes(hTYPE,hsq))+ 
  geom_violin(aes(fill=hTYPE),width=0.5,cex=0.1)+ 
  #geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_text(data = labels, aes(label = label1, x = group, y = 0.2, vjust = -5, hjust= -0.3), size = 4) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.2, vjust = -3, hjust= -0.3), size = 4) +
  labs(x = "",y = "Variance explained") + 
  #scale_fill_manual(values = c("#959595","#bc6c25"))+
  geom_boxplot(width=0.05,cex=0.5, fill="white",outlier.size = 1, outlier.alpha = 0.15)+
  theme(legend.position = 'none')


data=fread('04.GREML/res_plot/s5_2t/aframr/P3.comapre.hsq', col.names = c('pop','tt','phe','hTYPE','hsq','se','p'))
data$hsq=round(data$hsq, 6)
data=subset(data,data$pop=="AFR")

data$hTYPE<-factor(data$hTYPE,c('SGVs (marginal)', 'SVs (marginal)', 'SGVs (joint)', 'SVs (joint)', 'Sum (joint)'))

SGV_margin=subset(data,data$hTYPE=="SGVs (marginal)")
SV_margin=subset(data,data$hTYPE=="SVs (marginal)")
SGV_joint=subset(data,data$hTYPE=="SGVs (joint)")
SV_joint=subset(data,data$hTYPE=="SVs (joint)")
SUM_joint=subset(data,data$hTYPE=="Sum (joint)")

b=ggplot(data,aes(hTYPE,hsq))+ 
  geom_violin(aes(fill=hTYPE),width=0.5,cex=0.1)+ 
  #geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_text(data = labels, aes(label = label1, x = group, y = 0.2, vjust = -5, hjust= -0.3), size = 4) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.2, vjust = -3, hjust= -0.3), size = 4) +
  labs(x = "",y = "Variance explained") + 
  #scale_fill_manual(values = c("#959595","#bc6c25"))+
  geom_boxplot(width=0.05,cex=0.5, fill="white",outlier.size = 1, outlier.alpha = 0.15)+
  theme(legend.position = 'none')



data=fread('04.GREML/res_plot/s5_2t/aframr/P3.comapre.hsq', col.names = c('pop','tt','phe','hTYPE','hsq','se','p'))
data$hsq=round(data$hsq, 6)
data=subset(data,data$pop=="EUR")

data$hTYPE<-factor(data$hTYPE,c('SGVs (marginal)', 'SVs (marginal)', 'SGVs (joint)', 'SVs (joint)', 'Sum (joint)'))

SGV_margin=subset(data,data$hTYPE=="SGVs (marginal)")
SV_margin=subset(data,data$hTYPE=="SVs (marginal)")
SGV_joint=subset(data,data$hTYPE=="SGVs (joint)")
SV_joint=subset(data,data$hTYPE=="SVs (joint)")
SUM_joint=subset(data,data$hTYPE=="Sum (joint)")

c=ggplot(data,aes(hTYPE,hsq))+ 
  geom_violin(aes(fill=hTYPE),width=0.5,cex=0.1)+ 
  #geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_text(data = labels, aes(label = label1, x = group, y = 0.2, vjust = -5, hjust= -0.3), size = 4) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.2, vjust = -3, hjust= -0.3), size = 4) +
  labs(x = "",y = "Variance explained") + 
  #scale_fill_manual(values = c("#959595","#bc6c25"))+
  geom_boxplot(width=0.02,cex=0.5, fill="white",outlier.size = 1, outlier.alpha = 0.15)+
  theme(legend.position = 'none')


layout <- "
A
B
C
"
p=a + b + c + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size=16))

ggsave('SF_hsq_P3.pdf', p, width = 12, height = 10)


# for main_mg sv

mean_impu_mg_sv=data.frame()

xlist = fread('zInput_file_mg_sv.txt', header = F)

for(x in 1:lengths(xlist)){ 
  file_path = paste0('main_mg/', xlist[x])
  
  data=fread(file_path, 
             col.names = c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
                           "model", "q2", "cv_n", "RC", "cv_type", "iter", "mg_sv"), sep=',')
  
  data$mg_sv <- na.approx(data$mg_sv, rule=2)
  
  data_long<-melt(
    data,                       #待转换的数据集名称
    id.vars=c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
              "model", "q2", "cv_n", "RC", "cv_type", "iter"),  #要保留的主字段
    variable.name="hsq_type",         #转换后的分类字段名称（维度）
    value.name="estimate"#转换后的度量值名称
  )
  
  tmp = ddply(data_long, .(rmMsd, region, chr, lds_pool, model, q2, cv_n, RC, cv_type, hsq_type), function(.d)
    data.frame(mean_estimate=mean(.d$estimate)))
  
  mean_impu_mg_sv=rbind(mean_impu_mg_sv,tmp)

}


# for main_mg sgv

mean_impu_mg_sgv=data.frame()

xlist = fread('zInput_file_mg_sgv.txt', header = F)

for(x in 1:lengths(xlist)){ 
  file_path = paste0('main_mg/', xlist[x])
  
  data=fread(file_path, 
             col.names = c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
                           "model", "q2", "cv_n", "RC", "cv_type", "iter", "mg_sgv"), sep=',')
  
  data$mg_sgv <- na.approx(data$mg_sgv, rule=2)
  
  data_long<-melt(
    data,                       #待转换的数据集名称
    id.vars=c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
              "model", "q2", "cv_n", "RC", "cv_type", "iter"),  #要保留的主字段
    variable.name="hsq_type",         #转换后的分类字段名称（维度）
    value.name="estimate"#转换后的度量值名称
  )
  
  tmp = ddply(data_long, .(rmMsd, region, chr, lds_pool, model, q2, cv_n, RC, cv_type, hsq_type), function(.d)
    data.frame(mean_estimate=mean(.d$estimate)))
  
  mean_impu_mg_sgv=rbind(mean_impu_mg_sgv,tmp)
  
}



# for main_jt

mean_impu_jt=data.frame()

xlist = fread('zInput_file_jt.txt', header = F)

for(x in 1:lengths(xlist)){ 
  file_path = paste0('main_jt/', xlist[x])
  
  data=fread(file_path, 
             col.names = c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
                           "model", "q2", "cv_n", "RC", "cv_type", "iter", 
                           "jt_sgv", "jt_sv", "jt_sum"), sep=',')
  
  data$jt_sgv <- na.approx(data$jt_sgv, rule=2)
  data$jt_sv <- na.approx(data$jt_sv, rule=2)
  data$jt_sum <- na.approx(data$jt_sum, rule=2)
  
  data_long<-melt(
    data,                       #待转换的数据集名称
    id.vars=c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
              "model", "q2", "cv_n", "RC", "cv_type", "iter"),  #要保留的主字段
    variable.name="hsq_type",         #转换后的分类字段名称（维度）
    value.name="estimate"#转换后的度量值名称
  )
  
  tmp = ddply(data_long, .(rmMsd, region, chr, lds_pool, model, q2, cv_n, RC, cv_type, hsq_type), function(.d)
    data.frame(mean_estimate=mean(.d$estimate)))
  
  mean_impu_jt=rbind(mean_impu_jt,tmp)
  
}



## combine all dataframe

mean_impu_all=data.frame()

mean_impu_all=rbind(mean_impu_all,mean_impu_mg_sv)
mean_impu_all=rbind(mean_impu_all,mean_impu_mg_sgv)
mean_impu_all=rbind(mean_impu_all,mean_impu_jt)


fwrite(mean_impu_all, file = "mean_impu_all.csv", sep = ",")

###

# for main_jt diff en m2080

mean_impu_jt_diff_en_m2080=data.frame()

xlist = fread('zInput_file_jt_diff_en_m2080.txt', header = F)

for(x in 1:lengths(xlist)){ 
  file_path = paste0('main_jt/', xlist[x])
  
  data=fread(file_path, 
             col.names = c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
                           "model", "q2", "cv_n", "RC", "cv_type", "iter", 
                           "jt_sgv", "jt_sv", "jt_sum"), sep=',')
  
  data$jt_sgv <- na.approx(data$jt_sgv, rule=2)
  data$jt_sv <- na.approx(data$jt_sv, rule=2)
  data$jt_sum <- na.approx(data$jt_sum, rule=2)
  
  data_long<-melt(
    data,                       #待转换的数据集名称
    id.vars=c("FoS", "gp", "rmMsd", "rep","region", "chr", "lds_pool",
              "model", "q2", "cv_n", "RC", "cv_type", "iter"),  #要保留的主字段
    variable.name="hsq_type",         #转换后的分类字段名称（维度）
    value.name="estimate"#转换后的度量值名称
  )
  
  tmp = ddply(data_long, .(rmMsd, region, chr, lds_pool, model, q2, cv_n, RC, cv_type, hsq_type), function(.d)
    data.frame(mean_estimate=mean(.d$estimate)))
  
  mean_impu_jt_diff_en_m2080=rbind(mean_impu_jt_diff_en_m2080,tmp)
  
}



fwrite(mean_impu_jt_diff_en_m2080, file = "mean_impu_diff_en_m2080.csv", sep = ",")





## reload and plot

data=fread('mean_impu_all_with_diff_en_m2080_add.csv', 
           col.names = c("rmMsd", "Causal regions", "Chromosome", "LDS pool","Model", 
                         "Simulated hsq", "q2", "No. of causal variants", 
                         "Ratio of common and rare hsq", "comm_ratio", "Causal variants", 
                         "sv_ratio", "Components", "Estimates"), sep=',')


data$`Causal regions`<-factor(data$`Causal regions`,c('Whole regions', 'cCRE regions', 'Complex regions', 'Enriched cCRE regions', 'Enriched complex regions'))
data$Chromosome<-factor(data$Chromosome,c('chr2', 'chr10', 'chr20'))
data$Model<-factor(data$Model,c('Marginal analysis', 'Joint analysis', 'LD1M1', 'LD1M6', 'LD2M3', 'LD3M2', 'LD6M1'))
data$`Causal variants`<-factor(data$`Causal variants`,c('SV as causal', 'SV:SGV=1:1', 'SV:SGV=1:3', 'SV:SGV=1:6', 'SV:SGV=1:9', 'SGV as causal'))
data$Components<-factor(data$Components,c('SGVs (marginal)', 'SVs (marginal)', 'SGVs (joint)', 'SVs (joint)', 'Sum (joint)'))


# for main scenario

tmp=subset(data, 
           (data$Model=='Marginal analysis' | data$Model=='Joint analysis')
           & (data$`Chromosome`=='chr2' | data$`Chromosome`=='chr10' | data$`Chromosome`=='chr20')
           & data$`Simulated hsq`=='Total hsq=0.3'
           & data$`No. of causal variants`=='Causal N=60'
           & data$`Causal regions`=='Whole regions'
           & (data$`Causal variants`=='SGV as causal' | data$`Causal variants`=='SV as causal')
           & data$`Ratio of common and rare hsq`=='common:rare=8:2'
)



p=ggplot(tmp, aes(x=`Causal variants`, y=`Estimates`/(q2*comm_ratio), fill=`Components`, na.rm = TRUE)) +
  geom_point(aes(shape=Chromosome,), size=2, colour = "black") +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = c(brewer.pal(5,"Set1")), guide = guide_legend(override.aes = list(color = c(brewer.pal(5,"Set1")))))+
  facet_grid(.~`Model`)+
  labs(x = "",y = "Vg(imp) / Vg",title = "", subtitle = "", caption = "", tag = "")

p

tmp$y=tmp$`Estimates`/(tmp$q2*tmp$comm_ratio)

tmp_mean = ddply(tmp, .(Model, Components, `Causal variants`), function(.d)
  data.frame(mean_y=mean(.d$y)))



ggsave('Figure_R_whole_82_0.3_60.pdf', p, width = 6, height = 3.5)



### desect for the four propotions
col.names = c("rmMsd", 
              "Causal regions", 
              "Chromosome", 
              "LDS pool",
              "Model", 
              "Simulated hsq", 
              "q2", 
              "No. of causal variants", 
              "Ratio of common and rare hsq", 
              "comm_ratio", 
              "Causal variants", 
              "sv_ratio", 
              "Components", 
              "Estimates")

tmp=subset(data, 
           data$Model=='Joint analysis'
           & data$`Causal regions`=='Whole regions'
           & data$`No. of causal variants`=='Causal N=60'
           & data$`Ratio of common and rare hsq`=='common:rare=8:2'
           & (data$`Causal variants`=='SV as causal' | data$`Causal variants`=='SGV as causal')
           & (data$`Components`=='SVs (joint)' | data$`Components`=='SGVs (joint)')
)

tmp$comb=paste0(tmp$`sv_ratio`,':',tmp$`Components`)

z1=ggplot(tmp, aes(x=`Chromosome`, y=`Estimates`/(q2*comm_ratio), 
                fill=`comb`, color=`comb`, na.rm = TRUE)) +
  geom_point(aes(), size=2, position = position_jitterdodge(jitter.width = -1, dodge.width = 0))+
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`Simulated hsq`)+
  labs(x = "",y = "Proportion",title = "", subtitle = "", caption = "", tag = "")
z1


tmp=subset(data, 
           data$Model=='Joint analysis'
           & data$`Causal regions`=='Whole regions'
           & data$`Simulated hsq`=='Total hsq=0.3'
           & data$`Ratio of common and rare hsq`=='common:rare=8:2'
           & (data$`Causal variants`=='SV as causal' | data$`Causal variants`=='SGV as causal')
           & (data$`Components`=='SVs (joint)' | data$`Components`=='SGVs (joint)')
)

tmp$comb=paste0(tmp$`sv_ratio`,':',tmp$`Components`)

z2=ggplot(tmp, aes(x=`Chromosome`, y=`Estimates`/(q2*comm_ratio), 
                   fill=`comb`, color=`comb`, na.rm = TRUE)) +
  geom_point(aes(), size=2, position = position_jitterdodge(jitter.width = -1, dodge.width = 0))+
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`No. of causal variants`)+
  labs(x = "",y = "Proportion",title = "", subtitle = "", caption = "", tag = "")
z2


tmp=subset(data, 
           data$Model=='Joint analysis'
           & (data$`Causal regions`=='Whole regions'|data$`Causal regions`=='cCRE regions'|data$`Causal regions`=='Complex regions')
           #& data$`Chromosome`=='chr2'
           & data$`Simulated hsq`=='Total hsq=0.3'
           & data$`No. of causal variants`=='Causal N=60'
           & data$`Ratio of common and rare hsq`=='common:rare=8:2'
           & (data$`Causal variants`=='SV as causal' | data$`Causal variants`=='SGV as causal')
           & (data$`Components`=='SVs (joint)' | data$`Components`=='SGVs (joint)')
)

tmp$comb=paste0(tmp$`sv_ratio`,':',tmp$`Components`)

z3=ggplot(tmp, aes(x=`Chromosome`, y=`Estimates`/(q2*comm_ratio), 
                   fill=`comb`, color=`comb`, na.rm = TRUE)) +
  geom_point(aes(), size=2, position = position_jitterdodge(jitter.width = -1, dodge.width = 0))+
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`Causal regions`)+
  labs(x = "",y = "Proportion",title = "", subtitle = "", caption = "", tag = "")
z3


tmp=subset(data, 
           data$Model=='Joint analysis'
           & data$`Causal regions`=='Whole regions'
           & data$`Simulated hsq`=='Total hsq=0.3'
           & data$`No. of causal variants`=='Causal N=60'
           & (data$`Ratio of common and rare hsq`=='common:rare=9:1' | data$`Ratio of common and rare hsq`=='common:rare=8:2' | data$`Ratio of common and rare hsq`=='common:rare=7:3')
           & (data$`Causal variants`=='SV as causal' | data$`Causal variants`=='SGV as causal')
           & (data$`Components`=='SVs (joint)' | data$`Components`=='SGVs (joint)')
)

tmp$comb=paste0(tmp$`sv_ratio`,':',tmp$`Components`)

z4=ggplot(tmp, aes(x=`Chromosome`, y=`Estimates`/(q2*comm_ratio), 
                   fill=`comb`, color=`comb`, na.rm = TRUE)) +
  geom_point(aes(), size=2, position = position_jitterdodge(jitter.width = -1, dodge.width = 0))+
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`Ratio of common and rare hsq`)+
  labs(x = "",y = "Proportion",title = "", subtitle = "", caption = "", tag = "")+
  theme(plot.margin = margin(t = -20, r = 0, b = 0, l = 0, unit = "mm"))
z4


layout <- "
A
B
C
D
"
p=z3 + z1 + z2 + z4 + plot_layout(design = layout, heights = 0.1) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size=16))

ggsave('SF_four_pi_tmp.pdf', p, width = 8, height = 12)





## solve for final proportion

Four_prop=subset(data, 
                 data$Model=='Joint analysis'
                 & (data$`Ratio of common and rare hsq`!='Only common causal variants' & data$`Ratio of common and rare hsq`!='Only rare causal variants')
                 & (data$`Causal variants`=='SGV as causal' | data$`Causal variants`=='SV as causal')
                 & (data$`Components`=='SGVs (joint)' | data$`Components`=='SVs (joint)')
)


Four_prop$prop=Four_prop$`Estimates`/(Four_prop$q2*Four_prop$comm_ratio)

Four_prop$Estimates <- NULL
Four_prop$sv_ratio <- NULL

Four_prop_wide <- pivot_wider(Four_prop, 
                              names_from = c(`Causal variants`, Components), 
                              values_from = prop, 
                              names_sep = ":")


solve_linear_system <- function(a, b, c, d, z1, z2) {
  solutions <- vector("list", length(a))    
  for (i in seq_along(a)) {
    coeff_matrix <- matrix(c(a[i], b[i], c[i], d[i]), nrow = 2, byrow = TRUE)
    const_vector <- c(z1[i], z2[i])
    
    if(det(coeff_matrix) != 0) {
      solutions[[i]] <- solve(coeff_matrix, const_vector)
    } else {
      solutions[[i]] <- NA 
    }
  }
  
  do.call(rbind, solutions)
}

solutions <- solve_linear_system(Four_prop_wide$`SV as causal:SVs (joint)`, 
                                 Four_prop_wide$`SGV as causal:SVs (joint)`, 
                                 Four_prop_wide$`SV as causal:SGVs (joint)`,
                                 Four_prop_wide$`SGV as causal:SGVs (joint)`, 
                                 rep(0.011, times = 243),
                                 rep(0.185, times = 243))

solutions_df <- as.data.frame(solutions)

names(solutions_df) <- c("SV", "SGV")

solutions_df$Proportion <- solutions_df$SV / (solutions_df$SV + solutions_df$SGV)


Four_prop_wide$SV_estimates = solutions_df$SV
Four_prop_wide$SGV_estimates = solutions_df$SGV
Four_prop_wide$SV_proportion = solutions_df$Proportion


tmp=subset(Four_prop_wide, 
           Four_prop_wide$`Causal regions`=='Whole regions'
           & Four_prop_wide$`No. of causal variants`=='Causal N=60'
           & Four_prop_wide$`Ratio of common and rare hsq`=='common:rare=8:2'
)


a=ggplot(tmp, aes(x=`Chromosome`, y=`SV_proportion`, 
                           na.rm = TRUE)) +
  geom_point(aes(), size=2)+
  scale_y_continuous(limits = c(-0.1, 0.35), breaks = c(0, 0.1, 0.2)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`Simulated hsq`)+
  labs(x = "",y = "Predicted proportion",title = "", subtitle = "", caption = "", tag = "")
a


tmp=subset(Four_prop_wide, 
           Four_prop_wide$`Causal regions`=='Whole regions'
           & Four_prop_wide$`Simulated hsq`=='Total hsq=0.3'
           & Four_prop_wide$`Ratio of common and rare hsq`=='common:rare=8:2'
)


b=ggplot(tmp, aes(x=`Chromosome`, y=`SV_proportion`, 
                  na.rm = TRUE)) +
  geom_point(aes(), size=2)+
  scale_y_continuous(limits = c(-0.1, 0.35), breaks = c(0, 0.1, 0.2)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`No. of causal variants`)+
  labs(x = "",y = "Predicted proportion",title = "", subtitle = "", caption = "", tag = "")
b



tmp=subset(Four_prop_wide, 
           Four_prop_wide$`Simulated hsq`=='Total hsq=0.3'
           & Four_prop_wide$`No. of causal variants`=='Causal N=60'
           & Four_prop_wide$`Ratio of common and rare hsq`=='common:rare=8:2'
)

c=ggplot(tmp, aes(x=`Chromosome`, y=`SV_proportion`, 
                  na.rm = TRUE)) +
  geom_point(aes(), size=2)+
  scale_y_continuous(limits = c(-0.1, 0.35), breaks = c(0, 0.1, 0.2)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`Causal regions`)+
  labs(x = "",y = "Predicted proportion",title = "", subtitle = "", caption = "", tag = "")
#theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

c



tmp=subset(Four_prop_wide, 
           Four_prop_wide$`Causal regions`=='Whole regions'
           #& data$`Chromosome`=='chr2'
           & Four_prop_wide$`Simulated hsq`=='Total hsq=0.3'
           & Four_prop_wide$`No. of causal variants`=='Causal N=60'
           # Four_prop_wide$`Ratio of common and rare hsq`=='common:rare=8:2'
)

d=ggplot(tmp, aes(x=`Chromosome`, y=`SV_proportion`, 
                  na.rm = TRUE)) +
  geom_point(aes(), size=2)+
  scale_y_continuous(limits = c(-0.1, 0.35), breaks = c(0, 0.1, 0.2)) +
  scale_fill_manual(values = c(brewer.pal(5,"Set1")))+
  scale_color_manual(values = c(brewer.pal(5,"Set1")))+
  facet_grid(.~`Ratio of common and rare hsq`)+
  labs(x = "",y = "Predicted proportion",title = "", subtitle = "", caption = "", tag = "")
#theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

d


layout <- "
A
B
C
D
"
p=c + a + b + d + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size=16))

ggsave('SF_SV_prop.pdf', p, width = 8, height = 12)





## for corr for 1:1 1:3 1:6 1:9

data_cor_sv=subset(data,
                (data$`Causal variants`!='SGV as causal' & data$`Causal variants`!='SV as causal')
                & data$Components=='SVs (joint)'
                #& data$`Simulated hsq`!='Total hsq=0.2' # not included 0.2
                & (data$`Ratio of common and rare hsq`!='Only common causal variants' & data$`Ratio of common and rare hsq`!='Only rare causal variants')
                )
data_cor_sv$comm_SV_q2=data_cor_sv$q2*data_cor_sv$comm_ratio*data_cor_sv$sv_ratio
data_cor_sv$comm_SGV_q2=data_cor_sv$q2*data_cor_sv$comm_ratio*(1-data_cor_sv$sv_ratio)

data_cor_sv$predict = 0.392 * data_cor_sv$comm_SV_q2 + 0.013 * data_cor_sv$comm_SGV_q2



data_cor_sgv=subset(data,
                   (data$`Causal variants`!='SGV as causal' & data$`Causal variants`!='SV as causal')
                   & data$Components=='SGVs (joint)'
                   #& data$`Simulated hsq`!='Total hsq=0.2' # not included 0.2
                   & (data$`Ratio of common and rare hsq`!='Only common causal variants' & data$`Ratio of common and rare hsq`!='Only rare causal variants')
)
data_cor_sgv$comm_SV_q2=data_cor_sgv$q2*data_cor_sgv$comm_ratio*data_cor_sgv$sv_ratio
data_cor_sgv$comm_SGV_q2=data_cor_sgv$q2*data_cor_sgv$comm_ratio*(1-data_cor_sgv$sv_ratio)

data_cor_sgv$predict = 0.652 * data_cor_sgv$comm_SV_q2 + 0.949 * data_cor_sgv$comm_SGV_q2



data_cor=rbind(data_cor_sv,data_cor_sgv)


cor_stat = ddply(data_cor, .(`Causal regions`, `Ratio of common and rare hsq`), function(.d)
  data.frame(correlation=round(cor(.d$Estimates, .d$predict)^2,3)))

cor_stat$correlation=paste0('rsq = ', cor_stat$correlation)

tmp=subset(data_cor, 
           data_cor$`Ratio of common and rare hsq`=='common:rare=8:2'
           #&data_cor$`Causal regions`=='Whole regions'
           )

tmp_cor_stat=subset(cor_stat, 
                    cor_stat$`Ratio of common and rare hsq`=='common:rare=8:2'
           #&data_cor$`Causal regions`=='Whole regions'
)

e=ggplot(tmp) +
  #geom_abline() +
  geom_point(aes(x=tmp$Estimates, y=tmp$predict, color=`Causal variants`, fill=`Causal variants`, shape=`Components`),size=1.8, stroke = 0.3)+
  geom_text(data=tmp_cor_stat, aes(x = Inf, y = Inf, hjust = 2.7, vjust = 3,
                               label = correlation), size=2.8,
  )+
  facet_grid(.~`Causal regions`)+
  labs(x = "GREML estimates",y = "Predicted estimates",title = "", subtitle = "", caption = "", tag = "")+
  theme(axis.text.x = element_text())

e



tmp=subset(data_cor, 
           data_cor$`Causal regions`=='Whole regions'
           & (data_cor$`Ratio of common and rare hsq`!='Only common causal variants' & data_cor$`Ratio of common and rare hsq`!='Only rare causal variants')
) 

tmp_cor_stat=subset(cor_stat, 
                    cor_stat$`Causal regions`=='Whole regions'
                    & (cor_stat$`Ratio of common and rare hsq`!='Only common causal variants' & cor_stat$`Ratio of common and rare hsq`!='Only rare causal variants')
)

f=ggplot(tmp) +
  #geom_abline() +
  geom_point(aes(x=tmp$Estimates, y=tmp$predict, color=`Causal variants`, fill=`Causal variants`, shape=`Components`),size=1.8, stroke = 0.3)+
  geom_text(data=tmp_cor_stat, aes(x = Inf, y = Inf, hjust = 2.7, vjust = 3,
                                   label = correlation), size=2.8,
  )+
  facet_grid(.~`Ratio of common and rare hsq`)+
  labs(x = "GREML estimates",y = "Predicted estimates",title = "", subtitle = "", caption = "", tag = "")+
  theme(axis.text.x = element_text())

f


layout <- "
E
F
"

p= e + f + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size=16, color = "black"))

ggsave('Figure_R_predicted.pdf', p, width = 7, height = 7)






# LDS and LDF


man_theme <- ggthemr('fresh')

theme_set(
  man_theme$theme +
    theme(legend.title = element_text(size = 8))+
    #theme(text=element_text(family="Arial"))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1))+
    theme(axis.title.y = element_text(size = 8))+
    theme(axis.title.x = element_text(size = 8))+
    theme(axis.line = element_line(colour = "black"))+
    #theme(panel.grid = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.major.y = element_blank())+
    theme(strip.text.x = element_text(size = 8, colour = "Black"))+
    theme(legend.text = element_text(size = 8))+
    theme(legend.key.size = unit(0.8, "lines"))+
    theme(legend.title = element_blank())
)


### lds
data = read.table('11.LDA/plot/zSumm_lds.tsv.name', sep = '\t', header = T)
minor_ticks <- c(seq(1, 10, by = 1), seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000))
p1=ggplot(data, mapping = aes(x = type, y = ldscore)) +
  labs(x = "",y = "Cumulative LD score") + 
  geom_violin(aes(fill=type),width=1,cex=0.1,adjust = 2)+ 
  scale_fill_manual(values = c("#e07a5f","#f2cc8f","#81b29a","#457b9d"))+
  geom_boxplot(width=0.1,cex=0.5, fill="white",outlier.size = 1, outlier.alpha = 0.15)+
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, "10000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  ) +
  theme(legend.position = 'none')
p1
ggsave('LD_score.pdf', p, width = 5, height = 4)


#### ldf
data = read.table('11.LDA/plot/zSumm_ldf.tsv.name', sep = '\t', header = T)
minor_ticks <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000))
p2=ggplot(data, mapping = aes(x = type, y = nSNPs)) +
  labs(x = "",y = "Number of LD friends") + 
  geom_violin(aes(fill=type),width=1,cex=0.1,adjust = 2)+ 
  scale_fill_manual(values = c("#e07a5f","#f2cc8f","#81b29a","#457b9d"))+
  geom_boxplot(width=0.1,cex=0.5, fill="white",outlier.size = 1, outlier.alpha = 0.15)+
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, "10000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  ) +
  theme(legend.position = 'none')
p2
ggsave('LD_friend.pdf', p, width = 5, height = 4)

#

layout <- "
AB
"


p=p1 + p2  + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

ggsave('SF_LDA.pdf', p, width = 8, height = 5)



###

## SV length. sig vs non-sig

man_theme <- ggthemr('fresh')

theme_set(
  man_theme$theme +
    theme(legend.title = element_text(size = 12))+
    #theme(text=element_text(family="Arial"))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 12))+
    theme(axis.title.y = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 12))+
    theme(axis.line = element_line(colour = "black"))+
    theme(panel.grid = element_blank()) +
    #theme(panel.grid.major.x = element_blank()) +
    #theme(panel.grid.major.y = element_line(colour = "grey90"))+
    theme(strip.text.x = element_text(size = 12, colour = "Black"))+
    theme(legend.text = element_text(size = 12))+
    theme(legend.key.size = unit(0.8, "lines"))+
    #theme(legend.position = 'none')+
    theme(legend.title = element_blank())
)

colors_sig=c("#bc4749","#457b9d")

data=fread('05.GWAS/summ/bench_s5_2t/zAnnotSV_for_plot/zBench.Merged.EUR.maf0.01.hwe-6.chrAuto.fastGWA.type3.SV.uniq.annotSV.length',sep = '\t')

data$significant<-factor(data$significant,c('sig', 'nonSig'))

# length
minor_ticks <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
p=ggplot(data, aes(x=SV_type, y=SV_length, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.05, #outlier.shape = NA,
               notch = F, width=0.1, size=0.2, fill="white") +
  labs(x = "",y = "SV length (bp)",title = "", subtitle = "", caption = "", tag = "") +
  #theme(legend.position = 'none') +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  scale_fill_manual(values = colors_sig, labels = c("Traits-associated SVs", "Other SVs"))+
  scale_color_manual(values = c("black","black"), labels = c("Traits-associated SVs", "Other SVs"))+ stat_compare_means(aes(group = significant,label = paste0("p =", ..p.format..)))

ggsave('SF_bh_sig_sv_length.pdf', p, width = 7, height = 6)



###

# individual hsq, MAF, BEAT SV vs SGV

man_theme <- ggthemr('fresh')

theme_set(
  man_theme$theme +
    theme(legend.title = element_text(size = 12))+
    #theme(text=element_text(family="Arial"))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 12))+
    theme(axis.title.y = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 12))+
    theme(axis.line = element_line(colour = "black"))+
    theme(panel.grid = element_blank()) +
    #theme(panel.grid.major.x = element_blank()) +
    #theme(panel.grid.major.y = element_line(colour = "grey90"))+
    theme(strip.text.x = element_text(size = 12, colour = "Black"))+
    theme(legend.text = element_text(size = 12))+
    theme(legend.key.size = unit(0.8, "lines"))+
    #theme(legend.position = 'none')+
    theme(legend.title = element_blank())
)


# individual hsq
data=fread('05.GWAS/summ/bench_s5_2t/zBench.cont.clump.af.beta.indi-hsq.shared-phe',sep = '\t',
           col.names = c("af", "beta", "phe", "type", "hsq_rint", "hsq_rint_resi"))

data=subset(data,data$type=='SGV' | data$type=='SV')

data_SGV=subset(data,data$type=='SGV')
data_SV=subset(data,data$type=='SV')

labels <- data.frame(
  group = c('SGV', 'SV'),
  label0 = c(paste0('median: ',round(median(data_SGV$hsq_rint_resi), 5)), 
             paste0('median: ',round(median(data_SV$hsq_rint_resi), 5))),
  label1 = c(paste0('mean: ',round(mean(data_SGV$hsq_rint_resi), 5)), 
             paste0('mean: ',round(mean(data_SV$hsq_rint_resi), 5))),
  label2 = c(paste0('s.e.m: ',round(sd(data_SGV$hsq_rint_resi)/sqrt(length(data_SGV$hsq_rint_resi)), 6)), 
             paste0('s.e.m: ',round(sd(data_SV$hsq_rint_resi)/sqrt(length(data_SV$hsq_rint_resi)), 6)))
)


#minor_ticks <- c(seq(200, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
a=ggviolin(data, x = "type", y = "hsq_rint_resi", fill = "type",width = 0.8)+
  #p <- ggplot(data,aes(significant,SV_length))+ 
  geom_boxplot(width=0.1, cex=0.7, outlier.size = 0.1, outlier.alpha = 0.1, fill = "white")+
  geom_text(data = labels, aes(label = label0, x = group, y = 0.01, vjust = -5, hjust= -0.2), size = 2.5) +
  geom_text(data = labels, aes(label = label1, x = group, y = 0.01, vjust = -3, hjust= -0.2), size = 2.5) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.01, vjust = -1, hjust= -0.2), size = 2.5) +
  scale_fill_manual(values = c("#959595","#bc6c25"))+
  scale_y_log10(
    #  breaks = c(100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    #  labels = c(100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  labs(x = "", y = "log10(variance explained)")+
  theme(legend.position = 'none')

a


# individual af
data=fread('05.GWAS/summ/bench_s5_2t/zBench.cont.clump.af.beta.indi-hsq.shared-phe',sep = '\t',
           col.names = c("af", "beta", "phe", "type", "hsq_rint", "hsq_rint_resi"))

data=subset(data,data$type=='SGV' | data$type=='SV')

data_SGV=subset(data,data$type=='SGV')
data_SV=subset(data,data$type=='SV')

labels <- data.frame(
  group = c('SGV', 'SV'),
  label0 = c(paste0('median: ',round(median(data_SGV$af), 5)), 
             paste0('median: ',round(median(data_SV$af), 5))),
  label1 = c(paste0('mean: ',round(mean(data_SGV$af), 5)), 
             paste0('mean: ',round(mean(data_SV$af), 5)))
)


#minor_ticks <- c(seq(200, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
b=ggviolin(data, x = "type", y = "af", fill = "type", width = 0.5)+
  #p <- ggplot(data,aes(significant,SV_length))+ 
  geom_boxplot(width=0.1, cex=0.7, outlier.size = 0.1, outlier.alpha = 0.1, fill = "white")+
  geom_text(data = labels, aes(label = label0, x = group, y = 0.01, vjust = -13, hjust= -0.4), size = 2.5) +
  geom_text(data = labels, aes(label = label1, x = group, y = 0.01, vjust = -11, hjust= -0.4), size = 2.5) +
  scale_fill_manual(values = c("#959595","#bc6c25"))+
  labs(x = "", y = "Allele frequency")+
  theme(legend.position = 'none')

b

# individual beta
data=fread('05.GWAS/summ/bench_s5_2t/zBench.cont.clump.af.beta.indi-hsq.shared-phe',sep = '\t',
           col.names = c("af", "beta", "phe", "type", "hsq_rint", "hsq_rint_resi"))

data=subset(data,data$type=='SGV' | data$type=='SV')

data_SGV=subset(data,data$type=='SGV')
data_SV=subset(data,data$type=='SV')

labels <- data.frame(
  group = c('SGV', 'SV'),
  label0 = c(paste0('median: ',round(median(data_SGV$beta), 5)), 
             paste0('median: ',round(median(data_SV$beta), 5))),
  label1 = c(paste0('mean: ',round(mean(data_SGV$beta), 5)), 
             paste0('mean: ',round(mean(data_SV$beta), 5))),
  label2 = c(paste0('s.e.m: ',round(sd(data_SGV$beta)/sqrt(length(data_SGV$beta)), 5)), 
             paste0('s.e.m: ',round(sd(data_SV$beta)/sqrt(length(data_SV$beta)), 5)))
)


#minor_ticks <- c(seq(200, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
c=ggviolin(data, x = "type", y = "beta", fill = "type", width = 0.5)+
  #p <- ggplot(data,aes(significant,SV_length))+ 
  geom_boxplot(width=0.1, cex=0.7, outlier.size = 0.1, outlier.alpha = 0.1, fill = "white")+
  geom_text(data = labels, aes(label = label0, x = group, y = 0.01, vjust = -19, hjust= -0.2), size = 2.5) +
  geom_text(data = labels, aes(label = label1, x = group, y = 0.01, vjust = -17, hjust= -0.2), size = 2.5) +
  geom_text(data = labels, aes(label = label2, x = group, y = 0.01, vjust = -15, hjust= -0.2), size = 2.5) +
  scale_fill_manual(values = c("#959595","#bc6c25"))+
  scale_y_log10(
    #  breaks = c(100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    #  labels = c(100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  labs(x = "", y = "log10(Effect size)")+
  theme(legend.position = 'none')



# hsq, af dot
data=fread('05.GWAS/summ/bench_s5_2t/zBench.cont.clump.af.beta.indi-hsq.shared-phe.shuf10000',sep = '\t',
           col.names = c("af", "beta", "phe", "type", "hsq_rint", "hsq_rint_resi"))

data=subset(data,data$type=='SGV' | data$type=='SV')
data$type<-factor(data$type,c('SV', 'SGV'))

d=ggplot(data, aes(x=af, y=beta, fill=type)) +
  geom_point(size=1.3, colour = "black",shape = 21,stroke = 0) +
  labs(x = "Allele frequency",y = "Effect size") + 
  scale_fill_manual(values = c("#bc6c25","#959595"))
#facet_grid( .~model)+
#theme(legend.position = 'none')



layout <- "
AAACCC
DDDDDD
"


p=a + c + d + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))


ggsave('SF_individual_hsq_af_beta.pdf', p, width = 10, height = 8)








####

# VNTR pleo and nearest TSS and SS


man_theme <- ggthemr('fresh')

theme_set(
  man_theme$theme +
    theme(legend.title = element_text(size = 12))+
    #theme(text=element_text(family="Arial"))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 12))+
    theme(axis.title.y = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 12))+
    theme(axis.line = element_line(colour = "black"))+
    theme(panel.grid = element_blank()) +
    #theme(panel.grid.major.x = element_blank()) +
    #theme(panel.grid.major.y = element_line(colour = "grey90"))+
    theme(strip.text.x = element_text(size = 12, colour = "Black"))+
    theme(legend.text = element_text(size = 12))+
    theme(legend.key.size = unit(0.8, "lines"))+
    #theme(legend.position = 'none')+
    theme(legend.title = element_blank())
)


# pleiotropy

data=fread('06.VNTR/summ_s5_2t_bench/zBench.EUR.VNTR.sig.with-condi.annotSV.tsv.pheN.base')


colors1=c("#FFCF8F","#7C9FB0")
colors2=c("#F4A460","#3B5998")

a1=ggplot(data,aes(x =traits_N))+
  geom_histogram(aes(y=..count..), alpha=1, binwidth = 1, color="white")+
  labs(x = "Number of traits associated with VNTRs",y = "Number of VNTRs")+
  scale_fill_manual(values = colors1)+
  scale_color_manual(values = colors2)+
  scale_x_continuous(limits = c(0, 20))+
  theme(legend.position = 'none')
#geom_density(adjust = 1.8, alpha=0) # 密度曲线
a1

a2=ggplot(data,aes(x =traits_N))+
  geom_histogram(aes(y=..count..), alpha=1, binwidth = 0.35)+
  labs(x = "",y = "")+
  scale_fill_manual(values = colors1)+
  scale_color_manual(values = colors1)+
  scale_x_continuous(limits = c(20, 60), breaks = c(20, 40, 60))+
  theme(legend.position = c(1,0.5))
#geom_density(adjust = 1.8, alpha=0) # 密度曲线
a2

a=ggdraw(a1)+draw_plot(a2,x=0.15,y=0.12,scale = 0.6, width = 1,height = 1)


# Dist_nearest_
colors_sig=c("#bc4749","#457b9d")

data=fread('06.VNTR/summ_s5_2t_bench/for_anno_plot/sig.nonsig.vntr.anno.whole.tsv.ss.tss.rmNoneLine',sep = '\t')

data_long<-melt(
  data,                       #待转换的数据集名称
  id.vars=c("significant", "VNTR_ID", "ALE_N", "CHR", "REF_STR", "REF_END"),  #要保留的主字段
  variable.name="s_type",         #转换后的分类字段名称（维度）
  value.name="value"#转换后的度量值名称
)

data_long$significant<-factor(data_long$significant,c('sig', 'nonsig'))
minor_ticks <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000))
b=ggplot(data_long, aes(x=s_type, y=value, color=significant, fill=significant)) +
  geom_violin(position=position_dodge(1),width=0.5, size=0.2, adjust=1, color= NA)+
  geom_boxplot(position=position_dodge(1), outlier.size = 0.1, outlier.alpha = 0.05,
               notch = F, width=0.1, size=0.2, fill="white") +
  labs(x = "",y = "Distance (bp)",title = "", subtitle = "", caption = "", tag = "") +
  #theme(legend.position = 'none') +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 100000, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 1000, 10000, "100000", rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
  ) +
  scale_fill_manual(values = colors_sig, labels = c("Traits-associated VNTRs", "Other VNTRs"))+
  scale_color_manual(values = c("black","black"), labels = c("Traits-associated VNTRs", "Other VNTRs"))+ 
  stat_compare_means(aes(group = significant,label = paste0("p", ..p.format..)))


layout <- "
AB
"


p=a + b  + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

ggsave('SF_bh_VNTR_pleo_TSS_SS.pdf', p, width = 12, height = 5)



### COMP

man_theme <- ggthemr('pale')

## compare counts
data = read.table('16.COMPARE/plot/2t/summ_stat_long.txt', sep='\t', header = T)
data$TRAIT=factor(data$TRAIT, c('Vitamin D','Glucose','Direct Bilirubin','Total Bilirubin','Phosphate','Calcium','Urea','Alanine Aminotransferase','C-reactive Protein','Albumin','Total Protein','Aspartate Aminotransferase','Urate','SHBG','Gamma Glutamyltransferase','Cystatin C','HbA1c','Alkaline Phosphatase','Creatinine','IGF1','Basophil Count','Immature reticulocyte fraction','Neutrophil Count','Lymphocyte Count','High Light Scatter Reticulocyte Count','White blood cell (leukocyte) count ','Reticulocyte Count','Monocyte Count','Mean Sphered Cell Volume','Eosinophil Count','Mean Corpuscular Haemoglobin','Mean Reticulocyte Volume','Red blood cell (erythrocyte) count','Red Blood Cell Distribution Width','Platelet Distribution Width','Mean Corpuscular Volume','Mean Platelet Volume','Platelet Cout','Lipoprotein A','Apolipoprotein B','LDL Cholesterol','Cholesterol','Apolipoprotein A','HDL Cholesterol','Triglycerides','BMI','Height'))
data$type=factor(data$type,c('shared','unique to HGSVC2-HPRC panel','unique to 150K SRS panel'))

a=ggplot(data, aes( x = TRAIT, y=counts, fill=type))+
  labs(x="",y="Number of traits-associated loci", fill=" ",title="")+
  scale_fill_manual(values=c("#264a5f","#d15034","#f1ae2d")) +
  geom_bar(stat = "identity", position = "dodge",  width = 0.6) + 
  man_theme$theme +
  theme(legend.title = element_text(size = 4))+
  theme(axis.text.y = element_text(size = 7))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=6)) +
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 7))+
  theme(legend.title = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = c(0.8,0.85))
a

## compare P
colors3=c("#a47148","#f0c05a","#6c757d","#e2583e","#61a5c2")
data = read.table('16.COMPARE/plot/summ_stat_shared_closest_p.txt', sep = '\t', header = T)

minor_ticks <- c(seq(10, 100, by = 10), seq(200, 330, by = 100))
b=ggplot(data, aes(x=-log10(p_own), y=-log10(p_uoi))) +
  geom_abline() +
  geom_pointdensity(size = 0.5,adjust = 1) +
  #geom_point(size=1.8, colour = "black",shape = 21,stroke = 0.3) +
  labs(x = "-log10(P value) in this study",y = "-log10(P value) in 150K SRS study") + 
  scale_fill_manual(values = colors3)+
  facet_wrap( .~CATEGORY, scales = "free")+
  scale_x_log10(
    breaks = c(10, 100, 330, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 330, rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  ) +
  scale_y_log10(
    breaks = c(10, 100, 330, minor_ticks),  # Specify major and minor tick positions
    labels = c(10, 100, 330, rep("", length(minor_ticks))),  # Specify labels for major and minor ticks
    expand = c(0.03, 0),  # Ensure axis does not extend beyond data range
  ) +
  man_theme$theme + #scale_color_manual(values = colors)+
  theme(legend.title = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.title = element_blank())
#theme(legend.position = 'none')
b


x=-log10(data$p_own)
y=-log10(data$p_uoi)
x[is.infinite(x)] <- 0
y[is.infinite(y)] <- 0
correlation <- cor(x, y)


## compare beta
data = read.table('16.COMPARE/plot/summ_stat_shared.txt', sep = '\t', header = T)
data$beta_uoi[data$Freq_uoi > 0.5] <- -data$beta_uoi[data$Freq_uoi > 0.5]
data$Freq_uoi[data$Freq_uoi > 0.5] <- 1-data$Freq_uoi[data$Freq_uoi > 0.5]
data$fd=abs(data$Freq_own-data$Freq_uoi)


c=ggplot(data, aes(x=beta_own, y=beta_uoi,fill=CATEGORY)) +
  geom_abline() +
  #geom_pointdensity(size = 0.5,adjust = 1) +
  geom_point(size=1.8, colour = "black",shape = 21,stroke = 0.3) +
  labs(x = "Effect size in this study",y = "Effect size in 150K SRS study") + 
  scale_fill_manual(values = colors3)+
  #facet_wrap( .~CATEGORY, scales = "free")+
  man_theme$theme + #scale_color_manual(values = colors)+
  theme(legend.title = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.title = element_blank())
#theme(legend.position = 'none')
c
#ggsave('compare_to_UKB_150K_SV_shared_beta.pdf', p1, width = 7, height = 6)



## uniq to OWN region
data = read.table('16.COMPARE/plot/summ_stat_uniq_own_anno_regions.txt.anno.id.lenRefType_2', sep = '\t', header = T)
data$regions=factor(data$regions,c('complex regions', 'confident regions'))
d=ggplot(data, aes(x = regions,fill=regions)) +
  geom_bar(width = 0.8) +
  #geom_vline(xintercept = 0.01, linetype = "dashed", color = "orange")+
  labs(x='regions',y='Number of SVs', title = '')+
  scale_fill_manual(values = c("#e2583e","#f0c05a","#6c757d","#a47148"))+ 
  man_theme$theme +
  theme(axis.text.x = element_text(size = 8)) +
  theme(legend.title = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(panel.grid = element_blank())+
  theme(legend.position = 'none')
d


layout <- "
AAAAAA
BBBDDD
"
p=a + b + d + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size=10))


ggsave('SF_compare_SV.pdf', p, width = 9, height = 10)


# VNTR
data = read.table('16.COMPARE/plot/2t/summ_stat_VNTR.txt', sep='\t', header = T)
data$TRAIT=factor(data$TRAIT, c('Basophil Count','Glucose','Vitamin D','Forced vital capacity (FVC)','Lipoprotein A','Immature reticulocyte fraction','Phosphate','Diastolic blood pressure, automated reading','Total Bilirubin','Albumin','Direct Bilirubin','LDL Cholesterol','Systolic blood pressure, automated reading','Alanine Aminotransferase','Urea','BMI','Urate','FEV1/ FVC ratio Z-score','Apolipoprotein B','Calcium','C-reactive Protein','Triglycerides','Cholesterol','Total Protein','Aspartate Aminotransferase','SHBG','Lymphocyte Count','Creatinine','Gamma Glutamyltransferase','Apolipoprotein A','Cystatin C','HDL Cholesterol','High Light Scatter Reticulocyte Count','Heel bone mineral density (BMD) T-score, automated','Neutrophil Count','Reticulocyte Count','White blood cell (leukocyte) count ','Mean Sphered Cell Volume','Red blood cell (erythrocyte) count','Monocyte Count','Eosinophil Count','Alkaline Phosphatase','Platelet Cout','Platelet Distribution Width','Mean Corpuscular Haemoglobin','IGF1','Red Blood Cell Distribution Width','Mean Corpuscular Volume','Mean Platelet Volume','HbA1c','Height'))
data$type=factor(data$type,c('shared','unique to HGSVC2-HPRC panel','unique to 4,688 SRS panel'))

a=ggplot(data, aes( x = TRAIT, y=counts, fill=type))+
  labs(x="",y="Number of traits-associated loci", fill=" ",title="")+
  scale_fill_manual(values=c("#264a5f","#d15034","#f1ae2d")) +
  geom_bar(stat = "identity", position = "dodge",  width = 0.6) + 
  man_theme$theme +
  theme(legend.title = element_text(size = 4))+
  theme(axis.text.y = element_text(size = 7))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=6)) +
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 7))+
  theme(legend.title = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = c(0.8,0.85))
a



# shared P
colors3=c("#a47148","#f0c05a","#6c757d","#e2583e","#61a5c2")
data = read.table('16.COMPARE/plot/summ_stat_VNTR_shared_P.txt', sep = '\t', header = T)

minor_ticks <- c(seq(10, 100, by = 10), seq(200, 330, by = 100))
b=ggplot(data, aes(x=-log10(p_own), y=-log10(p_poru))) +
  geom_abline() +
  geom_pointdensity(size = 0.5,adjust = 1) +
  labs(x = "-log10(P value) in this study",y = "-log10(P value) in 4,688 SRS study") + 
  scale_fill_manual(values = colors3)+
  facet_wrap( .~CATEGORY, scales = "free")+
  scale_x_log10(
    breaks = c(10, 100, 330, minor_ticks),  
    labels = c(10, 100, 330, rep("", length(minor_ticks))),  
    expand = c(0.03, 0),  
  ) +
  scale_y_log10(
    breaks = c(10, 100, 330, minor_ticks),  
    labels = c(10, 100, 330, rep("", length(minor_ticks))),  
    expand = c(0.03, 0), 
  ) +
  man_theme$theme + 
  theme(legend.title = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.title = element_blank())
#theme(legend.position = 'none')


layout <- "
AAAAAA
BBBB##
"
p=a + b + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size=10))


ggsave('SF_compare_VNTR.pdf', p, width = 9, height = 10)


