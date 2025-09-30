library(ggsci)
library(gridExtra)
library(stringr)
library(data.table)
library(ggplot2)
library(ggVennDiagram)
library(ggpubr)
library(ggrepel)

# Significantly functional enrichment plot
setwd("~/Documents/2025/New_caledonia/")
data<-read.csv("Sif_enrichFunc_2.csv",header = TRUE)
names(data)
pd <- position_dodge(0.6)
ggplot(data,aes(Species,-log10(FDR),group=Type,colour=Type,size=Nr_Test))+ 
  theme_bw() + #背景变为白色
  geom_point(position=pd)+
  theme(axis.text.x=element_text(colour="black",family="Times",size=18), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=18,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=14), #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=14), #设置图例的总标题的字体属性
        legend.position = c(0.8,0.6),
        # panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  #geom_hline(yintercept=2, linetype="dashed")+
  geom_hline(yintercept=1.3, linetype="dashed")+
  ylab("-log(FDR)")+xlab("")+ 
  theme(strip.text = element_text(face="plain", family="Times", colour="black", size=14),
        strip.background = element_rect(color = "white"))+
  scale_color_manual(name = "Type", values = c("Circadian rhythm" = "#BB0021", "Energy metabolism"="#EFC000", 
                                               "Immune response"="#008280", "Neural transduction"="#631879", 
                                               "Oxygen response"="#92E9FF", "Vision perception"="#EE4C97", 
                                               "ribosome"="#808180"))
