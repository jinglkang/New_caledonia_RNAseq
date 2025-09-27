setwd("~/Documents/2025/New_caledonia/Reads_nb_matrix")
# DEGs nb
data<-read.table("DEGs_num_barplot.txt",header = TRUE,sep = "\t")
head(data)
ggplot(data, aes(x=Species, y=Num,fill=Tissue)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.6) + 
  theme_bw() + xlab("") + theme_classic()+
  theme(axis.text.x=element_text(colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.ticks.length=unit(.2, "cm"), # set tick length
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=15),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=15), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(), legend.position = c(0.15,0.88),  #不显示网格线
        panel.grid.minor = element_blank())+ylab("DEGs number") +  
  scale_fill_manual(name = "Tissue", values = c("brain" = "#925E9F", "gill"="#FDAF91"))

# Venn
dat <- read.table("DEGs_all_spe.txt",header=TRUE,sep="\t",fill=TRUE, na.strings = "")
library(VennDiagram)
Daru_brain <- dat$Daru_brain[!is.na(dat$Daru_brain)]
Daru_gill <- dat$Daru_gill[!is.na(dat$Daru_gill)]
Zlep_brain <- dat$Zlep_brain[!is.na(dat$Zlep_brain)]
Zlep_gill <- dat$Zlep_gill[!is.na(dat$Zlep_gill)]
venn.plot<-venn.diagram(list(Daru_brain=Daru_brain, Daru_gill=Daru_gill, 
                             Zlep_brain=Zlep_brain, Zlep_gill=Zlep_gill), 
                        fill=c("#925E9F","#FDAF91","#925E9F","#FDAF91"),filename = NULL, cex=2, cat.cex=2)
pdf(file="venn_DEGs.pdf")
grid.draw(venn.plot)
dev.off()
