# load the packages
library(clipr)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(openxlsx)
library(ggpubr)

# set work directory
setwd("C:/Users/Lenovo/Desktop/code/")

# load the taxa data
taxa <- read.csv("taxa.csv")[,c(1,9:14)]3

# rename the data
names(taxa) <- c("asv_id","Kingdom","Phylum","Class","Order","Family","Genus")

# load the EAF data which using 1,000 bootstrap resampling
data.pop <- read.csv("eaf_pro_1000.csv")[,-1]

# merge the taxa and EAF data to combined one data sheet
data.pop <- merge(data.pop, taxa, by="asv_id", all.x = TRUE, sort = F)

# filter the moss_control and moss_warming, repectively.
Moss_control <- filter(data.pop, treatment=="Control" )
Moss_Warming <- filter(data.pop, treatment=="Warming")

# To group microbes by Phylum and rank them by their abundance indicator (X50.) within each group, 
# making it easier to visualize or compare taxonomic composition and relative importance in subsequent analyses.
Moss_control=Moss_control[order(Moss_control$Phylum, Moss_control$X50., decreasing = T),]
Moss_Warming=Moss_Warming[order(Moss_Warming$Phylum, Moss_Warming$X50., decreasing = T),]

# To provide each OTU (row) with a unique rank/order identifier.
Moss_control$ID <- seq(1,length(Moss_control$Phylum),1)
Moss_Warming$ID <- seq(1,length(Moss_Warming$Phylum),1)

# To generate a ranked OTU plot for moss_control and moss_warming
# OTUs ordered by rank ID .
# Estimated 15N enrichment values with confidence intervals.
# Taxonomic grouping visualized by color-coded phyla.

# picture for moss_control
p_moss_control <- 
  ggplot(Moss_control, aes(x=X50.,y=ID)) +
  geom_errorbarh(aes(xmin=X2.5., xmax=X97.5.,colour=Phylum), alpha=0.3, height=0, position=position_dodge(0.5),size=1)+
  geom_point(position=position_dodge(0.5), alpha=0.5, size=2, shape=21, aes(colour=Phylum,fill = Phylum))+
  labs(x=expression(paste("Excess atom fraction" ~ of^15*N)),y="Ranked OTUs") +
  scale_x_continuous(limits = c(-0.3, 0.9),breaks = c(-0.3,0,0.3,0.6,0.9))+
  scale_color_manual(values=c("Actinobacteria"="#BFB1D0",
                              "Aquificae"="#C1C0BF",
                              "Bacteroidetes"="#DCD7C1",
                              "Cyanobacteria"="#6C91C2",
                              "Firmicutes"="#A7C0DE",
                              "Nitrospirae"="#1EB28D",
                              "Planctomycetes"="#EFD01D",
                              "Proteobacteria"="#CF433E",
                              "Tenericutes"="#0808B7",
                              "Spirochaetes"="#EF941D",
                              "unknown Bacteria"="#DBD8D7",
                              "Verrucomicrobia"="#FFEBAD"))+
  scale_fill_manual(values =c("Actinobacteria"="#BFB1D0",
                              "Aquificae"="#C1C0BF",
                              "Bacteroidetes"="#DCD7C1",
                              "Cyanobacteria"="#6C91C2",
                              "Firmicutes"="#A7C0DE",
                              "Nitrospirae"="#1EB28D",
                              "Planctomycetes"="#EFD01D",
                              "Proteobacteria"="#CF433E",
                              "Tenericutes"="#0808B7",
                              "Spirochaetes"="#EF941D",
                              "unknown Bacteria"="#DBD8D7",
                              "Verrucomicrobia"="#FFEBAD"))+
  geom_vline(aes(xintercept=0),linetype="dashed")+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",size=1),
        axis.line=element_line(colour="black",size=1))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(axis.text.x= element_text(size=12,color="black",vjust=0.5, hjust=0.5),
        element_line(colour="black",size=1),
        axis.title = element_text(size=12,color="black"))

# picture for moss_warming
p_moss_warming <- 
  ggplot(Moss_Warming, aes(x=X50.,y=ID)) +
  geom_errorbarh(aes(xmin=X2.5., xmax=X97.5.,colour=Phylum), alpha=0.3, height=0, position=position_dodge(0.5),size=1)+
  geom_point(position=position_dodge(0.5), alpha=0.5,size=2, shape=21, aes(colour=Phylum,fill = Phylum))+
  labs(x=expression(paste("Excess atom fraction" ~ of^15*N)),y="Ranked OTUs") +
  scale_x_continuous(limits = c(-0.3, 0.9),breaks = c(-0.3,0,0.3,0.6,0.9))+
  scale_color_manual(values=c("Actinobacteria"="#BFB1D0",
                              "Aquificae"="#C1C0BF",
                              "Bacteroidetes"="#DCD7C1",
                              "Cyanobacteria"="#6C91C2",
                              "Firmicutes"="#A7C0DE",
                              "Nitrospirae"="#1EB28D",
                              "Planctomycetes"="#EFD01D",
                              "Proteobacteria"="#CF433E",
                              "Tenericutes"="#0808B7",
                              "Spirochaetes"="#EF941D",
                              "unknown Bacteria"="#DBD8D7",
                              "Verrucomicrobia"="#FFEBAD"))+
  scale_fill_manual(values =c("Actinobacteria"="#BFB1D0",
                              "Aquificae"="#C1C0BF",
                              "Bacteroidetes"="#DCD7C1",
                              "Cyanobacteria"="#6C91C2",
                              "Firmicutes"="#A7C0DE",
                              "Nitrospirae"="#1EB28D",
                              "Planctomycetes"="#EFD01D",
                              "Proteobacteria"="#CF433E",
                              "Tenericutes"="#0808B7",
                              "Spirochaetes"="#EF941D",
                              "unknown Bacteria"="#DBD8D7",
                              "Verrucomicrobia"="#FFEBAD"))+
  geom_vline(aes(xintercept=0),linetype="dashed")+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",size=1),
        axis.line=element_line(colour="black",size=1))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(axis.text.x= element_text(size=12,color="black",vjust=0.5, hjust=0.5),
        element_line(colour="black",size=1),
        axis.title = element_text(size=12,color="black"))

# combine multiple ggplot objects into one figure.
ggarrange(p_moss_control,p_moss_warming,nrow = 1,ncol = 2,align = "hv")
# ggsave("EAF_All.pdf", dpi = 600,units = "cm", width=28,height=20)

