## * Genes (USMA) vs Chr size
# libraries
library(ggplot2)
library(data.table)
library(scales)
library(dplyr)
library(ggthemes)

rm(list = ls())

scientific <- function(x){
  ifelse(x==0, "0",
         parse(text=gsub("[+]", "",gsub("e", " %*% 10^", scientific_format()(x)))))
} # function for scientific notation

UstilagoGenome <- data.frame("Chr" = c("01", "02", "03", "04", "05", "06", "07", "08",
                                       "09", "10", "11", "12", "13", "14", "15", "16",
                                       "17", "18", "19", "20", "21", "22", "23"),
                             "Size" = c(2476501, 1879391, 1642070, 884984, 1393419, 1031380,
                                        957188, 813247, 733962, 692354, 690620, 650984, 606072,
                                        611467, 575481, 552767, 576628, 560724, 571809, 523884,
                                        470505, 403590, 344927))

UstilagoGenome$Chr <- as.numeric(UstilagoGenome$Chr)

UstilagoGenome$Chr <- paste("USMA_521_v2_", UstilagoGenome$Chr, sep = "")

df.gene <- fread("~/Dropbox (LaboratorioInternaci)/ShraredWithJorge/LabNotebook/USMA_Info/USMA_521_v2_Annotation.csv")

df.gene <- df.gene %>% group_by(Chromosome) %>% summarise(NumGenes = n()) %>% setDT()

df.gene <- df.gene[Chromosome != "USMA_521_v2_25"]
df.gene <- df.gene[Chromosome != "USMA_521_v2_26"]
df.gene <- df.gene[Chromosome != "USMA_521_v2_27"]

df.gene <- df.gene %>% left_join(select(UstilagoGenome, Chr, Size), 
                                 by = c("Chromosome" = "Chr"))
df.gene$GeneDensity <- df.gene$NumGenes / df.gene$Size

plot.1 <- ggplot(df.gene, aes(x = (Size/1000000), y = NumGenes, color = Size)) + 
  geom_point(aes(size = 2), shape = 18) +
  #scale_x_continuous(labels = scientific)+
  theme_classic()+
  labs(x = "Chromosome size  (Mb)", y = "Number of genes") +
  theme(axis.title = element_text(face = "bold", size = 15, color = "black"),
        axis.text = element_text(size = 13),
        legend.position = "none")
plot.1
ggsave(filename = "~/Dropbox (LaboratorioInternaci)/ShraredWithJorge/LabNotebook/2022/2022_Week_38/USMA_GeneNumber.png", plot = plot.1, width = 8, height = 6, units = "in", dpi = 300)

# gene density = NumGenes / Size(Mb)

df.gene$GeneDensity <- df.gene$NumGenes / (df.gene$Size/1000000)

df.gene$Chromosome <- gsub("USMA_521_v2_", "", df.gene$Chromosome)
df.gene$Chromosome <- as.numeric(df.gene$Chromosome)


plot.2 <- ggplot(df.gene, aes(x = Chromosome, y = GeneDensity)) + 
  geom_point(aes(size = 2), shape = 18) + 
  scale_y_continuous(limits = c(00,500)) +
  scale_x_continuous(breaks = seq(min(df.gene$Chromosome), max(df.gene$Chromosome), 1)) +
  labs(y = "Gene density (No. Genes/Mb)", x = "Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")

plot.2

ggsave(filename = "~/Dropbox (LaboratorioInternaci)/ShraredWithJorge/LabNotebook/2022/2022_Week_38/USMA_GeneDensity.png", plot = plot.2, width = 8, height = 6, units = "in", dpi = 300)
  


rm(list = ls())
