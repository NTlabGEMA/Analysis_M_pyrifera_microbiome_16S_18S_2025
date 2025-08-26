############################
library("reshape2")
library("stringr")
library("ggplot2")
library("RColorBrewer")
library(ggh4x)
library(ggalluvial)
library(microeco)
library(magrittr)

###Edit
# make all needless annotation to  sth like "k__" ... "s__"
mecoC<-clone(meco1)
mecoC
mecoC$tax_table %<>% tidy_taxonomy
mecoC$tax_table
#tmp <- mecoC$tax_table
tmp<-ASV_Tax
# search each row to add sth
for(i in 1:nrow(tmp)){
  if(any(grepl("__$", tmp[i, ]))){
    matchall <- grep("__$", tmp[i, ])
    tmp[i, matchall] <- paste0(tmp[i, matchall], "unclassified ", gsub(".__", "", tmp[i, matchall[1] - 1]))
  }
}
# then reassign tmp to the raw table
# mecoC$tax_table <- tmp
# mecoC$tax_table
# mecoC
T_table <- mecoC1$tax_table
O_table <- mecoC1$otu_table
S_table <- mecoC1$sample_table

S_table <- S_table %>%
  mutate(Site = recode(Site,
                           "Las Docas" = "Las Docas (N)",
                       "Algarrobo" = "Algarrobo (N)",
                       "Navidad" = "Navidad (C)",
                       "Topocalma" = "Topocalma (C)",
                       "Pargua" = "Pargua (SP)",
                       "Los Chonos" = "Los Chonos (SP)",
                       "Ilque" = "Ilque (IS)",
                       "San Antonio" = "San Antonio (IS)"))
# Verificar los nuevos niveles
levels(as.factor(S_table$Site))

# List of all replicates, ordered. The names should match the ones in the OTU table.
replicates_list <- as.character(colnames(O_table))

# List of the corresponding replicates groups to which the individual replicates belong to.
# The order should match the one of 'replicates_list'.
# If some samples were not replicated, simply re-copy here the 'replicates_list' just above.
# The order of this list determines the sample order in the bubble plot.
replicates_groups <- as.character(S_table$Site)
# Do these two list have the same length?
length(replicates_list) == length(replicates_groups)

# Taxonomic level that will be aggregated and displaeyed in the bubble plot.
# It should match a taxonomic level from the taxonomic annotation file.
# If you want to display OTUs/ASVs, simply write 'OTU'.
tax_aggr <- "Genus"

# Number of taxonomic bins to be displayed in the bubble plot.
tax_number <- 60

# Taxonomic level that will be used for colouring bubbles (i.e. categorical variable)
# It should match a taxonomic level from the taxonomic annotation file
# (and it should be a higher taxonomic level that the one from 'tax_aggr').
tax_col <- "Class"
tax_col1 <- "Family"

# Filename for saving the bubble plot (svg or pdf)
file_name <- "bubble_plot.svg"

# Dimension (in inch) of this file (first number is the width, the second number is the height)
plot_dim <- c(6,6)

# names(otu_tab)[1] <- "OTU"
# otu_tab[1:5,1:5]

#ASV_Tax<-cbind(OTU,T_table)
# names(T_table)[1] <- "OTU"
# 
# ASV_Tax[1:5,1:7]
# tax_tab<-ASV_Tax
# tax_tab[1:5,1:7]
# Delete all cells containing 'uncultured' or 'unkown'

for (col in 2:ncol(tax_tab)) {
  for (row in 1:nrow(tax_tab)) {
    if (grepl("uncultured",tax_tab[row,col],ignore.case = TRUE)) {
      tax_tab[row,col] <- ""
    }
    if (grepl("unknown",tax_tab[row,col],ignore.case = TRUE)) {
      tax_tab[row,col] <- ""
    }
  }
}

# Replace empty cells by 'NA'
tax_tab2 <- as.data.frame(apply(tax_tab, 2, function(x) gsub("^$|^ $", NA, x)))

# Remove columns containing only 'NA'
col_to_remove <- c()

for (col in 2:ncol(tax_tab2)) {
  x <- sum(is.na(tax_tab2[,col]))/nrow(tax_tab2)
  if (x == 1) {
    col_to_remove <- c(col_to_remove, col)
  }
}

if (length(col_to_remove) > 0) {
  tax_tab3 <- tax_tab2[,-col_to_remove]
} else {
  tax_tab3 <- tax_tab2
}

# Set taxonomic annotations as character variables
for (col in 2:ncol(tax_tab3)) {
  tax_tab3[,col] <- as.character(tax_tab3[,col])
}

# Fill all NAs

for (col in 2:ncol(tax_tab3)) {
  for (row in 1:nrow(tax_tab3)) {
    if (is.na(tax_tab3[row,col])) {
      if (!grepl("OTU", tax_tab3[row,col-1]) & !grepl("unassigned", tax_tab3[row,col-1])) {
        tax_tab3[row,col] <- paste0("unassigned ", tax_tab3[row,col-1])
      } else {
        tax_tab3[row,col] <- tax_tab3[row,col-1]
      }
    }
  }
}

#Compute the relative abundance of OTUs for each sample

head(otu_tab)
otu_counts <- colSums(otu_tab)
otu_tab2 <- otu_tab
otu_tab2 <- sweep(otu_tab, 2, otu_counts, `/`)
#otu_counts <- colSums(otu_tab2)
#otu_tab2[is.na(otu_tab2)] <- 0

colSums(otu_tab2)
dim(otu_tab2)

#Merge the OTU and taxonomic tables together
dim(tax_tab3)
m <- cbind(otu_tab2, tax_tab3)
head(m)
# Has the merged table the expected dimension?
dim(m)
colSums(m[,1:71])
#Aggregate the table to taxonomic level defined in the variable 'tax_aggr'

# First, we should save in a column the taxonomic information needed for computing the bubble plot
taxonomy1 <- c()
for (row in 1:nrow(m)) {
  taxonomy1 <- c(taxonomy1, paste0(m[row,names(m)==tax_col1], ";", m[row,names(m)==tax_aggr]))
}
dim(taxonomy1)
 taxonomy <- c()
 for (row in 1:nrow(m)) {
   taxonomy <- c(taxonomy, paste0(m[row,names(m)==tax_col], ";", m[row,names(m)==tax_aggr]))
 }
# Subset from the merged table the selected samples only
m2 <- m[,names(m) %in% replicates_list]
dim(m2)

# Aggregate 'm2' based on the selected taxonomic level
m3 <- aggregate(m2, by=list(taxonomy1), FUN=sum)
m3 <- aggregate(m2, by=list(taxonomy), FUN=sum)
dim(m3)
dim(m4)
m3[1:5,1:4]
write.table(m3,file="SuperFilterASVs_Fin.csv",sep=";" , dec=",") #Guardo la tabla

# Sort the taxonomic table

if (tax_number > nrow(m3)) {
  tax_number <- nrow(m3)
}

m3$average <- rowMeans(m3[,-1])
m3.sorted <- m3[order(-m3$average),]

# Aggregate the smaller taxonomic bins together
m3.sorted$selection <- rep("discarded", nrow(m3.sorted))
m3.sorted$selection[1:tax_number] <- "retained"
m3.sorted$Group.1[m3.sorted$selection == "discarded"] <- "Other;Other"
m3.sorted$average <- NULL
m3.sorted$selection <- NULL
m5 <- aggregate(m3.sorted[,-1], by=list(taxonomy=m3.sorted$Group.1), FUN=sum)
head(m5)


# What is the relative abundances of the taxonomic bins that were pooled together in the 'Other' bin?
m5[m5$taxonomy == "Other;Other", -1]
write.table(m4,file="SuperFilterASVs_Fin.csv",sep=";" , dec=",") #Guardo la tabla

getwd()

# If you find these numbers too big, you can simply increase the value of the 'tax_number' variable,
# Or alternatively choose a higher taxonomic level to display

n <- m5$taxonomy
m5.t <- as.data.frame(t(m5[,-1]))
colnames(m5.t) <- n
m5.t$sample <- rownames(m5.t)
rownames(m5.t) <- NULL

#Calculate the mean and the standard deviation for each sample

m5.t$replicate <- rep(NA, nrow(m5.t))
for (line in 1:(nrow(m5.t))){
  m5.t$replicate[line] <- replicates_groups[m5.t$sample[line] == replicates_list]
}

# Compute the mean
m5.t.mean <- aggregate(m5.t[,1:(ncol(m5.t)-2)],
                       by = list(m5.t$replicate),
                       FUN = "mean")
names(m5.t.mean)[1] <- "sample"

dim(m5.t.mean)                              


# Compute the standard deviation
m5.t.sd <- aggregate(m5.t[,1:(ncol(m5.t)-2)],
                     by = list(m5.t$replicate),
                     FUN = "sd")
names(m5.t.sd)[1] <- "sample"

dim(m5.t.sd) 

#Melt and merge the two dataframes
# Melt the dataframes
molten.mean <- melt(m5.t.mean, id.vars = "sample")
molten.mean$id <- paste0(molten.mean$sample, "-", molten.mean$variable)

molten.sd <- melt(m5.t.sd, id.vars = "sample")
molten.sd$id <- paste0(molten.sd$sample, "-", molten.sd$variable)

# Merge the dataframes
molten <- merge(molten.mean, molten.sd, by.x = "id", by.y = "id")
head(molten)

#Final rearragement of the dataframe
molten$id <- NULL
molten$sample.y <- NULL
molten$variable.y <- NULL
names(molten) <- c("sample", "taxonomy", "mean", "sd")

molten$tax_col <- str_split_fixed(molten$taxonomy, ";", 2)[,1]
molten$tax_bin <- str_split_fixed(molten$taxonomy, ";", 2)[,2]

# Reorder the taxonomic annotation for the plot
molten <- molten[order(molten$tax_col),]
tax_levels <- as.character(molten$tax_bin[!duplicated(molten$tax_bin)])
tax_levels <- tax_levels[tax_levels != "Other"]
tax_levels <- c(tax_levels, "Other")
molten$tax_bin <- factor(molten$tax_bin, levels = rev(tax_levels))

# Reorder the samples for the plot
correct.order_Site <- c("Las Docas (N)", "Algarrobo (N)", "Navidad (C)",
                        "Topocalma (C)", "Pargua (SP)", "Los Chonos (SP)",
                        "Ilque (IS)", "San Antonio (IS)")
molten$sample <-factor(molten$sample,
                       levels = correct.order_Site)

levels(molten$sample)

# Remove null values
molten2 <- molten[molten$mean > 0,]
head(molten2)

write.table(molten2,file="Bubble_Table_Fin.csv",sep=";" , dec=",") #Save the table
getwd()

molten2$tax_bin

Alphaproteobacteria<-c("Fretibacter","Hellea","Litoreibacter","Litorimonas","Pacificibacter","Sulfitobacter",
                       "unassigned Hyphomonadaceae","unassigned Rhodobacteraceae","unassigned Terasakiellaceae",
                       "Yoonia-Loktanella")
Gammaproteobacteria<-c("Allofrancisella","Arenicella","Cocleimonas","Colwellia","Granulosicoccus","Leucothrix",
                       "Paraglaciecola","Perspicuibacter","Profundimonas", "Psychrobium","Thalassotalea","unassigned Colwelliaceae")
Bacteroidia<-c("Algitalea","Aquimarina","Croceitalea","Dokdonia","Lacinutrix","Pibocella","Polaribacter","Portibacter","Reichenbachiella",
               "Rubidimonas","Tenacibaculum","unassigned Cyclobacteriaceae","unassigned Flavobacteriaceae","unassigned Saprospiraceae")
Planctomycetes<-c("Blastopirellula")
Cyanobacteriia<-c("Pleurocapsa PCC-7319")
Bdellovibrionia<-c("Peredibacter","unassigned Bacteriovoracaceae")
Acidimicrobiia<-c("unassigned Microtrichaceae")
Verrucomicrobiae<-c("Persicirhabdus","Roseibacillus","Rubritalea")
Lev_tax_bin<-c(Alphaproteobacteria,Gammaproteobacteria,Bacteroidia,Planctomycetes,Cyanobacteriia,Bdellovibrionia,Acidimicrobiia,Verrucomicrobiae)
Lev_tax_bin

correct.order_Site <- c("Las Docas (N)", "Algarrobo (N)", "Navidad (C)",
                        "Topocalma (C)", "Ilque (SP)", "San Antonio (SP)",
                        "Pargua (IS)", "Los Chonos (IS)")
molten2$sample <-factor(molten2$sample,
                        levels = correct.order_Site)

levels(molten2$sample)

# Alphaproteobacteria	#66C5CC
# Gammaproteobacteria	#A6DEE2
# Bacteroidia	#F89C74
# Planctomycetes	#F6CF71
# Cyanobacteriia	#87C55F
# Bdellovibrionia	#8BE0A4
# Acidimicrobiia	#DCB0F2
# Verrucomicrobiae	#FE88B1

Lev_tax_col<-c("Alphaproteobacteria","Gammaproteobacteria","Bacteroidia","Planctomycetes","Cyanobacteriia","Bdellovibrionia","Acidimicrobiia","Verrucomicrobiae")
molten2$tax_col = factor(molten2$tax_col, levels=Lev_tax_col)
levels(molten2$tax_col)


colores_Class<-c("#66C5CC","#A6DEE2","#F89C74","#F6CF71","#87C55F","#8BE0A4","#DCB0F2","#FE88B1")
Lev_color<-c("#66C5CC","#A6DEE2","#F89C74","#F6CF71","#87C55F","#8BE0A4","#DCB0F2","#FE88B1")
colores_Class = factor(colores_Class, levels=Lev_color)
levels(colores_Class)

colores_Class<-c("Alphaproteobacteria"="#66C5CC","Gammaproteobacteria"="#A6DEE2","Bacteroidia"="#F89C74","Planctomycetes"="#F6CF71","Cyanobacteriia"="#87C55F","Bdellovibrionia"="#8BE0A4","Acidimicrobiia"="#DCB0F2","Verrucomicrobiae"="#FE88B1")

# Fill all NAs
molten2 <- molten2[!is.na(molten2$tax_bin), ]

molten2$tax_bin = factor(molten2$tax_bin, levels=Lev_tax_bin)
levels(molten2$tax_bin)
molten2$tax_bin <- factor(molten2$tax_bin, levels = rev(Lev_tax_bin))

bubble_plot <- ggplot(molten2,aes(sample,tax_bin)) +
  #geom_point(aes(size=mean+sd), shape=16, color = "red") + 
  geom_point(aes(size=mean, fill=tax_col),shape=21,color="black") +
  scale_fill_manual(values = colores_Class) +
  #scale_color_manual(values = colores_Class) +
  #geom_point(alpha = 0.9)
  #theme_ipsum() +
  #theme(legend.position = "none")
  scale_size(range = c(1,8), name="Relative\nabundance") +
  theme_bw()+
  theme(panel.grid.major=element_line(linetype=0.5,color="black"),
        axis.text.x=element_text(angle=70,hjust=1,vjust=0),
        plot.margin = margin (15,15,15,20),
        panel.background = element_blank()
  ) +
  
  #theme(legend.position=<desired-position>) +
  labs(
    y = "Taxonomic bins",
    x = "",
    fill = "Taxonomic\nclade",
    size = "Relative\nabundance"
  )
bubble_plot

png(file=paste("Bubble_Genus",".png",sep=""), units="in", width=11, height=8.5, res=300)
bubble_plot
dev.off()


bubble_plot <- ggplot(molten2, aes(sample, tax_bin)) +
  geom_point(aes(size = mean, fill = tax_col), shape = 21, color = "black") +
  scale_fill_manual(values = colores_Class) +
  scale_size(
    range = c(1, 8), 
    name = "Relative\nabundance",
    breaks = c(0.001,0.05, 0.1, 0.2, 0.4, 0.6),  
    labels = c("1%","5%", "10%", "20%", "40%", "60%")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linetype = 0.5, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
    plot.margin = margin(5, 5, 5, 5),
    panel.background = element_blank(),
    aspect.ratio = 1.6
      ) +
  labs(
    y = "Taxonomic bins",
    x = "",
    fill = "Taxonomic\nclade",
    size = "Relative\nabundance"
  )


bubble_plot

png(file=paste("Bubble_Genus_Fin",".png",sep=""), units="in", width=11, height=8.5, res=350)
bubble_plot
dev.off()
