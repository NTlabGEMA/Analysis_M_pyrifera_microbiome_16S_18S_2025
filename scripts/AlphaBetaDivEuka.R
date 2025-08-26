library(usethis)
library(devtools)
library(SpiecEasi)
library(microeco)
library(file2meco)
library(tidyverse)
library(dplyr)
library(remotes)
library(rgexf)
library(Biostrings)
library(ggplot2)
library(GUniFrac)
library(phyloseq)


#### Create and visualize a microeco obj####

meco0<-microtable$new(sample_table = ASV_Env, otu_table = ASV_Ab, tax_table = ASV_Tax)

#Clono microeco obj
meco1<- clone(meco0)
meco1$tidy_dataset()
meco1
meco1$tax_table
meco1$sample_table
meco1$otu_table
print(meco1)

meco1$cal_abund()
t1 <- trans_abund$new(dataset = meco1, taxrank = "Class", ntaxa = 100)
t1$plot_bar(others_color = "grey70", facet = "Site", xtext_keep = FALSE, legend_text_italic = FALSE)

meco1$sample_sums() %>% range
meco1$sample_sums()
meco1$taxa_sums()
meco1$taxa_abund$Division

######Alpha Diver####
meco1$cal_alphadiv(PD = FALSE)
Alpha_18S<-meco1$alpha_diversity
write.table(Alpha_18S, "Alpha_18S.txt", row.names = T, quote = F, sep = '\t')

SampleTable<-meco1$sample_table
SampleTable

Alpha_18S$SampleID <- rownames(Alpha_18S)
SampleTable$SampleID <- rownames(SampleTable)
DataAlpha<- merge(Alpha_16S, SampleTable, by = "SampleID")

anova_result_Obs <- aov(Observed ~ Site, data = DataAlpha)
summary(anova_result_Obs)  # Resumen de la ANOVA
anova_result_Shan <- aov(Shannon ~ Site, data = DataAlpha)
summary(anova_result_Shan)  # Resumen de la ANOVA

anova_result_Obs <- aov(Observed ~ Genotype, data = DataAlpha)
summary(anova_result_Obs)  # Resumen de la ANOVA
anova_result_Shan <- aov(Shannon ~ Genotype, data = DataAlpha)
summary(anova_result_Shan)  # Resumen de la ANOVA

anova_result_Obs <- aov(Observed ~ Fenotype, data = DataAlpha)
summary(anova_result_Obs)  # Resumen de la ANOVA
anova_result_Shan <- aov(Shannon ~ Fenotype, data = DataAlpha)
summary(anova_result_Shan)  # Resumen de la ANOVA


meco1$sample_table
Mp_Alpha_S <- trans_alpha$new(dataset = meco1, group = "Site")
Mp_Alpha_S$data_stat

Mp_Alpha_S$cal_diff(method = "anova")
Mp_Alpha_S$res_diff


meco1$sample_table
Mp_Alpha <- trans_alpha$new(dataset = meco0, group = "Site")
Mp_Alpha$data_stat
Mp_Alpha$cal_diff(method = "anova")
# return t1$res_diff
head(Mp_Alpha$res_diff)

Mp_Alpha <- trans_alpha$new(dataset = meco1, group = "Genotype")
Mp_Alpha$data_stat
Mp_Alpha$cal_diff(method = "anova")
# return t1$res_diff
head(Mp_Alpha$res_diff)

######Beta Diver#####

meco1$sample_table

correct.order_GT <- c("North", "Center",
                      "South_Pacific","Inner_Sea")

correct.order_ST <- c("Algarrobo","Topocalma","Los Chonos",
                      "Ilque","San Antonio")

meco1$sample_table$Genotype<- factor(meco1$sample_table$Genotype,
                                     levels = correct.order_GT)
levels(meco1$sample_table$Genotype)
meco1$sample_table$Site<- factor(meco1$sample_table$Site,
                                 levels = correct.order_ST)

meco1$cal_betadiv()
t10 <- trans_beta$new(dataset = meco1, group = "Site", measure = "bray")

t10$cal_betadisper()
t10$res_betadisper

t11 <- trans_beta$new(dataset = meco1, group = "Genotype", measure = "bray")

t11$cal_betadisper()
t11$res_betadisper

t12 <- trans_beta$new(dataset = meco1, group = "Fenotype", measure = "bray")

t12$cal_betadisper()
t12$res_betadisper

t10$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t10$res_ordination)
pcoa1<-t10$res_ordination

#Docas → #F1A8A8
# Algarrobo → #F3C8A0
# Navidad → #E4D88C
# Topocalma → #A8A87C
# Ilque → #A3D3B4
# San Antonio → #A7C8D2
# Pargua → #B4B4F2
# Los Chonos → #A4A4D1

# Mejoramos la personalización del gráfico
p0 <- t10$plot_ordination(
  plot_color = "Site",
  plot_shape = "Site",
  #  plot_type = c("point", "ellipse", "centroid"),
  plot_type = c("point", "ellipse"),
  loading_arrow = TRUE,
  loading_text_italic = TRUE,
  loading_taxa_num = 5
)

# Ajustar manualmente las capas para resolver los problemas y cambiar el tamaño de los puntos
p1 <- p0 +
  theme_bw() +
  theme(
    text = element_text(size = 24, color = "black"), # Texto en negro y más grande
    axis.text = element_text(color = "black", size = 22), # Números de los ejes en negro
    axis.title = element_text(color = "black", size = 22, face = "bold"), # Títulos de ejes en negro y negritas
    axis.line = element_line(color = "black", linewidth = 0,5), # Líneas de ejes más oscuras
    axis.ticks = element_line(color = "black", linewidth = 1), # Marcas de los ejes más oscuras
    panel.border = element_rect(color = "black", linewidth = 1), # Marco del panel más grueso y negro
    legend.key = element_rect(fill = "transparent", color = NA), # Fondo transparente en la leyenda
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    panel.grid = element_blank(), # Eliminar cuadrícula
    panel.background = element_rect(fill = "white", color = NA) # Fondo limpio
  ) +
  # Personalizar las formas
  scale_shape_manual(
    values = c(1, 15, 17, 18, 18) # Círculos abiertos y formas sólidas
  ) +
  scale_color_manual(
    values = c("#E78B3A","#7A7A40","#5656AF","#46AF65","#4A96AB")
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(size = 5, fill = NA, stroke = 1.5) # Sin relleno en la leyenda
    ),
    color = guide_legend(
      override.aes = list(size = 5) # Ajustar tamaño en la leyenda
    )
  )

# Modificar manualmente las capas existentes en el gráfico
# Acceder a las capas generadas por plot_ordination
p1$layers <- lapply(p1$layers, function(layer) {
  if (inherits(layer$geom, "GeomPolygon")) {
    # Asegurarse de que las elipses no tengan relleno
    layer$aes_params$fill <- NA
    layer$aes_params$alpha <- 0
  }
  if (inherits(layer$geom, "GeomPoint")) {
    # Aumentar el espesor del borde de los puntos
    layer$aes_params$stroke <- 1.5 # Aumentar el espesor del borde
    layer$aes_params$size <- 6     # Mantener el tamaño de los puntos
  }
  layer
})

# Mostrar el gráfico
p1
png(file=paste("PCoA_Euk_sinSpider",".png",sep=""), units="in", width=11, height=8.5, res=300)
p1
dev.off()

#Crear Leyenda:
library(cowplot)
# Crear un dataframe con los datos para la leyenda
legend_data <- data.frame(
  Site = factor(c( "Algarrobo","Topocalma", "Los Chonos",
                   "Ilque","San Antonio"), 
                levels = c( "Algarrobo","Topocalma", "Los Chonos",
                            "Ilque","San Antonio")),
  x = 1:5,
  y = 1:5
)
legend_data
# Crea un gráfico "vacío" con solo la leyenda
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Site, shape = Site)) +
  # Aumentar significativamente el tamaño de los puntos y el grosor del contorno
  geom_point(size = 10, stroke = 2) +  # Aumentado de 5 a 10
  scale_color_manual(values = c("#E78B3A","#7A7A40","#5656AF","#46AF65","#4A96AB")) +
  scale_shape_manual(values = c( 1, 15, 17, 18, 18)) +
  theme_void() +
  theme(
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    # Aumentar el tamaño de las "keys" en la leyenda
    legend.key.size = unit(2, "cm"),  # Aumentado significativamente
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(1.5, "cm")
  ) +
  # Asegurar que los símbolos sean grandes en la leyenda
  guides(
    color = guide_legend(override.aes = list(size = 10)),
    shape = guide_legend(override.aes = list(size = 10))
  )

# Extraer solo la leyenda del gráfico
legend_only <- cowplot::get_legend(legend_plot)

# Crear un nuevo gráfico que contenga solo la leyenda
legend_final <- ggplot() + 
  theme_void() +
  annotation_custom(legend_only)
legend_final
# Guardar la leyenda como un archivo independiente
# Aumentar un poco el tamaño del archivo para dar más espacio
ggsave("legend_only.pdf", legend_final, width = 4, height = 5)
# También en formato PNG por si lo prefieres
ggsave("legend_only.png", legend_final, width = 4, height = 5, dpi = 300)

###PERMANOVA
meco1$sample_table<-meco0$sample_table

t1 <- trans_beta$new(dataset = meco1, group = "Fenotype", measure = "bray")
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova