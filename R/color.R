library(ggsci)
library(scales)

test_colour3 <- ggsci::pal_aaas(alpha = 1)(12)
show_col(colours = c(test_colour3), borders = NA)
test_colour3


test_colour4 <- ggsci::pal_aaas(alpha = 0.7)(12)
show_col(colours = c(test_colour4), borders = NA)
test_colour4



RColorBrewer::display.brewer.all()

test_color5 <- RColorBrewer::brewer.pal(n = 12,name = "Set3")
scales::show_col(test_color5)
test_color5



## ggplot2 theme

## color for different data
# scale_fill_manual(values = c("Clinical" = "#8DD3C7", "Metabolome" = "#FB8072",
#                              "Transcriptome" = '#80B1D3', 
#                              'Proteome' = '#FDB462')) 


#theme
#  theme(axis.title = element_text(size = 15),
# axis.text = element_text(size = 13),
# strip.text = element_text(size = 13))
# 
##ledgend position
##


library(scales)
scales::show_col(ggsci::rgb_gsea(n = 21))
scales::show_col(ggsci::pal_aaas()(21))
scales::show_col(ggsci::pal_d3()(21))
scales::show_col(ggsci::pal_futurama()(21))
scales::show_col(ggsci::pal_gsea()(21))
scales::show_col(ggsci::pal_igv()(21))
scales::show_col(ggsci::pal_jama()(21))
scales::show_col(ggsci::pal_jco()(21))
scales::show_col(ggsci::pal_lancet()(21))
scales::show_col(ggsci::pal_locuszoom()(21))
scales::show_col(ggsci::pal_material()(21))
scales::show_col(ggsci::pal_nejm()(21))
scales::show_col(ggsci::pal_npg()(21))
scales::show_col(ggsci::pal_rickandmorty()(21))
scales::show_col(ggsci::pal_simpsons()(21))
scales::show_col(ggsci::pal_startrek()(21))
scales::show_col(ggsci::pal_tron()(21))
scales::show_col(ggsci::pal_uchicago()(21))
scales::show_col(ggsci::pal_ucscgb()(21))
                 

library(wesanderson)
library(tidyverse)

wes_palettes %>% 
  lapply(length) %>% 
  unlist()

scales::show_col(colours = wes_palettes$BottleRocket1)
scales::show_col(colours = wes_palettes$BottleRocket2)
scales::show_col(colours = wes_palettes$Rushmore)
scales::show_col(colours = wes_palettes$Royal1)
scales::show_col(colours = wes_palettes$Royal2)
scales::show_col(colours = wes_palettes$IsleofDogs2)
scales::show_col(colours = wes_palettes$IsleofDogs1)
scales::show_col(colours = wes_palettes$GrandBudapest2)
scales::show_col(colours = wes_palettes$Cavalcanti1)

wes_palette(name = wes_palettes$BottleRocket1, n = 5)



scales::show_col(colours = wes_palettes$Cavalcanti1)


wes_palettes %>% 
  lapply(length) %>% 
  unlist()

wes_palette(name = "Cavalcanti1", n = 5, type = "discrete")



scales::show_col(ggsci::pal_d3()(10))

value <- 
  c("wearable" = ggsci::pal_d3()(10)[1],
    "exposome" = ggsci::pal_d3()(10)[2],
    "metabolome" = ggsci::pal_d3()(10)[3],
    "proteome" = ggsci::pal_d3()(10)[4],
    "microbiome" = ggsci::pal_d3()(10)[5],
    "gut_microbiome" = ggsci::pal_d3()(10)[6],
    "nasal_microbiome" = ggsci::pal_d3()(10)[7],
    "cytokine" = ggsci::pal_d3()(10)[8])
