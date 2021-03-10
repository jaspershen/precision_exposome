#####figure
library(MSnbase)
library(xcms)
library(mzR)
libray(metflow2)
sxtTools::setwd_project()
setwd("test/Positive/")
f.in <- list.files(path = ".",
                   pattern = '\\.(mz[X]{0,1}ML|cdf)',
                   recursive = TRUE)

sample_group <-
  unlist(lapply(stringr::str_split(string = f.in, pattern = "/"), function(x) {
    x[1]
  }))

sample_group[grep("\\.(mz[X]{0,1}ML|cdf)", sample_group)] <-
  "Group0"

pd <-
  data.frame(
    sample_name = sub(
      basename(f.in),
      pattern = ".mzXML",
      replacement = "",
      fixed = TRUE
    ),
    sample_group = sample_group,
    stringsAsFactors = FALSE
  )

## Define colors for different groups
group_colors <-
  paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_group))], "60")
# group_colors <-
#   grDevices::colorRampPalette(ggsci::pal_npg()(10))(100)[1:length(unique(sample_group))]
names(group_colors) <- unique(sample_group)

raw_data <- MSnbase::readMSData(
  files = f.in,
  pdata = new("NAnnotatedDataFrame", pd),
  mode = "onDisk",
  verbose = TRUE
)




cwp <- xcms::CentWaveParam(ppm = 15,
                           peakwidth = c(5, 30),)

xdata <-
  try(xcms::findChromPeaks(
    raw_data,
    param = cwp,
    BPPARAM = BiocParallel::SnowParam(workers = 4,
                                      progressbar = TRUE)
  ),
  silent = FALSE)

if (class(xdata) == "try-error") {
  stop("Error in xcms::findChromPeaks.\n")
}



------------------------------------------------------------#retention time correction
  #Alignment
  
  
  
  xdata2 <- try(xcms::adjustRtime(xdata,
                                  param = xcms::ObiwarpParam(binSize = 0.5)),
                silent = FALSE)

xdata2 <- xdata


cat(crayon::red(clisymbols::symbol$tick, "OK\n"))


###TIC

cat(crayon::green("Drawing TIC plot..."))
tic.plot <- xcms::chromatogram(
  object = xdata2,
  aggregationFun = "sum",
  BPPARAM =
    BiocParallel::SnowParam(workers = 4,
                            progressbar = TRUE)
)

plot <- metflow2::chromatogramPlot(object = tic.plot,
                                   title = "TIC",
                                   group.for.figure = "QC")


plot


plot +
  scale_color_manual(values = c("729" == alpha("white", alpha = 0), 
                                "A1" = alpha("white", alpha = 0),
                                "731" = "skyblue",
                                "A2" = "red"))




ggsave("731_A2.pdf", width = 8, height = 6)






###BPC

bpc.plot <- xcms::chromatogram(
  object = xdata2,
  aggregationFun = "max",
  BPPARAM =
    BiocParallel::SnowParam(workers = 4,
                            progressbar = TRUE)
)

plot <- metflow2::chromatogramPlot(object = bpc.plot,
                         title = "",
                         group.for.figure = "QC")

ggsave("plot.pdf", width = 8, height = 6)



######
peak_table <- readr::read_csv("Peak_table_for_cleaning_pos.csv")

peak_table <- 
  peak_table %>% 
  dplyr::select(name, 
                mz, rt,
                X729, 
                X731,
                A1, A2)

peak_table <- 
  peak_table %>% 
  as.data.frame() 
peak_table[which(is.na(peak_table), arr.ind = TRUE)] <- 1


peak_table <- 
  peak_table %>% 
  as_tibble()


peak_table2 <- 
  peak_table %>% 
  dplyr::select( X729, 
                 X731,
                 A1, A2) %>% 
  log(10) %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>%
  t() %>%
  as_tibble()



peak_table <- 
  peak_table %>% 
  select(name, mz, rt) %>% 
  data.frame(peak_table2)


peak_table <- 
  peak_table %>% 
  as_tibble()



  peak_table %>% 
  dplyr::mutate(value = A1 - X729) %>% 
    dplyr::filter(!is.na(value)) %>% 
  ggplot2::ggplot(aes(x = rt, y = mz, colour = value)) +
    scale_colour_gradient2(low = "skyblue", mid = "white", high = "red") +
  geom_point() +
    theme_bw() +
    labs(x = "Retention time (seconds)", y = "m/z") +
    theme(axis.title = element_text(size = 15), 
          axis.text = element_text(size = 13))
  
  
  peak_table %>% 
    dplyr::mutate(value = A2 - X731) %>% 
    dplyr::filter(!is.na(value)) %>% 
    ggplot2::ggplot(aes(x = rt, y = mz, colour = value)) +
    scale_colour_gradient2(low = "skyblue", mid = "white", high = "red") +
    geom_point() +
    theme_bw() +
    labs(x = "Retention time (seconds)", y = "m/z") +
    theme(axis.title = element_text(size = 15), 
          axis.text = element_text(size = 13))
  
  
  
  
  peak_table %>% 
    dplyr::mutate(value = A1 - X729) %>% 
    dplyr::filter(!is.na(value)) %>% 
    pull(value) %>% 
    `>`(0) %>% 
    sum() 
    
  
  peak_table %>% 
    dplyr::mutate(value = A1 - X729) %>% 
    dplyr::filter(!is.na(value)) %>% 
    nrow()
  
  
  
  peak_table %>% 
    dplyr::mutate(value = A2- X731) %>% 
    dplyr::filter(!is.na(value)) %>% 
    pull(value) %>% 
    `>`(0) %>% 
    sum() 
  
  
  peak_table %>% 
    dplyr::mutate(value = A2- X731) %>% 
    dplyr::filter(!is.na(value)) %>% 
    nrow()
  
  2153/4140
  
  
  
  













