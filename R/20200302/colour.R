library(scales)
library(ggsci)

library(RColorBrewer)


RColorBrewer::display.brewer.all()

color1 <- RColorBrewer::brewer.pal(name = "Set3", n = 10)
show_col(color1)
color1

# library(bbplot)
# 
# 
# ggplot(data = mtcars, aes(mpg, qsec, colour = factor(am))) +
#   geom_point() +
#   bbc_style()
# 
# 
# 
# line <- ggplot(line_df, aes(x = year, y = lifeExp)) +
#   geom_line(colour = "#007f7f", size = 1) +
#   geom_hline(yintercept = 0, size = 1, colour="#333333") +
#   bbc_style()


