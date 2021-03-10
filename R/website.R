tinyTools::setwd_project()
rmarkdown::render_site(input = "website_files/")
file.copy(from = "website_files/docs/", to = ".", recursive = TRUE, overwrite = TRUE)
