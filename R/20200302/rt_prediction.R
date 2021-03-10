sxtTools::setwd_project()
setwd("rt_prediction/")
load("../data/metabolite_table")
load("../data/metabolite_tags")

##m

data <- 
metabolite_tags %>% 
  filter(Note == "Confirmed and Published")


data$MetabID


InChIKey <- lapply(data$MetabID, function(x){
  metflow2::transID(query = x, 
                    from = "Chemical Name", 
                    to = "InChIKey", 
                    top = 1)[,2]  
}) %>% 
  unlist()

library(webchem)

inchi <- sapply(InChIKey, function(x){
  cs_inchikey_inchi(inchikey = x)  
})


smiles <- sapply(inchi, function(x){
  cs_inchi_smiles(inchi = x)  
})


smiles <- unname(smiles)



library(rcdk)


descNames <-
  unique(unlist(sapply(get.desc.categories(), get.desc.names)))


##get all the md descriptors for the standards.
standard.mols = vector(mode = "list", length = length(smiles))
for (i in 1:length(standard.mols)) {
  cat(i, " ")
  temp <- try(parse.smiles(smiles[i]))
  if (class(temp) == "try-error") {
    cat('Error\n')
    temp <- NA
  }
  standard.mols[[i]] <- temp[[1]]
}

save(standard.mols, file = "standard.mols")

standard.md <-
  pbapply::pblapply(standard.mols, function(x) {
    md <- try(eval.desc(x, descNames))
    if (class(md) == "try-error") {
      md <- matrix(rep(NA, 286), nrow = 1)
    }
    md
  })

col_name = colnames(standard.md[[1]])

standard.md <-
  lapply(standard.md, function(x) {
    colnames(x) = col_name
    x
  })



standard.md <-
  standard.md %>%
  do.call(rbind, .)


##remove MD are all NA
remove.idx <-
  which(apply(standard.md, 2, function(x) {
    all(is.na(x))
  })) %>%
  unname()

if(length(remove.idx) > 0){
  standard.md <-
    standard.md[, -remove.idx]  
}


standard.md <- 
  t(standard.md)


colnames(standard.md) <- data$MetabID


head(standard.md)


sum(is.na(standard.md))



train.y <- data$rt
train.x <- t(standard.md)


names(train.y) <- rownames(train.x)


plot(train.y, train.x[,2])


#use the marker we used from MetDNA
marker.name <- c('XLogP', "WTPT.4", "WTPT.5", "ALogp2", "BCUTp.1l")


# idx <- match(marker.name, colnames(train.x))

# train.x <- train.x[, idx]

para <- NULL
ntree1 <- seq(300, 1000, by = 200)
mtry1 <- seq(1, length(marker.name), by = 1)
for (i in 1:length(ntree1)) {
  para <- rbind(para, cbind(ntree1[i], mtry1))
}
colnames(para) <- c("ntree", "mtry")
mse <- NULL
rf.reg <- list()
cat("\n")
cat("Find the optimal parameters\n")
for (i in 1:nrow(para)) {
  cat(i)
  cat(" ")
  temp.ntree <- para[i, 1]
  temp.mtry <- para[i, 2]
  rf.reg[[i]] <-
    randomForest::randomForest(
      x = train.x,
      y = train.y,
      ntree = temp.ntree,
      mtry = temp.mtry,
      replace = TRUE,
      importance = TRUE,
      proximity = TRUE
    )
  mse[i] <- mean(rf.reg[[i]]$mse)
}
cat("\n")
result <- data.frame(para, mse, stringsAsFactors = FALSE)
temp.idx <- which.min(result$mse)
ntree <- result$ntree[temp.idx]
mtry <- result$mtry[temp.idx]
##


rf.reg <- randomForest::randomForest(
  x = train.x,
  y = train.y,
  ntree = ntree,
  mtry = mtry,
  replace = TRUE,
  importance = TRUE,
  proximity = TRUE
)


predict.rt <- NULL
for(i in 1:length(train.y)){
  rf.reg <- randomForest::randomForest(
    x = train.x[-i,],
    y = train.y[-i],
    ntree = ntree,
    mtry = mtry,
    replace = TRUE,
    importance = TRUE,
    proximity = TRUE
  )
  
predict.rt[i] <-
  predict(object = rf.reg, newdata = train.x[i, , drop = FALSE])
}



data.frame(predict.rt, train.y, stringsAsFactors = FALSE) %>% 
  filter(train.y < 20) %>% 
  ggplot(aes(x = predict.rt, y = train.y)) +
  geom_abline(slope = 1, intercept = 0) +  
  geom_point() +
  theme_bw() +
  scale_x_continuous(limits = c(3,15)) +
  scale_y_continuous(limits = c(3,15)) +
  labs(x = "Predicted RT", y = "Actual RT")








