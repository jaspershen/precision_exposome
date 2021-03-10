##consturct the kegg network
# sxtTools::setwd_project()
# setwd("data_analysis_2020_02_13/metabolome_data_analysis/marker_for_labor/module_analysis/")
# load("edge_info")
# load("node_info")
# load("reaction_info")
#
# edge_info[[1]]
#
# library(igraph)
#
# kegg_network <-
#   do.call(rbind, reaction_info)
#
#
# kegg_network <-
# kegg_network %>%
#   dplyr::select(substrateName, productName) %>%
#   distinct(substrateName, productName)
#
# ##remove the G compound
# kegg_network <-
#   kegg_network %>%
#   dplyr::filter(!stringr::str_detect(substrateName, 'G')) %>%
#   dplyr::filter(!stringr::str_detect(productName, 'G')) %>%
#   dplyr::filter(!stringr::str_detect(substrateName, 'D')) %>%
#   dplyr::filter(!stringr::str_detect(productName, 'D'))
#
#
# kegg_network <-
#   kegg_network %>%
#   dplyr::mutate(substrateName = stringr::str_replace(substrateName, pattern = " cpd", ""),
#                 productName = stringr::str_replace(productName, pattern = " cpd", "")
#                 )
#
# id <-
#   c(kegg_network$substrateName,kegg_network$productName)
#
# ##remove duplicated rows
# id2 <- paste(kegg_network$productName, kegg_network$substrateName, sep = ";")
#
# remove_idx <- vector(mode = "list", length = length(id2))
#
# for(i in 1:nrow(kegg_network)){
#   cat(i, " ")
#   temp <- paste(kegg_network[i,1], kegg_network[i,2], sep = ";")
#   idx <- match(temp, id2)
#   remove_idx[[i]] <- c(i, idx)
# }
#
#
# remove_idx <-
#   remove_idx %>%
#   do.call(rbind, .)
#
# remove_idx <-
# remove_idx %>%
#   as.data.frame() %>%
#   dplyr::filter_all(all_vars(!is.na(.))) %>%
#   dplyr::filter(V2 > V1)
#
# kegg_network[remove_idx[100,1],]
# kegg_network[remove_idx[100,2],]
#
#
# kegg_network <- kegg_network[-remove_idx$V2,]
#
# nrow(kegg_network)
#
# length(unique(c(kegg_network$substrateName, kegg_network$productName)))
#
# ##remove metabolites which are not in KEGG compound database
# kegg_network <-
#   kegg_network %>%
#   dplyr::filter(substrateName %in% kegg.compound$ID,
#                 productName %in% kegg.compound$ID)
#
# ##get the kegg compound database
# load("kegg.compound.rda")
#
# kegg_id <- unique(c(kegg_network$substrateName, kegg_network$productName)) %>%
#   data.frame(id = ., stringsAsFactors = FALSE)
#
# which(is.na(match(kegg_id$id, kegg.compound$ID)))
#
# node_info <-
#   kegg_id %>%
#   dplyr::left_join(kegg.compound, by = c("id" = "ID"))
#
#
# kegg_network <-
# igraph::graph_from_data_frame(d = kegg_network, directed = FALSE, vertices = node_info)
#
#
# save(kegg_network, file = "kegg_network")


setGeneric(name = "getModule",
           function(node.id,
                    metabolic.network,
                    threads = 3,
                    max.reaction = 3,
                    progressbar = TRUE,
                    calculate.impact = TRUE) {
unique.id <- unique(node.id)
unique.id <-
  unique.id[unique.id %in% V(metabolic.network)$name]
total.id.number <- length(unique.id)

##find the hidden metabolite which can connect two input metabolites in up to
## three reaction
temp.fun <- function(name,
                     metabolic.network,
                     unique.id,
                     max.reaction) {
  options(warn = -1)
  library(
    igraph,
    quietly = TRUE,
    logical.return = FALSE,
    warn.conflicts = FALSE
  )
  idx <- match(name, unique.id)
  to <- unique.id[-c(1:idx)]
  
  distance <-
    igraph::distances(graph = metabolic.network, v = name, to = to)[1,]
  distance[is.infinite(distance)] <- 5
  to <- to[which(distance <= max.reaction)]
  if (length(to) == 0)
    return(NULL)
  result <-
    igraph::shortest_paths(graph = metabolic.network,
                           from = name,
                           to = to)[[1]]
  # len <- unlist(lapply(result, function(x) length(x)))
  # result <- result[len<=4]
  result <- lapply(result, names)
  result <- unique(unlist(result))
  result <- setdiff(result, unique.id)
  rm(list = c("idx", "to", "distance"))
  gc()
  return(result)
}

if (progressbar)
  cat("Find hidden metabolites.\n")
hidden.id <-
  BiocParallel::bplapply(
    unique.id[-length(unique.id)],
    FUN = temp.fun,
    BPPARAM = BiocParallel::SnowParam(workers = threads,
                                      progressbar = progressbar),
    metabolic.network = metabolic.network,
    unique.id = unique.id,
    max.reaction = max.reaction
  )
             
hidden.id <- unlist(hidden.id) %>% 
  unique()

subnetwork <-
  igraph::subgraph(graph = metabolic.network,
                   v = unique(c(unique.id, hidden.id)))

wt <- igraph::cluster_walktrap(subnetwork)

node.name <- wt$names
membership <- wt$membership
membership1 <- table(membership)
             
group <- lapply(names(membership1), function(x) {
  node.name[which(membership == as.numeric(x))]
})
             
##remove group which less than 3 nodes
group <- group[unlist(lapply(group, length)) > 3]

###the group whose node more than 100 should be re-subgraph
idx <- which(unlist(lapply(group, length)) > 100)
             
while (length(idx) > 0) {
  large.group <- group[idx]
  group <- group[-idx]
  new.group <- lapply(large.group, function(x) {
    temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
    temp.wt <- igraph::cluster_walktrap(temp.graph)
    temp.node.name <- temp.wt$names
    temp.membership <- temp.wt$membership
    temp.membership1 <- table(temp.membership)
    
    temp.group <-
      lapply(names(temp.membership1), function(x) {
        temp.node.name[which(temp.membership == as.numeric(x))]
      })
    
    ##remove group which less than 3 nodes
    temp.group <-
      temp.group[unlist(lapply(temp.group, length)) > 3]
    temp.group
  })
  
  if (length(new.group) == 1) {
    new.group <- new.group[[1]]
  } else{
    new.group1 <- new.group[[1]]
    for (i in 2:length(new.group)) {
      new.group1 <- c(new.group1, new.group[[i]])
    }
    new.group <- new.group1
  }
  group <- c(group, new.group)
  idx <- which(unlist(lapply(group, length)) > 100)
}
             
##remove group which edge number is less than node number
remove.idx <- which(unlist(lapply(group, function(x) {
  temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
  node.number <- length(igraph::V(temp.graph))
  edge.number <- length(igraph::E(temp.graph))
  diff.number <- node.number - edge.number
  degree <- igraph::degree(graph = temp.graph, v = x)
  if (diff.number > 0 & all(degree < 3)) {
    return(TRUE)
  } else{
    FALSE
  }
})))
             
group <- group[-remove.idx]
             
##remove the hidden nodes whose degree is less than 2
if (progressbar)
  cat("Remove the hidden metabolites whose degrees are less than two.\n")
group <- lapply(group, function(x) {
  temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
  degree <- igraph::degree(graph = temp.graph, v = x)
  
  class <- sapply(names(degree), function(x) {
    ifelse(is.element(x, unique.id), "Detetcted", "Hidden")
  })
  
  remove.idx <- which(degree == 1 & class == "Hidden")
               
  while (length(remove.idx) > 0) {
    x <- x[-remove.idx]
    temp.graph <-
      igraph::subgraph(graph = metabolic.network, v = x)
    
    degree <- igraph::degree(graph = temp.graph, v = x)
    
    class <- sapply(names(degree), function(x) {
      ifelse(is.element(x, unique.id), "Detetcted", "Hidden")
    })
    remove.idx <-
      which(degree == 1 & class == "Hidden")
    
  }
  x
})
             
group <- group[unlist(lapply(group, length)) > 3]

# pbapply::pboptions(type = "timer", style = 1)
if (progressbar) {
  cat("\n")
  cat("Calculate activity scores of dysregulated modules.\n")
}
             
temp.fun <-
  function(x,
           unique.id,
           hidden.id,
           total.id.number,
           metabolic.network = metabolic.network) {
    library(
      igraph,
      quietly = TRUE,
      logical.return = FALSE,
      warn.conflicts = FALSE
    )
    options(warn = -1)
    
    temp.graph <-
      igraph::subgraph(graph = metabolic.network, v = x)
    input.number <- sum(x %in% unique.id)
    hidden.number <- sum(x %in% hidden.id)
    
    m <- length(igraph::E(metabolic.network))
    e <- length(igraph::E(temp.graph))
    node.name <- names(igraph::V(temp.graph))
    
    value <- sapply(node.name, function(node) {
      temp.value <- sapply(node.name, function(x) {
        temp1 <- igraph::degree(graph = metabolic.network, v = x)
        temp2 <-
          igraph::degree(graph = metabolic.network, v = node)
        temp1 * temp2 / (4 * m ^ 2)
      })
      sum(temp.value)
    })
    
    modularity <- e / m - sum(value)
    # modularity
    q <- modularity * sqrt(total.id.number / length(x))
    a <- input.number * q / length(x)
    a
  }
             
             
system.time(
  activity.score <-
    BiocParallel::bplapply(
      group,
      temp.fun,
      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                        progressbar = progressbar),
      unique.id = unique.id,
      hidden.id = hidden.id,
      total.id.number = total.id.number,
      metabolic.network = metabolic.network
    )
)
             
activity.score <- unlist(activity.score)


if (calculate.impact) {
  ##Module impact
  if (progressbar)
    cat("Calculate impacts of dysregulated modules.\n")
  impact <- lapply(group, function(temp.group) {
    temp.graph <- igraph::subgraph(graph = metabolic.network,
                                   v = temp.group)
    centrality <-
      getCentrality(graph = temp.graph, type = "d")
    temp.idx <- which(temp.group %in% unique.id)
    impact <-
      sum(centrality[temp.idx]) / sum(centrality)
    impact
  })
  impact <- unlist(impact)
  module <- group
  result <-
    list(module, activity.score, hidden.id, impact)
  names(result) <-
    c("module", "activity.score", "hidden.id", "impact")
} else{
  module <- group
  result <- list(module, activity.score, hidden.id)
  names(result) <-
    c("module", "activity.score", "hidden.id")
}
result <- result
           })


#-----------------------------------------------------------------------------
# title getCentrality
# description Get the centrality of all nodes in nwtwork.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param graph The graph.
# param type degree or betweenness.
# return  The centrality of all nodes.


setGeneric(
  name = "getCentrality",
  def = function(graph,
                 type = c("degree", "betweenness")) {
    type = match.arg(type)
    node <- igraph::V(graph)$name
    if (type == "degree") {
      node.degree <- igraph::degree(graph = graph, v = node)
      return(node.degree)
    }
    
    node.betweenness <- lapply(node, function(x) {
      i <- match(x, node)
      if (i == length(node))
        return(NULL)
      from <- x
      to <- node[(i + 1):length(node)]
      temp  <- igraph::shortest_paths(graph = graph,
                                      from = from,
                                      to = to)[[1]]
      temp <- lapply(temp, function(x) {
        if (length(x) <= 2)
          return(NULL)
        names(x)[-c(1, length(x))]
      })
      unlist(temp)
    })
    
    node.betweenness <- unlist(node.betweenness)
    node.betweenness <- table(node.betweenness)
    node0 <- setdiff(node, names(node.betweenness))
    if (length(node0) == 0)
      return(node.betweenness)
    add <- rep(0, length(node0))
    names(add) <- node0
    node.betweenness <- c(node.betweenness, add)
    node.betweenness <-
      node.betweenness[match(node, names(node.betweenness))]
    return(node.betweenness)
  }
)



setGeneric(
  name = "keggMatch",
  def = function(mz,
                 rt,
                 metabolite,
                 mz.tol = 25,
                 rt.tol = 30,
                 polarity = c("positive", "negative")) {
    
options(warn = -1)
polarity <- match.arg(polarity)
metabolite <- metabolite[metabolite$Mode == polarity, ]
metabolite <-
  metabolite[metabolite$Accurate.mass <= max(mz), ]
    
#----------------------------------------------------------------
metabolite.mz <-
  sapply(metabolite$Accurate.mass, function(mass) {
    ifelse(mass >= 400, mass, 400)
  })
    
result <- mapply(function(x, y) {
  # mz.error <- abs(x - metabolite$Accurate.mass)*10^6/metabolite$Accurate.mass
  mz.error <-
    abs(x - metabolite$Accurate.mass) * 10 ^ 6 / metabolite.mz
  # rt.error <- abs(y - metabolite$RT) * 100 / metabolite$RT
  index <-
    which(mz.error <= mz.tol)
  if (length(index) == 0)
    return(NULL)
  temp <- data.frame(metabolite[index, , drop = FALSE],
                     mz.error = mz.error[index],
                     # rt.error[index],
                     stringsAsFactors = FALSE)
  list(temp)
},
x = mz,
y = rt)
    
result <- lapply(result, function(x) {
  if (is.null(x))
    return(NULL)
  if (length(x) == 1)
    return(x[[1]])
  return(x)
})

result <- result
  }
)


setGeneric(
  name = "getNullDistribution",
  def = function(all.mz.pos,
                 all.mz.neg,
                 metabolic.network = kegg.rpair2,
                 size.pos,
                 size.neg,
                 times = 20,
                 threads = 5,
                 mz.tol = 25,
                 progressbar = FALSE,
                 kegg.adduct.pos,
                 kegg.adduct.neg) {
acti.score <- vector(mode = "list", length = times)
    
for (i in seq_len(times)) {
  # if(i %% 10 == 0 | i == 1) cat(i); cat(" ")
  cat(i, "/", times)
  cat("\n")
  # if (polarity == "positive" | polarity == "both") {
  idx.pos <- sample(1:length(all.mz.pos), size.pos) %>% sort()
  temp.mz.pos <- as.numeric(all.mz.pos[idx.pos])
    
  match.result.pos <- keggMatch(
    mz = temp.mz.pos,
    rt = rep(100, length(temp.mz.pos)),
    metabolite = kegg.adduct.pos,
    mz.tol = mz.tol,
    rt.tol = 1000000,
    polarity = "positive"
  )

  idx.neg <- sample(1:length(all.mz.neg), size.neg) %>% sort()
  temp.mz.neg <- as.numeric(all.mz.neg[idx.neg])
  
  match.result.neg <- keggMatch(
    mz = temp.mz.neg,
    rt = rep(100, length(temp.mz.neg)),
    metabolite = kegg.adduct.neg,
    mz.tol = mz.tol,
    rt.tol = 1000000,
    polarity = "negative"
  )

  unique.id <-
    c(
      lapply(match.result.pos, function(x) {
        if (is.null(x))
          return(NULL)
        x$ID
      }) %>%
        unlist() %>%
        unique(),
      lapply(match.result.neg, function(x) {
        if (is.null(x))
          return(NULL)
        x$ID
      }) %>%
        unlist() %>%
        unique()
    ) %>%
    unique()
  
  ac.score <- getModule(
    node.id = unique.id,
    metabolic.network = kegg_network,
    threads = 3,
    progressbar = FALSE,
    max.reaction = 3,
    calculate.impact = FALSE
  )[[2]]

  acti.score[[i]] <- unlist(ac.score)
    }
             acti.score <- acti.score
           })
