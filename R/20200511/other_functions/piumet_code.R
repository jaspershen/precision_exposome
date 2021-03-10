###----------------------------------------------------------------------------
readPIUMet <-
  function(path = ".",
           marker_table,
           text = TRUE,
           layout = "kk",
           size_range = c(3, 8),
           width_range = c(0.2, 0.5),
           marker_name = NULL) {
    
    annotation_result <-
      read.table(
        file.path(
          path,
          "peaks_putative_metabolites_w10.0_b2.0_mu0.0005_R1.txt"
        ),
        sep = "\t",
        header = TRUE
      )
    
    if (nrow(annotation_result) == 0) {
      cat("No result.\n")
      return(NULL)
    }
    
    annotation_result <-
      annotation_result %>%
      dplyr::mutate(mz =
                      stringr::str_replace(mz.peak, "m/z=", "") %>%
                      as.numeric() %>%
                      round(4)) %>%
      dplyr::mutate(mz2 = as.character(mz))
    
    if (!is.null(marker_name)) {
      marker <-
        variable_info %>%
        dplyr::filter(name %in% marker_name) %>%
        dplyr::mutate(polarity = case_when(
          stringr::str_detect(name, "POS") ~ "positive",
          stringr::str_detect(name, "NEG") ~ "negative"
        )) %>%
        dplyr::mutate(
          mz2 = case_when(
            polarity == "positive" ~ as.character(round(mz, 4) - 1),
            polarity == "negative" ~ as.character(round(mz, 4) + 1)
          )
        )
    } else{
      marker_table <-
        marker_table %>% 
        dplyr::mutate(
          mz2 = case_when(
            polarity == "positive" ~ as.character(round(mz, 4) - 1),
            polarity == "negative" ~ as.character(round(mz, 4) + 1)
          )
        )
    }
    
    
    
    annotation_result <-
      annotation_result %>%
      dplyr::left_join(marker_table, by = "mz2") %>%
      dplyr::select(-c(mz.y)) %>%
      dplyr::rename(mz = mz.x)
    
    edge_attr <-
      read.table(
        file.path(path, "result_edge_frequency_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(edge = V1)
    
    edge_data <-
      read.table(
        file.path(path, "result_union_net_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(from = V1, to = V2) %>%
      dplyr::mutate(edge = paste(from, "(pp)", to, sep = " ")) %>%
      dplyr::left_join(edge_attr, by = "edge")
    
    node_data <-
      read.table(
        file.path(path, "result_node_frequency_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(node = V1,
                    node_class = V3,
                    HMDB_ID = V4) %>%
      dplyr::left_join(annotation_result[, c("name", "Metabolite.Name", "super.class")],
                       by = c("node" = "Metabolite.Name"))
    
    node <-
      node_data$node %>%
      stringr::str_replace("m/z=", "") %>%
      as.numeric() %>%
      round(4) %>%
      as.character()
    
    node <- marker_table$name[match(node, marker_table$mz2)]
    
    rename <-
      data.frame(name1 = node_data$node[!is.na(node)],
                 name2 = node[!is.na(node)])
    
    node_data$node <-
      sapply(node_data$node, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    edge_data$from <-
      sapply(edge_data$from, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    edge_data$to <-
      sapply(edge_data$to, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    
    edge_data$edge <-
      paste(edge_data$from, "(pp)", edge_data$to, sep = " ")
    
    node_data <-
      node_data %>%
      dplyr::select(-name) %>%
      dplyr::distinct()
    
    node_data$node_class[grep("Metabolite", node_data$node_class)] <-
      "Metabolite"
    
    graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(Degree = tidygraph::centrality_degree(mode = 'all'))
    
    fill <-
      c(
        "m/z Peak" = ggsci::pal_aaas()(10)[9],
        "Metabolite" = ggsci::pal_aaas()(10)[1],
        "Protein" = ggsci::pal_aaas(alpha = 0.1)(10)[2]
        # "Protein" = "tomato"
      )
    
    col <-
      c(
        "m/z Peak" = ggsci::pal_aaas()(10)[9],
        "Metabolite" = ggsci::pal_aaas()(10)[1],
        "Protein" = ggsci::pal_aaas()(10)[2]
      )
    
    shape = c("m/z Peak" = 24,
              "Metabolite" = 22,
              # "Metabolite_others" = 22,
              "Protein" = 21)
    
    require(ggraph)
    if (text) {
      plot <-
        ggraph(graph,
               layout = layout) +
        geom_edge_link(
          aes(edge_width = V3),
          alpha = 1,
          color = "black",
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(
            size = V2,
            fill = node_class,
            shape = node_class
          ),
          alpha = 1,
          show.legend = TRUE
        ) +
        scale_shape_manual(values = shape) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggraph::geom_node_text(aes(label = node,
                                   color = node_class),
                               repel = TRUE,
                               size = 3) +
        ggraph::scale_edge_width(range = width_range) +
        scale_size_continuous(range = size_range) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = col) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
    } else{
      plot <-
        ggraph(graph,
               layout = layout) +
        geom_edge_link(
          aes(edge_width = V3),
          alpha = 1,
          color = "black",
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(
            size = V2,
            fill = node_class,
            shape = node_class
          ),
          alpha = 1,
          show.legend = TRUE
        ) +
        scale_shape_manual(values = shape) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggraph::scale_edge_width(range = width_range) +
        scale_size_continuous(range = size_range) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = col) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
    }
    
    
    
    output_path <- file.path(path, "Result")
    dir.create(output_path)
    save(edge_data, file = file.path(output_path, "edge_data"))
    save(node_data, file = file.path(output_path, "node_data"))
    save(graph, file = file.path(output_path, "graph"))
    save(annotation_result, file = file.path(output_path, "annotation_result"))
    
    ggsave(
      plot,
      filename = file.path(output_path, "graph_plog.pdf"),
      width = 7,
      height = 7
    )
    
    plot
  }