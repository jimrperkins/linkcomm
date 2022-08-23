weighted_edges_to_symm_matrix <- function(weighted_edges) {
  unique_elements <- unique(c(weighted_edges[,1], weighted_edges[,2]))
  n_elements <- length(unique_elements)
  adj_matrix<-matrix(0, n_elements, n_elements)
  row.names(adj_matrix) <- colnames(adj_matrix) <- unique_elements
  adj_matrix[as.matrix(weighted_edges)[,1:2]] <- weighted_edges[,3]
  adj_matrix[as.matrix(weighted_edges)[,2:1]] <- weighted_edges[,3]
  return(adj_matrix)
}

weighted_edges_to_sparse_symm_matrix <- function(weighted_edges) {
  unique_nodes <- unique(c(weighted_edges[,1], weighted_edges[,2]))
  nodes_to_numbers <- data.frame(pair=unique_nodes, numid=seq_along(unique_nodes))
  weighted_edges$i <- nodes_to_numbers[match(weighted_edges[,1], nodes_to_numbers$pair), 2]
  weighted_edges$j <- nodes_to_numbers[match(weighted_edges[,2], nodes_to_numbers$pair), 2]
  # Sorting step needed to fill matrix above diagonal - i needs to be less than j
  weighted_edges[weighted_edges$i > weighted_edges$j,] <- weighted_edges[weighted_edges$i > weighted_edges$j, c(2,1,3,5,4)]
  sparse_matrix <- Matrix::sparseMatrix(i=weighted_edges$i, j=weighted_edges$j, x=weighted_edges[,3], 
    symmetric=TRUE, dimnames=list(nodes_to_numbers$pair, nodes_to_numbers$pair))
  #save(sparse_matrix, file="sparse_mx.RData")
  return(sparse_matrix)
}

calc_tanimoto_precomp <- function(ai, aj, crossprods) {
  a_i_j <- crossprods[ai, aj]
  a_i2 <- crossprods[ai, ai]
  a_j2 <- crossprods[aj, aj]
  a_i_j / (a_i2 + a_j2 - a_i_j)
}

sort_mat_by_row <- function(mat) {
  d <- data.table::as.data.table(mat)
  d[, row := .I]
  d <- data.table::melt(d, id.vars = "row") #wide to long format
  data.table::setkey(d, row, value) #sort
  d[, variable := rep(1:ncol(mat) ,nrow(mat))] #decreasing orde
  as.matrix(dcast(d, row ~ variable)[, row := NULL])
}

calc_edge_dens <- function(e,n) {(e*(e-n+1))/((n-2)*(n-1))}

get_link_densities <- function(hcedges, edges, clusnums, numcl, heights){
  tot_edges <- nrow(edges)
  tot_clust <- length(clusnums)
  height_subsets <- lapply(
    split(hcedges$merge, rep(1:tot_clust, clusnums)), matrix, ncol=2
  )
  all_comms <- list()
  i_glob <- 1
  pdens <- vector(length=tot_clust)
  all_scores <- vector(length=tot_clust)
  ldens <- 0
  for(r in 1:tot_clust) {
    height_subset <- height_subsets[[r]]
    hs_rows <- nrow(height_subset)
    ld_scores <- vector(length=hs_rows)
    comm_list <- list()
    for(i in 1:hs_rows) {
      neg_score <- 0
      pair <- height_subset[i,]
      #print(pair)
      if(pair[1] < 0 && pair[2] < 0) {
        comm <- edges[abs(pair),]
      } else if (pair[1] > 0 && pair[2] < 0) {
        comm <- rbind(edges[abs(pair[2]),],
                   all_comms[[pair[1]]])
        neg_score <- all_scores[pair[1]]
      } else if (pair[1] < 0 && pair[2] > 0) {
        comm <- rbind(edges[abs(pair[1]),],
                      all_comms[[pair[2]]])
        neg_score <- all_scores[pair[2]]
      } else if (pair[1] > 0 && pair[2] > 0) {
        comm <- rbind(all_comms[[pair[1]]],
                      all_comms[[pair[2]]])
        neg_score <- all_scores[pair[1]] + all_scores[pair[2]]
      }
      all_comms[[i_glob]] <- comm
      e <- nrow(comm) 
      n <- length(unique(as.vector(comm)))
      raw_ld_score <- calc_edge_dens(e,n)
      all_scores[i_glob] <- raw_ld_score
      ld_scores[i] <-  raw_ld_score - neg_score
      i_glob <- i_glob+1
    }
    ldens <- ldens + sum(ld_scores)
    pdens[r] <-  (2/tot_edges)*ldens
  }
  pdmax <- heights[which(pdens == max(pdens))][1]

  all_cl <- list()
  i_glob <- 1
  for(r in 1:which(pdens == max(pdens))[1]) {
    height_subset <- height_subsets[[r]]
    hs_rows <- nrow(height_subset)
    for(i in 1:hs_rows) {
      pair <- height_subset[i,]
      # print(pair)
      cl <- vector()
      if(pair[1] < 0 && pair[2] < 0) {
        cl <- pair
      } else if (pair[1] < 0 && pair[2] > 0) {
        cl <- c(pair[1], all_cl[[pair[2]]])
        all_cl[[pair[2]]] <- NA
      } else if (pair[1] > 0 && pair[2] < 0) {
        cl <- c(pair[2], all_cl[[pair[1]]])
        all_cl[[pair[1]]] <- NA
      } else if (pair[1] > 0 && pair[2] > 0) {
        cl <- c(all_cl[[pair[1]]], all_cl[[pair[2]]])
        all_cl[[pair[1]]] <- NA
        all_cl[[pair[2]]] <- NA
      }
      all_cl[[i_glob]] <- cl
      i_glob <- i_glob + 1
    }
  }
  # csize <- sum(length(unique(unlist(all_cl))))
  all_cl <- all_cl[! is.na(all_cl)]
  all_cl <- all_cl[sapply(all_cl, length) > 2]
  all_cl <- lapply(all_cl, function(x) sort(abs(x)))
  return(list(
    pdens = pdens,
    heights = heights[,1],
    pdmax = pdmax,
    csize = length(all_cl),
    linkcomm_clusters = all_cl)
  )
}

getLinkCommunities_fast <- function(network, hcmethod = "complete") 
  {
      all_nodes <- c(network[,1], network[,2])
      all_terms <- sort(unique(unlist(network[,c(1,2)])))
      assoc_term <- matrix(0, length(all_terms), length(all_terms), dimnames = list(all_terms, all_terms))
      assoc_term[as.matrix(network[,c(1,2)])] <- network[,3]
      assoc_term[as.matrix(network[,c(2,1)])] <- network[,3]
      tm <- tanimoto_edge_clustering(assoc_term, hcmethod=hcmethod)
      hc <- tm$hc
      edges <- tm$pair
      lc <- hclust_to_lc_obj(hc, edges)
      return(lc)
  }

tanimoto_edge_clustering <- function(sim_matrix, hcmethod) {
  k_ij_pairs <- apply(sim_matrix, 1, function(k_all_nodes_vec) {
    k_nodes_vec <- k_all_nodes_vec[k_all_nodes_vec > 0]
    if(length(k_nodes_vec) < 2) return()
    t(combn(sort(names(k_nodes_vec)), 2))
  })
  all_linked_pairs <- do.call(rbind, k_ij_pairs)
  unique_pairs <- unique(all_linked_pairs)

#  if(! is.null(opt$rows_per_submatrix)) {
#    crossprods_matrix <- crossprod_sections(sim_matrix, opt$rows_per_submatrix)
#  } else {
    crossprods_matrix <- crossprod(sim_matrix, sim_matrix)
#  }

  tanimoto_unique_pairs <- apply(unique_pairs, 1, function(x) calc_tanimoto_precomp(x[1], x[2], crossprods_matrix))
  tanimoto_matrix <- weighted_edges_to_symm_matrix(cbind(unique_pairs, tanimoto_unique_pairs))

  lengths_kij <- sapply(k_ij_pairs, length)/2
  k_names_vec <- rep(names(k_ij_pairs), lengths_kij)
  k_ij_pairs_df <- do.call("rbind", k_ij_pairs)
  col1 <- cbind(k_names_vec, k_ij_pairs_df[,1])
  col2 <- cbind(k_names_vec, k_ij_pairs_df[,2])
  col1 <- sort_mat_by_row(col1)
  col2 <- sort_mat_by_row(col2)
  col1 <- paste(col1[,1], col1[,2], sep="_")
  col2 <- paste(col2[,1], col2[,2], sep="_")
  k_pairs <- cbind(col1, col2)
  tanimoto_pairs <- cbind(k_pairs,
    apply(k_ij_pairs_df, 1, function(x) tanimoto_matrix[x[1],x[2]])
  )
  tanimoto_pairs <- data.frame(tanimoto_pairs)
  tanimoto_pairs[,3] <- as.numeric(tanimoto_pairs[,3])
  tanimoto_pairs[,3] <- round(tanimoto_pairs[,3], digits = 6)

  tanimoto_sparse <- weighted_edges_to_sparse_symm_matrix(tanimoto_pairs)

  pairs <- as.data.frame(t(as.data.frame(strsplit(row.names(tanimoto_sparse),"_"))))
  hc <- sparseAHC::sparseAHC(tanimoto_sparse, hcmethod, TRUE)
  return(list(hc=hc, pairs=pairs))
}


hclust_to_lc_obj <- function(hcedges, el) {
  verbose <- TRUE
  intel <- linkcomm::integer.edgelist(el) # Edges with numerical node IDs.
  edges <- intel$edges
  node.names <- names(intel$nodes)
  hcedges$order <- rev(hcedges$order)
  hh <- unique(round(hcedges$height, digits = 5)) # Round to 5 digits to prevent numerical instability affecting community formation.
  countClusters <- function(x,ht){return(length(which(ht==x)))}

  clusnums <- sapply(hh, countClusters, ht = round(hcedges$height, digits = 5)) # Number of clusters at each height.
  numcl <- length(clusnums)
  len <- nrow(el) # Number of edges.
  nnodes <- length(unique(c(as.character(el[,1]),as.character(el[,2])))) # Number of nodes.
  # Make this less horrid
  directed <- FALSE
  bipartite <- FALSE

  ldlist <- get_link_densities(hcedges, edges, clusnums, numcl, hh)
  pdens <- c(0,ldlist$pdens)
  heights <- c(0,hh)
  pdmax <- ldlist$pdmax
  csize <- ldlist$csize
  clus <- ldlist$linkcomm_clusters

  if(csize == 0){
    stop("\nno clusters were found in this network; maybe try a larger network\n")
  }

  if(verbose){
    cat("\n   Maximum partition density = ",max(pdens),"\n")
  }

  ecn <- data.frame()
  ee <- data.frame()
  lclus <- length(clus)
  for(i in 1:lclus){
    if(verbose){
      mes<-paste(c("   Finishing up...2/4... ",floor((i/lclus)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    ee <- rbind(ee,cbind(el[clus[[i]],],i))
    nodes <- node.names[unique(c(edges[clus[[i]],]))]
    both <- cbind(nodes,rep(i,length(nodes)))
    ecn <- rbind(ecn,both)
  }
  colnames(ecn) <- c("node","cluster")
  colnames(ee) <- c("node1","node2","cluster")


  # Extract the node-size of each edge cluster and order largest to smallest.
  ss <- NULL
  unn <- unique(ecn[,2])
  lun <- length(unn)
  for(i in 1:length(unn)){
    if(verbose){
      mes<-paste(c("   Finishing up...3/4... ",floor((i/lun)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    ss[i] <- length(which(ecn[,2]==unn[i]))
  }
  names(ss) <- unn
  ss <- sort(ss,decreasing=T)

  # Extract the number of edge clusters that each node belongs to.
  unn <- unique(ecn[,1])

  iecn <- as.integer(as.factor(ecn[,1]))
  iunn <- unique(iecn)
  lunn <- length(iunn)
  nrows <- nrow(ecn)

  oo <- rep(0,lunn)

  oo <- .C("getNumClusters", as.integer(iunn), as.integer(iecn), counts = as.integer(oo), as.integer(lunn), as.integer(nrows), as.logical(verbose))$counts

  names(oo) <- unn

  if(verbose){cat("\n")}

  pdplot <- cbind(heights,pdens)

  # Add nodeclusters of size 0.
  missnames <- setdiff(node.names,names(oo))
  m <- rep(0,length(missnames))
  names(m) <- missnames
  oo <- append(oo,m)

  all <- list()

  all$numbers <- c(len,nnodes,length(clus)) # Number of edges, nodes, and clusters.
  all$hclust <- hcedges # Return the 'hclust' object. To plot the dendrogram: 'plot(lcobj$hclust,hang=-1)'
  all$pdmax <- pdmax # Partition density maximum height.
  all$pdens <- pdplot # Add data for plotting Partition Density as a function of dendrogram height.
  all$nodeclusters <- ecn # n*2 character matrix of node names and the cluster ID they belong to.
  all$clusters <- clus # Clusters of edge IDs arranged as a list of lists.
  all$edges <- ee # Edges and the clusters they belong to, arranged so we can easily put them into an edge attribute file for Cytoscape.
  all$numclusters <- sort(oo,decreasing=TRUE) # The number of clusters that each node belongs to (named vector where the names are node names).
  all$clustsizes <- ss # Cluster sizes sorted largest to smallest (named vector where names are cluster IDs).
  all$igraph <- igraph::graph.edgelist(as.matrix(el), directed = directed) # igraph graph.
  all$edgelist <- el # Edge list.
  all$directed <- directed # Logical indicating if graph is directed or not.
  all$bipartite <- bipartite # Logical indicating if graph is bipartite or not.

  class(all) <- "linkcomm"
  return(all)
}


