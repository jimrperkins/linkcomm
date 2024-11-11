
sort_mat_by_row <- function(mat) {
  d <- data.table::as.data.table(mat)
  d[, row := .I]
  d <- data.table::melt(d, id.vars = "row") #wide to long format
  data.table::setkey(d, row, value) #sort
  d[, variable := rep(1:ncol(mat) ,nrow(mat))] #decreasing orde
  as.matrix(dcast(d, row ~ variable)[, row := NULL])
}


weighted_edges_to_symm_matrix2 <- function(weighted_edges) {
  # unique_elements <- unique(c(weighted_edges[,1], weighted_edges[,2]))
  unique_elements <- unique(as.character(t(weighted_edges[,1:2])))
  n_elements <- length(unique_elements)
  adj_matrix<-matrix(0, n_elements, n_elements)
  row.names(adj_matrix) <- colnames(adj_matrix) <- unique_elements
  adj_matrix[as.matrix(weighted_edges)[,1:2]] <- weighted_edges[,3]
  adj_matrix[as.matrix(weighted_edges)[,2:1]] <- weighted_edges[,3]
  return(adj_matrix)
}

calc_tanimoto_precomp2 <- function(aijs) {  
  node_n <- ncol(aijs)
  ai2s <- matrix(diag(aijs), nrow=node_n, ncol=node_n, byrow=TRUE)
  aj2s <- matrix(diag(aijs), nrow=node_n, ncol=node_n)
  aijs / (ai2s + aj2s - aijs)
}

getLinkCommunities_fast <- function(network, hcmethod = "average", verbose = TRUE, plot = TRUE, diagnorm = FALSE, minedges = 0) {
print("start")
      if(ncol(network) > 2) {
        network[,1] <- as.character(network[,1])
        network[,2] <- as.character(network[,2])
        assoc_term <- weighted_edges_to_symm_matrix2(network)
        unconnected_edges <- which(rowSums(assoc_term > 0) == 1) |> names() # Remove unconnected pairs
        network <- network[! (network[,1] %in% unconnected_edges & network[,2] %in% unconnected_edges), ]
        assoc_term <- weighted_edges_to_symm_matrix2(network)
        diag(assoc_term) <- 0 # is this really necessary?
        if(diagnorm) {
          print("in diagnorm")
          # diag(assoc_term) <- NA
          diag(assoc_term) <- apply(assoc_term, 1, function(x) mean(x[x>0]))
        } else {
          diag(assoc_term) <- 1
        }
        print("head assoc")
        print(head(assoc_term))
        orig_el <- integer.edgelist(network[,1:2])
        row.names(assoc_term) <- colnames(assoc_term) <-  orig_el$nodes
        tm <- tanimoto_edge_clustering(assoc_term, orig_el, hcmethod=hcmethod, verbose=verbose)
      } else if (ncol(network) == 2) {
        print("jaccard area")
        network[,1] <- as.character(network[,1])
        network[,2] <- as.character(network[,2])
        unique_nodes <- gtools::mixedsort(unique(c(network[,1], network[,2])))

        assoc_term <- matrix(0, length(unique_nodes), length(unique_nodes), 
          dimnames= list(unique_nodes, unique_nodes))
        assoc_term[as.matrix(network)] <- 1
        assoc_term[as.matrix(network)[,2:1]] <- 1

        tm <- tanimoto_edge_clustering(assoc_term, hcmethod=hcmethod, verbose=verbose)
      } else {
          error("strange input network format")
      }

      hc <- tm$hc
      # write.table(data.frame(hc$merge,hc$height), file="hctree", row.names=FALSE, col.names=FALSE, sep = "\t")
      edges <- tm$pair
      # write.table(edges, file="edges", row.names=FALSE, col.names=FALSE, sep="_")
      lc <- hclust_to_lc_obj(hc, network[,1:2], verbose=verbose, plot=plot, minedges=minedges)
      
      return(lc)
}

tanimoto_edge_clustering <- function(sim_matrix, el, hcmethod, verbose=verbose) {
  print("Start edge clust")
  print(proc.time())
  print(pryr::mem_used())
  crossprods_matrix <- crossprod(sim_matrix, sim_matrix)
  tanimoto_matrix <- calc_tanimoto_precomp2(crossprods_matrix)
  diag(sim_matrix) <- 0
  print("Precalculation of Tanimoto Values Done")
  print(proc.time())
  print(pryr::mem_used())

  kij_list <- apply(sim_matrix, 1, function(x) { if(sum(x>0) < 2) return(); combn(which(x > 0), 2) |> t()})
  tanim_vals <- lapply(kij_list, function(p) tanimoto_matrix[p])

  ij_df <- do.call(rbind, kij_list)
  k_vector <- rep(1:length(kij_list), times=sapply(kij_list, length)/2)
  print("About to sort")
  print(proc.time())
  print(pryr::mem_used())

  m1 <- sort_mat_by_row(data.frame(k_vector, ij_df[,1]))
  m2 <- sort_mat_by_row(data.frame(k_vector, ij_df[,2]))
  sorted_edges <- sort_mat_by_row(el$edges)
    print("Sorted")
  print(proc.time())
  print(pryr::mem_used())

  sparse_x_vals <- match(do.call(paste, data.frame(m1)),
    do.call(paste, data.frame(sorted_edges)))

  sparse_y_vals <- match(do.call(paste, data.frame(m2)),
    do.call(paste, data.frame(sorted_edges)))
    print("Matched")
  print(proc.time())
  edge_n <- nrow(el$edges)
  #  save(sparse_x_vals, sparse_y_vals, tanim_vals, edge_n, file="test.RData")
  tanimoto_sparse <- Matrix::Matrix(0, nrow = edge_n ,ncol = edge_n,sparse = T)
  # tanimoto_sparse <- matrix(0, nrow = edge_n, ncol = edge_n)
  tanimoto_sparse[cbind(sparse_x_vals, sparse_y_vals)] <- unlist(tanim_vals)

  # dij3 <- with(df, structure(Distance,
  #                            Size = length(nams),
  #                            Labels = nams,
  #                            Diag = FALSE,
  #                            Upper = FALSE,
  #                            method = "user",
  #                            class = "dist"))

  tanimoto_sparse[cbind(sparse_y_vals, sparse_x_vals)] <- unlist(tanim_vals)
  # save(list = ls(all.names = TRUE), file = "environment.RData")

  if(verbose) {
     tanimoto_full <- as.matrix(tanimoto_sparse)
     names_matrix <- apply(el$edges, 2, function(x) names(el$nodes)[x])
     names_matrix <- apply(names_matrix, 1, sort)
     tanimoto_names <- apply(names_matrix, 2, paste0, collapse="_")
     tanimoto_fromto <- data.frame(tanimoto_names[sparse_x_vals], tanimoto_names[sparse_y_vals])
     tanimoto_fromto <-  apply(tanimoto_fromto, 1, sort) |> t()
     write.table(data.frame(tanimoto_fromto, unlist(tanim_vals)), quote=FALSE,
      file="edge_scores.txt", sep="\t", row.names=FALSE, col.names=FALSE)
  }

    print("Sparse matrix made")
  print(proc.time())
  print(pryr::mem_used())

   #  hc <- fastcluster::hclust(1- as.dist(tanimoto_sparse), method="single")
    print("After Fast Cluster")
  print(proc.time())
  print(pryr::mem_used())
   hc <- sparseAHC::sparseAHC(tanimoto_sparse, linkage=hcmethod)
   # The other version
   #hc <- hclust(1-as.dist(tanimoto_full), method=hcmethod)

    print("Sparse cluster")
  print(proc.time())
  print(pryr::mem_used())
  pairs <- data.frame(sorted_edges)
  colnames(pairs) <- c("V1","V2")
  print("At the end of edge clust")
  print(proc.time())
  # q()
  return(list(hc=hc, pairs=pairs))
}





hclust_to_lc_obj <- function(hcedges, el, verbose, plot, minedges) {
  print("Start comm detect")
  print(proc.time())

  #  hcedges <- tm$hcedges
  #  el <- tm$pair
  # write.table(cbind(hcedges$merge, hcedges$height), file=paste0(fnam, "_hclust.txt"))

  intel <- linkcomm::integer.edgelist(el) # Edges with numerical node IDs.
  edges <- intel$edges
  node.names <- names(intel$nodes)
  #hcedges$order <- rev(hcedges$order)
  # hh <- unique(round(hcedges$height, digits = 5)) # Round to 5 digits to prevent numerical instability affecting community formation.
  hh <- unique(hcedges$height)
  countClusters <- function(x,ht){return(length(which(ht==x)))}

  # clusnums <- sapply(hh, countClusters, ht = round(hcedges$height, digits = 5)) # Number of clusters at each height.
  clusnums <- sapply(hh, countClusters, ht = hcedges$height) # Number of clusters at each height.
  numcl <- length(clusnums)
  len <- nrow(el) # Number of edges.
  nnodes <- length(unique(c(as.character(el[,1]),as.character(el[,2])))) # Number of nodes.
  # Make this less horrid
  directed <- FALSE
  bipartite <- FALSE



  ldlist <- get_link_densities(hcedges, edges, clusnums, numcl, hh, minedges)
  #write.table(ldlist$pdens, file='densities', sep="\t", row.names = FALSE, col.names = FALSE)


  pdens <- c(0,ldlist$pdens)
  heights <- c(0,hh)
  pdmax <- ldlist$pdmax
  csize <- ldlist$csize
  clus <- ldlist$linkcomm_clusters

  print("Got link density stuff")
  print(proc.time())

  if(csize == 0){
    stop("\nno clusters were found in this network; maybe try a larger network\n")
  }

  if(verbose){
    cat("\n   Maximum partition density = ",max(pdens),"\n")
  }

  ecn <- data.frame()
  ee <- data.frame()
  lclus <- length(clus)
  print("finishing up 2")
  for(i in 1:lclus){
    if(verbose){
      mes<-paste(c("   Finishing up...2/4... ",floor((i/lclus)*100),"%"),collapse="")
      # cat(mes,"\r")
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
  print("finishing up 3")
  for(i in 1:length(unn)){
    if(verbose){
      mes<-paste(c("   Finishing up...3/4... ",floor((i/lun)*100),"%"),collapse="")
      # cat(mes,"\r")
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
    oo_orig <- oo

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
  all$edgelist <- as.matrix(el) # Edge list.
  all$directed <- directed # Logical indicating if graph is directed or not.
  all$bipartite <- bipartite # Logical indicating if graph is bipartite or not.

  class(all) <- "linkcomm"

  if(plot){
    if(verbose){
      cat("   Plotting...\n")
      }
    if(len < 1500){ # Will be slow to plot dendrograms for large networks.
      if(len < 500){
        all <- plot(all, type="summary", droptrivial = removetrivial, verbose = verbose)
      }else{ # Slow to reverse order of large dendrograms.
        all <- plot(all, type="summary", right = FALSE, droptrivial = removetrivial, verbose = verbose)
        }
    }else if(len <= edglim){
      oldpar <- par(no.readonly = TRUE)
      par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
      plot(hcedges,hang=-1,labels=FALSE)
      abline(pdmax,0,col='red',lty=2)
      plot(pdens,heights,type='n',xlab='Partition Density',ylab='Height')
      lines(pdens,heights,col='blue',lwd=2)
      abline(pdmax,0,col='red',lty=2)
      par(oldpar)
    }else{
      plot(heights,pdens,type='n',xlab='Height',ylab='Partition Density')
      lines(heights,pdens,col='blue',lwd=2)
      abline(v = pdmax,col='red',lwd=2)
      }
    }
  print("At the end of link density stuff")
  print(proc.time())
  return(all)
}

calc_edge_dens <- function(e,n) {(e*(e-n+1))/((n-2)*(n-1))}

get_link_densities <- function(hcedges, edges, clusnums, numcl, heights, minedges){
  tot_edges <- nrow(edges)
  tot_nodes <- max(edges)
  tot_clust <- length(clusnums)

  height_subsets <- lapply(
    split(hcedges$merge, rep(1:tot_clust, clusnums)), matrix, ncol=2
  )
  all_comms <- list()
  i_glob <- 1
  pdens <- vector(length=tot_clust)
  all_scores <- vector(length=tot_clust)
  ldens <- 0
  save(list = ls(all.names = TRUE), file = "environment.RData")

  # Here we are going through each cluster branching height in the tree,
  # looking to calculate a partition density (pdens) value at each height.
  # The pdens value is calculated based on the numbers of nodes and edges among the members
  # of each cluster at a given height.
  for(r in 1:tot_clust) {
    #cat(r)
    # cat(c(i_glob, r, "\t"), file = paste0(fnam, "_pdens.txt"), append=TRUE)
    print(paste("--", hcedges[r], heights[r]))
    height_subset <- height_subsets[[r]]
    # This is to catch a weird bug in the sparse clustering method that means sometimes the merge table
    # contains decimal values (this shouldn't happen as the merge table refers to the individual items
    # or clusters to be merged.)
    #if(any(abs(height_subset) < 1)) {
    if(any(height_subset %% 1 != 0) & any(abs(height_subset) < 1)) {
      print("The merge table contains decimal values, so we assume we are at the top of the tree")
      final_ld_score <- calc_edge_dens(tot_edges, tot_nodes)
      pdens[r] <-  (2/tot_edges)*final_ld_score
      next #
    }
    # cat(c("height:", unlist(height_subset, "\n")), file = paste0(fnam, "_pdens.txt"), append=TRUE)
    hs_rows <- nrow(height_subset)
    ld_scores <- vector(length=hs_rows)
    # For each merge at that cluster branching height, we get a pair of merged branches
    # We then work out if the branches were singletons or clusters being merged, based on
        # the sign (negative singleton, positive cluster). We then add either new singletons to the
    # comm variable using edges, or clusters, or both, depending on the merge.
    # Thus, the comm variable include all pairs of nodes from the original graph that correspond
    # To that cluster.
    # We then use

    # print()
    total_n <- 0
    for(i in 1:hs_rows) {
      #cat(i, "of", hs_rows, "\n")
      neg_score <- 0
      pair <- height_subset[i,]
      if(pair[1] < 0 && pair[2] < 0) {
        comm <- edges[abs(pair),]
      } else if (pair[1] > 0 && pair[2] < 0) {
        comm <- rbind(edges[abs(pair[2]),], all_comms[[pair[1]]])
        neg_score <- all_scores[pair[1]]
      } else if (pair[1] < 0 && pair[2] > 0) {
        comm <- rbind(edges[abs(pair[1]),], all_comms[[pair[2]]])
        neg_score <- all_scores[pair[2]]
      } else if (pair[1] > 0 && pair[2] > 0) {
        comm <- rbind(all_comms[[pair[1]]], all_comms[[pair[2]]])
        neg_score <- all_scores[pair[1]] + all_scores[pair[2]]
      }
      all_comms[[i_glob]] <- comm
      e <- nrow(comm)
      n <- length(unique(as.vector(comm)))
      total_n <- total_n + n
      raw_ld_score <- calc_edge_dens(e,n)
      all_scores[i_glob] <- raw_ld_score
      ld_scores[i] <-  raw_ld_score - neg_score
            i_glob <- i_glob+1
    }
    ldens <- ldens + sum(ld_scores)
    pdens[r] <-  (2/tot_edges)*ldens
    # cat(pdens[r], "\n", file = paste0(fnam, "_pdens.txt"), append = TRUE )

  }
  # Find the height that corresponds to the maximum density
  pdmax <- heights[which(pdens == max(pdens))][1]

  print("FORMING COMMUNITIES")
  # print(pdens)
  # write.table(pdens, file="pdens", row.names=FALSE, col.names=FALSE, sep="\t")
  # write.table(max(pdens), file="maxpdens", row.names=FALSE, col.names=FALSE, sep="\t")
  all_cl <- list()
  i_glob <- 1
  for(r in 1:(which(pdens == max(pdens))[1])) {
    # print(r)
    height_subset <- height_subsets[[r]]
    hs_rows <- nrow(height_subset)
    for(i in 1:hs_rows) {
      pair <- height_subset[i,]
      #print("The pair is", pair)
      cat("the pair is ->","\t",pair)
if(-2 %in% unlist(pair)) #  print("test")
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
      # print(cl)
      all_cl[[i_glob]] <- cl
      i_glob <- i_glob + 1
    }
  }
  # csize <- sum(length(unique(unlist(all_cl))))
  all_cl <- all_cl[! is.na(all_cl)]
  all_cl <- all_cl[sapply(all_cl, length) >= minedges ] 
  all_cl <- lapply(all_cl, function(x) sort(abs(x)))
  return(list(
    pdens = pdens,
    heights = heights,
    pdmax = pdmax,
    csize = length(all_cl),
    linkcomm_clusters = all_cl)
  )
}


