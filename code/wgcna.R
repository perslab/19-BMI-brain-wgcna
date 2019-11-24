#' @title Single cell WGCNA modules
#' @author Jon Thompson rkm916 at ku dot dk
#' @importFrom Seurat 2.3.3

######################################################################
######################## INITIAL PACKAGES ############################
######################################################################

suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("optparse"))

######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  
  make_option("--raw.data_path", type="character",
              help = "Provide full path to gene * sample expression matrix in delimited file format (may be compressed). First column should be 'gene'."),
  make_option("--meta.data_path", type="character",
              help = "Provide full path to sample * variable metadata table in deliminted file format (may be compressed). First column should be sample names."),
  make_option("--project_dir", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir. [default %default]"),
  make_option("--scratch_dir", type="character", default="/scratch/tmp-wgcna/",
              help = "Directory for temporary files, [default %default]"),
  make_option("--data_prefix", type="character", default="sc_data",
              help = "Dataset prefix for the project, [default %default]"),
  make_option("--data_type", type="character", default="sc",
              help = "Expression data type, one of 'sc' or 'bulk', [default %default]"),
  make_option("--run_prefix", type="character", default="run1",
              help = "Run prefix to distinguish runs e.g. with different parameters, [default %default]"),
  make_option("--autosave", type="logical", default=T,
              help = "Autosave session images at regular intervals to resume later? [default %default]."),
  make_option("--resume", type="character", default=NULL,
              help = "Resume from a checkpoint? Must have same path and data_prefix as provided. Options are 'checkpoint_1' - '5' [default %default]"),
  make_option("--metadata_subset_col", type="character",
              help = "Specify a seurat@meta.data$... column to use for subsetting the data."),
  make_option("--metadata_corr_col", type="character", default='NULL',
              help = "Specify seurat_obj@meta.data$... column(s) for which to compute correlations with gene modules. Takes a character with a vector in single (double) quotes of seurat_obj@meta.data column names in double (single) quotes, without whitspace, e.g. 'nUMI' or 'c('Sex','Age')'. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution.  [default %default]"),
  make_option("--metadata_corr_filter_vals", type="character", default='NULL',
              help = "Specify one or more values within the seurat@meta.data$... column(s). Retain only modules which are significantly (anti-) correlated (at present, there is no threshold value for the correlation). Takes a character with a vector of meta.data column names without whitespace, e.g. 'Female' or 'c('fasted', 'HFD')'. Case-insensitive [default %default]."),
  make_option("--regress_out", type="character", default='c("nUMI", "percent.mito", "percent.ribo")',
              help="Provide arguments to Seurat's ScaleData function in the form of a vector in quotes, defaults to c('nUMI', 'percent.mito', 'percent.ribo') [default %default]"),
  make_option("--min.cells", type="integer", default=20L,
              help="What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, [default %default]."),
  make_option("--minCellClusterSize", type="integer", default=50L,
              help="What is the minimum number of cells in a subset to continue? Integer, [default %default]."),
  make_option("--map_genes_to_ensembl", type="logical", default=TRUE,
              help="Required if applying MAGMA gene set test [default %default]."),
  make_option("--genes_use", type="character", default="PCA",
              help="One of 'all', 'var.genes' or 'PCA' for genes with significant loading on at least one significant PC. 'All' is not recommended. [default %default]"), 
  make_option("--pca_genes", type="character", default="var.genes",
              help="'all' or 'var.genes'. 'All' is computationally expensive but allows for selecting genes based on PC loading p-values rather than magnitudes [default %default]"), 
  make_option("--n_genes_use", type="integer", default=5000L,
              help="If using PCA with var.genes or jackstrawnReplicate=0 and therefore selecting genes based on PC loadings, how many top loading genes to use [default %default]"), 
  make_option("--corFnc", type="character", default="bicor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function (https://en.wikipedia.org/wiki/Biweight_midcorrelation). [default %default]"), 
  make_option("--networkType", type="character", default = "signed",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); ''c('signed hybrid')'' (quoted vector) sets negative correlations to zero. [default %default]"),
  make_option("--hclustMethod", type="character", default="average",
              help = "Hierarchical clustering agglomeration method. One of 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC). See hclust() documentation for further information. [default %default]"),
  make_option("--minClusterSize", type="character", default="15L",
              help = "Minimum genes needed to form a module, or an initial cluster before the Partitioning Around Medoids-like step. WGCNA authors recommend decreasing the minimum cluster size when using higher settings of deepSplit. Takes a character with a vector, without whitespace, of integer values to try 'c(15,20,25)' [default %default]"),
  make_option("--deepSplit", type="character", default="3L",
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes a character with a vector, without whitespace, of integer values between 0-4 to try, e.g. 'c(1,2,3)', [default %default]"),
  make_option("--moduleMergeCutHeight", type="character", default="c(0.15)",
              help = "Cut-off level for 1-cor(eigen-gene, eigen-gene) for merging modules. Takes a character with a vector, without whitespace, of double values to try, e.g. 'c(0.1, 0.2)', [default %default]"),
  make_option("--pamStage", type="character", default="c(TRUE)",
              help = "cutreeHybrid. Perform additional Partition Around Medroids step? Users favoring specificity over sensitivity (that is, tight clusters with few mis-classifications) should set pamStage = FALSE, or at least specify the maximum allowable object-cluster dissimilarity to be zero (in which case objects are assigned to clusters only if the object-cluster dissimilarity is smaller than the radius of the cluster). Takes a character with a vector, without whitespace, of logicals to try, e.g. 'c(TRUE,FALSE)', [default %default]"),
  make_option("--kM_reassign", type="logical", default="TRUE",
              help = "Following hierarchical clustering, do additional k-means clustering? [default %default]"),
  make_option("--kM_signif_filter", type="logical", default="TRUE",
              help = "Do a t test to filter out insignificant genes from modules? If fuzzyModMembership=='kME', carries out a correlation t-test; if kIM, does a two-sample t-test between IntraModular and ExtraModular connectivities"),
  make_option("--jackstrawnReplicate", type="integer", default=500L,
              help = "Number of times to re-run PCA after permuting a small proportion of genes to perform empirical significance tests, i.e. the `JackStraw` procedure (see `pca_genes` above), [default %default]."),
  make_option("--TOMnReplicate", type="integer", default=100L,
              help = "Number of times to resample the dataset when finding the consensus TOM [default %default]"),
  make_option("--fuzzyModMembership", type="character", default="kIM",
              help="Which 'fuzzy' measure of gene membership in a module should be used? Options are 'kME' (correlation between gene expression and module PC1 expression) and 'kIM' (sum of edges between a gene and genes in the module, normalized by number of genes in the module; i.e. average distance in the TOM between a gene and genes in the module [default %default]."),  
  make_option("--PPI_filter", type="logical", default=T,
              help="Valiate gene modules using Protein-Protein Interactions?"),
  make_option("--data_organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus', [default %default]"),
  make_option("--list_genesets_path", type="character", default = NULL,
              help = "Path to a list of geneset vectors saved as a RDS or RData object. If map_genes_to_ensembl is TRUE, if necessary lists will be remapped to ensembl [default %default]"),
  make_option("--magma_gwas_dir", type="character", default = NULL,
              help = "MAGMA input GWAS data directory, e.g. '/projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/', '/projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/'. NULL skips the magma GWAS step [default %default]"),
  make_option("--gwas_filter_traits", type="character", default = NULL,
              help = "Filter out modules not significantly correlated with matching gwas studies within the magma_gwas_dir. Takes a character with a vector, without whitespace, of character names to match within the filename of the GWAS , e.g. ''c('body_BMI_Locke2015')'' or ''c('BMI', 'T1D', 'T2D')''. Case-insensitive. [default %default]"),
  make_option("--RAM_Gb_max", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)

######################################################################
################### GET COMMAND LINE OPTIONS #########################
######################################################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 

opt <- parse_args(OptionParser(option_list=option_list))

resume <- opt$resume

data_prefix <- opt$data_prefix 

scratch_dir <- opt$scratch_dir

run_prefix <- opt$run_prefix

# load saved image?
if (!is.null(resume)) {
  message(paste0("loading from ", resume))
  tryCatch({load(file=sprintf("%s%s_%s_%s_image.RData.gz", scratch_dir, data_prefix, run_prefix, resume))},
           error = function(x) {stop(paste0(resume, " session image file not found in ", scratch_dir))})
  
}

opt <- parse_args(OptionParser(option_list=option_list))

raw.data_path <- opt$raw.data_path 

meta.data_path <- opt$meta.data_path 

project_dir <- opt$project_dir


data_type = opt$data_type

autosave <- opt$autosave

metadata_subset_col <- opt$metadata_subset_col

metadata_corr_col <- opt$metadata_corr_col
if (!is.null(metadata_corr_col)) metadata_corr_col <- eval(parse(text=metadata_corr_col))

metadata_corr_filter_vals <- opt$metadata_corr_filter_vals
if (!is.null(metadata_corr_filter_vals)) metadata_corr_filter_vals <- eval(parse(text=metadata_corr_filter_vals))

regress_out <- opt$regress_out
if (!is.null(regress_out)) regress_out <- eval(parse(text=regress_out))

min.cells <- opt$min.cells

minCellClusterSize <- opt$minCellClusterSize

map_genes_to_ensembl <- opt$map_genes_to_ensembl

genes_use <- opt$genes_use

pca_genes <- opt$pca_genes

n_genes_use <- opt$n_genes_use

corFnc <- opt$corFnc 

networkType <- opt$networkType
if (grepl("hybrid", networkType)) {
  networkType <- eval(parse(text=opt$networkType))
}

hclustMethod <- opt$hclustMethod

minClusterSize <- eval(parse(text=opt$minClusterSize))

deepSplit <- eval(parse(text=opt$deepSplit))

moduleMergeCutHeight <- eval(parse(text=opt$moduleMergeCutHeight))

pamStage <- eval(parse(text=opt$pamStage))

kM_reassign <- opt$kM_reassign

kM_signif_filter <- opt$kM_signif_filter

jackstrawnReplicate <- opt$jackstrawnReplicate

TOMnReplicate <- opt$TOMnReplicate

fuzzyModMembership <- opt$fuzzyModMembership

PPI_filter <- opt$PPI_filter

list_genesets_path <- opt$list_genesets_path

magma_gwas_dir <- opt$magma_gwas_dir 

gwas_filter_traits <- opt$gwas_filter_traits
if (!is.null(gwas_filter_traits)) gwas_filter_traits <- eval(parse(text=gwas_filter_traits))

data_organism <- opt$data_organism

RAM_Gb_max <- opt$RAM_Gb_max

############################################################################################################################################################
############################################################## DEFINE FUNCTIONS ############################################################################
############################################################################################################################################################

FilterGenes <- function(seurat_obj_sub, min.cells) {
  
  if (min.cells > 0 & !is.null(dim(seurat_obj_sub@raw.data))) { 
    num.cells <- rowSums(seurat_obj_sub@raw.data > 0)
    genes.use <- names(x = num.cells[which(x = num.cells >= min.cells)])
    if (length(genes.use)>0) {
      seurat_obj_sub@raw.data <- seurat_obj_sub@raw.data[genes.use, ]
    } else {
      seurat_obj_sub <- NULL
    }
  }
  
  return(seurat_obj_sub)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

wrapJackStraw = function(seurat_obj_sub, n_cores, jackstrawnReplicate, pvalThreshold, n_genes_use=5000) {

  prop.freq <- max(0.016, round(4/length(seurat_obj_sub@var.genes),3)) # to ensure we have at least 3 samples so the algorithm works well
  # see https://github.com/satijalab/seurat/issues/5

  
  pcs.compute = ncol(seurat_obj_sub@dr$pca@gene.loadings)
  
  if (jackstrawnReplicate > 0) {
    
    seurat_obj_sub <- JackStraw(object = seurat_obj_sub,
                                num.pc = pcs.compute,
                                num.replicate = jackstrawnReplicate, 
                                display.progress = T,
                                do.par = T,
                                num.cores = n_cores,
                                prop.freq = prop.freq) # https://github.com/satijalab/seurat/issues/5
    

    pAll <- GetDimReduction(seurat_obj_sub, reduction.type = "pca", slot = "jackstraw")@emperical.p.value

    pAll <- pAll[, 1:pcs.compute, drop = FALSE]
    pAll <- as.data.frame(pAll)
    pAll$Contig <- rownames(x = pAll)
    pAll.l <- melt(data = pAll, id.vars = "Contig")
    colnames(x = pAll.l) <- c("Contig", "PC", "Value")
    
    score.df <- NULL
    
    for (i in (1:pcs.compute)) {
      pc.score <- suppressWarnings(prop.test( 
        x = c(length(x = which(x = pAll[, i] <= pvalThreshold)), floor(x = nrow(x = pAll) * pvalThreshold)),
        n = c(nrow(pAll), nrow(pAll)))$p.val)
      if (length(x = which(x = pAll[, i] <= pvalThreshold)) == 0) {
        pc.score <- 1
      }
      if (is.null(x = score.df)) {
        score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
      } else {
        score.df <- rbind(score.df, data.frame(PC = paste0("PC",i), Score = pc.score))
      }
    }
    
    PC_select_idx <- which(score.df$Score < pvalThreshold)
    
    if (nrow(pAll) == length(seurat_obj_sub@var.genes) | length(PC_select_idx)<5) {
      
      seurat_obj_sub <- ProjectPCA(seurat_obj_sub, 
                                   do.print = F, 
                                   pcs.print = NULL, 
                                   pcs.store = pcs.compute, 
                                   genes.print = NULL, 
                                   replace.pc=F, 
                                   do.center=T)
      
      loadings <- abs(seurat_obj_sub@dr$pca@gene.loadings.full[,PC_select_idx, drop=F])
      max_loadings <- apply(loadings, 1, function(x) max(x))
      names_genes_use <- names(max_loadings[order(max_loadings, decreasing = T)])[1:5000]
      
    } else if (nrow(pAll) > length(seurat_obj_sub@var.genes) & length(PC_select_idx)>=5) {
      
      pAll[,sapply(pAll, function(x) class(x)!="numeric")] <- NULL # remove the column of gene names
      row_min <- apply(pAll[,PC_select_idx, drop=F], MARGIN = 1, FUN = function(x) min(x))
      names_genes_use <- rownames(pAll)[row_min < pvalThreshold]
      
      if (length(names_genes_use) < 1000) names_genes_use <- rownames(pAll)[row_min < pvalThreshold*2]
    }
    
  } else if (jackstrawnReplicate == 0) {
    
    seurat_obj_sub <- ProjectPCA(seurat_obj_sub, 
                                 do.print = F, 
                                 pcs.print = NULL, 
                                 pcs.store = pcs.compute, 
                                 genes.print = NULL, 
                                 replace.pc=F, 
                                 do.center=T)
    
    loadings <- abs(seurat_obj_sub@dr$pca@gene.loadings.full[,1:(min(13,pcs.compute))])
    max_loadings <- apply(loadings, MARGIN=1, FUN=function(x) max(x))
    names_genes_use <- names(max_loadings[order(max_loadings, decreasing = T)])[1:min(n_genes_use, length(max_loadings))]
    
  }
  
  datExpr <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] %>% t() #%>% as.matrix() 
  colnames(datExpr) <- rownames(seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,])
  
  return(datExpr)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

seurat_to_datExpr = function(seurat_obj_sub, idx_genes_use) {
  datExpr <- seurat_obj_sub@scale.data[idx_genes_use,] %>% t()
  colnames(datExpr) <- rownames(seurat_obj_sub@scale.data[idx_genes_use,])
  return(datExpr)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

bootstrap <- function(datExpr,
                      nPermutations,
                      replace,
                      fraction,
                      randomSeed)
  #' @Usage: Resample a dataset.
  #' @args:
  #'       datExpr: Dataset with samples in the rows, is coerced to matrix
  #'       nPermutations: runs
  #'       replace: sample with replacement?
  #'       fraction: sample what fraction of the total each iteration?
  #'       randomSeed: initial random seed
  #' @return:
  #'       result: list of resampled datasets in matrix format. If replace = T, first item is the unpermuted dataset.
  #' @author: Jonatan Thompson jjt3f2188@gmail.com
  #' @date: 180222

{
  startRunIndex = 1
  endRunIndex = if (replace==T) nPermutations+1 else nPermutations 
  result = vector("list", length= if (replace==T) nPermutations+1 else nPermutations);
  nSamples = nrow(datExpr);
  nGenes = ncol(datExpr);
  
  try(
    for (run in startRunIndex:endRunIndex)
      
    {
      
      if (run == startRunIndex & replace == T) {
        useSamples = c(1:nSamples) # The first run just returns the full dataset
      } else if (run>startRunIndex | replace == F) 
      {  useSamples = sample(nSamples, as.integer(nSamples * fraction), replace = replace)
      } 
      
      
      samExpr = as.matrix(datExpr[useSamples, ]);
      
      result[[run]]$data <- samExpr
    })
  
  return(result)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

cutreeHybrid_for_vec <- function(comb, geneTree, dissTOM, maxPamDist, useMedoids) {
  # Utility function for more easily parallellising the cutreeHybrid function
  tree = cutreeHybrid(dendro = geneTree, 
                      cutHeight = NULL,
                      minClusterSize = comb[[1]], 
                      distM=as.matrix(dissTOM),
                      deepSplit=comb[[2]], 
                      pamStage=comb[[3]],
                      pamRespectsDendro=comb[[3]],
                      maxPamDist = maxPamDist,
                      useMedoids = useMedoids) 
  # Gandal et al 2018:  cutHeight = 0.999
  return(tree)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

mergeCloseModskIM = function(datExpr,
                             colors,
                             kIMs,
                             dissTOM = NULL,
                             moduleMergeCutHeight,
                             verbose=2,
                             cellType) {
  
  #' @usage: compute correlations between gene module cell embeddings. 
  #'         merge modules with a positive correlation > 1-moduleMergeCutHeight
  #'
  #' @args: datExpr: a cell * gene expression matrix
  #'       colors: vector module assignments with gene names
  #'       kIMs: generalised IntraModular connectivity, dataframe of gene row and module column. Should not include the grey module
  #'       dissTOM: only needed if iterate=T, in order to recompute kIMs
  #'       cellModEmbed_mat: embedding matrix of cell rows and module columbs. Passing the matrix saves time; otherwise will be computed
  #'       moduleMergeCutHeight: max distance (1-correlation coefficient) for merging modules
  #'       verbose: verbosity of WGCNA::cor function, 1-5
  #
  #' @returns: list with entries colors=colors_merged, kIMs = kIMs_merged
  
  
  corr_clust <- character(length = length(unique(colors)))
  colors_original <- colors
  # Remove 'grey' (unassigned) module 
  #if (any(colnames(kIMs) == "grey")) kIMs[["grey"]] <- NULL
  
  if (length(unique(colors)) > 1) { # if there is more than one non-grey module
    while (TRUE) {
      
      message(paste0(cellType, ": computing cell-module embedding matrix"))
      cellModEmbed_mat <- cellModEmbed(datExpr=datExpr, 
                                       colors=colors, 
                                       latentGeneType = "IM",
                                       cellType=NULL,
                                       kMs=kIMs)
      
      #cellModEmbed_mat <- cellModEmbed_mat[,-grep("^grey$", colnames(cellModEmbed_mat))]
      # Cluster modules using the Pearson correlation between module cell embeddings
      mod_corr <- WGCNA::cor(x=cellModEmbed_mat, method=c("pearson"), verbose=verbose)
      mod_corr[mod_corr<0] <- 0 #only keep positive correlations
      corr_dist <- as.dist(m = (1-mod_corr), diag = F, upper = F) # convert to (1-corr) distance matrix
      corr_dendro <- hclust(d = corr_dist, method = "average")
      
      # use a simple cut to determine clusters
      corr_clust = cutreeStatic(dendro = corr_dendro, 
                                cutHeight = moduleMergeCutHeight, 
                                minSize=1)
      names(corr_clust) =  corr_dendro$labels
      
      # At this point if corr_clust has as many unique modules as the original partition, the while loop will exit
      if(length(unique(corr_clust)) == length(unique(colors))) {  
        message(paste0(cellType, " done"))
        break
      }
      
      # if not we merge modules
      n_merge <-  length(unique(colors))  - length(unique(corr_clust)) 
      if (verbose>0) message(paste0(cellType, ": ", n_merge, " modules to be merged with others"))
      
      merge_idx <- sapply(unique(corr_clust), function(x) sum(corr_clust==x)>1)
      
      # merge modules
      for (clust in unique(corr_clust)[merge_idx]) { # loop over numeric cluster assignments 
        new_color = names(corr_clust)[corr_clust==clust][1] # use the colors name of the first module in the cluster
        for (color in names(corr_clust[corr_clust==clust])) { # loop over all module colors in the cluster
          colors[colors==color] <- new_color # give them all the new colour
        }
      }
      
      # compute merged module kIMs 
      message(paste0("Computing merged module kIMs for ", cellType))
      # Compute new kIMs
      kIMs <- kIM_eachMod_norm(dissTOM = dissTOM, 
                               colors = colors,
                               verbose = 1,
                               excludeGrey = F)
      
    }
    names(colors) = names(colors_original)
  }
  return(list(colors=colors, kIMs = kIMs))
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parMatchColors <- function(list_colors) {
  # Inherits n_combs from the parent environment, same for every subsetName, no need to pass as an argument
  
  if (length(list_colors)>1) {
    # Match colors between iterations by reference to the set of colors with the highest number of unique colors
    max_colors_idx = which.max(lapply(list_colors, function(x) length(table(x))))
    
    diffParams_colors = matrix(0, nrow=length(list_colors[[1]]), ncol=n_combs)
    diffParams_colors[, max_colors_idx] = list_colors[[max_colors_idx]] # use the set of colors with the highest number of modules as 'reference'
    
    for (j in (1:n_combs)[-max_colors_idx]) {
      if (length(table(list_colors[[j]])) > 1) {
        diffParams_colors[,j] <- matchLabels(source = list_colors[[j]], # traverse gene colour assignments
                                             reference = diffParams_colors[,max_colors_idx],
                                             pThreshold = 5e-2,
                                             na.rm = TRUE,
                                             extraLabels = standardColors())#paste0("Extra_", 0:500)) # labels for modules in source that cannot be matched to any modules in reference
        
      } else if (length(table(list_colors[[j]])) == 1) {
        diffParams_colors[,j] <- "grey"
      }
    }
    # Split the matrix of matched colors back into a list 
    list_colors_matched <- split(diffParams_colors, rep(1:ncol(diffParams_colors), each = nrow(diffParams_colors))) 
    
  } else if (length(list_colors)==1)  {
    diffParams_colors <- as.matrix(list_colors[[1]], nrow=length(list_colors[[1]]), ncol=1)
    list_colors_matched <- list_colors
  } 
  
  
  
  return(list_colors_matched)
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

pkMs_fnc = function(kMs, colors) {
  pkMs <- NULL
  if (!is.null(kMs) & length(unique(colors))>1) {
    for (i in 1:length(colors)) {
      pkMs[i] <- if (colors[i]=="grey") 0 else kMs[colors[i]] [i,]
    }
    names(pkMs) <- names(colors)
  } else {
    pkMs = rep(0, length(colors))
    names(pkMs) <- names(colors)
  }
  return(pkMs)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kM_reassign_fnc = function(cols,
                           fuzzyModMembership,
                           dissTOM = NULL,
                           datExpr = NULL,
                           verbose=2,
                           corFnc=NULL,
                           max_iter=3,
                           reassign_threshold=1.2,
                           cellType) {
  
  #' @Usage: iteratively reassign genes to modules based on kIM / kME until the number to reassign <= stop_condition
  #' @Args: 
  #'   fuzzyModMembership: "kME" or "kIM"
  #'   dissTOM: list of TOM distance matrices. Only required if fuzzyModMembership = "kIM"
  #'   datExpr: cell x gene expression matrix. Only required if fuzzyModMembership = "kME"
  #'   cols: vector of module assignments with gene names
  #'   corFnc: WGCNA correlation function; only needed if fuzzyModMembership = "kME"
  #'   verbose: 0-5, controls message output
  #'   max_iter: integer; max iterations.
  #'   reassign_threshold: threshold ratio of competing to current module kME for reassignment to occur 
  #' @Returns:
  #'   list with entries "cols", "kMs" and "log". 
  
  print(paste0(cellType, ": Reassigning genes to modules with higher ", fuzzyModMembership))
  
  # initialise
  tryCatch({
    cols_original <- cols_new <- cols
    reassign_total  <- reassign_t_1 <- logical(length(cols))
    reassign_t <- ! logical(length(cols))
    min_change = 5 # the minimum change from one iteration to the next required to continue
    t = 1
    MEs <- NULL
    kMs <- NULL 
    log=NULL 
    
    if (!is.null(cols) & length(unique(cols))>1) {
      while(TRUE) {
        if (fuzzyModMembership == "kME") {
          message(paste0(cellType, ": Computing Module Eigengenes"))
          MEs <- moduleEigengenes_uv(expr=datExpr, # TODO do we need to make it into a dataframe with names?..
                                     colors = cols,
                                     excludeGrey=T)$eigengenes
          
          message(paste0(cellType, ": Computing ", fuzzyModMembership, "s"))
          kMs <- signedKME(as.matrix(datExpr),
                           MEs,
                           outputColumnName = "",
                           corFnc = corFnc)
          
        }  else if (fuzzyModMembership == "kIM") {
          
          kMs <- kIM_eachMod_norm(dissTOM=dissTOM, 
                                  colors=cols, 
                                  verbose=verbose,
                                  excludeGrey=T)
        }
        
        message(paste0(cellType, ": Computing primary "), fuzzyModMembership, "s")
        
        pkMs <- pkMs_fnc(kMs=kMs, colors=cols)
        
        maxkMs <- max.col(kMs, ties.method = "random") #  integer vector of length nrow(kMs)
        # Reassign if there is a kM value to another module which is more than reassign_threshold times the current pkM value
        cols_new <- mapply(function(i, j, pkM) if (kMs[i,j] > reassign_threshold*pkM) colnames(kMs)[j] else cols[i], i = 1:length(maxkMs), j = maxkMs, pkM=pkMs, SIMPLIFY=T) # get new module assignment vector
        cols_new[cols=="grey"] <- "grey" # Do not reallocate genes previously allocated to grey, since they are likely to get reallocated as "grey" is not a real cohesive module
        names(cols_new) <- names(cols)
        
        reassign_t <- cols_new != cols
        
        if ((t>1 & sum(reassign_t_1) - sum(reassign_t) < min_change) | sum(reassign_t) < min_change | t >= max_iter) break #else message(paste0(cellType, ", iteration ", t, ": reassigning ",  sum(reassign_t), " genes to new modules based on their ", fuzzyModMembership))
        
        cols <- cols_new
        
        reassign_total <- reassign_total | reassign_t
        reassign_t_1 <- reassign_t
        
        t = t+1
      }
      
      log <- data.frame("gene" = names(cols)[reassign_total], 
                        "original_module" = cols_original[reassign_total], 
                        "new_module" = cols_new[reassign_total], stringsAsFactors=F, row.names=NULL)
      
      if (verbose > 0) print(paste0(cellType, ": a total of ", sum(reassign_total), " genes reassigned to new modules"))
      
    } else {
      warning(paste0(cellType, ": one or no modules, nothing to reassign"))
    }
    return(list("colors" = cols_new, "kMs" = kMs, "log" = log))
  }, 
  error = function(c) {
    warning(paste0("kM_reassign_fnc failed for ", cellType, " with the following error: ", c))
    return(list("colors" = cols_original, 
                "kMs" = NULL,
                "log" = NULL))}
  )
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kIM_eachMod_norm = function(dissTOM, colors, verbose=2, excludeGrey=T, do.par=F, n_cores=3) {
  # compute intramodularConnectivity for every gene with regard to every module
  # Args:
  #   dissTOM: distance matrix from WGCNA in sparse matrix form
  #   colors: vector of module assignment 'color' labels with gene names
  #   verbose: controls messages
  #   excludeGrey: compute "grey" kME?
  # Value:
  #   kIMs: gene x module dataframe of normalised intramodular gene-module connectivity scores,
  #         which are just the mean distance between the gene and every gene in that module
  
  unique_colors = sort(names(table(colors)))#[-which(sort(names(table(colors))) == "grey")]
  if (excludeGrey ==T) unique_colors <- unique_colors[!unique_colors=="grey"]
  # Convert the distance matrix to a proximity matrix. Nb: It is square with 0 diagonal
  mTOM <- as.matrix(1-dissTOM)
  
  rm(dissTOM)
  
  colnames(mTOM) <- rownames(mTOM) <- names(colors)
  
  # make a gene * unique_colors matrix for results: each gene's degree w.r.t each color (module)
  results = matrix(nrow=nrow(mTOM), ncol=length(unique_colors))
  
  for (j in 1:length(unique_colors)) { # modules in columns
    
    idx_mTOM_col <- match(names(colors)[colors %in% unique_colors[[j]]], colnames(mTOM))
    idx_mTOM_col <- idx_mTOM_col[!is.na(idx_mTOM_col)]
    norm_term <- max(1,length(idx_mTOM_col)-1)
    #norm_term <- (sum(colors == unique_colors[[j]])-1)
    mTOM_mod_sub <- mTOM[, idx_mTOM_col, drop=F]
    cl <- NULL
    if (do.par) {
      invisible(gc()); invisible(R.utils::gcDLLs())
      cl <- try(makeCluster(spec=n_cores, type="FORK", timeout=30))
      if (!"try-error" %in% class(cl)) {
        if (verbose>2) message(paste0("Computing kIM for ", unique_colors[j], " using ", n_cores, " cores"))
        results[,j] <- parSapply(cl=cl, X = 1:nrow(mTOM_mod_sub), FUN=function(i) {
          sum(mTOM_mod_sub[i, ]) / norm_term
        })
        stopCluster(cl)
      }
    }
    if ("try-error" %in% class(cl) | !do.par) {
      do.par <- F
      if (verbose>2) message(paste0("Computing kIM for ", unique_colors[j], " serially"))
      results[,j] <- sapply( X = 1:nrow(mTOM_mod_sub), FUN=function(i) {
        sum(mTOM_mod_sub[i, ]) / norm_term
      })
    }
    rm(mTOM_mod_sub)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    if (verbose>2) message(paste0("Done computing kIM for ", unique_colors[j]))
  }
  
  if (verbose>0) message("Done computing kIMs for set of colors")
  
  # Output a dataframe with genes in rows and modules (colors) as columns
  
  results <- as.data.frame(results, row.names= rownames(mTOM))
  colnames(results) = unique_colors
  return(results)
}




############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

count_grey_in_list_of_vec = function(list_colors) {
  vec_n_grey <- sapply(list_colors, function(x) sum(x=="grey"), simplify=T)
  return(vec_n_grey)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

load_obj <- function(f) {
  # Usage: Loads (gzip compressed) file from .RData, .RDS, .loom, .csv, .txt, .tab, .delim
  #
  # Args: 
  #   f: path to file
  # returns: 
  #   RObject
  
  compressed = F
  
  if (grepl(pattern = "\\.gz|\\.gzip", x=f))  {
    compressed <- T
    f = paste0("gzfile('",f,"')")
  }
  
  if (grepl(pattern = "\\.RDS", x = f, ignore.case = T)) {
    out <- readRDS(file=if(compressed) eval(parse(text=f)) else f)
  } else if (grepl(pattern="\\.RData|\\.Rda", x=f, ignore.case = T)) { 
    env <- new.env()
    nm <- load(f, env)[1]
    out <- env[[nm]]
  } else if (grepl(pattern="\\.loom", x=f)) {
    out <- connect(filename=f, mode = "r+")
  } else if (grepl(pattern = "\\.csv", x=f)) {
    out <- read.csv(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } else if (grepl(pattern = "\\.tab|\\.tsv", x=f)) {
    out <- read.table(file=if(compressed) eval(parse(text=f)) else f, sep="\t", stringsAsFactors = F, quote="", header=T) 
  } else if (grepl(pattern = "\\.txt", x=f)) {
    out <- read.delim(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } 
  out
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

plotLabel_for_vec <- function(comb) {
  # Utility function for more easily parallelising making labels
  label = paste0("MMS=", comb[[1]], ",DS=", comb[[2]],",PAM=",comb[[3]], ",CUT=",comb[[4]],sep="") 
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

factorToIndicator <- function(vec) {
  if (!class(vec) %in% c("character", "factor")) return(vec)
  vec <- as.character(x = vec)
  mat <- matrix(data=0, nrow=length(vec), ncol=length(unique(vec)))
  for (l in 1:length(unique(vec))) {
    mat[,l] <- ifelse(vec==unique(vec)[l], yes=1, no=0)
  }
  #class(mat) <- "numeric"
  dim(mat) <- c(length(vec), length(unique(vec)))
  colnames(mat) <- unique(vec)
  return(mat)
}
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

name_for_vec = function(to_be_named, given_names, dimension=NULL){
  # dimension = 1 names the rows, 2 names columns 
  
  if (!is.null(dim(to_be_named))) { # it's a matrix or data frame
    if (dimension == 1) {
      rownames(to_be_named) <- given_names
    } else if (dimension ==2) {
      colnames(to_be_named) <- given_names
    }
  } else if (is.null(dim(to_be_named))) {
    names(to_be_named) <- given_names
  } 
  return(to_be_named)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

PPI_outer_for_vec = function(colors, 
                             pkMs, 
                             STRINGdb_species, 
                             PPI_pkM_threshold, 
                             pvalThreshold) {
  
  # Rather than for parallelising over modules within a set of modules, parallelise over several sets of colors (produced by comparing parameters)
  # It calls PPI_innver_for_vec as a subroutine
  # Args:
  #   colors
  #   pkMs
  #   STRINGdb_species
  #   PPI_pkM_threshold: numeric scalar
  #   pvalThreshold: numeric scalar
  #   
  # Returns: 
  #   colors_PPI: colors where modules that did not satisfy the PPI threshold are set to grey
  
  #suppressPackageStartupMessages(library(STRINGdb)) 
  
  module_PPI <- NULL
  module_PPI_signif <- NULL
  unique_colors <- NULL
  unique_colors_PPI <- NULL
  colors_PPI <- NULL
  
  unique_colors <- unique(colors[pkMs>PPI_pkM_threshold])
  
  if (length(unique_colors)>0){ # grey isn't included since the pkMs will not pass the threshold
    
    # generate a STRINGdb object instance
    string_db <- STRINGdb$new(version="10", 
                              species = STRINGdb_species, 
                              score_threshold=0, 
                              input_directory="") 
    
    # For each module, check STRING_db enrichment
    sapply(unique_colors, function(x) PPI_inner_for_vec(color=x, 
                                                        unique_colors = unique_colors, 
                                                        colors = colors, 
                                                        string_db = string_db), simplify=T) %>% t() -> module_PPI 
    
    module_PPI <- data.frame("colors"= unique_colors, module_PPI, stringsAsFactors=F, row.names = NULL)
    module_PPI$q.value <- as.numeric(module_PPI$q.value)
    module_PPI$expected.interactions <- as.numeric(module_PPI$expected.interactions)
    
  } else {
    module_PPI <- data.frame(colors= "grey", q.value=1, expected.interactions=0)
  }
  # FILTER MODULES ON PPI ENRICHMENT  
  
  if (!is.null(module_PPI)) {
    if (sum(module_PPI$q.value < pvalThreshold)>0) {
      module_PPI_signif <- module_PPI[module_PPI$q.value < pvalThreshold,]
      unique_colors_PPI = unique_colors[module_PPI$q.value < pvalThreshold]
      genes_PPI_idx <- colors %in% unique_colors_PPI
      colors_PPI <- colors
      colors_PPI[!genes_PPI_idx] <- "grey"
    } else {
      colors_PPI <- rep("grey", length(colors))
    }
  } else {
    module_PPI_signif <- NULL
    colors_PPI <- rep("grey", length(colors))
  }  
  
  return(list("colors_PPI" = colors_PPI, "module_PPI" = module_PPI, "module_PPI_signif" = module_PPI_signif))
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

PPI_inner_for_vec <- function(color, unique_colors, colors, string_db) {
  # Usage: 
  #   Parallelise over modules within a set of modules
  #   utility function to parallelise the PPI STRINGdb call
  # args: 
  #   color should be one of the unique_colors
  # value:
  #   module_PPI: a list with named entries 'p-value' and 'expected interactions'
  
  ppi <- data.frame(gene = names(colors[colors==color])) # extract the genes with the corresponding color to dataframe
  module_PPI <- list('q.value' = 1, 'expected.interactions' = NA)
  
  tryCatch({
    
    example1_mapped <- string_db$map(ppi, 'gene', removeUnmappedRows = TRUE ) # Check the dataframe genes' PPI. May produce error: we couldn't map to STRING 100% of your identifiers
    # (ctd.) .. Error in if (hitList[i] %in% V(ppi_network)$name) { : argument is of length zero
    hits <- example1_mapped$STRING_id
    module_PPI['q.value'] = p.adjust(string_db$get_ppi_enrichment(hits)$enrichment,method = 'fdr',n = length(names(colors)))
    module_PPI['expected.interactions'] = string_db$get_ppi_enrichment(hits)$lambda
    
  }, error = function(c) {
    module_PPI <- list('q.value'=1, 'expected.interactions' = NA) # probably unnecessary since defined above but just in case it has been changed in place
  }) 
  
  return(module_PPI)
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

mapMMtoHs = function(modulekM,
                     colors = NULL,
                     log_dir,
                     data_prefix,
                     run_prefix,
                     mapping_orthology) {
  
  if (!is.null(modulekM$genes)) modulekM$genes <- NULL
  
  log_not_mapped_filepath = paste0(log_dir,"genes_orthology_not_mapped_",data_prefix,"_", run_prefix, ".tab")
  
  mapping = data.frame(ensembl.mouse=row.names(modulekM))
  # orthology mapping
  mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse, mapping_orthology$mmusculus_homolog_ensembl_gene) ]
  #mapping$ensembl.human[mapping$ensembl.human == ""] = NA
  df_not_mapped = mapping[is.na(mapping$ensembl.human),]
  append = !file.exists(log_not_mapped_filepath)
  col.names = !append 
  write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t", col.names = col.names, row.names=F, append=append)
  
  modulekM$ensembl = mapping$ensembl.human
  
  # Also name colors
  if (!is.null(colors)) {
    colors_ens <- colors[!is.na(modulekM$ensembl)]
    names(colors_ens) <- na.omit(mapping$ensembl.human)
  }
  
  modulekM <- na.omit(modulekM) # does this not produce a na.omit out object?
  
  ### 180508_v.18_dev2
  #tmp = within(modulekM, rm("symbol","ensembl"))
  tmp = within(modulekM, rm("ensembl"))
  ###
  
  # Average duplicated gene IDs
  modulekM_ens <-aggregate(tmp, by=list(modulekM$ensembl),FUN=mean, na.rm=TRUE)
  rownames(modulekM_ens) = modulekM_ens$Group.1
  modulekM_ens = within(modulekM_ens, rm("Group.1"))
  
  if (!is.null(colors)) colors_ens <- colors_ens[names(colors_ens) %in% rownames(modulekM_ens)] else NULL 
  
  modulekM_ens$genes <- NULL
  
  prop.mapped <- round(sum(!is.na(mapping$ensembl.human))/nrow(mapping),2)
  
  return(list(kM=modulekM_ens, colors=colors_ens, prop.mapped=prop.mapped))
  
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

cor_magma_pval <- function(gwas_pvals, 
                           data, 
                           indices) {
  x <- gwas_pvals
  y <- data[indices] # allows boot to select sample 
  cor = cor.test(x,y,method="spearman", exact=F)
  return(cor$estimate)
} 

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kM_magma <- function(cellType,
                     modulekM,
                     gwas,
                     test_type = c("full_kME_spearman", "module_genes_spearman", "module_genes_t"),
                     genes_background_ensembl_hs,
                     colors_hs = NULL) {
  # Usage: Subroutine to calculate spearman's correlation between gene module membership and GWAS gene significance
  # Args: 
  # Returns:
  
  tryCatch({
    unique_colors = colnames(modulekM)
    table.kM.cor.p = table.kM.cor.r = table.kM.cor.emp.p <- matrix(NA,nrow=length(unique_colors),ncol=length(gwas)) 
    rownames(table.kM.cor.r) = rownames(table.kM.cor.p) = rownames(table.kM.cor.emp.p) = unique_colors
    colnames(table.kM.cor.r) = colnames(table.kM.cor.p) = colnames(table.kM.cor.emp.p) = names(gwas) 
    
    for (col in unique_colors) {
      for (j in 1:length(gwas)) {
        genes = if (test_type=="full_kME_spearman") intersect(rownames(modulekM), gwas[[j]]$gene_name) else if (test_type %in% c("module_genes_spearman", "module_genes_t")) intersect(names(colors_hs)[colors_hs==col], gwas[[j]]$gene_name)  
        
        if (test_type %in% c("full_kME_spearman", "module_genes_spearman")) {
          x = -log10(gwas[[j]]$P[match(genes, gwas[[j]]$gene_name)])
          y = modulekM[match(genes,rownames(modulekM)), col]
          
          cor = cor.test(x,y,method="spearman", exact=F)
          table.kM.cor.r[col,j] <- cor$estimate
          table.kM.cor.p[col,j] <- cor$p.value
          
          # Generate 1000 (10000) permutation distribution samples and their associated correlation coefficients
          # For a 0.05 significance threshold this gives an error of the p-value of sqrt(p(1-p)/k) ~ 0.00689202437 (0.0021794497)
          # see https://stats.stackexchange.com/questions/80025/required-number-of-permutations-for-a-permutation-based-p-value
          boot_out <- boot(data = y,
                           statistic=cor_magma_pval,
                           R=R,
                           sim = "permutation",
                           gwas_pvals=x,
                           parallel = "no")
          
          # compute the empirical probability of the p-value - should correspond to the p-value, ideally..
          table.kM.cor.emp.p[col,j] <- ecdf(boot_out$t)(boot_out$t0)
          
        } else if (test_type == "module_genes_t") {
          
          x = -log10(gwas[[j]]$P[match(genes, gwas[[j]]$gene_name)])
          shared_genes_tmp <- intersect(genes_background_ensembl_hs, gwas[[j]]$gene_name)
          #shared_genes_tmp <- shared_genes_tmp[!shared_genes_tmp %in% genes]
          y = -log10(gwas[[j]]$P[match(shared_genes_tmp, gwas[[j]]$gene_name)])
          
          test = tryCatch({t.test(x=x,y=y, alternative = "g")}, error = function(c) {return(list(p.value=1))})
          
          table.kM.cor.p[col,j] <- test$p.value
          table.kM.cor.r[col,j] <- 1e-5 # to avoid 0, which will lead to problems when taking sign of the correlation to adjust p-values
          table.kM.cor.emp.p[col,j] <- 1        
        }
      }
    }
    
    rownames(table.kM.cor.r) <- rownames(table.kM.cor.p) <- rownames(table.kM.cor.emp.p) <- paste0(cellType, "__", rownames(table.kM.cor.p))
    
    return(list('p.val'= table.kM.cor.p, 'corrCoef' = table.kM.cor.r, 'emp.p.val' = table.kM.cor.emp.p))}, 
    error = function(c) {warning(paste0("kM_magma failed for ", cellType))})
}



############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

cellModEmbed <- function(datExpr, 
                         colors=NULL, 
                         latentGeneType,
                         cellType=NULL,
                         kMs=NULL,
                         excludeGrey=T,
                         dissTOM=NULL) {
  # datExpr should be cell x gene matrix or data.frame with ALL genes and ALL cells in the analysis
  # colors is a character vector with gene names
  # the genes should be ordered identically
  
  # prepare celltype character for naming columns
  tryCatch({
    if (is.null(cellType)) cellType_prefix <- "" else cellType_prefix <- paste0(cellType, "__")
    if (latentGeneType == "ME") { 
      
      colors_full <- rep("grey", times = ncol(datExpr))
      names(colors_full) <- colnames(datExpr)
      
      for (col in unique(colors)) {
        colors_full[names(colors_full) %in% names(colors[colors==col])] <- col
      }
      
      embed_mat <- moduleEigengenes_uv(expr=datExpr, 
                                       colors = colors_full, 
                                       impute = TRUE, 
                                       nPC = 1, 
                                       align = "along average", 
                                       excludeGrey = excludeGrey, 
                                       subHubs = TRUE,
                                       trapErrors = FALSE, 
                                       returnValidOnly = trapErrors, 
                                       scale = TRUE,
                                       verbose = 0, 
                                       indent = 0)$eigengenes
      
      # list_datExpr <- lapply(unique(colours)[unique(colours)!="grey"], function(x) datExpr[,match(names(colours)[colours==x], colnames(datExpr))])
      # embed_mat <- as.matrix(sapply(list_datExpr, function(x) prcomp_irlba(t(x), n=1, retx=T)$rotation, simplify = T))
      
      colnames(embed_mat) <- paste0(cellType_prefix, gsub("ME", "", colnames(embed_mat)))
      
    } else if (latentGeneType == "IM") {
      
      list_datExpr <- lapply(colnames(kMs), function(x) datExpr[ , match(names(colors[colors==x]), colnames(datExpr))])
      
      list_kMs <- lapply(colnames(kMs), function(module) kMs[match(names(colors)[colors==module], rownames(kMs)),module])
      
      embed_mat <- as.matrix(mapply(function(x,y) x %*% as.matrix(y),
                                    x = list_datExpr,
                                    y = list_kMs,
                                    SIMPLIFY=T))
      # datExpr_1 = datExpr[, match(rownames(kMs), colnames(datExpr), nomatch=0) [match(rownames(kMs), colnames(datExpr),  nomatch=0)>0]]
      # embed_mat <- as.matrix(datExpr_1) %*% as.matrix(kMs)
      colnames(embed_mat) <- paste0(cellType_prefix, colnames(kMs)) 
      rownames(embed_mat) <- rownames(datExpr)
    } 
    return(embed_mat)
  }, error = function(c) {warning(paste0("cellModEmbed failed for ", cellType))})
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

moduleEigengenes_uv <- function (expr, colors, 
                                 impute = TRUE, 
                                 nPC = 1, 
                                 align = "along average", 
                                 excludeGrey = FALSE, 
                                 grey = if (is.numeric(colors)) 0 else "grey", 
                                 subHubs = TRUE, 
                                 trapErrors = FALSE, 
                                 returnValidOnly = trapErrors, 
                                 softPower = 8, 
                                 scale = TRUE, 
                                 verbose = 0, 
                                 indent = 0) 
{
  # moduleEigengenes function modified to also return left (u) eigenvectors of length n_gene_in_module
  
  spaces = indentSpaces(indent)
  if (verbose == 1) 
    printFlush(paste(spaces, "moduleEigengenes: Calculating", 
                     nlevels(as.factor(colors)), "module eigengenes in given set."))
  if (is.null(expr)) {
    stop("moduleEigengenes: Error: expr is NULL. ")
  }
  if (is.null(colors)) {
    stop("moduleEigengenes: Error: colors is NULL. ")
  }
  if (is.null(dim(expr)) || length(dim(expr)) != 2) 
    stop("moduleEigengenes: Error: expr must be two-dimensional.")
  if (dim(expr)[2] != length(colors)) 
    stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
  if (is.factor(colors)) {
    nl = nlevels(colors)
    nlDrop = nlevels(colors[, drop = TRUE])
    if (nl > nlDrop) 
      stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
                 "Use colors[, drop=TRUE] to get rid of them."))
  }
  if (softPower < 0) 
    stop("softPower must be non-negative")
  alignRecognizedValues = c("", "along average")
  if (!is.element(align, alignRecognizedValues)) {
    printFlush(paste("ModulePrincipalComponents: Error:", 
                     "parameter align has an unrecognised value:", align, 
                     "; Recognized values are ", alignRecognizedValues))
    stop()
  }
  maxVarExplained = 10
  if (nPC > maxVarExplained) 
    warning(paste("Given nPC is too large. Will use value", 
                  maxVarExplained))
  nVarExplained = min(nPC, maxVarExplained)
  modlevels = levels(factor(colors))
  if (excludeGrey) 
    if (sum(as.character(modlevels) != as.character(grey)) > 0) {
      modlevels = modlevels[as.character(modlevels) != 
                              as.character(grey)]
    } else {
      stop(paste("Color levels are empty. Possible reason: the only color is grey", 
                 "and grey module is excluded from the calculation."))
    }
  PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
                                ncol = length(modlevels)))
  
  PrinComps_l = vector(mode="list", length= length(modlevels))  
  
  # PrinComps_l = data.frame(matrix(NA, nrow = dim(expr)[[2]], 
  #                               ncol = length(modlevels)))
  
  averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
                               ncol = length(modlevels)))
  
  averExpr_l = vector(mode="list", length= length(modlevels)) 
  # averExpr_l = data.frame(matrix(NA, nrow = dim(expr)[[2]], 
  #                              ncol = length(modlevels)))
  
  
  varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
  validMEs = rep(TRUE, length(modlevels))
  validAEs = rep(FALSE, length(modlevels))
  isPC = rep(TRUE, length(modlevels))
  isHub = rep(FALSE, length(modlevels))
  validColors = colors
  names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
                           sep = "")
  names(PrinComps_l) = modlevels
  
  names(averExpr) = paste("AE", modlevels, sep = "")
  names(averExpr_l) = paste("AE", modlevels, sep = "")
  
  for (i in c(1:length(modlevels))) {
    if (verbose > 1) 
      printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", 
                       modlevels[i]))
    modulename = modlevels[i]
    restrict1 = as.character(colors) == as.character(modulename)
    if (verbose > 2) 
      printFlush(paste(spaces, " ...", sum(restrict1), 
                       "genes"))
    datModule = as.matrix(t(expr[, restrict1])) # datModule is a genes_mod * cells = n * p matrix
    n = dim(datModule)[1]
    p = dim(datModule)[2]
    
    svd_1 = try({
      if (nrow(datModule) > 1 && impute) {
        seedSaved = FALSE
        if (exists(".Random.seed")) {
          saved.seed = .Random.seed
          seedSaved = TRUE
        }
        if (any(is.na(datModule))) {
          if (verbose > 5) 
            printFlush(paste(spaces, " ...imputing missing data"))
          datModule = impute.knn(datModule, k = min(10, 
                                                    nrow(datModule) - 1))
          try({
            if (!is.null(datModule$data)) 
              datModule = datModule$data
          }, silent = TRUE)
        }
        if (seedSaved) 
          .Random.seed <<- saved.seed
      }
      if (verbose > 5) 
        printFlush(paste(spaces, " ...scaling"))
      if (scale) 
        datModule = t(scale(t(datModule)))
      
      if (verbose > 5) 
        printFlush(paste(spaces, " ...calculating SVD"))
      
      svd_1 = svd(datModule, nu = min(n, p, nPC), nv = min(n,p, nPC))
      
      if (verbose > 5) 
        printFlush(paste(spaces, " ...calculating PVE"))
      veMat = cor(svd_1$v[, c(1:min(n, p, nVarExplained))], 
                  t(datModule), use = "p")
      varExpl[c(1:min(n, p, nVarExplained)), i] = rowMeans(veMat^2, 
                                                           na.rm = TRUE)
      svd_1
    }, silent = TRUE)
    
    pc <- svd_1$v[, 1]
    pc_l <- svd_1$u[, 1]
    
    ###############################
    
    if (class(pc) == "try-error") {
      if ((!subHubs) && (!trapErrors)) 
        stop(pc)
      if (subHubs) {
        if (verbose > 0) {
          printFlush(paste(spaces, " ..principal component calculation for module", 
                           modulename, "failed with the following error:"))
          printFlush(paste(spaces, "     ", pc, spaces, 
                           " ..hub genes will be used instead of principal components."))
        }
        isPC[i] = FALSE
        pc = try({
          scaledExpr = scale(t(datModule))
          covEx = cov(scaledExpr, use = "p")
          covEx[!is.finite(covEx)] = 0
          modAdj = abs(covEx)^softPower
          kIM = (rowMeans(modAdj, na.rm = TRUE))^3
          if (max(kIM, na.rm = TRUE) > 1) 
            kIM = kIM - 1
          kIM[is.na(kIM)] = 0
          hub = which.max(kIM)
          alignSign = sign(covEx[, hub])
          alignSign[is.na(alignSign)] = 0
          isHub[i] = TRUE
          pcxMat = scaledExpr * matrix(kIM * alignSign, 
                                       nrow = nrow(scaledExpr), ncol = ncol(scaledExpr), 
                                       byrow = TRUE)/sum(kIM)
          pcx = rowMeans(pcxMat, na.rm = TRUE)
          varExpl[1, i] = mean(cor(pcx, t(datModule), 
                                   use = "p")^2, na.rm = TRUE)
          pcx
        }, silent = TRUE)
      }
    }
    if (class(pc) == "try-error") {
      if (!trapErrors) 
        stop(pc)
      if (verbose > 0) {
        printFlush(paste(spaces, " ..ME calculation of module", 
                         modulename, "failed with the following error:"))
        printFlush(paste(spaces, "     ", pc, spaces, 
                         " ..the offending module has been removed."))
      }
      warning(paste("Eigengene calculation of module", 
                    modulename, "failed with the following error \n     ", 
                    pc, "The offending module has been removed.\n"))
      validMEs[i] = FALSE
      isPC[i] = FALSE
      isHub[i] = FALSE
      validColors[restrict1] = grey
    } else {
      PrinComps[, i] = pc
      ae = try({
        if (isPC[i]) 
          scaledExpr = scale(t(datModule))
        averExpr[, i] = rowMeans(scaledExpr, na.rm = TRUE)
        if (align == "along average") {
          if (verbose > 4) 
            printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
          corAve = cor(averExpr[, i], PrinComps[, i], 
                       use = "p")
          if (!is.finite(corAve)) 
            corAve = 0
          if (corAve < 0) 
            PrinComps[, i] = -PrinComps[, i]
        }
        0
      }, silent = TRUE)
      if (class(ae) == "try-error") {
        if (!trapErrors) 
          stop(ae)
        if (verbose > 0) {
          printFlush(paste(spaces, " ..Average expression calculation of module", 
                           modulename, "failed with the following error:"))
          printFlush(paste(spaces, "     ", ae, spaces, 
                           " ..the returned average expression vector will be invalid."))
        }
        warning(paste("Average expression calculation of module", 
                      modulename, "failed with the following error \n     ", 
                      ae, "The returned average expression vector will be invalid.\n"))
      }
      validAEs[i] = !(class(ae) == "try-error")
      
      ########## ADDED ##########
      # Also align left singular components with average gene expression
      names(pc_l) <- rownames(datModule)
      PrinComps_l[[i]] = pc_l
      
      # if (isPC[i]) {
      #   scaledExpr_l = scale(t(datModule)) # a cell * gene_mod_i matrix
      # }
      averExpr_l[[i]] = colMeans(scaledExpr, na.rm = TRUE) # for each gene averaging over cells
      if (align == "along average") {
        if (verbose > 4) 
          printFlush(paste(spaces, " .. aligning reverse module eigengene with average gene expression."))
        corAve_l = cor(averExpr_l[[i]], PrinComps_l[[i]], 
                       use = "p")
        if (!is.finite(corAve_l)) 
          corAve_l = 0
        if (corAve_l < 0) 
          PrinComps_l[[i]] = -PrinComps_l[[i]]
      }
      
      #################
    }
  }
  allOK = (sum(!validMEs) == 0)
  if (returnValidOnly && sum(!validMEs) > 0) {
    PrinComps = PrinComps[, validMEs]
    #### ADDED #####
    PrinComps_l = PrinComps[validMEs]
    ################
    averExpr = averExpr[, validMEs]
    #### ADDED #####
    averExpr_l = averExpr_l[validMEs]
    ################
    varExpl = varExpl[, validMEs]
    validMEs = rep(TRUE, times = ncol(PrinComps))
    isPC = isPC[validMEs]
    isHub = isHub[validMEs]
    validAEs = validAEs[validMEs]
  }
  allPC = (sum(!isPC) == 0)
  allAEOK = (sum(!validAEs) == 0)
  list(eigengenes = PrinComps, 
       ##### ADDED ####
       u = PrinComps_l,
       #####
       averageExpr = averExpr, 
       ##### ADDED ####
       averageExpr_l = averExpr_l, 
       ################
       
       varExplained = varExpl, 
       nPC = nPC, 
       validMEs = validMEs, 
       validColors = validColors, 
       allOK = allOK, 
       allPC = allPC, 
       isPC = isPC, 
       isHub = isHub, 
       validAEs = validAEs, 
       allAEOK = allAEOK)
}

######################################################################
########################### GENE REMAP ###############################
######################################################################

gene_map <- function(df,
                     idx_gene_column = NULL,
                     mapping,
                     from="hgnc", 
                     to="ensembl",
                     replace = F,
                     na.rm = T) {
  # args:
  #   df: data.frame with a column or rownames containing genes to remap
  #   idx_gene_column: df gene column. NULL implies rownames (implement later)
  #   mapping: data.frame or matrix with columns corresponding to 'from' and 'to' arguments (using grep matching)
  #   replace: boolean; T to replace original gene names, F to add a new column to the data.frame
  # value:
  #   df with new gene names in place of or in addition to the original gene names
  
  # usage: 
  # df_test <- gene_map(df=load_obj("/projects/jonatan/tmp-mousebrain/tables/mousebrain_Vascular_ClusterName_1_PER3_kIM.csv"),
  #                     idx_gene_column =1,
  #                     mapping=load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"),
  #                     from="ensembl", 
  #                     to="gene_name_optimal",
  #                     replace=T,
  #                     na.rm=T)
  
  orig_class <- class(df)
  
  df <- as.data.frame(df)
  
  if (is.null(idx_gene_column) & replace==T) warning("Duplicate genes will be averaged and merged to keep row.names unique")
  
  if (!from %in% colnames(mapping)) from <- grep(pattern=from, x = colnames(mapping), ignore.case=T, value = T)
  if (!to %in% colnames(mapping)) to <- grep(pattern=to, x = colnames(mapping), ignore.case=T, value = T)
  
  stopifnot(length(from)>0 & length(to)>0)
  
  genes_from <- if (is.null(idx_gene_column)) rownames(df) else df[[idx_gene_column]]
  
  genes_to <- mapping[[to]][match(genes_from, mapping[[from]])]
  
  
  if (replace) { # remove NAs
    if (is.null(idx_gene_column)) {
      # average identical gene names to ensure unique row names 
      df_aggr <- aggregate(df, by= list(genes_to), FUN=mean, na.rm=T)
      rownames(df_aggr) <- df_aggr[["Group.1"]]
      df <- within(df_aggr, rm("Group.1"))
    } else {
      if (na.rm) {
        df <- df[!is.na(genes_to),]
        genes_to <- genes_to[!is.na(genes_to)]
      }
      df[[idx_gene_column]] <- genes_to
      colnames(df)[idx_gene_column] <- to
    }
  }  else {
    if (na.rm) {
      df <- df[!is.na(genes_to),]
      df[[to]] <- genes_to[!is.na(genes_to)]
    }
    df[[to]] <- genes_to  
  }
  
  if (orig_class != "data.frame") {
    class(df) <- orig_class
  }
  
  return(df)
}

#############

signedKME = function (datExpr, datME, outputColumnName = "kME", corFnc = "cor", 
                      corOptions = "use = 'p'") 
{
  datExpr = data.frame(datExpr)
  datME = data.frame(datME)
  output = list()
  if (dim(as.matrix(datME))[[1]] != dim(as.matrix(datExpr))[[1]]) 
    stop("Number of samples (rows) in 'datExpr' and 'datME' must be the same.")
  varianceZeroIndicatordatExpr = as.vector(apply(as.matrix(datExpr), 
                                                 2, var, na.rm = TRUE)) == 0
  varianceZeroIndicatordatME = as.vector(apply(as.matrix(datME), 
                                               2, var, na.rm = TRUE)) == 0
  if (sum(varianceZeroIndicatordatExpr, na.rm = TRUE) > 0) 
    warning("Some genes are constant. Hint: consider removing constant columns from datExpr.")
  if (sum(varianceZeroIndicatordatME, na.rm = TRUE) > 0) 
    warning(paste("Some module eigengenes are constant, which is suspicious.\n", 
                  "    Hint: consider removing constant columns from datME."))
  no.presentdatExpr = as.vector(apply(!is.na(as.matrix(datExpr)), 
                                      2, sum))
  if (min(no.presentdatExpr) < 4) 
    warning(paste("Some gene expressions have fewer than 4 observations.\n", 
                  "    Hint: consider removing genes with too many missing values or collect more arrays."))
  corExpr = parse(text = paste("data.frame(", corFnc, "(datExpr, datME ", 
                               prepComma(corOptions), "))"))
  output = eval(corExpr)
  output[no.presentdatExpr < 4, ] = NA
  names(output) = paste(outputColumnName, substring(names(datME), 
                                                    first = 3, last = 100), sep = "")
  
  dimnames(output)[[1]] = names(datExpr)
  output
}

#############

detectCores_plus <- function(Gb_max=250, 
                             additional_Gb=1) {
  # args: 
  #  Gb_max: ceiling on session memory usage in Gb, assuming that each worker duplicates the session memory
  #  additional_Gb: max additional memory requirement for new (temporary) objects created within a parallel session 
  # value:
  #   n_cores (integer)
  obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x))))) / 1024^3)
  max(1, min(detectCores(), Gb_max %/% (obj_size_Gb + additional_Gb))-1) 
}

############################################################################################################################################################
########################################################## safeParallel ####################################################################################
############################################################################################################################################################

safeParallel = function(fun, args, simplify=F, MARGIN=NULL, n_cores=NULL, Gb_max=NULL, outfile=NULL,  ...) {
  # summary: calls the appropriate parallel computing function, with load balancing, 
  #          and falls back on vectorised equivalent if makeCluster hangs or the parallelised computation fails.
  #          simplify defaults to FALSE
  # args:
  #   fun: function to run in parallel. The (unevaluated) object
  #   args: a named list of arguments for fun to iterate over
  #   n_cores: n_cores. If NULL, a safe estimate is made based on size of objects in the global environment 
  #        and the length of the iterables in args
  #   outfile: passed to the parallelising function
  #   .. : additional arguments for fun to evaluate, but not iterate over
  # value:
  #   list (or if args contains "SIMPLIFY" = T, a vector or matrix); in case of failure, NULL
  
  if (length(args)==1) names(args) <- "X"
  
  list_out <- NULL
  
  if (is.null(n_cores)) {
    if (is.null(Gb_max)) Gb_max=200
    additional_Gb = max(as.numeric(sapply(args, FUN = function(x) object.size(x), simplify = T)))/1024^3
    obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
    n_cores <- min(max(sapply(args, length)), min(detectCores()%/%3, Gb_max %/% (obj_size_Gb + additional_Gb))-1)
  }
  
  if (n_cores >= 2) {
    cl <-  if (!is.null(outfile)) try(makeCluster(spec=max(1,n_cores), type="FORK", timeout=30, outfile = outfile)) else try(makeCluster(spec=max(1,n_cores), type="FORK", timeout=30))
  } else {
    cl <- "failed"
    class(cl) <- "try-error"
  }
  #cl <- "hi"
  #class(cl) <- "try-error"
  if (!"try-error" %in% class(cl)) {
    if (length(args)>1) {
      fnc <- "clusterMap"
      args[[".scheduling"]] = c("dynamic")
      args[["SIMPLIFY"]] <- simplify
    } else {
      fnc <- "parLapplyLB"
      if (simplify) {
        fnc <- "parSapplyLB"
        args[["simplify"]] <- simplify
      }
      if (!is.null(MARGIN)) {
        fnc <- "parApplyLB"
        args[["MARGIN"]] <- MARGIN
        args[["simplify"]] <- NULL
      }
    }
    args[["fun"]] <- fun
    args[["cl"]] = cl
    
  } else if ("try-error" %in% class(cl)) {
    
    if (length(args)>1) {
      
      fnc = "mapply"
      args[["SIMPLIFY"]] <- simplify
      
    } else {
      
      fnc <- "lapply"
      if (simplify) {
        fnc <- "sapply"
        args[["simplify"]] = simplify
      }
      if (!is.null(MARGIN)) {
        fnc <- "apply"
        args[["MARGIN"]] <- MARGIN
        args[["simplify"]] <- NULL
      }
    }
    
    args[["FUN"]] <- fun
    
  }
  
  list_out <- tryCatch({ 
    do.call(what=fnc, args = args)
  }, error = function(err) {
    invisible(gc())
    
    if (fnc=="clusterMap") {
      fnc <- "mapply" 
    } else if (fnc=="parLapplyLB") {
      fnc <- "lapply"
    } else if (fnc=="parSapplyLB") {
      fnc <- "sapply"
      args[["SIMPLIFY"]] <- NULL
      args[["simplify"]] <- simplify
    } else if (fnc=="parApplyLB") {
      fnc <- "apply"
    }
    
    args[["cl"]] <- args[["fun"]] <- args[[".scheduling"]] <- NULL
    args[["FUN"]] <- fun
    do.call(what=fnc, args = args)
  })
  
  if (!"try-error" %in% class(cl)) try(stopCluster(cl))
  invisible(gc())
  
  list_out 
  
}


######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(raw.data_path) & is.null(resume)) stop("Input data path not found and no previous session to resume")

# if no output directory provided, use here to infer 

if (is.null(project_dir)) {
  #pos <- tail(gregexpr("/", raw.data_path)[[1]],2)[1]
  project_dir = here()#substr(raw.data_path, 1, pos)
}

# if specified output directory doesn't exist, create it 
if (!file.exists(project_dir)) {
  dir.create(project_dir) 
  message("Project directory not found, new one created")
}

plots_dir = paste0(project_dir,"plots/")
if (!file.exists(plots_dir)) dir.create(plots_dir) 

tables_dir = paste0(project_dir,"tables/")
if (!file.exists(tables_dir)) dir.create(tables_dir)

RObjects_dir = paste0(project_dir,"RObjects/")
if (!file.exists(RObjects_dir)) dir.create(RObjects_dir)

log_dir = paste0(project_dir,"log/")
if (!file.exists(log_dir)) dir.create(log_dir)

#flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)
t_start <- as.character(Sys.time())

######################################################################
######################## HARDWIRED PARAMETERS ########################
######################################################################

############################## VARIOUS ###############################

pvalThreshold = 5e-2
options(stringsAsFactors = F)
############################## ADJACENCY #############################

#selectCols = NULL
# similarity = NULL # a (signed) similarity matrix: square, symmetric matrix with entries between -1 and 1.
type = networkType # a_{ij} = cor^Beta / "signed" # a_{ij} = ((cor+1)/2)^Beta / "signed hybrid" # a_{ij} = cor^Beta if cor>0; 0 otherwise / "distance" # Note that in some functions the parameter is called 'type', in others 'networktype'. We save the value under 'networkType'
# Contentious https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/TechnicalReports/signedTOM.pdf and http://www.peterlangfelder.com/signed-or-unsigned-which-network-type-is-preferable/.
# "signed" throws away negative correlations, while "unsigned" treats them equally (at least for the adjancy between two nodes) but this can be remediated at
# least for the part of the TOM value that takes into account neighbours by using signed TOM.
# The crux of the matter is that the TOM network has to have positive weights (?) meaning that at some point we'll have to jettison information
corOptions = if (corFnc == "cor") list(use = 'p') else NULL # "use = 'p', method = 'spearman'"  or list(use = 'p', method = 'spearman') to obtain Spearman correlation
#weights = NULL # 	 optional observation weights for datExpr to be used in correlation calculation. A matrix of the same dimensions as datExpr, containing non-negative weights. Only used with Pearson correlation.
#distFnc = "dist" #	character string specifying the function to be used to calculate co-expression similarity for distance networks. Defaults to the function dist. Any function returning non-negative values can be used. # TWEAKPARAM
#distOptions = "method = 'euclidean'" #	character string or a list specifying additional arguments to be passed to the function given by distFnc. For example, when the function dist is used, the argument method can be used to specify various ways of computing the distance.
#weightArgNames = c("weights.x", "weights.y") # character vector of length 2 giving the names of the arguments to corFnc that represent weights for variable x and y. Only used if weights are non-NULL.

############################## BICOR #################################

#x = NULL
#y = NULL
#robustX = TRUE
#robustY = TRUE
use = 'pairwise.complete.obs' #TWEAKPARAM
maxPOutliers = 0.05 # Suggested by package author in https://support.bioconductor.org/p/65124/, https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
quick = 0
pearsonFallback = "individual" # in case of columns with zero MAD
#cosine = FALSE
#cosineX = cosine
#cosineY = cosine

############################## BLOCKWISEMODULES ######################

#blocks = NULL # optional specification of blocks in which hierarchical clustering and module detection should be performed. If given, must be a numeric vector with one entry per gene of multiExpr giving the number of the block to which the corresponding gene belongs.
maxBlockSize = 5e3 # integer giving maximum block size for module detection. Should be higher than the number of genes to avoid block-wise module detection. Ignored if blocks above is non-NULL. Otherwise, if the number of genes in datExpr exceeds maxBlockSize, genes will be pre-clustered into blocks whose size should not exceed maxBlockSize.
blockSizePenaltyPower = 10 # number specifying how strongly blocks should be penalized for exceeding the maximum size. Set to a lrge number or Inf if not exceeding maximum block size is very important.
#nPreclusteringCenters = as.integer(min(ncol(datExpr0)/20, 100*ncol(datExpr0)/maxBlockSize))

loadTOM = FALSE # load TOM from previously saved file?

# Network construction arguments: correlation options

corType = if (corFnc == "cor") "pearson" else "bicor" # */!\* This same parameter is called 'corFnc' in adjancency().
# /!\: "bicor" gives: Error in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs",  : unused arguments (weights.x = NULL, weights.y = NULL) # TWEAKPARAM
quickCor = 0 # as defined for bicor. Real number between 0 and 1 that controls the handling of missing data in the calculation of correlations. 1 speeds up computation at the risk of large errors.
cosineCorrelation = FALSE # Logical: should the cosine version of the correlation calculation be used? The cosine calculation differs from the standard one in that it does not subtract the mean.

# Adjacency function options
#networkType = # Set value to same as in adjacency()
replaceMissingAdjacencies = TRUE # Should missing values in the calculation of adjacency be set to 0?
suppressTOMForZeroAdjacencies = FALSE

# Saving or returning TOM
getTOMs = NULL
saveTOMs = FALSE # should the networks (topological overlaps) be saved for each run?
#saveTOMFileBase = "TOM" # RObject name for a saved TOM

# Basic tree cut options
#detectCutHeight = 0.995 # Default value # TWEAKPARAM
#minModuleSize = minClusterSize # TWEAKPARAM

#useBranchEigennodeDissim = FALSE # should branch eigennode (eigengene) dissimilarity be considered when merging branches in Dynamic Tree Cut?
#minBranchEigennodeDissim = mergeCutHeight # Minimum consensus branch eigennode (eigengene) dissimilarity for branches to be considerd separate. The branch eigennode dissimilarity in individual sets is simly 1-correlation of the eigennodes; the consensus is defined as quantile with probability consensusQuantile.

#stabilityLabels = NULL # 	Optional matrix of cluster labels that are to be used for calculating branch dissimilarity based on split stability. The number of rows must equal the number of genes in multiExpr; the number of columns (clusterings) is arbitrary. See branchSplitFromStabilityLabels for details.
#stabilityCriterion = c("Individual fraction", "Common fraction") # One of c("Individual fraction", "Common fraction"), indicating which method for assessing stability similarity of two branches should be used. We recommend "Individual fraction" which appears to perform better; the "Common fraction" method is provided for backward compatibility since it was the (only) method available prior to WGCNA version 1.60.
#minStabilityDissim = NULL # Minimum stability dissimilarity criterion for two branches to be considered separate. Should be a number between 0 (essentially no dissimilarity required) and 1 (perfect dissimilarity or distinguishability based on stabilityLabels). See branchSplitFromStabilityLabels for details.

# Gene reassignment, module trimming, and module "significance" criteria

reassignThreshold = 1e-6
minCoreKME = 0.5 # a number between 0 and 1. If a detected module does not have at least minModuleKMESize genes with eigengene connectivity at least minCoreKME, the module is disbanded (its genes are unlabeled and returned to the pool of genes waiting for module detection). 0.5 is default.
#minCoreKMESize = minModuleSize/3 # see minCoreKME above. # TWEAKPARAM
minKMEtoStay = 0.2 # 0.3 # 	genes whose eigengene connectivity to their module eigengene is lower than minKMEtoStay are removed from the module. # TWEAKPARAM
kME_reassign_threshold = 1.25
# Module merging options

#mergeCutHeight = 0.2 # Dendrogram cut height for module merging. Note that the parameter moduleMergeCutHeight, used in mergeCloseModules, = mergeCutHeight. # TWEAKPARAM
# Value of 0.1 is used in Gandal,.., Geschwind et al 2018. Default 0.15.

############################## BOOTSTRAP #############################

# home-made - defined further down
#nPermutations = if (test_run==FALSE) 100 else 3 # TODO set to 100 as in Gandal,...,Geschwind et al 2018
startRunIndex = 1
replace = F # Sample with replacement, or 2/3 without replacement? Galdal et al use FALSE. Using replacement increases variance of the sampling mean. TWEAKPARAM
fraction = if (replace) 1.0 else 0.66 # TWEAKPARAM.

############################## BOOT #############################

#Number of times to repeat the calculation on permuted data (for kME MAGMA)
R = 1e4

############################ CONSENSUSKME ############################

signed = if (networkType == "signed" | networkType == "signed hybrid") TRUE else if (networkType == "unsigned") FALSE # In signed networks (TRUE), negative kME values are not considered significant and the corresponding p-values will be one-sided.
# In unsigned networks (FALSE), negative kME values are considered significant and the corresponding p-values will be two-sided.
useModules = NULL # Optional specification of module labels to which the analysis should be restricted. This could be useful if there are many modules, most of which are not interesting.
metaAnalysisWeights = NULL

# Correlation options

#corAndPvalueFnc = if (corFnc == "cor") corAndPvalue else bicorAndPvalue # TODO check that this works
corOptions = NULL # if (corFnc == "cor") list(use = 'p')
corComponent = if (corFnc == "cor") "cor" else "bicor"
getQvalues = TRUE # should q-values (estimates of FDR) be calculated?
useRankPvalue = TRUE # Logical: should the rankPvalue function be used to obtain alternative meta-analysis statistics?
rankPvalueOptions = list(calculateQvalue = getQvalues, pValueMethod = "scale")
setNames = NULL # names for the input sets. If not given, will be taken from names(multiExpr). If those are NULL as well, the names will be "Set_1", "Set_2", ....
excludeGrey = TRUE

########### CONSENSUSTOM / BLOCKWISECONSENSUSMODULES #################

# Finds ConsensusModules by computing adjacenies and TOMs, aligning and 'merging' the TOM matrices.
# See https://rdrr.io/cran/WGCNA/man/consensusTOM.html

# TOM precalculation arguments, if available

individualTOMInfo = NULL
useIndivTOMSubset = NULL # If individualTOMInfo is given, this argument allows to only select a subset of the individual set networks contained in individualTOMInfo. It should be a numeric vector giving the indices of the individual sets to be used. Note that this argument is NOT applied to multiExpr.

# Save individual TOMs?

saveIndividualTOMs = T 
individualTOMFileNames = "individualTOM-Set%s-Block%b.RData"

# Consensus calculation options: network calibration

networkCalibration =  "single quantile" #full quantile" #   "none" # network calibration method. One of "single quantile", "full quantile", "none" (or a unique abbreviation of one of them).
# Full quantile normalization, implemented in normalize.quantiles, adjusts the TOM matrices such that all quantiles equal each other (and equal to the quantiles of the component-wise average of the individual TOM matrices).

# Simple quantile calibration options

calibrationQuantile = 0.95 # if networkCalibration is "single quantile", topological overlaps (or adjacencies if TOMs are not computed) will be scaled such that their calibrationQuantile quantiles will agree.

sampleForCalibration = TRUE # if TRUE, calibration quantiles will be determined from a sample of network similarities. Note that using all data can double the memory footprint of the function and the function may fail.
sampleForCalibrationFactor = 1000 # determines the number of samples for calibration: the number is 1/calibrationQuantile * sampleForCalibrationFactor. Should be set well above 1 to ensure accuracy of the sampled quantile.
getNetworkCalibrationSamples = FALSE # logical: should samples used for TOM calibration be saved for future analysis? This option is only available when sampleForCalibration is TRUE.

# Consensus definition

consensusQuantile = 0.20 # The desired quantile to use in the consensus similarity calculation. Lower is conservative. Gandal,..,Geschwind (2018) use 0.2
# Set to quartile recommended for different datasets in http://www.genetics.ucla.edu/courses/statgene/networks/files/Langfelder-Thursday-ConsensusModules.pdf.
useMean = FALSE # should the consensus be determined from a (possibly weighted) mean across the data sets rather than a quantile?
#setWeights = NULL # Optional vector (one component per input set) of weights to be used for weighted mean consensus. Only used when useMean above is TRUE.

# Saving the consensus TOM

saveConsensusTOMs = TRUE # should the consensus topological overlap matrices for each block be saved and returned?
#consensusTOMFilePattern = "consensusTOM-block.%b.RData" # 	character string containing the file namefiles containing the consensus topological overlaps. The tag %b will be replaced by the block number.

# Internal handling of TOMs

useDiskCache = FALSE # should calculated network similarities in individual sets be temporarily saved to disk? Saving to disk is somewhat slower than keeping all data in memory, but for large blocks and/or many sets the memory footprint may be too big.
chunkSize = NULL # network similarities are saved in smaller chunks of size chunkSize.
cacheBase = ".blockConsModsCache" # character string containing the desired name for the cache files. The actual file names will consists of cacheBase and a suffix to make the file names unique.
cacheDir = "." # character string containing the desired path for the cache files.

# Alternative consensus TOM input from a previous calculation

consensusTOMInfo = NULL # optional list summarizing consensus TOM, output of consensusTOM. It contains information about pre-calculated consensus TOM. Supplying this argument replaces TOM calculation, so none of the individual or consensus TOM calculation arguments are taken into account.

# Basic tree cut options
checkMinModuleSize = TRUE # should sanity checks be performed on minModuleSize?

# Gene reassignment and trimming from a module, and module "significance" criteria

reassignThresholdPS = 1e-4 # per-set p-value ratio threshold for reassigning genes between modules.
trimmingConsensusQuantile = consensusQuantile # a number between 0 and 1 specifying the consensus quantile used for kME calculation that determines module trimming according to the arguments below.

equalizeQuantilesForModuleMerging = FALSE # specific to blockwiseConsensusModules
quantileSummaryForModuleMerging = "mean"  # specific to blockwiseConsensusModules

#Module merging options

#equalizeQuantilesForModuleMerging = FALSE # If TRUE, the quantiles of the eigengene correlation matrices (interpreted as a single vectors of non-redundant components) will be equalized across the input data sets. Note that although this seems like a reasonable option, it should be considered experimental and not necessarily recommended.
#quantileSummaryForModuleMerging = "median" # "mean" # One of "mean" or "median". If quantile equalization of the module eigengene networks is performed, the resulting "normal" quantiles will be given by this function of the corresponding quantiles across the input data sets.

mergeConsensusQuantile = consensusQuantile # consensus quantile for module merging. See mergeCloseModules for details.

# Output options

numericLabels = FALSE # should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)?

# detectCutHeight	= cutHeight # dendrogram cut height for module detection
minHeight = 0.1 # TODO if we want to use blockwiseConsensusModules

########################### CUTREEHYBRID ############################
# See https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/BranchCutting/
# TODO: read up on this clustering algorithm

#minClusterSize = 10 # min(20, ncol(datExpr0)/2 ) # 40 is the value in Gandal et al, 2018 # Set value in notebook after evaluation
#deepSplit =  2 # integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive. 

# cutHeight = NULL # Maximum joining heights that will be considered. For method=="tree" it defaults to 0.99. For method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram. #TODO research this algorithm
# /!\ the parameter name is aliased wth the cutHeight parameter of the mergeCloseModules function which instead denotes maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging. If you need to set a non-default value, define it here as 'cutreeHybridCutHeight'

# Basic tree cut options - DEFINE THESE IN THE SCRIPT
method = "hybrid"
#distM = "dissTOM" # Only used for method "hybrid". The distance matrix used as input to hclust. If not given and method == "hybrid", the function will issue a warning and default to method = "tree".

# Advanced options /!\ TODO: check if we need to set these to non-default values
#maxCoreScatter = NULL #	 Only used for method "hybrid". Maximum scatter of the core for a branch to be a cluster, given as the fraction of cutHeight relative to the 5th percentile of joining heights.
#minGap = NULL # Only used for method "hybrid". Minimum cluster gap given as the fraction of the difference between cutHeight and the 5th percentile of joining heights.
#maxAbsCoreScatter = NULL # Only used for method "hybrid". Maximum scatter of the core for a branch to be a cluster given as absolute heights. If given, overrides maxCoreScatter.
#minAbsGap = NULL # 	Only used for method "hybrid". Minimum cluster gap given as absolute height difference. If given, overrides minGap.

#minSplitHeight = NULL # Minimum split height given as the fraction of the difference between cutHeight and the 5th percentile of joining heights. Branches merging below this height will automatically be merged. Defaults to zero but is used only if minAbsSplitHeight below is NULL.
#minAbsSplitHeight = NULL # Minimum split height given as an absolute height. Branches merging below this height will automatically be merged. If not given (default), will be determined from minSplitHeight above.

# External (user-supplied) measure of branch split
#externalBranchSplitFnc = NULL # Optional function to evaluate split (dissimilarity) between two branches. [...] This argument is only used for method "hybrid".
#minExternalSplit = NULL #	 Thresholds to decide whether two branches should be merged. It should be a numeric vector of the same length as the number of functions in externalBranchSplitFnc above. Only used for method "hybrid".
#externalSplitOptions = list() # Further arguments to function externalBranchSplitFnc. [If only one external function is specified in externalBranchSplitFnc above, externalSplitOptions can be a named list of arguments or a list with one component that is in turn the named list of further arguments for externalBranchSplitFnc[[1]]. [...] Only used for method "hybrid".
#externalSplitFncNeedsDistance = NULL # 	Optional specification of whether the external branch split functions need the distance matrix as one of their arguments. Either NULL or a logical vector with one element per branch split function that specifies whether the corresponding branch split function expects the distance matrix as one of its arguments. The default NULL is interpreted as a vector of TRUE. When dealing with a large number of objects, setting this argument to FALSE whenever possible can prevent unnecessary memory utilization.
#assumeSimpleExternalSpecification = TRUE # 	Logical: when minExternalSplit above is a scalar (has length 1), should the function assume a simple specification of externalBranchSplitFnc and externalSplitOptions? If TRUE, externalBranchSplitFnc is taken as the function specification and externalSplitOptions the named list of options. This is suitable for simple direct calls of this function. If FALSE, externalBranchSplitFnc is assumed to be a list with a single component which specifies the function, and externalSplitOptions is a list with one component that is in turn the named list of further arguments for externalBranchSplitFnc[[1]].

# Partitioning Around Medoids (PAM) stage options - these are set in the notebook script after evaluating the effects of different values
#pamStage = T # Only used for method "hybrid". If TRUE, the second (PAM-like) stage will be performed. Default = TRUE. Increases specificity at the cost of sensitivity. Gandal et al 2018 set this to FALSE.
pamRespectsDendro = T # Only used if pamStage = TRUE. Logical, only used for method "hybrid". If TRUE, the PAM stage will respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters that belong to the same branch that the objects or small clusters being assigned belong to. 
useMedoids = F # If pamStage = T, compute distance using medoids (most centrally located objects in each cluster). 
# If not, uses average pairwise distance (recommended in the function documentation)
maxPamDist = 0 #	Only used for method "hybrid" and only if labelUnlabeled==TRUE. 
# maximum object distance to closest cluster that will result in the object assigned to that cluster. 
# Defaults to cutHeight. If zero,  objects are assigned to clusters only if the object–cluster dissimilarity is smaller than the radius of the cluster (see Langfelder et al, 2007, supplement)
# https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/BranchCutting/Supplement.pdf
respectSmallClusters = TRUE #	 Only used for method "hybrid" and only if labelUnlabeled==TRUE. 
# If TRUE, branches that failed to be clusters in stage 1 only because of insufficient size will be assigned together in stage 2. 
# If FALSE, all (lowest-level) objects will be assigned individually. Seems best to respect local structure

############################ GENERAL WGCNA ###########################

# The following setting is important for WGCNA
options(stringsAsFactors = FALSE);
# Multiple threads possible in RStudio server but not in desktop
#enableWGCNAThreads(nThreads = n_cores)#
randomSeed = 12345

options(pairwise.complete.obs = T) # How to handle missing values in correlations
checkMissingData = TRUE # should data be checked for excessive missing values in genes and samples, and for genes with zero variance?
# Problem: in blockwiseConsensusModules, our consensus TOM seems to lose genes, which makes it impossible to find eigengenes

trapErrors = FALSE
verbose = 5 
indent = 0

############################## HCLUST ################################

# Takes the place of flashClust as of 2014

### 180611
#hclustMethod = "complete" # Gandal et al 2018 use "average".
###
# We call it hclustMethod because method is aliased with argument for cutreedynamic

################ HIERARCHICALCONSENSUSCALCULATION ####################

# Not implemented

################# HIERARCHICALCONSENSUSMODULES #######################

# Not implemented

###################### HIERARCHICALCONSENSUSTOM ######################

# Not implemented

######################### MERGECLOSEMODULES ##########################

# See also pquantile()

# Optional starting eigengenes
# MEs = NULL

# mergeCutHeight = 0.15 # defined in blockwise modules. value of 1-cor of eigengene correlation at which to merge modules. This is the default value

# Optional restriction to a subset of all sets
# useSets = NULL

# Input handling options
checkDataFormat = TRUE # If TRUE, the function will check exprData and MEs for correct multi-set structure. If single set data is given, it will be converted into a format usable for the function. If FALSE, incorrect structure of input data will trigger an error.
unassdColor = if (is.numeric(colors)) 0 else "grey"

# Options for eigengene network construction
useAbs = FALSE # 	Specifies whether absolute value of correlation or plain correlation (of module eigengenes) should be used in calculating module dissimilarity. Default: FALSE

# Options for constructing the consensus
equalizeQuantiles = TRUE # Should quantiles of the eigengene dissimilarity matrix be equalized ("quantile normalized")? The default is FALSE for reproducibility of old code; when there are many eigengenes (e.g., at least 50), better results may be achieved if quantile equalization is used. # TWEAKPARAM
quantileSummary = "mean" # One of "mean" or "median". Controls how a reference dissimilarity is computed from the input ones (using mean or median, respectively). # TWEAKPARAM

#TODO How does this actually work?

# Merging options
#moduleMergeCutHeight = mergeCutHeight # Maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging. /!\ This parameter is actually called cutHeight and hence is aliased with 'cutHeight' in cutreeHybrid.  Renamed to avoid confusion. # TWEAKPARAM
iterate = TRUE # Controls whether the merging procedure should be repeated until there is no change. If FALSE, only one iteration will be executed

# Output options
relabel = FALSE # Controls whether, after merging, color labels should be ordered by module size.
colorSeq = NULL # Color labels to be used for relabeling. Defaults to the standard color order used in this package if colors are not numeric, and to integers starting from 1 if colors is numeric.
getNewMEs = TRUE # 	Controls whether module eigengenes of merged modules should be calculated and returned.
getNewUnassdME = TRUE # When doing module eigengene manipulations, the function does not normally calculate the eigengene of the 'module' of unassigned ('grey') genes. Setting this option to TRUE will force the calculation of the unassigned eigengene in the returned newMEs, but not in the returned oldMEs.

############################ MODULEEIGENGENES ########################

# NB: If you want to call blockwiseModules where method="complete" rather than "average", you need to re-define the function!

impute = TRUE # should WGCNA call impute.knn to impute NA values for module eigengene calculation? Note that '0's will NOT be imputed and that this is only for calculating the eigengenes (which is sensitive to NAs)
nPC_wgcna = 10 # Number of principal components and variance explained entries to be calculated. Note that only the first principal component is returned; the rest are used only for the calculation of proportion of variance explained. The number of returned variance explained entries is currently min(nPC, 10)
align = "along average" # Controls whether eigengenes, whose orientation is undetermined, should be aligned with average expression (align = "along average", the default) or left as they are (align = ""). Any other value will trigger an error.
excludeGrey = FALSE # Should the improper module consisting of 'grey' genes be excluded from the eigengenes?
grey = if (is.numeric(colors)) 0 else "grey" # Value of colors designating the improper module. Note that if colors is a factor of numbers, the default value will be incorrect.
subHubs = TRUE #Controls whether hub genes should be substituted for missing eigengenes. If TRUE, each missing eigengene (i.e., eigengene whose calculation failed and the error was trapped) will be replaced by a weighted average of the most connected hub genes in the corresponding module. If this calculation fails, or if subHubs==FALSE, the value of trapErrors will determine whether the offending module will be removed or whether the function will issue an error and stop.
returnValidOnly = trapErrors  # logical; controls whether the returned data frame of module eigengenes contains columns corresponding only to modules whose eigengenes or hub genes could be calculated correctly (TRUE), or whether the data frame should have columns for each of the input color labels (FALSE).
#softPower = 6 # This is defined based on plotting the fit with scale-free topology
scale = TRUE # logical; can be used to turn off scaling of the expression data before calculating the singular value decomposition. The scaling should only be turned off if the data has been scaled previously, in which case the function can run a bit faster. Note however that the function first imputes, then scales the expression data in each module. If the expression contain missing data, scaling outside of the function and letting the function impute missing data may lead to slightly different results than if the data is scaled within the function. # TWEAKPARAM

######################### MODULEPRESERVATION #########################

# Compares test sets with a module assignments and test sets optionally with modules assigned.
# Computes the *preservation* of modules through statistics that compare adjacencies and overall module statistics.
# Also computes the *quality* of modules **in the reference set**. The quality statistics are calculated with respect
# to genes in common wiht the test set, hence the function returns *a set* of quality statistics for each reference-test pair.
# This may be counter-intuitive but it is convenient for returning the quality and preservation statistics together, even though
# the former only concern the reference network.

# multiData	  # expression data or adjacency data in the multi-set format (see checkSets). A vector of lists, one per set.
# Each set must contain a component data that contains the expression or adjacency data. If expression data are used,
# rows correspond to samples and columns to genes or probes. In case of adjacencies, each data matrix should be a symmetric
# matrix with entries between 0 and 1 and unit diagonal. Each component of the outermost list should be named.
# multiColor	# a list in which every component is a vector giving the module labels of genes in multiExpr.
# The components must be named using the same names that are used in multiExpr; these names are used top match labels to
# expression data sets.
dataIsExpr = TRUE # logical: if TRUE, multiData will be interpreted as expression data; if FALSE, multiData will be interpreted as adjacencies.
#referenceNetworks = c(2:length(multiExpr))
# A vector giving the indices of expression data to be used as reference networks.
# Reference networks must have their module labels given in multiColor.
#testNetworks = list(rep(1,length(referenceNetworks))) #	a list with one component per each entry in referenceNetworks above, giving the test networks in which
# to evaluate module preservation for the corresponding reference network. If not given, preservation will
# be evaluated in all networks (except each reference network). If referenceNetworks is of length 1,
# testNetworks can also be a vector (instead of a list containing the single vector).
# NB: We
#nPermutations = nPermutations # TODO nPermutations #	specifies the number of permutations that will be calculated in the permutation test. TODO: How costly is it?
includekMEallInSummary = FALSE  # logical: should cor.kMEall be included in the calculated summary statistics?
# Because kMEall takes into account all genes in the network, this statistic measures preservation
# of the full network with respect to the eigengene of the module.
# This may be undesirable, hence the default is FALSE.
restrictSummaryForGeneralNetworks = TRUE  # logical: should the summary statistics for general (not correlation) networks be restricted
# (density to meanAdj, connectivity to cor.kIM and cor.Adj)? The default TRUE corresponds to published work.
calculateQvalue = F # logical: should q-values (local FDR estimates) be calculated? Package qvalue must be installed for this calculation. Fails.
#  Note that q-values may not be meaningful when the number of modules is small and/or most modules are preserved.
maxGoldModuleSize = 1000 # maximum size of the "gold" module, i.e., the random sample of all network genes. TODO what is this?
maxModuleSize = 1000 # alternatively, max(levels(multiColor)) ? maximum module size used for calculations. Larger modules will be reduced
# by randomly sampling maxModuleSize genes.
ccTupletSize = 2 # 	tuplet size for co-clustering calculations. TODO ?
calculateCor.kIMall = FALSE # logical: should cor.kMEall be calculated? This option is only valid for adjacency input.
# If FALSE, cor.kIMall will not be calculated, potentially saving significant amount of time
# if the input adjacencies are large and contain many modules.
calculateClusterCoeff = FALSE # logical: should statistics based on the clustering coefficient be calculated? TOTO
# While these statistics may be interesting, the calculations are also computationally expensive.
useInterpolation = FALSE  # logical: should permutation statistics be calculated by interpolating an artificial set of evenly spaced modules?
# This option may potentially speed up the calculations, but it restricts calculations to density measures.
checkData = FALSE #TRUE  # logical: should data be checked for excessive number of missing entries? See goodSamplesGenesMS for details.
#greyName = NULL # label used for unassigned genes. Traditionally such genes are labeled by grey color or numeric label 0.
# These values are the default when multiColor contains character or numeric vectors, respectively.
savePermutedStatistics = TRUE
loadPermutedStatistics = FALSE # If interpolation is used (see useInterpolation above), the function can optionally generate
# diagnostic plots that can be used to assess whether the interpolation makes sense.
permutedStatisticsFile = if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData"
plotInterpolation = TRUE # file name to save the interpolation plots into.
interpolationPlotFile = "modulePreservationInterpolationPlots.pdf"
discardInvalidOutput = TRUE # should output columns containing no valid data be discarded? This option may be useful when input dataIsExpr is FALSE and some of the output statistics cannot be calculated. This option causes such statistics to be dropped from output.
parallelCalculation = FALSE  # Note that parallel calculations are turned off by default and will lead to somewhat DIFFERENT results
# than serial calculations because the random seed is set differently. For the calculation to actually run in parallel mode,
# a call to enableWGCNAThreads must be made before this function is called.
# Better control this locally

# out:
# The function returns a nested list of preservation statistics.
# At the top level, the list components are:
#
# quality
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of quality statistics. All logarithms are in base 10.
# **relevant for evalating clustering parameters**
#
# preservation
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of density and connectivity preservation statistics. All logarithms are in base 10.
#
# accuracy
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of cross-tabulation statistics. All logarithms are in base 10.
#
# referenceSeparability
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of module separability in the reference network. All logarithms are in base 10.
# **relevant for evalating clustering parameters**
#
# testSeparability
# observed values, Z scores, p-values, Bonferoni-corrected p-values, and (optionally) q-values of module separability in the test network. All logarithms are in base 10.
# permutationDetails	results of individual permutations, useful for diagnostics

############# SAMPLEDHIERARCHICALCONSENSUSMODULES ####################

# Not implemented

############################ SEURAT ##################################

#min.cells = 5 # Filter out genes with few cells 
do.center = T#if (corFnc == "bicor") FALSE else TRUE # for ScaleData() # update: just used scale_data==F
nPC_seurat = 120 # for RunPCA() and ProjectPCA
maxit = 1000 # for RunPCA() IRLBA package - default is 1000 https://cran.r-project.org/web/packages/irlba/irlba.pdf
fastpath = T # for RunPCA()

############################### STRINGdb ##################################

PPI_pkME_threshold = 10e-2
#pvalThreshold = 5e-2

############################### TOM ##################################

TOMType = "signed" # Takes into account the sign of the adjacency between neighbours
TOMDenom = "mean" # "min" gives the standard TOM described in Zhang and Horvath (2005), and "mean" in which the min
# function in the denominator is replaced by mean.
# The "mean" may produce better results but at this time should be considered experimental.

######################################################################
######################## REMAINING PARAMETERS ########################
######################################################################

# Entrezgene to ensemb_gene_id mapping file path for MAGMA

mapping_hs_filepath = dir(pattern="gene_annotation_hsapiens.txt.gz", path = here("data"), full.names = T) # columns: ensembl_gene_id, entrezgene, hgnc_symbol

if (data_organism == "mmusculus") {
  # For mapping symbol to ensembl
  mapping_mm_filepath = dir(pattern="Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz", path=here("data"), full.names = T)
  # Synonyms
  mapping_mm_synonyms_filepath = dir(pattern="Mus_musculus.gene_info_symbol2ensembl.gz", path= here("data"), full.names=T)
  # Mouse to human ortholog mapping  
  mapping_hs_mm_filepath =  dir(pattern="gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz", path= here("data"), full.names=T)
  
}

set.seed(randomSeed)

options(stringsAsFactors = F)

if (is.null(resume)) {
  
  ######################################################################
  ############################## PACKAGES ##############################
  ######################################################################
  
  # load packages into namespace
  
  suppressPackageStartupMessages(library("data.table"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  #suppressPackageStartupMessages(library("hdf5r", lib.loc ="~/R/x86_64-pc-linux-gnu-library/"))
  suppressPackageStartupMessages(library("Seurat"))
  
  message("Packages loaded")
  
  ######################################################################
  ############################ CHECK INPUT #############################
  ######################################################################
  
  if (min.cells < 0) stop("min.cells must be a non-negative integer") 
  
  if (!(sapply(c("all", "var.genes", "PCA"), function(x) grepl(x, genes_use, ignore.case=T)) %>% any())) {
    stop("genes_use must be one of 'all', 'var.genes' or 'PCA'")
  }
  
  if (!pca_genes %in% c('all', 'var.genes')) {
    warning("Invalid pca_genes argument, reverted to var.genes")
    pca_genes <- var.genes
  }
  
  if (!corFnc %in% c("cor", "bicor")) stop("corFnc must be one of 'cor' for Pearson's or 'bicor' for biweighted midcorrelation")
  
  if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid' (not 'signed_hybrid')")
  
  if (length(c(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight))>4) warning("Comparing different parameters increases processing time")
  
  if (min(minClusterSize) < 5) stop("minClusterSize must be a vector of integers over 5")
  
  if (max(deepSplit) > 4 | min(deepSplit) < 0) stop("deepSplit must be a vector of integers between 0 and 4")
  
  if (min(moduleMergeCutHeight) < 0 | max(moduleMergeCutHeight) > 1) stop("moduleMergeCutHeight must be a vector of doubles between 0 and 1, recommended range is between 0.1 and 0.2")
  
  if (!fuzzyModMembership %in% c("kME", "kIM")) {
    warning("Invalid fuzzeModuleMembership value -  reset to kIM (default)")
    fuzzyModMembership <- "kIM"
  }
  
  if (! (jackstrawnReplicate >= 0)) stop("jackstrawnReplicate must be 0 or higher")
  
  if (! (TOMnReplicate >= 0)) stop("TOMnReplicate must be non-negative")
  

  if (!map_genes_to_ensembl) {
    if (!is.null(magma_gwas_dir)) {
      magma_gwas_dir <- gwas_filter_traits<-NULL
      message("map_genes_to_ensembl must be TRUE to perform MAGMA gene set analysis")
    }
  }
  if (!data_type %in% c("sc", "bulk")) stop("data_type must be 'sc' or 'bulk'")
  
  if (!RAM_Gb_max >= 1 | !RAM_Gb_max <= 1000) stop("verify RAM_Gb_max value")
  
  ######################################################################
  ######################## CREATE ERROR LOG ############################
  ######################################################################
  
  path_cellClusters_dropped_log <- paste0(log_dir, data_prefix, "_", run_prefix, "_cellClusters_dropped.txt")
  
  ######################################################################
  ################# LOAD DATA AND CREATE SEURAT OBJECT ######################
  ######################################################################
  
  message("Loading expression data..")
  
  dt_raw <- fread(raw.data_path)
  dt_metadata <- fread(meta.data_path)
  
  dt_raw[,-1] %>% as.matrix %>% Matrix(., sparse = TRUE) -> spMat_raw
  
  rownames(spMat_raw) <- dt_raw[[1]]
  rm(dt_raw)
  
  df_metadata <- as.data.frame(dt_metadata)
  
  rownames(df_metadata) <- colnames(spMat_raw)
  
  seurat_obj <- CreateSeuratObject(raw.data = spMat_raw, 
                                   project=data_prefix, 
                                   meta.data = df_metadata)
    
  ######################################################################
  #################### SET SEURAT IDENT AND SUBSET NAMES ###############
  ######################################################################
  
  seurat_obj <- SetAllIdent(object = seurat_obj,id = metadata_subset_col)
  ident <- seurat_obj@ident
  sNames_0 <- names(table(seurat_obj@ident))
  
  ######################################################################
  ################### COMPUTE PERCENT MITO, PERCENT RIBO ###############
  ######################################################################
  
  message("Computing percent.mito and percent.ribo")
  
  if (data_type=="sc") {
    idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(x = seurat_obj@raw.data), ignore.case=T)
    idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = seurat_obj@raw.data), ignore.case=T) 
    percent.mito <- Matrix::colSums(seurat_obj@raw.data[idx_mito.genes, ])/Matrix::colSums(seurat_obj@raw.data)
    percent.ribo <- Matrix::colSums(seurat_obj@raw.data[idx_ribo.genes, ])/Matrix::colSums(seurat_obj@raw.data)
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.ribo, col.name = "percent.ribo")
  }
  
  ######################################################################
  ################# IF SYMBOL, REMAP TO ENSEMBL ID #####################
  ######################################################################
  
  if (map_genes_to_ensembl) {
    if (!any(grepl("ENSG|ENSMUSG",rownames(seurat_obj@raw.data)))) {
      
      if (is.null(data_organism)) {
        
        stop("Error: seurat object gene names are not in ensembl ID format. To remap from hgnc to ensembl, please provide a data_organism 'mmusculus' or 'hsapiens'")
        
      } else if (!is.null(data_organism)) {
        
        message("Mapping genes from hgnc to ensembl")
        
        if (data_organism == "mmusculus") { 
          
          # Step 1: direct mapping
          mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
          mapping = data.frame(symbol=row.names(seurat_obj@raw.data), ensembl = mapping_direct$ensembl_gene_id[ match(toupper(gsub("-|_", ".", row.names(seurat_obj@raw.data))), 
                                                                                                                      toupper(gsub("-|_", ".", mapping_direct$gene_name_optimal))) ])
          
          # Step 2: map remaining using synonyms
          mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
          mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( toupper(gsub("-|_", ".", mapping$symbol[which(is.na(mapping$ensembl)) ])) , 
                                                                                              toupper(gsub("-|_", ".", mapping_synonyms$symbol))) ]
          rm(mapping_direct, mapping_synonyms)
          # save mapping file for reference and later use
          write.csv(mapping, file=sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism), quote = F, row.names = F)
          
        } else if (data_organism == "hsapiens") {
          
          mapping_direct = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T, stringsAsFactors = F) # columns: ensembl_gene_id, entrezgene, hgnc_symbol
          mapping_direct$hgnc_symbol <- gsub("-|_", ".", mapping_direct$hgnc_symbol)
          # Step 1: direct mapping
          mapping = data.frame(symbol=row.names(seurat_obj@raw.data), ensembl = mapping_direct$ensembl_gene_id[ match(toupper(gsub("-|_", ".", row.names(seurat_obj@raw.data))), 
                                                                                                                      toupper(gsub("-|_", ".", mapping_direct$hgnc_symbol))) ])
          rm(mapping_direct)
          # save mapping file for reference and later use
          write.csv(mapping, file=sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism), quote = F, row.names = F)
          
        }
        
        # Make a log of unmapped genes 
        log_not_mapped_filepath = sprintf("%s%s_%s_not_mapped_to_ensembl.tab", log_dir, data_prefix, run_prefix)
        log_duplicate_filepath = sprintf("%s%s_%s_duplicate_ensembl_id_genenames.tab", log_dir, data_prefix, run_prefix)
        
        df_not_mapped = mapping[is.na(mapping$ensembl),]
        write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
        
        # Make a log of duplicate genes
        idx_duplicate_genes <- duplicated(mapping$ensembl)
        df_duplicate <- mapping[idx_duplicate_genes,]
        write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)
        
        # Filter out unmapped and duplicate genes from Seurat object
        seurat_obj@raw.data <- seurat_obj@raw.data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
        
        # rename Seurat object rows where mapping was successful to ensembl ID
        rownames(seurat_obj@raw.data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
        
      } 
      
    } else if (any(grepl("ENS", rownames(seurat_obj@raw.data)))) {
      
      # Still do the mapping from ensembl to symbol to allows us to easily find mitochrondrial and ribosomomal genes to regress out confounding variance
      
      message("Mapping genes from ensembl to hcgn")
      
      if (any(grepl("ENSG", rownames(seurat_obj@raw.data)))) {
        
        message("Homo sapiens ensembl id gene names detected")
        data_organism <- "hsapiens"
        
        mapping_direct = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T, stringsAsFactors = F) # columns: ensembl_gene_id, entrezgene, hgnc_symbol
        
        # Step 1: direct mapping
        mapping = data.frame(ensembl=row.names(seurat_obj@raw.data), symbol = mapping_direct$hgnc_symbol[match(toupper(row.names(seurat_obj@raw.data)), toupper(mapping_direct$ensembl_gene_id)) ])
        rm(mapping_direct)
        # save mapping file for reference and later use
        write.csv(mapping, file=sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism), quote = F, row.names = F)
        
      } else if (any(grepl("ENSMUSG", rownames(seurat_obj@raw.data)))) {
        
        message("mus musculus ensembl id gene names detected")
        data_organism <- "mmusculus"
        
        # Step 1: direct mapping
        mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
        mapping = data.frame(ensembl=row.names(seurat_obj@raw.data), symbol = mapping_direct$gene_name_optimal[ match(toupper(row.names(seurat_obj@raw.data)), toupper(mapping_direct$ensembl_gene_id)) ])
        
        # Step 2: map remaining using synonyms
        mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
        mapping$symbol[ which(is.na(mapping$symbol)) ] = mapping_synonyms$symbol[ match( toupper(mapping$ensembl[which(is.na(mapping$symbol)) ]) , toupper(mapping_synonyms$ensembl)) ]
        rm(mapping_direct, mapping_synonyms)
        # save mapping file for reference and later use
        write.csv(mapping, file=sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism), quote = F, row.names = F)
      }
      
      # Make a log of unmapped genes 
      log_not_mapped_filepath = sprintf("%s%s_%s_not_mapped_to_ensembl.tab", log_dir, data_prefix, run_prefix)
      log_duplicate_filepath = sprintf("%s%s_%s_duplicate_ensembl_id_genenames.tab", log_dir, data_prefix, run_prefix)
      
      df_not_mapped = mapping[is.na(mapping$symbol),]
      write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
      
      # Make a log of duplicate genes
      idx_duplicate_genes <- duplicated(mapping$symbol)
      df_duplicate <- mapping[idx_duplicate_genes,]
      write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)
    }
  }
  
  ######################################################################
  ######## DO SEURAT PROCESSING ON FULL EXPRESSION MATRIX ##############
  ######################################################################
  
  message("Normalizing, scaling and saving full expression matrix")
  
  # Normalise
  seurat_obj <- NormalizeData(object = seurat_obj, display.progress = T)
  
  # Scale and regress out confounders
  vars.to.regress = if (!is.null(regress_out)) regress_out[regress_out %in% names(seurat_obj@meta.data)] else NULL
  
  seurat_obj <- ScaleData(object = seurat_obj,
                          vars.to.regress = vars.to.regress,
                          model.use="linear",
                          do.par=T,
                          num.cores = min(10, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = as.numeric(object.size(seurat_obj@data))/1024^3)-1),
                          do.scale=T,
                          do.center=T,
                          display.progress = T)
  
  # Make sure we close socket workers
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  scale_data <- seurat_obj@scale.data
  
  # Save scale and regressed whole expression matrix with ensembl rownames for later use
  saveRDS(scale_data, file = sprintf("%s%s_%s_scale_regr_data_ensembl.RDS.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")

  rm(scale_data)
  
  ######################################################################
  ######## EXTRACT METADATA AND CONVERT FACTORS TO MODEL MATRIX ########
  ######################################################################
  
  metadat_names <- colnames(seurat_obj@meta.data)
  
  # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
  if (!is.null(metadata_corr_col)) {
    if (any(colnames(seurat_obj@meta.data) %in% metadata_corr_col)) {
      metadata <- matrix(NA, nrow=nrow(seurat_obj@meta.data), ncol=1)
      include <- seurat_obj@meta.data[,colnames(seurat_obj@meta.data) %in% metadata_corr_col, drop=F]
      for (i in 1:ncol(include)) {
        if (class(include[,i]) %in% c("factor", "character")) {
          metadata <- cbind(metadata, factorToIndicator(include[,i, drop=T]))        
        } else {
          metadata <- cbind(metadata, include[,i, drop=F])
        }
      }
      #rownames(metadata) <- rownames(seurat_obj@meta.data)
      metadata <- metadata[,-1, drop=F]
      if (!is.null(metadata_corr_filter_vals)) metadata <- metadata[, toupper(colnames(metadata)) %in% toupper(metadata_corr_filter_vals), drop = F]
      metadata <- as.data.frame(metadata)
      
      # Filter out any metadata columns where all the values are identical
      metadata <- metadata[apply(metadata, MARGIN=2, FUN = function(x) length(unique(x))>1)]
      if (ncol(metadata) == 0) metadata <- NULL    
      rownames(metadata) = rownames(seurat_obj@meta.data)
    } else metadata <- NULL
  } else metadata <- NULL
  
  
  ######################################################################
  ######## DO SEURAT PROCESSING ON SUBSETTED EXPRESSION MATRICES #######
  ######################################################################
  
  message("Subsetting the dataset")
  
  subsets <- lapply(sNames_0, function(name) SubsetData(seurat_obj,
                                                        ident.use = name, 
                                                        do.scale = F, 
                                                        do.center = F,
                                                        subset.raw = T,
                                                        do.clean = T))
  
  
  # args = list("X"=subsets)
  # fun =  function(obj) {AddMetaData(obj, metadata = seurat_obj@meta.data[match(colnames(obj@raw.data), rownames(seurat_obj@meta.data)),])}
  # subsets <- safeParallel(fun = fun, 
  #                         args=args, 
  #                         seurat_obj=seurat_obj)
  
  names(subsets) = sNames_0
  
  # Save stats for outputting at the end
  n_cells_subsets = sapply(subsets, function(seurat_obj) ncol(seurat_obj@raw.data), simplify = T) %>% unlist(.,use.names = F)
  n_genes_subsets = sapply(subsets, function(seurat_obj) nrow(seurat_obj@raw.data), simplify=T) %>% unlist(.,use.names = F)
  
  # Free up space 
  rm(seurat_obj)
  
  ### Filter out cell types that have too few cells
  # We do this do avoid downstream problems with Seurat or WGCNA. 
  # E.g. Seurat will ScaleData will fail if regressing out variables when there are only 2 cells in the data.
  
  if (data_type=="sc") {
    message(paste0("Filtering out cell clusters with fewer than ", minCellClusterSize, " cells"))
    idx_cellcluster_ok <- n_cells_subsets>=minCellClusterSize
  } else if (data_type=="bulk") { # no need to filter out any subsets
    idx_cellcluster_ok = !logical(length = length(subsets))
  }
  
  if (!all(idx_cellcluster_ok)) {
    log_entry <- paste0(sNames_0[!idx_cellcluster_ok], ": had <", minCellClusterSize," cells and was therefore dropped")
    warning(log_entry)
    cat(text = log_entry, file =  path_cellClusters_dropped_log, append=T, sep = "\n")
    
  }
  # Filter
  subsets <- subsets[idx_cellcluster_ok]
  sNames_1 <- sNames_0[idx_cellcluster_ok]
  names(subsets) <- sNames_1 
  
  # Filter genes expressed in fewer than min.cells in a subset as these will also lead to spurious associations and computational difficulties
  if (data_type=="sc") {
    
    message(paste0("Filtering out genes expressed in fewer than ", min.cells, " cells"))
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "FilterGenes.txt")
    args = list("X"=subsets)
    fun = function(x) FilterGenes(x, min.cells = min.cells)
    subsets <- safeParallel(fun=fun, 
                            args=args, 
                            outfile=outfile, 
                            min.cells=min.cells)
  }
  
  
  # Filter
  idx_cellcluster_ok <- sapply(subsets, function(subset) {!is.null(rownames(subset@raw.data))}, simplify=T)
  if (!all(idx_cellcluster_ok)) {
    log_entry <-paste0(sNames_1[!idx_cellcluster_ok], ": had no genes left after filtering out genes expressed in fewer than ", min.cells, " cells, therefore dropped")
    warning(log_entry)
    cat(log_entry, file = path_cellClusters_dropped_log, append=T, sep = "\n")
    
  }
  # filter
  subsets <- subsets[idx_cellcluster_ok]
  sNames_2 <- sNames_1[idx_cellcluster_ok]
  
  if (data_type=="sc") {
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "NormalizeData.txt")
    args = list("X"=subsets)
    fun = function(seurat_obj) { NormalizeData(object = seurat_obj) }
    subsets <- safeParallel(fun=fun, 
                            args=args, 
                            outfile=outfile)
    
  }
  
  message("Scaling and regressing data subsets")
  
  vars.to.regress = if (!is.null(regress_out)) regress_out[regress_out %in% metadat_names] else NULL  
  
  # Scale data and regress out confounders
  subsets <- mapply(function(seurat_obj, name) {
    tryCatch({
      ScaleData(object = seurat_obj,
                vars.to.regress = if (data_type=="sc" & length(vars.to.regress)>0) vars.to.regress else NULL,
                model.use="linear",
                do.par=T, 
                num.cores = min(10, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = as.numeric(object.size(seurat_obj))/1024^3)-1),
                do.scale=T,
                do.center=do.center)}, 
      error = function(err) {
        message(paste0(name,": ScaleData failed with the error: ", err))
        seurat_obj
      })
  },
  seurat_obj = subsets,
  name = names(subsets),
  SIMPLIFY=F)
  
  message("Finding variable genes")
  
  args = list("X"=subsets)
  
  fun <- function(seurat_obj) {
    FindVariableGenes(object = seurat_obj,
                      x.low.cutoff = 0.0125,
                      x.high.cutoff = 8,
                      y.cutoff = 1,
                      do.plot=F,
                      display.progress=F)}
  
  subsets <- safeParallel(fun=fun, 
                          args=args)
  
  if (grepl("PCA", genes_use, ignore.case = T)) {
    
    message("Performing Principal Component Analysis..")
    start_time <- Sys.time()
    
    args <- list("seurat_obj" = subsets, "name" = names(subsets))
    
    outfile <- paste0(log_dir, data_prefix, "_", run_prefix, "_PCA_log.txt")
    
    fun = function(seurat_obj, name) {
      seurat_tmp <- tryCatch({ 
        RunPCA(object = seurat_obj,
               pc.genes = if (pca_genes == 'all') rownames(seurat_obj@data) else if (pca_genes == "var.genes") seurat_obj@var.genes,
               pcs.compute = min(nPC_seurat, min(if (pca_genes == 'all') nrow(seurat_obj@data) else length(seurat_obj@var.genes) %/% 2, ncol(seurat_obj@data) %/% 2)),
               use.imputed = F, 
               weight.by.var = F,
               do.print = F,
               seed.use = randomSeed,
               maxit = maxit, # set to 500 as default
               fastpath = fastpath, 
               verbose=T) }, 
        error = function(err) {
          message(paste0(name, ": RunPCA's IRLBA algorithm failed with the error: ", err))
          message("Trying RunPCA with var.genes, half the number of components, double max iterations and fastpath == F")
          tryCatch({
            RunPCA(object = seurat_obj,
                   pc.genes = seurat_obj@var.genes,
                   pcs.compute = min(nPC_seurat, min(length(seurat_obj@var.genes) %/% 2, ncol(seurat_obj@data) %/% 2))%/%2,
                   use.imputed = F, 
                   weight.by.var = F,
                   do.print = F,
                   seed.use = randomSeed,
                   maxit = maxit*2, # set to 500 as default
                   fastpath = F,
                   verbose=T)
          }, error = function(err1) {
            message(paste0(name, ": RunPCA's IRLBA algorithm failed again with error: ", err1))
            message("Returning the original Seurat object with empty dr$pca slot")
            return(seurat_obj)}) 
        })
      if (pca_genes == "all") {
        seurat_tmp@dr$pca@gene.loadings.full <- seurat_tmp@dr$pca@gene.loadings
      }
      return(seurat_tmp)
    }
    
    subsets <- safeParallel(fun=fun, 
                            args=args, 
                            outfile=outfile, 
                            pca_genes=pca_genes, 
                            nPC_seurat=nPC_seurat,
                            randomSeed=randomSeed,
                            maxit=maxit,
                            fastpath=fastpath)
    
    end_time <- Sys.time()
    
    message(sprintf("PCA done, time elapsed: %s seconds", round(end_time - start_time,2)))
    
    # Select significant PCs using empirical p-value based on JackStraw resampling
    # Source: https://rdrr.io/cran/Seurat/src/R/plotting.R
    # score the PCs using JackStraw resampling to get an empirical null distribution to get p-values for the PCs based on the p-values of gene loadings
    
    if (jackstrawnReplicate > 0 & genes_use == "PCA") message(sprintf("Performing JackStraw with %s replications to select genes that load on significant PCs", jackstrawnReplicate))
    
    list_datExpr <- mapply(function(seurat_obj, name) {
      tryCatch({
        wrapJackStraw(seurat_obj_sub = seurat_obj, 
                      n_cores = min(length(subsets), detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = as.numeric(max(sapply(subsets, object.size)))/1024^3)-1), 
                      jackstrawnReplicate = jackstrawnReplicate, 
                      pvalThreshold = pvalThreshold,
                      n_genes_use=n_genes_use)
      }, error = function(err) {
        message(paste0(name, ": PCA gene selection (with n Jackstraw resampling == ", jackstrawnReplicate, " failed with error: ", err))
        message("using var.genes instead")
        datExpr <- seurat_to_datExpr(seurat_obj_sub = seurat_obj, idx_genes_use = rownames(seurat_obj@scale.data) %in% seurat_obj@var.genes)
      })
    },
    seurat_obj = subsets, 
    name = names(subsets),
    SIMPLIFY=F)
    
    # stopCluster(cl)
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
  } else if (genes_use == "all") {
    list_datExpr <- lapply(subsets, function(x) 
      seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rep(TRUE, nrow(x@scale.data))))
  } else if (genes_use == "var.genes") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rownames(x@scale.data) %in% x@var.genes))
  } 
  
  # Clear up
  rm(subsets) #, builtInParams)
  invisible(gc())
  
  ######################################################################
  ########################### CHECKPOINT ###############################
  ######################################################################
  
  # Save or load session image 
  resume = "checkpoint_1"
  if (autosave) save.image(file=sprintf("%s%s_%s_checkpoint_1_image.RData.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")
  
} 

if (resume == "checkpoint_1") {
  
  ######################################################################
  ###################### (UN) LOAD PACKAGES ############################
  ######################################################################
  
  # Unload packages
  try(detach("package:Seurat", unload=TRUE))
  
  # Free up DLLs
  invisible(R.utils::gcDLLs())
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("WGCNA"))
  
  try(disableWGCNAThreads())
  
  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  message("Computing soft threshold powers to maximise the fit of a scale free topology to the adjacency matrix")
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  softPower <- 8 # Set a default value as fall back
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_compute_softPowers.txt")
  args = list("datExpr" = list_datExpr, "subsetName"=sNames_2)
  fun = function(datExpr, subsetName) {
    powers = c(1:30)
    sft = pickSoftThreshold(data=datExpr,
                            powerVector = powers,
                            blockSize = min(maxBlockSize , ncol(datExpr)), #try to prevent crashing
                            corFnc = corFnc,
                            corOptions =  corOptions,
                            networkType = networkType,
                            verbose = 1)
    pdf(sprintf("%s%s_%s_%s_pickSoftThresholdSFTFit.pdf", plots_dir, data_prefix, run_prefix, subsetName),width=10,height=5)
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",
         type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of 0.9
    abline(h=0.90,col="red")
    dev.off()
    # Mean connectivity as a function of the soft-thresholding power
    pdf(sprintf("%s%s_%s_%s_pickSoftThresholdMeanCon.pdf", plots_dir, data_prefix, run_prefix, subsetName),width=10,height=5)
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",
         ylab="Mean Connectivity",
         type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
    # Select softPower: of lower .95 percentile connectivity, if several softPowers achieve 0.9 R.sq, 
    # take smallest softPower; else take best fitting one
    
    fitIndices <- as.data.frame(sft$fitIndices)
    fitIndices %>% dplyr::filter(median.k. <= quantile(median.k.,0.95, na.rm=T)) -> fitIndices_filter 
    
    if (sum(fitIndices_filter$SFT.R.sq >= 0.93) > 1) {
      fitIndices_filter_2 <- fitIndices_filter[fitIndices_filter$SFT.R.sq>=0.93,] 
      fitIndices_filter_2[which.min(fitIndices_filter_2$Power), c(1,2,6)] -> df_sft
    } else {
      (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)))[1,c(1,2,6)]  -> df_sft
    }
    return(df_sft)
  }
  
  list_sft = safeParallel(fun=fun, 
                          args=args, 
                          outfile=outfile, 
                          maxBlockSize=maxBlockSize, 
                          corFnc=corFnc, 
                          corOptions=corOptions, 
                          networkType=networkType,
                          plots_dir=plots_dir, 
                          data_prefix=data_prefix, 
                          run_prefix=run_prefix)
  
  names(list_sft) <- sNames_2
  
  # TOM files will be saved here by consensusTOM / blockwiseConsensusModules
  setwd(scratch_dir)
  
  if (TOMnReplicate > 0) {
    
    ######################################################################
    ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
    ######################################################################
    
    message(sprintf("Resampling the expression data with fraction %s with replace = %s and %s replicates", fraction, replace, TOMnReplicate))
    
    message(sprintf("Computing the consensus Topological Overlap Matrix with %s permutations", TOMnReplicate))
    
    invisible(gc())
    
    fun <- function(datExpr,name) {
      # tryCatch({
        
        multiExpr <- bootstrap(datExpr=datExpr, 
                               nPermutations = TOMnReplicate,
                               replace = replace,
                               fraction = fraction,
                               randomSeed = randomSeed)
        
        consensusTOM(multiExpr = multiExpr, 
                     checkMissingData = checkMissingData,
                     maxBlockSize = maxBlockSize, 
                     blockSizePenaltyPower = blockSizePenaltyPower, 
                     randomSeed = randomSeed,
                     corType = corType,
                     maxPOutliers = maxPOutliers,
                     quickCor = quickCor,
                     pearsonFallback = pearsonFallback,
                     cosineCorrelation = cosineCorrelation,
                     replaceMissingAdjacencies = replaceMissingAdjacencies,
                     power = list_sft[[name]][["Power"]],
                     networkType = networkType,
                     TOMDenom = TOMDenom,
                     saveIndividualTOMs = saveIndividualTOMs,
                     individualTOMFileNames = paste0(data_prefix, "_", run_prefix, "_", name, "_individualTOM-Set%s-Block%b.RData"),
                     networkCalibration = networkCalibration,
                     sampleForCalibration = sampleForCalibration,
                     sampleForCalibrationFactor = sampleForCalibrationFactor,
                     getNetworkCalibrationSamples = getNetworkCalibrationSamples,
                     consensusQuantile = consensusQuantile,
                     useMean = useMean,
                     saveConsensusTOMs = saveConsensusTOMs,
                     consensusTOMFilePattern = paste0(data_prefix,"_", run_prefix, "_", name,"_consensusTOM-block.%b.RData"),
                     returnTOMs = F,
                     useDiskCache = T,
                     cacheDir = scratch_dir,
                     cacheBase = ".blockConsModsCache",
                     verbose = verbose,
                     indent = indent)
        
      # }, error = function(err) {
      #   adjacency = adjacency(datExpr=list_datExpr[[name]], 
      #                         type=type, 
      #                         power = list_sft[[name]]$Power, 
      #                         corFnc = corFnc, 
      #                         corOptions = corOptions)
      #   
      #   consTOMDS = TOMsimilarity(adjMat=adjacency,
      #                             TOMType=TOMType,
      #                             TOMDenom=TOMDenom,
      #                             verbose=verbose,
      #                             indent = indent)
      #   colnames(consTOMDS) <- rownames(consTOMDS) <- colnames(list_datExpr[[name]])
      #   save(consTOMDS, file=sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, run_prefix, name)) # Save TOM the way consensusTOM would have done
      # })
    }
    args <- list("datExpr"=list_datExpr, "name"=sNames_2)
    outfile = outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_consensusTOM.txt")
    
    additional_Gb = max(as.numeric(sapply(args, FUN = function(x) {object.size(x)*TOMnReplicate}, simplify = T)))/1024^3 
    obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
    n_cores <- min(max(sapply(args, length)), min(detectCores()%/%3, RAM_Gb_max %/% (obj_size_Gb + additional_Gb))-1)
    
    list_consensus <- safeParallel(fun=fun, 
                   args=args,
                   outfile=outfile,
                   n_cores=n_cores,
                   nPermutations=TOMnReplicate,
                   replace=replace,
                   fraction=fraction,
                   randomSeed=randomSeed,
                   checkMissingData = checkMissingData,
                   maxBlockSize = maxBlockSize, 
                   blockSizePenaltyPower = blockSizePenaltyPower, 
                   randomSeed = randomSeed,
                   corType = corType,
                   corFnc = corFnc, 
                   corOptions = corOptions,
                   maxPOutliers = maxPOutliers,
                   quickCor = quickCor,
                   pearsonFallback = pearsonFallback,
                   cosineCorrelation = cosineCorrelation,
                   replaceMissingAdjacencies = replaceMissingAdjacencies,
                   list_sft = list_sft,
                   networkType = networkType,
                   type = networkType,
                   TOMType=TOMType,
                   TOMDenom = TOMDenom,
                   saveIndividualTOMs = saveIndividualTOMs,
                   data_prefix=data_prefix, 
                   run_prefix=run_prefix,
                   networkCalibration = networkCalibration,
                   sampleForCalibration = sampleForCalibration,
                   sampleForCalibrationFactor = sampleForCalibrationFactor,
                   getNetworkCalibrationSamples = getNetworkCalibrationSamples,
                   consensusQuantile = consensusQuantile,
                   useMean = useMean,
                   saveConsensusTOMs = saveConsensusTOMs,
                   returnTOMs = F,
                   useDiskCache = T,
                   cacheDir = scratch_dir,
                   cacheBase = ".blockConsModsCache",
                   verbose = verbose,
                   indent = indent,
                   list_datExpr=list_datExpr,
                   scratch_dir = scratch_dir)
    
  } else if (TOMnReplicate==0) {
    
    message("Computing the Topological Overlap Matrix")
    
    ######################################################################
    ##### COMPUTE THE ADJACENCY MATRIX AND TOM WITHOUT RESAMPLING ########
    ######################################################################
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_TOM_for_par.txt")
    
    fun <- function(datExpr,name) { 
      adjacency = adjacency(datExpr=datExpr, 
                            type=type, 
                            power = list_sft[[name]]$Power, 
                            corFnc = corFnc, 
                            corOptions = corOptions)
      consTOMDS = TOMsimilarity(adjMat=adjacency,
                                TOMType=TOMType,
                                TOMDenom=TOMDenom,
                                verbose=verbose,
                                indent = indent)
      colnames(consTOMDS) <- rownames(consTOMDS) <- colnames(datExpr)
      save(consTOMDS, file=sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, run_prefix, name)) # Save TOM the way consensusTOM would have done
    }
    args = list("datExpr"=list_datExpr, "name"=sNames_2)
    
    safeParallel(fun=fun, args=args, outfile=outfile, type=type, list_sft=list_sft, corFnc=corFnc, corOptions=corOptions, TOMType=TOMType, TOMDenom=TOMDenom, verbose=verbose, indent=indent, data_prefix=data_prefix, run_prefix=run_prefix, scratch_dir=scratch_dir)
  }
  
  #rm(list_datExpr)
  
  #if (TOMnReplicate > 0) rm(list_multiExpr) #TODO: we'll need it for plotting permuted
  
  # Load TOMs into memory and convert to distance matrices
  fun = function(name) {
    load_obj(sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, run_prefix, name))
    
  }
  args = list("name"=sNames_2)
  list_consTOM = safeParallel(fun=fun, 
                              args=args, 
                              scratch_dir=scratch_dir, 
                              data_prefix=data_prefix, 
                              run_prefix=run_prefix, 
                              load_obj=load_obj)
  
  # Filter datExpr to retain genes that were kept by goodSamplesGenesMS, called by consensusTOM 
  if (TOMnReplicate > 1) { 
    fun = function(datExpr, consensus) {
      datExpr[,consensus$goodSamplesAndGenes$goodGenes]
    }
    args=list(datExpr=list_datExpr, consensus=list_consensus)
    list_datExpr_gg <- safeParallel(fun=fun, args=args)
  } else {
    list_datExpr_gg <- list_datExpr
  }
  rm(list_datExpr)
  
  # Convert proximity TOMs to distance matrices 
  args= list("X"=list_consTOM)
  fun = function(consTOM) { 
    1-as.dist(consTOM) 
  }
  list_dissTOM <- safeParallel(fun=fun, args=args)  
  rm(list_consTOM)
  # Save list_dissTOM to harddisk
  names(list_dissTOM) <- sNames_2
  
  saveRDS(list_dissTOM, file=sprintf("%s%s_%s_list_dissTOM.rds.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")
  
  ######################################################################
  ######################### CLEAR UP TOMS IN MEMORY AND HD #############
  ######################################################################
  
  rm(list_dissTOM)
  
  for (subsetName in sNames_2) {
    if (file.exists(sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, run_prefix, subsetName))) {
      try(file.remove(sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, run_prefix, subsetName)))
    }
    # Delete individual TOMs from disk
    indivTOM_paths <- dir(path=scratch_dir, pattern=paste0(data_prefix, "_", run_prefix, "_", subsetName, "_individualTOM"), full.names = T)
    for (indivTOM in indivTOM_paths) {
      try(file.remove(indivTOM))
    }
  }
  
  ######################################################################
  ############################ CHECKPOINT ##############################
  ######################################################################
  
  resume="checkpoint_2"
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_2_image.RData.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")
  
  
} 

if (resume == "checkpoint_2") {
  
  ######################################################################
  #################### (UN) LOAD PACKAGES #############################
  ######################################################################
  
  # Free up DLLs
  invisible(R.utils::gcDLLs())
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("WGCNA"))

  ######################################################################
  ##################### LOAD AND FILTER DISSTOMS #######################
  ######################################################################
  # TODO: Could filter earlier
  
  list_dissTOM_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_", run_prefix, "_list_dissTOM"), full.names = T)
  list_dissTOM <- load_obj(list_dissTOM_path)
  
  # Remove NULL
  idx_cellcluster_ok <- !sapply(list_dissTOM, is.null)
  
  if (!all(idx_cellcluster_ok)) {
    log_entry <-paste0(sNames_2[!idx_cellcluster_ok], ": Topological Overlap Matrix step failed, dropped from analysis")
    warning(log_entry)
    cat(log_entry, file = path_cellClusters_dropped_log, append=T, sep = "\n")
    
  }
  # Edit objects
  list_dissTOM <- list_dissTOM[idx_cellcluster_ok]
  list_datExpr_gg <- list_datExpr_gg[idx_cellcluster_ok]  
  sNames_3 <- sNames_2[idx_cellcluster_ok]
  
  ######################################################################
  ####################### CLUSTER ON THE TOM ###########################
  ######################################################################
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_parHclust.txt")
  fun <- function(dissTOM) {hclust(d=dissTOM, method=hclustMethod)}
  args = list("X"=list_dissTOM)
  
  list_geneTree <- safeParallel(fun=fun, args=args, outfile = outfile, hclustMethod=hclustMethod)
  names(list_geneTree) = sNames_3 # used for PlotDendro
  
  ######################################################################
  ######################### COMPUTE MODULES ############################
  ######################################################################
  
  message("Computing modules on the Topological Overlap Matrix")
  
  dims = sapply(list(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight), length) # list of number of values per parameter
  n_combs <- prod(dims) # Total number of combinations of parameter values
  list_comb = vector("list", length = n_combs)
  
  # Make a vector of different parameter combinations
  k=1
  for (size in 1:dims[1]) {
    for (ds in 1:dims[2]) { # Gandal et al 2018 use c(50,100, 200)
      for (pam in 1:dims[3]) {
        for (dthresh in 1:dims[4]) {
          list_comb[[k]] <- c(minClusterSize[size],deepSplit[ds], pamStage[pam], moduleMergeCutHeight[dthresh])
          k = k+1
        }
      }
    }
  }
  
  list_list_comb <- lapply(1:length(sNames_3), function(i) return(list_comb)) # just multiply..
  
  attr(x=list_list_comb, which = "names") <- sNames_3 # for some reason, just setting names(list_list_comb) <- sNames_3 filled the list with NULLs...
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_cutreeHybrid_for_vec.txt")
  
  fun <- function(geneTree,dissTOM) {
    lapply(list_comb, function(comb) {
      tryCatch({
        cutreeHybrid_for_vec(comb=comb, 
                             geneTree=geneTree, 
                             dissTOM = dissTOM,
                             maxPamDist = maxPamDist, 
                             useMedoids = useMedoids)
      }, error= function(err) {
        message(paste0(comb), ": cutreeHybrid failed")
        return(NULL)
      })
    }
    )}
  
  args = list("geneTree"=list_geneTree, "dissTOM"=list_dissTOM)
  
  list_list_cutree <- safeParallel(fun=fun, args=args, outfile=outfile, list_comb=list_comb, maxPamDist=maxPamDist, useMedoids=useMedoids)
  
  # free up memory
  if (fuzzyModMembership=="kME") rm(list_dissTOM)
  
  # Produce labels for plotting the modules found by different parameters
  list_plot_label <- lapply(list_comb, function(comb) plotLabel_for_vec(comb=comb)) # list of labels for plot
  # Make a copy of the labels for each subset
  list_list_plot_label = list()
  list_list_plot_label <- lapply(sNames_3, function(name) list_list_plot_label[[name]] = list_plot_label)
  
  # Remove NULLs in nested lists where cutreeHybrid failed
  list_vec_idx.combcutreeNULL <- lapply(list_list_cutree, function(list_cutree) {
    sapply(list_cutree, is.null)
  })
  
  names(list_vec_idx.combcutreeNULL) = sNames_3
  
  # filter
  if (any(sapply(list_vec_idx.combcutreeNULL, function(idx) idx==T))) {
    
    list_list_cutree <- mapply(function(list_cutree, vec_idx.combcutreeNULL) {
      list_cutree[!vec_idx.combcutreeNULL]
    }, list_cutree=list_list_cutree, vec_idx.combcutreeNULL = list_vec_idx.combcutreeNULL, SIMPLIFY=F)
    
    list_list_comb <- mapply(function(list_comb, vec_idx.combcutreeNULL) {
      list_comb[!vec_idx.combcutreeNULL]
    }, list_comb=list_list_comb, vec_idx.combcutreeNULL = list_vec_idx.combcutreeNULL, SIMPLIFY=F)
    
    list_list_plot_label <- mapply(function(list_plot_label, vec_idx.combcutreeNULL) {
      list_plot_label[!vec_idx.combcutreeNULL]
    }, list_plot_label = list_list_plot_label, vec_idx.combcutreeNULL = list_vec_idx.combcutreeNULL, SIMPLIFY=F)
  }
  
  names(list_list_cutree) <- names(list_list_comb) <- names(list_list_plot_label) <- sNames_3
  
  # Filter at the cell cluster level
  idx_cellcluster_ok <- sapply(list_list_cutree, function(x) length(x)>0)
  
  if (!all(idx_cellcluster_ok)) {
    log_entry <-paste0(sNames_3[!idx_cellcluster_ok], " : cutreeHybrid failed, removed from analysis")
    warning(log_entry)
    cat(log_entry, file = path_cellClusters_dropped_log, append=T, sep = "\n")
    
  }
  list_list_cutree <- list_list_cutree[idx_cellcluster_ok]# Filter(f=length, x= list_list_cutree)
  list_list_comb <-list_list_comb[idx_cellcluster_ok] #Filter(f=length, x=list_list_comb)
  list_list_plot_label <-list_list_plot_label[idx_cellcluster_ok]# Filter(f=length, x=list_list_plot_label)
  list_datExpr_gg <- list_datExpr_gg[idx_cellcluster_ok]#list_datExpr_gg[names(list_datExpr_gg) %in% sNames_4]
  
  sNames_4 <- sNames_3[idx_cellcluster_ok]# names(list_list_cutree) 
  ######################################################################
  ######################### MERGE CLOSE MODULES ########################
  ######################################################################
  
  message(paste0("Merging modules with a distance below ", moduleMergeCutHeight))
  
  if (fuzzyModMembership=="kME") {
    
    fun = function(list_comb,list_cutree, datExpr, cellType) {
      mapply(function(comb, cutree) {
        colors <- rep("grey", times=ncol(datExpr))
        MEs <- NULL 
        out <- list("colors"= colors, "MEs" = MEs)
        if (length(unique(cutree$labels))>1) {
          if (length(unique(cutree$labels))>2) {
            merged <- try({mergeCloseModules(exprData=as.matrix(datExpr), 
                                             colors = cutree$labels, 
                                             impute =T,
                                             corFnc = corFnc,
                                             corOptions = corOptions,
                                             cutHeight = comb[4],
                                             iterate = T,
                                             getNewMEs = F,
                                             getNewUnassdME = F)
            })
            if (class(merged)=="try-error") {
              colors = rep("grey", times=ncol(datExpr))
            } else {
              colors = labels2colors(merged$colors)
            }
          } else {
            warning(paste0(cellType, ": Only one proper module found, nothing to merge"))
            colors = labels2colors(cutree$labels) 
          }
          MEs <- try(moduleEigengenes_uv(expr = as.matrix(datExpr),
                                         colors=colors,
                                         excludeGrey = T))
          if (class(MEs)=="try-error") MEs <- NULL
        } else {
          warning(paste0(cellType, ": No modules found"))
        }
        out <- list("colors" = colors, "MEs"= MEs)
        out
      },
      cutree = list_cutree,
      comb = list_comb,
      SIMPLIFY=F)}
    
    args = list(list_comb=list_list_comb,
                list_cutree = list_list_cutree,
                datExpr = list_datExpr_gg,
                cellType = names(list_datExpr_gg))
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_mergeCloseModules.txt")
    
    list_list_merged <- safeParallel(fun=fun, args=args, outfile=outfile, corFnc = corFnc,corOptions=corOptions)
    
    names(list_list_merged) = sNames_4
    
    # get colors
    fun = function(list_merged,datExpr,cellType) {
      list_colors = lapply(list_merged, function(merged) {
        colors <- merged$colors
        names(colors) = colnames(datExpr)
        colors
      }) # list of merged colors
    } 
    
    args=list("list_merged"=list_list_merged, 
              datExpr=list_datExpr_gg, 
              cellType=sNames_4)
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "getColors.txt")
    
    # Extract the colors from the list returned by mergeCloseModules
    list_list_colors <- safeParallel(fun=fun, args=args, outfile=outfile)
    names(list_list_colors) <- sNames_4
    
    # Extract the Module Eigengenes from the list returned by mergeCloseModules
    fun = function(list_merged) {lapply(list_merged, function(merged) {
      if (!is.null(merged$MEs)) merged$MEs$eigengenes else NULL})}
    args = list("X"=list_list_merged)
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "getMEs.txt")
    list_list_MEs <- safeParallel(fun=fun, args=args, outfile=outfile)
    names(list_list_MEs) <- sNames_4
    
    # Compute kMEs 
    fun = function(list_MEs, datExpr) {
      lapply(list_MEs, function(MEs) {
        kMEs <- NULL
        if (!is.null(MEs)) {
          kMEs <- signedKME(datExpr=as.matrix(datExpr),
                            datME=MEs,
                            outputColumnName="",
                            corFnc = corFnc )
          kMEs
        } else NULL
      })}
    args = list(list_MEs=list_list_MEs,  
                datExpr=list_datExpr_gg)
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "compute_kMEs.txt")
    
    list_list_kMs <- safeParallel(fun=fun, args=args, outfile=outfile)
    
    
  } else if (fuzzyModMembership=="kIM") { 
    
    fun = function(x) lapply(x, function(y) labels2colors(y$labels))
    args = list("X"=list_list_cutree)
    list_list_colors <- safeParallel(fun=fun, args=args)
    
    fun =  function(list_colors,datExpr) lapply(list_colors, function(colors) name_for_vec(to_be_named=colors, given_names = colnames(datExpr), dimension=NULL))
    args = list("list_colors"=list_list_colors, "datExpr"=list_datExpr_gg)
    list_list_colors <- safeParallel(fun=fun, args=args)
    
    
    fun = function(list_colors,list_plot_label) name_for_vec(to_be_named=list_colors, given_names = as.character(list_plot_label), dimension=NULL)
    args = list(list_colors = list_list_colors,
                list_plot_label = list_list_plot_label)
    list_list_colors <- safeParallel(fun=fun, args=args)
    
    
    names(list_list_colors) <- sNames_4
    
    # Merge close modules using correlation between IM embeddings
    fun = function(dissTOM, list_colors) lapply(list_colors, function(colors) kIM_eachMod_norm(dissTOM = dissTOM, 
                                                                                               colors = colors,
                                                                                               verbose=verbose,
                                                                                               excludeGrey=F,
                                                                                               do.par=F))
    args = list(dissTOM = list_dissTOM,
                list_colors = list_list_colors)
    list_list_kMs <- safeParallel(fun = fun, args=args)
    
    
    
    fun = function(list_colors,list_kMs,datExpr,dissTOM, cellType) mapply( function(colors,kMs) mergeCloseModskIM(datExpr=datExpr,
                                                                                                                  colors = colors,
                                                                                                                  kIMs  = kMs,
                                                                                                                  dissTOM = dissTOM,
                                                                                                                  moduleMergeCutHeight=moduleMergeCutHeight,
                                                                                                                  verbose=verbose,
                                                                                                                  cellType = cellType),
                                                                           colors = list_colors,
                                                                           kMs = list_kMs,
                                                                           SIMPLIFY=F)
    
    args = list(list_colors = list_list_colors,
                list_kMs = list_list_kMs,
                datExpr = list_datExpr_gg,
                dissTOM = list_dissTOM,
                cellType = names(list_list_colors))
    list_list_merged <- safeParallel(fun=fun, args=args)
    
    
    fun = function(cellType) {tryCatch({lapply(list_list_merged[[cellType]], function(y) y$colors)}, 
                                       error = function(c) {
                                         warning(paste0(cellType, " failed"))
                                       })}
    args = list("X"=names(list_list_merged))
    list_list_colors <- safeParallel(fun=fun, args=args)
    
    fun = function(x) lapply(x, function(y) y[['colors']])
    args= list("X"=list_list_merged)
    list_list_colors <- safeParallel(fun=fun, args=args)
    
    fun = function(x) lapply(x, function(y) y[['kIMs']])
    args= list("X"=list_list_merged)
    list_list_kMs <- safeParallel(fun=fun, args=args)
    
    list_list_MEs <- NULL
    
    
  }
  names(list_list_colors) <- names(list_list_kMs) <- sNames_4
  
  ################################################  ######################
  ###################### REASSIGN GENES BASED ON kM ####################
  ######################################################################
  # see https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6
  
  if (kM_reassign) {
    message(paste0("Reassigning genes based on ", fuzzyModMembership))
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_par_kM_reassign.txt")
    
    if (fuzzyModMembership=="kIM") {
      fun = function(list_colors,
                     dissTOM,
                     datExpr,
                     cellType) lapply(list_colors, function(cols) kM_reassign_fnc(cols = cols,
                                                                                  fuzzyModMembership = fuzzyModMembership,
                                                                                  dissTOM = if (fuzzyModMembership=="kIM") dissTOM else NULL,
                                                                                  datExpr = if (fuzzyModMembership=="kME") datExpr else NULL,
                                                                                  corFnc = if (fuzzyModMembership=="kME") corFnc else NULL,
                                                                                  verbose=verbose,
                                                                                  max_iter = 3,
                                                                                  cellType = cellType)) 
      args = list(list_colors = list_list_colors,
                  dissTOM = list_dissTOM,
                  datExpr = list_datExpr_gg,
                  cellType = names(list_list_colors))
      
    } else if (fuzzyModMembership=="kME"){ # to avoid storing list_dissTOM
      fun = function(list_colors,
                     datExpr,
                     cellType) lapply(list_colors, function(cols) kM_reassign_fnc(cols = cols,
                                                                                  fuzzyModMembership = fuzzyModMembership,
                                                                                  #dissTOM = if (fuzzyModMembership=="kIM") dissTOM else NULL,
                                                                                  datExpr = if (fuzzyModMembership=="kME") datExpr else NULL,
                                                                                  corFnc = if (fuzzyModMembership=="kME") corFnc else NULL,
                                                                                  verbose=verbose,
                                                                                  max_iter = 3,
                                                                                  cellType = cellType))
      args = list(list_colors = list_list_colors,
                  datExpr = list_datExpr_gg,
                  cellType = names(list_list_colors))
      
    }
    
    
    list_list_reassign = safeParallel(fun=fun, args=args, outfile=outfile, fuzzyModMembership=fuzzyModMembership, corFnc=corFnc, verbose=verbose)
    
    
    # Extract new colors and kMs from the list of lists of lists returned by the vectorised kM_reassign 
    list_list_colors_reassign <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$colors))
    list_list_kMs_reassign <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$kMs))
    list_list_reassign_log <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$log))
    
    #invisible(gc()); invisible(R.utils::gcDLLs())
  } else {
    list_list_colors_reassign <- list_list_colors
    list_list_kMs_reassign <- list_list_kMs
    list_list_reassign_log <- NULL
  }  
  
  ######################################################################
  #################### EXTRACT PRIMARY kMEs / kIMs #####################
  ######################################################################
  # Extract the primary kMs - i.e. those of each gene w.r.t. the module to which it belongs
  # use these for only submitting relatively 'core' genes for PPI validation
  
  message(paste0("Computing primary "), fuzzyModMembership, "s")
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_par_primary_kMs_for_PPI.txt")
  fun = function(list_kMs,
                 list_colors) mapply(function(kMs, colors) {
                   if (length(unique(colors))>1) {
                     pkMs_fnc(kMs=kMs, colors=colors)  
                   } else {
                     pkMs<-numeric(length=length(colors))
                     names(pkMs) <- names(colors)
                     pkMs}
                 }, 
                 kMs = list_kMs, 
                 colors=list_colors, SIMPLIFY=F)
  
  args=list(list_kMs = list_list_kMs_reassign, 
            list_colors = list_list_colors_reassign)
  
  list_list_pkMs <- safeParallel(fun=fun, args=args, outfile=outfile)
  
  
  
  # The pkM vectors already have gene names. Now name each vector, and each list of vectors
  fun = function(list_pkMs,list_plot_label) name_for_vec(to_be_named = list_pkMs, 
                                                         given_names = as.character(list_plot_label), 
                                                         dimension = NULL)
  args = list(list_pkMs=list_list_pkMs, 
              list_plot_label=list_list_plot_label)
  list_list_pkMs <- safeParallel(fun=fun, args=args)
  
  names(list_list_pkMs) <- sNames_4
  
  
  ######################################################################
  ##### FILTER OUT GENES WITH INSIGNIFICANT GENE-MOD CORRELATION #######
  ######################################################################
  # Module expression as defined by ME - Module Eigengene or IM - IntraModular connectivity embedding
  
  if (kM_signif_filter) {
    
    message("Computing gene-module correlation t-test p-values")
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_par_t.test.txt")
    
    if (fuzzyModMembership=="kME") {
      
      fun = function(list_pkMs, list_colors, list_comb, cellType)  {
        mapply(function(pkMs, colors, comb) {
          tryCatch({
            if (length(unique(colors))>1){
              vec_t <- ifelse(colors!="grey", (pkMs*sqrt(length(pkMs)-2))/sqrt(1-pkMs^2), 0) # compute t.stats, for a vector of rho
              vec_p.val <- ifelse(colors!="grey", 
                                  stats::pt(q=abs(vec_t),  
                                  lower.tail = F, 
                                  df = length(pkMs)-2) + stats::pt(q=-(abs(vec_t)),  lower.tail = T, df = length(pkMs)-2), 1) # compute p.values. We are using both tails
          
              vec_q.val <- p.adjust(p = vec_p.val, method = "fdr")
              vec_idx_signif <- vec_q.val < pvalThreshold # compute boolean significance vector. Set "grey" to TRUE so we don't count them
              vec_idx_signif[colors=="grey"] <- TRUE
              out <- data.frame("p.val" = vec_p.val, "q.val" = vec_q.val, "signif" = vec_idx_signif)
              out } 
            else {
              out <- NULL
            }
          }, 
          error = function(err) {
            message(paste0(cellType, ", ", comb, ": ", "Failed to compute gene p-values with error: ", err))  
            out <- data.frame("p.val" = rep(0, length(colors)), "q.val" = rep(0, length(colors)), "signif" = rep(TRUE, length(colors)))
            out
          })
        },
        pkMs=list_pkMs, 
        colors=list_colors, 
        comb = list_comb,
        SIMPLIFY=F)}
      
      args = list(list_pkMs = list_list_pkMs, 
                  list_colors = list_list_colors_reassign, 
                  list_comb = list_list_comb,
                  cellType = sNames_4)
      
      list_list_geneMod_t.test <- safeParallel(fun=fun, args=args, pvalThreshold=pvalThreshold) 
      
      
    } else if (fuzzyModMembership=="kIM") {
      
      fun = function(datExpr, list_colors, cellType, list_kMs, dissTOM, list_comb) { 
        mapply(function(colors, kMs, comb) {
          tryCatch({
            if (length(unique(colors))>1) {
              message(paste0(cellType, ": Computing module kIM embeddings"))
              cellModEmbed(datExpr=datExpr,
                           colors=colors,
                           latentGeneType="IM",
                           cellType=NULL, # to avoid having it in the embedding matrix column names
                           kMs = kMs,
                           dissTOM = dissTOM)
            } else {
              NULL
            }
          }, error = function(err) {
            message(paste0(cellType, ", ", comb, ": ", "failed to compute embedding matrix with error: ", err))
            NULL
          })
        },
        colors=list_colors,
        kMs = list_kMs,
        comb= list_comb,
        SIMPLIFY=F)
      }
      
      args  = list(datExpr = list_datExpr_gg,
                   list_colors = list_list_colors_reassign,
                   cellType = sNames_4,
                   list_kMs = list_list_kMs_reassign,
                   list_comb=list_list_comb,
                   dissTOM = list_dissTOM)
      
      list_list_embed_mat <- safeParallel(fun=fun, args=args)
      
      
      fun = function(datExpr, list_embed_mat, list_colors, list_comb, cellType) {
        mapply(function(embed_mat, colors, comb){
          tryCatch({
            if (length(unique(colors))>1) {
              message(paste0(cellType, ": Computing gene correlations with module kIM embeddings"))
              vec_p.val <- numeric(length=length(colors))
              names(vec_p.val) <- names(colors)
              vec_p.val[colors=="grey"] <- 0
              
              for (color in names(table(colors))[!grepl("^grey$",names(table(colors)))]) {
                vec_p.val[colors==color] <- apply(X = datExpr[,colors==color], MARGIN = 2, FUN = function(datExpr_col) cor.test(datExpr_col, embed_mat[, color,drop=F], 
                                                                                                                                alternative = "two.sided", 
                                                                                                                                method = "pearson", 
                                                                                                                                conf.level = 1-pvalThreshold)$p.value)
              }
              vec_q.val <- p.adjust(p = vec_p.val, method = "fdr")
              vec_idx_signif <- vec_q.val < pvalThreshold
              vec_idx_signif[colors=="grey"] <- TRUE# compute boolean significance vector. Set "grey" to TRUE so we don't count them
              out <- data.frame("p.val" = vec_p.val, "q.val" = vec_q.val, "signif" = vec_idx_signif)
              out } else out <- NULL
          }, error= function(err) {
            message(paste0(cellType, ", ", comb, ": ", "failed to compute t.test with error: ", err))
            out <- data.frame("p.val" = rep(0, length(colors)), "q.val" = rep(0, length(colors)), "signif" = rep(TRUE, length(colors)))
            out
          })
        },
        embed_mat=list_embed_mat,
        colors=list_colors,
        comb = list_comb,
        SIMPLIFY=F)}
      
      args = list(datExpr = list_datExpr_gg,
                  list_embed_mat = list_list_embed_mat,
                  list_colors = list_list_colors_reassign,
                  list_comb = list_list_comb,
                  cellType = sNames_4)
      
      list_list_geneMod_t.test <- safeParallel(fun=fun, args=args, outfile=outfile, pvalThreshold=pvalThreshold)
      
      
      rm(list_list_embed_mat)
    } 
    #stopCluster(cl)
  } else {
    list_list_geneMod_t.test <- NULL 
  }  
  
  if (fuzzyModMembership=="kIM") rm(list_dissTOM)
  
  ########################################################################################
  #################### REASSIGN NON-SIGNIFICANT GENES TO "GREY" ##########################
  ########################################################################################
  
  list_list_colors_t.test <- list_list_colors_reassign
  list_list_pkMs_t.test <- list_list_pkMs
  
  if(kM_signif_filter) if (!all(sapply(list_list_geneMod_t.test, 
                                       function(list_geneMod_t.test) {
                                         sapply(list_geneMod_t.test, function(geneMod_t.test) {
                                           if (!is.null(geneMod_t.test)) all(geneMod_t.test$signif) else T}
                                           , simplify=T)}, simplify=T))) {
    message("Filtering out genes with insignificant gene-module correlation p-values in expression profile")
    
    # filter colors
    fun = function(list_colors_reassign, list_geneMod_t.test) {
      mapply(function(colors_reassign, geneMod_t.test) {
        
        if (!is.null(geneMod_t.test)) {
          colors_t.test <- colors_reassign
          colors_t.test[!geneMod_t.test$signif] <- rep("grey", times=sum(!geneMod_t.test$signif))
          return(colors_t.test)
        } else {
          colors_reassign
        }
      }, colors_reassign = list_colors_reassign, 
      geneMod_t.test = list_geneMod_t.test, 
      SIMPLIFY=F)
    }
    args = list(list_colors_reassign = list_list_colors_reassign, 
                list_geneMod_t.test = list_list_geneMod_t.test)
    list_list_colors_t.test <- safeParallel(fun=fun, args=args)
    
    fun = function(list_pkMs, list_geneMod_t.test) {
      
      mapply(function(pkMs, geneMod_t.test) {
        
        if (!is.null(geneMod_t.test)) {
          
          pkMs_t.test <- pkMs
          pkMs_t.test[!geneMod_t.test$signif] <- rep(0, times=sum(!geneMod_t.test$signif))
          return(pkMs_t.test)
        } else {
          pkMs 
        }
      }, 
      pkMs = list_pkMs, 
      geneMod_t.test = list_geneMod_t.test, 
      SIMPLIFY=F)
    }
    args = list(list_pkMs = list_list_pkMs, 
                list_geneMod_t.test = list_list_geneMod_t.test)
    list_list_pkMs_t.test <- safeParallel(fun=fun, args=args)
    
  } 
  
  ######################################################################
  ################### Filter out very small modules ####################
  ######################################################################
  
  if(kM_signif_filter) if (!all(sapply(list_list_geneMod_t.test, 
                                       function(list_geneMod_t.test) {
                                         sapply(list_geneMod_t.test, function(geneMod_t.test) {
                                           if (!is.null(geneMod_t.test)) all(geneMod_t.test$signif) else T}
                                           , simplify=T)}, simplify=T))) {
    
    message(paste0("Filtering out modules with fewer than ", minClusterSize, " genes left after t.testing"))
    
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_par_minClusterSize_filter.txt")
    
    # filter pkMs
    fun = function(list_pkMs, list_colors) {
      mapply(function(pkMs, colors) {
        mods_too_small <- names(table(colors))[table(colors) < minClusterSize]
        pkMs[colors %in% mods_too_small] <- rep(0, times=sum(colors %in% mods_too_small))
        return(pkMs)
      }, pkMs = list_pkMs, 
      colors = list_colors, 
      SIMPLIFY=F)
    }
    args= list(list_pkMs = list_list_pkMs_t.test, 
               list_colors = list_list_colors_t.test)
    list_list_pkMs_t.test <- safeParallel(fun=fun, args=args, outfile=outfile)
    
    
    # filter colors after pkMs as we needed colors to filter pkMs
    fun = function(list_colors) {
      lapply(list_colors, function(colors) {
        mods_too_small <- names(table(colors))[table(colors) < minClusterSize]
        colors[colors %in% mods_too_small] <- rep("grey", times=sum(colors %in% mods_too_small))
        colors
      })
    }
    args= list("X"=list_list_colors_t.test)
    list_list_colors_t.test <- safeParallel(fun=fun, args=args, outfile=outfile)
    
  } 
  
  ######################################################################
  ############## Match colors between parameter settings ###############
  ######################################################################
  # This step is purely for plotting. 
  # For each celltype, we will anyway take only modules found using one set of parameters
  
  message("Matching module color labels between parameter settings")
  
  # Align / match the colors so similar modules found with different sets of parameters have the same name
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_parMatchColors.txt")
  fun = function(x) parMatchColors(list_colors = x)
  args = list("X"=list_list_colors_t.test)
  list_list_colors_matched <- safeParallel(fun=fun, args=args, outfile=outfile)
  
  # Rename highest level list entries (cell clusters)
  names(list_list_colors_matched) = sNames_4
  
  # Name each entry of the vectors of color assignments with the corresponding genes
  fun =  function(list_colors_matched,datExpr) lapply(list_colors_matched, 
                                                      function(colors_matched) name_for_vec(to_be_named=colors_matched, 
                                                                                            given_names=colnames(datExpr), 
                                                                                            dimension = NULL))
  args = list(list_colors_matched=list_list_colors_matched, 
              datExpr=list_datExpr_gg)
  list_list_colors_matched <- safeParallel(fun=fun, args=args, outfile=outfile)
  
  ######################################################################
  ######################### REMOVE ALL-GREY RUNS #######################
  ######################################################################
  
  message("Removing all-grey runs")
  
  # count the grey
  list_vec_n_grey <- lapply(X=list_list_colors_matched, 
                            FUN = function(list_colors) sapply(X = list_colors, 
                                                               FUN=function(colors) sum(as.numeric(colors=="grey")), simplify = T))
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  list_logical_params_ok <- mapply(function(vec_n_grey,list_colors_matched) {
    as.logical(mapply(function(n_grey,colors_matched) {
      n_grey!=length(colors_matched)
    }, 
    n_grey=vec_n_grey, 
    colors_matched=list_colors_matched, SIMPLIFY=F))
  }, 
  vec_n_grey=list_vec_n_grey, 
  list_colors_matched=list_list_colors_matched, 
  SIMPLIFY=F)
  
  idx_cellcluster_ok <- sapply(list_logical_params_ok, any, simplify = T)
  
  # First, for each run, remove results of any parametrisation for which all the genes were assigned to grey
  # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_list_plot_label_ok <- mapply(function(list_plot_label,logical_params_ok) list_plot_label[logical_params_ok], 
                                    list_plot_label = list_list_plot_label, 
                                    logical_params_ok = list_logical_params_ok, 
                                    SIMPLIFY = F)
  
  list_list_cutree_ok <- mapply(function(list_cutree,logical_params_ok) list_cutree[logical_params_ok], 
                                list_cutree = list_list_cutree,
                                logical_params_ok = list_logical_params_ok,
                                SIMPLIFY = F)
  
  list_list_colors_matched_ok <- mapply(function(list_colors_matched,logical_params_ok) {
    list_colors_matched[logical_params_ok]
  }, 
  list_colors_matched=list_list_colors_matched, 
  logical_params_ok=list_logical_params_ok, 
  SIMPLIFY = F)
  
  list_list_pkMs_ok <- mapply(function(x,y) x[y], x=list_list_pkMs_t.test, y=list_logical_params_ok, SIMPLIFY = F)
  list_list_reassign_log_ok <- if (kM_reassign) mapply(function(x,y) x[y], x=list_list_reassign_log, y=list_logical_params_ok, SIMPLIFY = F)  else NULL
  
  # Make a warning and report if any whole subsets produced no modules
  if (!all(idx_cellcluster_ok)) {
    log_entry <- paste0(sNames_4[idx_cellcluster_ok], " had no modules, dropped from the analysis")
    warning(log_entry)
    cat(log_entry, file = path_cellClusters_dropped_log, append=T, sep = "\n")
    
  }
  # filter
  list_list_plot_label_ok <- list_list_plot_label[idx_cellcluster_ok]
  list_list_cutree_ok <- list_list_cutree[idx_cellcluster_ok]
  list_datExpr_ok <- list_datExpr_gg[idx_cellcluster_ok] # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_geneTree_ok <- list_geneTree[idx_cellcluster_ok]
  list_list_colors_matched_ok <- list_list_colors_matched_ok[idx_cellcluster_ok]
  list_list_pkMs_ok <- list_list_pkMs_ok[idx_cellcluster_ok]
  list_list_reassign_log_ok <- if (kM_reassign) list_list_reassign_log_ok[idx_cellcluster_ok] else NULL 
  
  sNames_5 <- sNames_4[idx_cellcluster_ok]
  
  # Assign gene names to each color vector
  list_list_colors_matched_ok <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), 
                                                                                         dimension = NULL)), 
                                        x = list_list_colors_matched_ok, 
                                        y = list_datExpr_ok, SIMPLIFY=F)
  
  rm(list_list_plot_label, list_list_cutree, list_datExpr_gg, list_geneTree, list_list_colors_matched, list_list_pkMs_t.test, list_list_reassign_log, list_list_MEs)#, list_idx_mapped_gg)
  
  invisible(gc())
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_3"
  message("Reached checkpoint 3, saving session image")
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_3_image.RData.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")

} 

if (resume == "checkpoint_3") {
  
  ##########################################################################
  ############################ (UN)LOAD PACKAGES ############################
  ##########################################################################
  
  # Free up DLLs
  #try(detach("package:STRINGdb", unload=TRUE))
  invisible(R.utils::gcDLLs())
  # Load packages
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("STRINGdb"))
  
  ######################################################################
  #################### CHECK MODULES FOR PPI ENRICHMENT ################
  ######################################################################
  # in:
  #   list_list_colors_matched_ok
  #   list_list_pkMs
  #
  # out:
  #   list_list_colors_matched_ok_PPI
  
  message("Checking modules for significant Protein-Protein Interactions through STRINGdb")
  
  if (data_organism == "hsapiens") STRINGdb_species <- 9606 else if (data_organism == "mmusculus") STRINGdb_species <- 10090
  
  if (fuzzyModMembership=="kME") {
    # Value already set in parameters; just make it into a vector of identical values, one entry per celltype 
    list_PPI_pkM_threshold <- as.list(rep(PPI_pkME_threshold, times = length(list_list_pkMs_ok)))
  } else if (fuzzyModMembership=="kIM") {
    # If using pkIMs to filter out genes with low module membership from the lists submitted to the PPI analysis, for each
    # cellType subset (upper level of the list_list_pkMs), compute the 0.05 quantile of gene primary module membership values (pkIMs)
    # and set that as a threshold for that datatype. Code below returns a numeric vector of quantile values.
    list_PPI_pkM_threshold <- lapply(list_list_pkMs_ok, function(x) quantile(unlist(x, use.names = F), probs = c(0.05), names=F))
  }
  
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_PPI_outer_for_vec.txt")
  fun = function(list_colors,list_pkMs,PPI_pkME_threshold) {
    mapply(function(colors,pkMs) PPI_outer_for_vec(colors = colors,
                                                   pkMs = pkMs,
                                                   STRINGdb_species = STRINGdb_species,
                                                   PPI_pkM_threshold = PPI_pkME_threshold,
                                                   pvalThreshold = pvalThreshold),
           colors = list_colors,
           pkMs = list_pkMs,
           SIMPLIFY=F)
  }
  args = list(list_colors = list_list_colors_matched_ok,
              list_pkMs = list_list_pkMs_ok,
              PPI_pkME_threshold = list_PPI_pkM_threshold)
  list_list_PPI <- safeParallel(fun=fun, args=args, outfile=outfile, STRINGdb_species=STRINGdb_species)
  
  
  list_list_colors_PPI <- if (PPI_filter) lapply(list_list_PPI, function(x) lapply(x, function(y) y$colors_PPI)) else list_list_colors_matched_ok
  list_list_module_PPI <- lapply(list_list_PPI, function(x) lapply(x, function(y) data.frame(y$module_PPI, stringsAsFactors = F)))
  list_list_module_PPI_signif <- lapply(list_list_PPI, function(x) lapply(x, function(y) data.frame(y$module_PPI_signif, stringsAsFactors = F)))
  
  names(list_list_colors_PPI) <- names(list_list_module_PPI) <- names(list_list_module_PPI_signif) <- list_list_colors_matched_ok 
  
  # Name the PPI enrichment results as they will be output in this state; 
  # we will however soon remove a layer from the list_list_colors_PP double nested list of vector so no need to name them
  list_list_module_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension=NULL),
                                 x=list_list_module_PPI,
                                 y=list_list_colors_matched_ok, 
                                 SIMPLIFY = F)
  
  list_list_module_PPI_signif <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension=NULL),
                                        x=list_list_module_PPI_signif,
                                        y=list_list_colors_matched_ok, 
                                        SIMPLIFY = F)
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_4"
  message("Reached checkpoint 4, saving session image")
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_4_image.RData.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")

  
} 

if (resume == "checkpoint_4") {
  
  ######################################################################
  ######################### LOAD PACKAGES ##############################
  ######################################################################
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("WGCNA"))
  if (!is.null(list_genesets_path)) suppressPackageStartupMessages(library("liger"))
  if (!is.null(magma_gwas_dir)) suppressPackageStartupMessages(library("boot"))
  
  ######################################################################
  ################ ORDER PARAMETER SETS BY PPI ENRICHMENT ##############
  ######################################################################
  
  message("Selecting parameters with the highest number of genes assigned to modules significantly enriched for Protein-Protein Interactions")
  
  # Count how many genes were assigned to grey under each parametrisation and order the parametrisations
  list_PPI_vec_n_grey <- lapply(list_list_colors_PPI, count_grey_in_list_of_vec)
  
  # Order all the outputs by how many genes were assigned to a (non-grey) module
  if (kM_reassign) {
    list_list_reassign_log_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_reassign_log_ok, y = list_PPI_vec_n_grey, SIMPLIFY=F) 
  } else {
    list_list_reassign_log_order <- NULL
  }
  list_list_plot_label_ok_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_plot_label_ok, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_colors_matched_ok_order <- mapply(function(x,y) x[order(y, decreasing=F)], x  =  list_list_colors_matched_ok , y = list_PPI_vec_n_grey, SIMPLIFY = F )
  list_list_colors_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_colors_PPI, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_module_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_module_PPI, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_module_PPI_signif_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_module_PPI_signif, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_geneMod_t.test_order <- if (kM_signif_filter) mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_geneMod_t.test[names(list_list_geneMod_t.test) %in% sNames_5], y = list_PPI_vec_n_grey, SIMPLIFY=F) else NULL
  
  ######################################################################
  ##### FOR EACH SUBSET SELECT PARAMETERS WITH BEST PPI ENRICHMENT #####
  ######################################################################
  
  # Eliminate a layer of nesting by selecting only the best parametrisation per celltype
  list_plot_label_final <- lapply(list_list_plot_label_ok_order, function(x) x[[1]])
  list_reassign_log <- if (kM_reassign) lapply(list_list_reassign_log_order, function(x) x[[1]]) else NULL
  list_colors_PPI <- lapply(list_list_colors_PPI_order, function(x) x[[1]])
  list_colors_all <- lapply(list_list_colors_matched_ok_order, function(x) x[[1]])
  list_module_PPI <- lapply(list_list_module_PPI_order, function(x) x[[1]])
  list_module_PPI_signif <- lapply(list_list_module_PPI_signif_order, function(x) x[[1]])
  list_geneMod_t.test <- if (kM_signif_filter)  lapply(list_list_geneMod_t.test_order, function(x) x[[1]]) else NULL 
  
  # Name by cell clusters
  names(list_plot_label_final)  <-  names(list_colors_PPI) <- names(list_colors_all) <- names(list_module_PPI) <- names(list_module_PPI_signif)<- sNames_5
  if (kM_reassign) names(list_reassign_log) <- sNames_5
  if (kM_signif_filter)  names(list_geneMod_t.test) <- sNames_5
  
  # Make list of list of final parameters
  param_names = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  list_list_cutree_params_final <- lapply(list_PPI_vec_n_grey, function(x) list_comb[order(x,decreasing = F)][[1]])
  list_list_cutree_params_final <- lapply(list_list_cutree_params_final, function(x) name_for_vec(to_be_named = x, given_names = param_names, dimension = NULL)) 
  
  ######################################################################
  ################### DELETE RUNS WITH ONLY GREY #######################
  ######################################################################
  
  PPI_vec_n_grey <- count_grey_in_list_of_vec(list_colors_PPI)
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  logical_subsets_PPI_ok <- as.logical(mapply(function(x,y) x != length(y), x = PPI_vec_n_grey, y = list_colors_PPI, SIMPLIFY = F))
  
  sNames_PPI <- sNames_5[logical_subsets_PPI_ok]
  list_colors_PPI <- list_colors_PPI[logical_subsets_PPI_ok]
  list_module_PPI_ok <- list_module_PPI[logical_subsets_PPI_ok]
  list_module_PPI_signif_ok <- list_module_PPI_signif[logical_subsets_PPI_ok]
  list_datExpr_PPI <- list_datExpr_ok[logical_subsets_PPI_ok]
  
  list_colors_all_PPI_ok <- list_colors_all[logical_subsets_PPI_ok]
  list_list_cutree_params_final_PPI <- list_list_cutree_params_final[logical_subsets_PPI_ok]
  list_plot_label_final_PPI  <- list_plot_label_final[logical_subsets_PPI_ok]
  list_geneTree_PPI <- list_geneTree_ok[logical_subsets_PPI_ok]
  
  ######################################################################
  ##################### MAKE MODULE COLORS UNIQUE ######################
  ######################################################################
  
  message("Making module colors unique across cell clusters")
  
  # Get nested list of modules
  list_mods <- lapply(list_colors_PPI, function(x) names(table(x)))
  
  list_mods <- lapply(list_mods, function(mods) {
    if (any(grepl("^grey$", mods))) {
      mods <- mods[mods!="grey"] 
    } else {
      mods
    }
  })
  
  mods <- unlist(list_mods)
  names(mods) <- NULL
  
  all_cols_nogrey_uniq <- unique(gsub("\\d+", "", colors()[-grep("grey|gray", colors())]))
  all_cols_nogrey <- colors()[-grep("grey|gray", colors())]
  
  # Replace colors with new unique colors
  if (length(mods) <= length(all_cols_nogrey_uniq) ) { # if there are enough unique colors without adding numbers
    mods_uniq <- all_cols_nogrey_uniq[sample(x=1:length(all_cols_nogrey_uniq), size=length(mods), replace=F)]
  } else if (length(mods) > length(all_cols_nogrey_uniq) & length(mods) < length(all_cols_nogrey) ) { # if there aren't enough unique colors unless they have numbers added
    mods_uniq <- all_cols_nogrey[sample(x=1:length(all_cols_nogrey), size=length(mods), replace=F)]
  } else if (length(mods) > length(all_cols_nogrey)) { # if there aren't enough unique colors in R
    mods_uniq <- paste0(all_cols_nogrey_uniq, "_", 1:(length(mods)))
  }
  
  list_mods_uniq <- vector(mode="list", length=length(list_mods))
  
  k=1
  for (j in 1:length(list_mods)) {
    for (i in 1:length(list_mods[[j]])) {
      list_mods_uniq[[j]][i] <- mods_uniq[k]
      k = k+1
    }
  }
  
  names(list_mods_uniq) <- names(list_colors_PPI)
  
  # Update color vectors
  list_colors_PPI_uniq <- mapply(function(x,y,z) z[match(x, y)],
                                 x = list_colors_PPI, 
                                 y = list_mods, 
                                 z = list_mods_uniq, SIMPLIFY = F)
  
  # Update also the old module assignment colors, but only change those that were validated by PPI (and therefore among those replaced with unique colors)
  list_colors_uniq <-  mapply(function(x,y,z) ifelse(x==y, yes=z, no = x),
                              x = list_colors_all_PPI_ok, 
                              y = list_colors_PPI, 
                              z = list_colors_PPI_uniq,
                              SIMPLIFY = F)
  
  # Replace NAs (i.e. no match in non-grey modules) with "grey"
  list_colors_uniq <- lapply(list_colors_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  list_colors_PPI_uniq <- lapply(list_colors_PPI_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  
  # Give gene names to vectors
  list_colors_uniq <- mapply(function(x,y) name_for_vec(to_be_named= x, given_names = names(y), dimension=NULL), x = list_colors_uniq, y = list_colors_all_PPI_ok, SIMPLIFY=F)
  list_colors_PPI_uniq <- mapply(function(x,y) name_for_vec(to_be_named= x, given_names = names(y), dimension=NULL), x = list_colors_PPI_uniq, y = list_colors_PPI, SIMPLIFY=F)
  
  # Update list_module_PPI_ok colors
  list_module_PPI_uniq <-  mapply(function(colors, colors_uniq, module_PPI) {
    module_PPI$colors <- colors_uniq[match(module_PPI$colors, colors)]
    module_PPI
  },
  colors = list_colors_PPI, 
  colors_uniq = list_colors_PPI_uniq,
  module_PPI = list_module_PPI_ok,
  SIMPLIFY = F)
  
  list_module_PPI_signif_uniq <-  mapply(function(colors, colors_uniq, module_PPI_signif) {
    module_PPI_signif$colors <- colors_uniq[match(module_PPI_signif$colors, colors)]
    module_PPI_signif
  },
  colors = list_colors_PPI, 
  colors_uniq = list_colors_PPI_uniq,
  module_PPI_signif = list_module_PPI_signif_ok,
  SIMPLIFY = F)
  
  names(list_module_PPI_uniq) <- names(list_module_PPI_signif_uniq) <- sNames_PPI
  
  ######################################################################
  ########## RECOMPUTE MEs, kMs, pkMs AFTER PPI FILTER #########
  ######################################################################
  # TODO: Surely we can just remove the grey ones?
  # Yes but then we would have to remap colours, which would be time consuming
  
  message("Computing kMs after making colors unique and computing Protein-Protein Interactions")
  
  if (fuzzyModMembership=="kIM"){
    list_dissTOM_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_", run_prefix, "_list_dissTOM"), full.names = T)
    list_dissTOM <- load_obj(list_dissTOM_path)
    list_dissTOM_PPI <- list_dissTOM[match(sNames_PPI,names(list_dissTOM))]
    rm(list_dissTOM)
  }
  
  if (fuzzyModMembership == "kME") {
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_kMEs_PPI.txt")
    fun = function(x,y,cellType) {
      out <- try({
        moduleEigengenes_uv(expr = as.data.frame(x, col.names=col.names(x)),
                            colors=y,
                            excludeGrey=T)#,
        
      })
      if (class(out) == "try-error"){
        warning(paste0(cellType, ": moduleEigengenes failed with the error: ", out))
        out <- list(eigengenes=NULL, u = NULL)
      } else {
        out
      }
    }
    args = list(x = list_datExpr_PPI, 
                y = list_colors_PPI_uniq,
                cellType = names(list_datExpr_PPI))
    
    list_ModuleEigengenes_out_PPI <- safeParallel(fun=fun, args=args, outfile=outfile)
    
    list_MEs_PPI <- lapply(list_ModuleEigengenes_out_PPI, function(x) x$eigengenes)
    
    list_u_PPI <- lapply(list_ModuleEigengenes_out_PPI, function(x) x$u) # u is actually a list, since the left eigenvectors have different lengths depending on the number of genes in the module
    
    names(list_MEs_PPI) <- names(list_u_PPI) <- sNames_PPI
    
    fun = function(x,y) {
      if(!is.null(y)) {
        signedKME(datExpr = as.matrix(x),
                  datME = y,
                  outputColumnName = "",
                  corFnc = corFnc)
      } else {
        NULL  
      }
    }
    args = list(x=list_datExpr_PPI, 
                y=list_MEs_PPI)
    list_kMs_PPI <- safeParallel(fun=fun, args=args, outfile=outfile)
    
    
    
    # Remove 'ME' from eigengenes
    list_MEs_PPI <- lapply(list_MEs_PPI, function(x) {
      if (!is.null(x)) name_for_vec(to_be_named = x, given_names = gsub(pattern="ME" , replacement="", x = colnames(x), ignore.case = F), dimension = 2) else NULL
    })
    
    # add cell row names
    list_MEs_PPI <- mapply(function(x,y) {
      if (!is.null(x)) name_for_vec(to_be_named = x, given_names = rownames(y), dimension=1) else NULL
    }, 
    x = list_MEs_PPI,
    y = list_datExpr_PPI,
    SIMPLIFY=F)
    
  } else if (fuzzyModMembership == "kIM"){
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_kIMs_PPI.txt")
    fun = function(x,y) kIM_eachMod_norm(dissTOM = x, 
                                         colors = y,
                                         excludeGrey=T,
                                         do.par=F)
    args = list(x = list_dissTOM_PPI, 
                y = list_colors_PPI_uniq)
    list_kMs_PPI <- safeParallel(fun=fun, args=args, outfile=outfile)
    
    
    
    names(list_kMs_PPI) <- sNames_PPI
    
    list_MEs_PPI <- NULL
    
    # Output list of list of kM gene weights, one vector per module, i.e. u is a list!
    fun = function(kMs, colors) {
      
      u <- lapply(colnames(kMs), function(module) {
        out <- kMs[match(names(colors)[colors==module], rownames(kMs)),module]
        names(out) = rownames(kMs)[colors==module]
        return(out)
      })
      names(u) <- colnames(kMs)
      return(u)
    }
    args = list(kMs = list_kMs_PPI,
                colors = list_colors_PPI_uniq)
    list_u_PPI <- safeParallel(fun=fun, args=args, outfile=outfile)
    
    
    
    invisible(gc()); invisible(R.utils::gcDLLs())
  } 
  
  if (fuzzyModMembership == "kIM") rm(list_dissTOM_PPI)
  
  ##########################################################################
  ###################### COMPUTE GENE SET ENRICHMENT #######################
  ##########################################################################
  
  n_modules_GSEA_enriched <- rep(NA, times = length(sNames_0))
  names(n_modules_GSEA_enriched) <- sNames_0
  
  if (!is.null(list_genesets_path)) {
    
    message("Performing Gene Set Enrichment Analysis")
    
    list_genesets <- load_obj(f=list_genesets_path)
    
    if (map_genes_to_ensembl) {
      if (!all(sapply(list_genesets, function(geneset) sum(grepl(pattern="ENSG|ENSMUSG", x=geneset))==length(geneset)))) {
        # retrieve gene hgnc symbols in mapping file
        mapping <-  read.csv(file=sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism), stringsAsFactors = F) 
        if (!is.null(names(list_genesets))) names_tmp <- names(list_genesets)
        list_genesets <- lapply(list_genesets, function(geneset) {
          out <- gene_map(data.frame("genes"=geneset),
                          idx_gene_column = 1,
                          mapping=mapping,
                          from="hgnc|symbol|optimal", 
                          to="ensembl",
                          replace = T,
                          na.rm = T)[[1]]
          if (!is.null(names(geneset))) names(out) <- names(geneset)
          return(out)
        })
        names(list_genesets) <- names_tmp; rm(names_tmp)
      }
    }
    
    # convert kMs from dataframes to lists
    list_list_kMs_PPI <- lapply(list_kMs_PPI, function(x) as.list(x)) 
    
    # give gene names to each vector of kMs
    list_list_kMs_PPI <- mapply(function(a,b) lapply(a, function(x) name_for_vec(to_be_named = x, given_names = rownames(b), dimension = NULL)),
                                a = list_list_kMs_PPI,
                                b = list_kMs_PPI,
                                SIMPLIFY = F)
    
    names(list_list_kMs_PPI) = sNames_PPI 
    
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "log_GSEA.txt")
    fun = function(list_values) lapply(list_values, 
                                       function(values) {
                                         out<-iterative.bulk.gsea(values = values, set.lists = list_genesets)
                                         #colnames(out) <- c("p.val", "edge.score", "edge.value", "scaled.score")
                                         return(out)})
    args =  list("X"=list_list_kMs_PPI)
    list_list_GSEA <- safeParallel(fun=fun, args=args, outfile=outfile, list_genesets=list_genesets)
    
    # dissolve the top level celltype list so we have a list of module enrichment vectors
    unlist(list_list_GSEA, recursive = F, use.names = T)  -> list_GSEA 
    
    celltype_col <- rep(sNames_PPI, times=sapply(list_list_GSEA, length))
    mods_col <- gsub("^.*\\.", "", names(list_GSEA))
    GSEA_cols <- Reduce( f = rbind, list_GSEA)
    
    GSEA_df  <- data.frame(celltype=celltype_col, module=mods_col, GSEA_cols, row.names=NULL, stringsAsFactors = F)
    
    # adjust p.vals for multiple testing
    p.adjust(GSEA_df[['p.val']], method = "fdr") %>% -log10(.) -> GSEA_df$q.val 
    
    # find row idx of enriched modules
    idx_row_GSEA <- GSEA_df$q.val > -log10(pvalThreshold)
    
    # count number of enriched module per celltype for summary stats
    n_modules_GSEA_enriched[names(n_modules_GSEA_enriched) %in% sNames_PPI] <- sapply(sNames_PPI, function(name) sum(GSEA_df$celltype[idx_row_GSEA]==name))
  }  
  ##########################################################################
  ############# COMPUTE MAGMA COMMON VARIANTS GWAS ENRICHMENT ##############
  ##########################################################################
  
  # MAGMA docs: http://ctg.cncr.nl/software/MAGMA/doc/manual_v1.06.pdf
  # TODO: Need to add non-protein coding genes to MAGMA
  
  if (!is.null(magma_gwas_dir)) {
    
    message("Scoring modules for enrichment with genes linked by GWAS to phenotypes of interest")
    
    magma_test_type = "module_genes_t" # TODO: make this an argument?
    
    if (data_organism == "mmusculus") {
      
      # map mouse to human gene orthologs  
      mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
      
      outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "mapMMtoHs_par.txt")
      # cl <- makeCluster(min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = 2)-1), type="FORK", outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "mapMMtoHs_par.txt"))
      fun = function(modulekM,colors_PPI_uniq) mapMMtoHs(modulekM = modulekM,
                                                         colors = colors_PPI_uniq,
                                                         log_dir = log_dir, 
                                                         data_prefix = data_prefix, 
                                                         run_prefix = run_prefix,
                                                         mapping_orthology = mapping_orthology)
      args = list(modulekM = list_kMs_PPI,
                  colors_PPI_uniq =list_colors_PPI_uniq)
      list_MMtoHsmapping <- safeParallel(fun=fun, args=args, outfile=outfile, log_dir=log_dir, data_prefix=data_prefix, run_prefix=run_prefix, mapping_orthology=mapping_orthology)
      
      
      list_kMs_hs <- lapply(list_MMtoHsmapping, function(x) x$kM) 
      list_colors_hs <- lapply(list_MMtoHsmapping, function(x) x$colors)  
      
      vec_MMtoHsmapping_prop.mapped <- sapply(list_MMtoHsmapping, function(x) x$prop.mapped, simplify = T)
      
      names(list_kMs_hs) <- names(list_colors_hs) <- names(list_kMs_PPI)
      
      # Also map all 10x genes to homo sapiens
      
      genes_background_ensembl_mm <- read.csv(sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism))[["ensembl"]]
      genes_background_ensembl_hs <- mapping_orthology$ensembl_gene_id[match(genes_background_ensembl_mm, mapping_orthology$mmusculus_homolog_ensembl_gene)]
      genes_background_ensembl_hs <- as.character(na.omit(genes_background_ensembl_hs)) 
      
    } else if (data_organism == "hsapiens") {
      
      list_kMs_hs <- list_kMs_PPI
      list_colors_hs <- list_colors_PPI_uniq
      names(list_kMs_hs) <- names(list_colors_hs) <- names(list_kMs_PPI)
      genes_background_ensembl_hs <- read.csv(sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism))[["ensembl"]]
      
    }
    
    # Load gene mapping and annotation files
    
    # columns: ensembl_gene_id, entrezgene, hgnc_symbol. Use this for mapping entrezgene to ensembl
    mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T, stringsAsFactors = F)
    
    # Load MAGMA genes and remap to Ensembl gene IDs
    
    d = dir(path=magma_gwas_dir, pattern="[.]genes.out", recursive = T)
    gwas = vector(mode="list")
    for(i in 1:length(d)) {
      gwas[[i]] = read.table(paste(magma_gwas_dir, d[[i]],sep=""),head=T, check.names = FALSE, stringsAsFactors = F)
    }
    names(gwas) = gsub(".genes.out", "", d)
    
    # Remap from human Entrez to human Ensembl gene IDs
    for (i in 1:length(gwas)) {
      idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
      mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
      gwas[[i]]$gene_name = mapping$ensembl
    }
    
    
    # invisible(gc()); invisible(R.utils::gcDLLs())
    # cl <- makeCluster(min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = 2)-1), type="FORK", outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "kM_magma_par.txt"))
    outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "kM_magma_par.txt")
    fun = function(cellType,kMs_hs,colors_hs) kM_magma(cellType =cellType, 
                                                       modulekM = kMs_hs,
                                                       gwas = gwas,
                                                       test_type = magma_test_type,
                                                       colors_hs = colors_hs,
                                                       genes_background_ensembl_hs = genes_background_ensembl_hs)
    args = list(cellType = names(list_kMs_hs),
                kMs_hs = list_kMs_hs,
                colors_hs = list_colors_hs)
    magma_results <- safeParallel(fun=fun, args=args, outfile = outfile, gwas=gwas, magma_test_type=magma_test_type, genes_background_ensembl_hs=genes_background_ensembl_hs)
    
    
    # Prepare tables for results
    list_magma.r <- list_magma.p <- list_magma.emp.p <- vector(mode="list", length = length(list_kMs_hs))
    
    # Extract the coefficient, analytical p-value and empirical p-value dataframes for each module, put them into lists
    
    for (i in 1:length(magma_results)) { 
      list_magma.r[[i]] <- magma_results[[i]][['corrCoef']]
      list_magma.p[[i]] <- magma_results[[i]][['p.val']]
      list_magma.emp.p[[i]] <- magma_results[[i]][['emp.p.val']]
    }
    
    # rowbind each list into a single table
    magma.r.all <- Reduce(f=rbind, x = list_magma.r)
    magma.p.all <- Reduce(f=rbind, x = list_magma.p)
    magma.emp.p.all <- Reduce(f=rbind, x = list_magma.emp.p)
    
    # We will plot the analytical and empirical p-values side-by-side later as a diagnostic
    
    # adjust analytical p-values for multiple testing
    magma.p.fdr.all = p.adjust(magma.p.all, method="fdr")
    dim(magma.p.fdr.all) = dim(magma.p.all);  dimnames(magma.p.fdr.all) = dimnames(magma.p.all)
    
    # Transform analytical p-values to retain those associated with positive correlation coefficients only
    magma.p.fdr.signed = as.data.frame(magma.p.fdr.all*sign(magma.r.all))
    magma.p.fdr.signed[magma.p.fdr.signed<0]=1 #Only look for positive enrichment
    magma.p.fdr.log = -log10(magma.p.fdr.signed)
    
    # Convert to dataframes
    magma.p.fdr.log <- magma.p.fdr.log %>% as.data.frame
    magma.p.all <- magma.p.all %>% as.data.frame
    magma.emp.p.all <- magma.emp.p.all %>% as.data.frame
    magma.r.all <- magma.r.all %>% as.data.frame
    
    celltype <- gsub("__.+","",rownames(magma.p.fdr.log))    
    module <- gsub(".+__","",rownames(magma.p.fdr.log))
    
    magma.p.fdr.log <- cbind(celltype, module, magma.p.fdr.log)
    magma.p.all <- cbind(celltype, module, magma.p.all)
    magma.emp.p.all <- cbind(celltype, module, magma.emp.p.all)
    magma.r.all <- cbind(celltype, module, magma.r.all)
    
  } else if (is.null(magma_gwas_dir)) {
    
    magma.r.all <- NULL
    magma.p.all <- NULL
    magma.emp.p.all <- NULL 
    magma.p.fdr.log <- NULL
    vec_MMtoHsmapping_prop.mapped <- NULL
  }
  
  rm(gwas)
  
  ##########################################################################
  ############## FILTER MODULES ON GENES WITH GWAS ENRICHMENT ##############
  ##########################################################################
  
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) {
    
    # Check if the gwas_filter_traits strings provided by the user can be found in the magma_gwas_dir files
    
    if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(magma.p.fdr.log), ignore.case=T)), simplify = T) %>% any) {
      
      # Identify columns to keep
      idx_col_keep <- sapply(gwas_filter_traits, function(x) grep(paste0(x,"|module|celltype"), colnames(magma.p.fdr.log), ignore.case=T), simplify = T) %>% Reduce(union,.)# %>% as.logical 
      
      # fdr-corrected log-transformed analytical p-values
      magma.p.fdr.log.sig <- magma.p.fdr.log[,idx_col_keep]
      
      # correlation coefficients
      magma.r.sig <- magma.r.all[,idx_col_keep]
      
      # Identify rows to keep
      idx_row_gwas <- apply(magma.p.fdr.log.sig[,-grep("module|celltype", colnames(magma.p.fdr.log.sig))], MARGIN = 1, max) > -log10(pvalThreshold)
      
      # rows enriched for MAGMA GWAS or GSEA
      idx_row_keep <- idx_row_gwas | (if (!is.null(list_genesets_path)) idx_row_GSEA[match(GSEA_df$module, magma.p.fdr.log$module)] else rep(TRUE, length(idx_row_gwas)) )
      
      if (sum(idx_row_keep)>0) {
        
        magma.p.fdr.log.sig <- magma.p.fdr.log.sig[idx_row_keep,]
        magma.r.sig <- magma.r.sig[idx_row_keep,]
        
      } else if (sum(idx_row_keep)==0) { # Warn if there are no significant enrichments
        
        gwas_filter_traits <- NULL # this will activate the logic in the next if statement
        warning("After filtering for GWAS enrichment no modules remain! Skipping GWAS filtering step")
        
      }
    } else {
      
      gwas_filter_traits <- NULL # this will activate the logic in the next if statement 
      warning("None of gwas_filter_traits strings found in magma_gwas_dir files! Skipping GWAS filtering step")
      
    }
  }  
  
  if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    magma.p.fdr.log.sig <- NULL # magma.p.fdr.log
    magma.r.sig <- NULL #  magma.r.all
  }
  
  # count number of enriched module per celltype for summary stats
  n_modules_gwas_enriched <- rep(NA, times = length(sNames_0))  
  names(n_modules_gwas_enriched) <- sNames_0
  n_modules_gwas_enriched[names(n_modules_gwas_enriched) %in% sNames_PPI] <- if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) sapply(sNames_PPI, function(x) sum(magma.p.fdr.log$celltype[idx_row_gwas]==x)) else rep(NA, times=length(sNames_PPI))
  
  ######################################################################
  ####### FILTER MODULES ON GWAS SIGNIFICANCE ALSO IN OTHER FILES ######
  ######################################################################
  
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) { # NB: if no significant gwas correlations found in previous section, gwas_filter_traits <- NULL
    
    # get cell types with modules with gwas enrichment
    sNames_gwas <- sNames_PPI[match(names(table(magma.p.fdr.log.sig$celltype)[table(magma.p.fdr.log.sig$celltype)>0]),sNames_PPI)]
    
    list_module_gwas <- vector(mode="list", length=length(sNames_gwas))
    
    names(list_module_gwas) <- sNames_gwas
    
    # Get a list of vectors with modules per celltype
    for (name in sNames_gwas) {
      idx <- grep(pattern = name, x = magma.p.fdr.log.sig$celltype)
      list_module_gwas[[name]] <- as.character(sapply(idx, function(j) magma.p.fdr.log.sig$module[j], simplify = T))
    }
    
    # Filter out cell types with no significant modules
    list_colors_gwas <- list_colors_PPI_uniq[match(sNames_gwas,sNames_PPI)]
    
    # relabel non-GWAS modules 'grey'
    list_colors_gwas <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                               x = list_colors_gwas, 
                               y = list_module_gwas,
                               SIMPLIFY = F)
    
    # give gene names to color assignment vectors
    list_colors_gwas <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                               x = list_colors_gwas,
                               y = list_colors_PPI_uniq[match(sNames_gwas,sNames_PPI), drop=F],
                               SIMPLIFY = F)
    
    # remove any cell clusters without gwas enrichment: expression matrices
    list_datExpr_gwas <- list_datExpr_PPI[match(sNames_gwas,sNames_PPI), drop=F] 
    
    # filter kMs
    list_kMs_gwas <- list_kMs_PPI[match(sNames_gwas,sNames_PPI), drop = F]
    list_kMs_gwas <- mapply(function(x,y) x[colnames(x) %in% y],
                            x = list_kMs_gwas,
                            y = list_module_gwas,
                            SIMPLIFY=F)
    
    # Filter left eigenvectors / kIMs
    list_u_gwas <- list_u_PPI[match(sNames_gwas,sNames_PPI), drop = F]
    
    list_u_gwas <- mapply(function(x,y) x[names(x) %in% y],
                          x = list_u_gwas,
                          y = list_module_gwas,
                          SIMPLIFY=F)
    
    # filter MEs
    if (fuzzyModMembership=="kME") {
      
      list_MEs_gwas <- list_MEs_PPI[match(sNames_gwas,sNames_PPI), drop = F]
      list_MEs_gwas <- mapply(function(x,y) x[colnames(x) %in% y], 
                              x = list_MEs_gwas,
                              y = list_module_gwas,
                              SIMPLIFY=F)
      
    } else {
      list_MEs_gwas <- NULL
      list_u_gwas <- NULL
    }
    
  } else if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    
    sNames_gwas <- sNames_PPI
    list_module_gwas <- lapply(list_kMs_PPI, function(x) colnames(x)) 
    names(list_module_gwas) <- sNames_gwas
    list_colors_gwas <- list_colors_PPI_uniq
    list_datExpr_gwas <- list_datExpr_PPI
    list_kMs_gwas <- list_kMs_PPI 
    list_MEs_gwas <- list_MEs_PPI
    list_u_gwas <- list_u_PPI
  }
  
  #####################################################################
  ####### COMPUTE MODULE - METADATA CORRELATION IN EACH CELL CLUSTER ###
  ######################################################################
  
  if (fuzzyModMembership=="kIM") {
    list_dissTOM_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_", run_prefix, "_list_dissTOM"), full.names = T)
    list_dissTOM <- load_obj(list_dissTOM_path)
    list_dissTOM_gwas <- list_dissTOM[match(sNames_gwas, names(list_dissTOM))]
    rm(list_dissTOM)  
  }
  
  if (!is.null(metadata_corr_col)) {
    
    if (!is.null(metadata)) {
      
      message("Computing module-metadata correlation in each celltype")
      
      list_metadata <- lapply(list_datExpr_gwas, function(x) metadata[match(rownames(x), rownames(metadata)), , drop=F]) # get list of cell * metadata. The space between the commas is intentional!
      
      if (fuzzyModMembership=="kME") {
        # Compute correlation between metadata (columns) and eigengenes (columns). 
        list_mod_metadata_corr_rho <- mapply(function(x,y) cor(x=as.matrix(x), 
                                                               y=as.matrix(y), 
                                                               method = c("pearson"), 
                                                               use = 'pairwise.complete.obs'), 
                                             x=list_metadata, 
                                             y=list_MEs_gwas,  
                                             SIMPLIFY = F) 
        # Remove 'ME' from eigengene names
        for (j in 1:length(list_mod_metadata_corr_rho)) {
          colnames(list_mod_metadata_corr_rho[[j]]) <- gsub("^ME", "", colnames(list_mod_metadata_corr_rho[[j]]), ignore.case=F)
        }
        
      } else if (fuzzyModMembership == "kIM") {
        
        # compute the equivalent to eigengenes (i.e. embeddings) but using kIMs as eigenvectors on which to project each cell's gene-length vector
        list_embed_mat <- mapply(function(a,b,c,d,e) cellModEmbed(datExpr=a, 
                                                                  colors=b, 
                                                                  latentGeneType="IM",
                                                                  cellType=c,
                                                                  kMs=d,
                                                                  dissTOM=e),
                                 a = list_datExpr_gwas,
                                 b = list_colors_gwas,
                                 c = names(list_datExpr_gwas),
                                 d = list_kMs_gwas,
                                 e = list_dissTOM_gwas,
                                 SIMPLIFY=F)
        
        # Get correlations
        list_mod_metadata_corr_rho <- mapply(function(x,y) cor(x=as.matrix(x), 
                                                               y=as.matrix(y), 
                                                               method = c("pearson"), 
                                                               use = 'pairwise.complete.obs'), 
                                             x=list_metadata, 
                                             y=list_embed_mat,  
                                             SIMPLIFY = F) 
        
        # Clear space
        rm(list_embed_mat)
        
        list_mod_metadata_corr_rho <- lapply(list_mod_metadata_corr_rho, function(x) name_for_vec(to_be_named = x, given_names = colnames(metadata), dimension = 1)) 
        
      }
      
      
      # Name the rows of the correlation matrix with the metadata columns
      
      # Compute p values
      list_mod_metadata_corr_pval <- mapply(function(x,y) WGCNA::corPvalueStudent(x, 
                                                                                  n = nrow(y)), 
                                            x=list_mod_metadata_corr_rho, 
                                            y=list_metadata, 
                                            SIMPLIFY=F)
      
      # Convert p values into one big matrix in order to adjust p-values for the number of modules tested across *all* celltypes
      corr_pval <- Reduce(x=list_mod_metadata_corr_pval, f = cbind)
      
      # set NAs to 1
      corr_pval[is.na(corr_pval)] <- 1
      
      # Compute the false discovery rates
      corr_fdr <- apply(corr_pval, MARGIN=1, FUN = function(x) p.adjust(x, method = "fdr")) %>% as.matrix # apply outputs the vectors as columns
      if (ncol(corr_pval)>1) corr_fdr <- t(corr_fdr)
      corr_fdr.log <- -log10(corr_fdr)     
      colnames(corr_fdr.log) <- colnames(corr_pval)
      
      # Split single dataframe back into a list of celltypes
      list_mod_metadata_corr_fdr <- list_mod_metadata_corr_fdr.log <- vector(mode = "list", length = length(list_mod_metadata_corr_pval))
      names(list_mod_metadata_corr_fdr) <- names(list_mod_metadata_corr_fdr.log) <- names(list_mod_metadata_corr_rho)
      
      k = 0
      for (j in 1:length(list_mod_metadata_corr_fdr.log)) {
        list_mod_metadata_corr_fdr[[j]] <- as.data.frame(corr_fdr[,(k+1):(k+ncol(list_mod_metadata_corr_pval[[j]])), drop=F])
        list_mod_metadata_corr_fdr.log[[j]] <- as.data.frame(corr_fdr.log[,(k+1):(k+ncol(list_mod_metadata_corr_pval[[j]])), drop=F])
        k <- k + ncol(list_mod_metadata_corr_pval[[j]])
      }
      
      # Also make a single dataframe for saving 
      dim(corr_fdr.log) <- c(ncol(metadata), length(unlist(list_module_gwas)))
      rownames(corr_fdr.log) <- rownames(corr_pval)
      colnames(corr_fdr.log) <- colnames(corr_pval)
      corr_fdr.log <- as.data.frame(corr_fdr.log)
      
      # Also single rho df
      corr_rho <- Reduce(x=list_mod_metadata_corr_rho, f=cbind)
      rownames(corr_rho) <- rownames(corr_pval)
      colnames(corr_rho) <- colnames(corr_pval)
      
    } 
    
  } else {
    metadata = NULL
  }
  
  if (fuzzyModMembership=="kIM" ) rm(list_dissTOM_gwas)
  
  ##########################################################################
  #### FILTER COLORS VECS, GENE AND KME LISTS FOR METADATA CORRELATIONS ####
  ##########################################################################
  
  if (!is.null(metadata_corr_col) & !is.null(metadata_corr_filter_vals)) {
    
    if (!is.null(metadata)) {
      
      # Get a list of logical vectors indicating significant correlations
      list_idx_module_meta_sig <- lapply(list_mod_metadata_corr_fdr.log, function(x) apply(x, 2, function(y) any(y>-log10(pvalThreshold)))) # logical
      
      if (sum(sapply(list_idx_module_meta_sig, sum))==0) {
        
        warning("After filtering for metadata correlations no modules remain! Skipping metadata filtering step")
        sNames_meta <- sNames_gwas
        list_colors_meta <- list_colors_gwas
        list_module_meta <- list_module_gwas
        list_kMs_meta <- list_kMs_gwas
        list_MEs_meta <- list_MEs_gwas
        list_u_meta <- list_u_gwas 
        list_datExpr_meta <- list_datExpr_gwas
        
        metadata_corr_filter_vals = NULL
        metadata_corr_col = NULL        
        metadata=NULL
        
      } else {
        
        # Only keep significantly correlated modules within each celltype
        list_module_meta <- mapply(function(x,y) colnames(x[,y, drop=F]),x=list_mod_metadata_corr_fdr.log, y=list_idx_module_meta_sig, SIMPLIFY=F) 
        list_module_meta <- lapply(list_module_meta, function(x) gsub("^.*__", "", x)) %>% Filter(f=length)
        
        sNames_meta <- sNames_gwas[match(names(list_module_meta),sNames_gwas)] # ordered correctly
        
        # Keep the order
        list_module_meta <- list_module_meta[match(sNames_meta, names(list_module_meta))]
        
        list_datExpr_meta <- list_datExpr_gwas[match(sNames_meta, names(list_module_meta))]
        
        # reassign genes of filtered out modules to grey and remove any empty cell clusters
        list_colors_meta <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                                   x = list_colors_gwas[match(sNames_meta,names(list_colors_gwas))],
                                   y = list_module_meta,
                                   SIMPLIFY = F)
        
        # give gene names to color assignment vectors
        list_colors_meta <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                                   x = list_colors_meta,
                                   y = list_colors_gwas[match(sNames_meta,names(list_colors_gwas))],
                                   SIMPLIFY = F)
        
        # Filter kM lists
        list_kMs_meta <- list_kMs_gwas[match(sNames_meta,names(list_kMs_gwas))]
        list_kMs_meta <- mapply(function(x,y) x[,match(y,colnames(x)), drop=F],x=list_kMs_meta, y = list_module_meta, SIMPLIFY=F) %>% Filter(f=length)
        
        # Filter left eigenvectors / kIMs 
        list_u_meta <- list_u_gwas[match(sNames_meta,names(list_u_gwas))]
        list_u_meta <- mapply(function(x,y) x[match(y,names(x))],
                              x = list_u_meta, 
                              y = list_module_meta,
                              SIMPLIFY=F) %>% Filter(f=length)
        
        # Filter ME lists
        if (fuzzyModMembership=="kME") {
          list_MEs_meta <- list_MEs_gwas[match(sNames_meta,names(list_MEs_gwas))]
          list_MEs_meta <- mapply(function(x,y) x[match(y,colnames(x))],
                                  x = list_MEs_meta, 
                                  y = list_module_meta,
                                  SIMPLIFY=F) %>% Filter(f=length)
          
          list_u_meta <- list_u_gwas[match(sNames_meta,names(list_u_gwas))]
          list_u_meta <- mapply(function(x,y) x[match(y, names(x))],
                                x = list_u_meta, 
                                y = list_module_meta,
                                SIMPLIFY=F) %>% Filter(f=length)
          
          
        }  else {
          list_MEs_meta <- list_MEs_gwas 
          list_u_meta <- list_u_gwas
        }
      }
      
    } else {
      
      sNames_meta <- sNames_gwas
      list_colors_meta <- list_colors_gwas
      list_module_meta <- list_module_gwas
      list_kMs_meta <- list_kMs_gwas
      list_MEs_meta <- list_MEs_gwas
      list_u_meta <- list_u_gwas
      list_datExpr_meta <- list_datExpr_gwas
    }
    
  } else if (is.null(metadata_corr_col) | is.null(metadata_corr_filter_vals)) {
    
    sNames_meta <- sNames_gwas
    list_colors_meta <- list_colors_gwas
    list_module_meta <- list_module_gwas
    list_kMs_meta <- list_kMs_gwas 
    list_MEs_meta <- list_MEs_gwas
    list_u_meta <- list_u_gwas
    list_datExpr_meta <- list_datExpr_gwas
  }
  
  # count number of enriched module per celltype for summary stats
  n_modules_meta_enriched <- rep(NA, times=length(sNames_0))
  names(n_modules_meta_enriched) <-sNames_0
  if (!is.null(metadata_corr_col) & !is.null(metadata_corr_filter_vals)) n_modules_meta_enriched[names(n_modules_meta_enriched) %in% sNames_0] <- sapply(list_idx_module_meta_sig, sum) 
  
  ######################################################################
  ############################ COMPUTE PRIMARY KMs #####################
  ######################################################################
  # Extract the primary kMs - i.e. those of each gene w.r.t. the module to which it belongs
  
  message(paste0("Computing primary "), fuzzyModMembership, "s")
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "list_pkMs_meta.txt")
  fun = function(kMs, colors) pkMs_fnc(kMs=kMs, colors=colors)
  args = list("kMs" = list_kMs_meta, "colors"=list_colors_meta)
  list_pkMs_meta <- safeParallel(fun=fun, args=args, outfile=outfile)
  
  names(list_pkMs_meta) <- sNames_meta
  
  ######################################################################
  ##################### COMPUTE CELL x EIGENGENE MATRIX ################
  ######################################################################
  
  message("Computing all cell embeddings on all modules, across celltypes")
  
  scale_data_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_", run_prefix, "_scale_regr_data_ensembl"), full.names = T)
  load_obj(scale_data_path[1]) %>% t -> datExpr 
  
  invisible(gc())
  
  latentGeneType <- if (fuzzyModMembership == "kME") "ME" else if (fuzzyModMembership=="kIM") "IM"
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "make_cell_embed_mat.txt")
  fun = function(x,y,z) cellModEmbed(datExpr = datExpr, 
                                     colors = x,
                                     latentGeneType = latentGeneType,
                                     cellType = y,
                                     kMs = if (latentGeneType== "IM") z else NULL)#,
  
  args = list(x = list_colors_meta, 
              y = names(list_colors_meta),
              z = list_kMs_meta)#,
  
  list_cellModEmbed_mat <- safeParallel(fun=fun, args=args, outfile=outfile, datExpr = datExpr, latentGeneType = latentGeneType)
  
  
  list_cellModEmbed_mat %>% Reduce(function(mat1, mat2) cbind(mat1, mat2), .) -> cellModEmbed_mat
  
  rownames(cellModEmbed_mat) <- NULL
  
  cellModEmbed_mat_annot <- cbind(rep(x = data_prefix,  times=nrow(cellModEmbed_mat)), ident, rownames(datExpr), cellModEmbed_mat)
  colnames(cellModEmbed_mat_annot) <- c("data", "ident", "cell_id", colnames(cellModEmbed_mat)) 
  
  rm(datExpr)
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_5"
  message("Reached checkpoint 5, saving session image")
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_5_image.RData.gz", scratch_dir, data_prefix, run_prefix), compress = "gzip")
  
} 

if (resume == "checkpoint_5") {
  
  ##########################################################################
  ############################ (UN)LOAD PACKAGES ############################
  ##########################################################################
  
  invisible(R.utils::gcDLLs())
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("readr"))
  
  ##########################################################################
  ######### PREPARE GENES LISTS AND DATAFRAME WITH MODULES, GENES ##########
  ##########################################################################
  # retrieve gene hgnc symbols in mapping file
  mapping <- if (map_genes_to_ensembl) read.csv(file=sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, run_prefix, data_organism), stringsAsFactors = F) else NULL
  
  outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "write_outputs.txt")
  
  message("Preparing outputs")
  
  # prepare nested lists of module genes
  fun = function(a,b) lapply(b, function(x) names(a)[a==x])
  args = list(a=list_colors_meta, b=list_module_meta)
  list_list_module_meta_genes = safeParallel(fun=fun,args=args, outfile=outfile)
  fun = function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL)
  args=  list(x=list_list_module_meta_genes, y=list_module_meta)
  list_list_module_meta_genes <- safeParallel(fun=fun,args=args,outfile=outfile)
  
  # make a copy for pkMs
  list_list_module_meta_pkMs <- list_list_module_meta_genes # just a template
  list_list_module_gene_loadings <- if (fuzzyModMembership == "kME") list_list_module_meta_genes else NULL # just a template
  # order the gene lists by pkM and get pkM value
  
  # iterate over celltypes
  for (i in 1:length(list_list_module_meta_genes)) {
    # iterate over modules
    for (j in 1:length(list_list_module_meta_genes[[i]])) {
      pkMs = list_pkMs_meta[[i]]
      if (fuzzyModMembership == "kME") mod_u <- list_u_meta[[i]][[j]]
      mod_genes <- list_list_module_meta_genes[[i]][[j]]
      mod_pkMs_sorted <- sort(pkMs[match(mod_genes,names(pkMs)), drop=F], decreasing=T)
      list_list_module_meta_genes[[i]][[j]] <- names(mod_pkMs_sorted)  # sort the genes
      list_list_module_meta_pkMs[[i]][[j]] <- mod_pkMs_sorted
      if (fuzzyModMembership == "kME") list_list_module_gene_loadings[[i]][[j]] <- mod_u[match(names(mod_pkMs_sorted), names(mod_u)), drop=F]
    }
  }
  
  message("Preparing module genes dataframe")
  
  cell_cluster <- rep(sNames_meta, times=unlist(lapply(list_list_module_meta_genes, FUN=function(x) sum(sapply(x, function(y) length(y), simplify=T)))))
  module <- unlist(lapply(list_list_module_meta_genes, function(x) rep(names(x), sapply(x, function(y) length(y), simplify = T))), use.names = F)
  ensembl <- if (map_genes_to_ensembl) unlist(list_list_module_meta_genes, recursive = T, use.names = F) else NULL 
  hgnc <- if (!map_genes_to_ensembl) unlist(list_list_module_meta_genes, recursive = T, use.names = F) else NULL 
  pkMs <- unlist(list_list_module_meta_pkMs, recursive = T, use.names = F)
  if (fuzzyModMembership == "kME") gene_loadings <- unlist(list_list_module_gene_loadings, recursive = T, use.names = F)
  data <- rep(data_prefix, times= length(pkMs))
  run <- rep(run_prefix, times= length(pkMs))    
  
  if (map_genes_to_ensembl) {
    fun = function(x) lapply(x, function(y) mapping$symbol[match(y, mapping$ensembl)])
    args = list("X"=list_list_module_meta_genes)
    list_list_module_meta_genes_hgnc = safeParallel(fun=fun,args=args, outfile=outfile, mapping=mapping)
    
    fun = function(x,y) name_for_vec(to_be_named = x, 
                                     given_names = y, 
                                     dimension = NULL)
    args = list(x=list_list_module_meta_genes_hgnc, 
                y=list_module_meta)
    list_list_module_meta_genes_hgnc = safeParallel(fun=fun, args=args, outfile=outfile)
    hgnc <- unlist(list_list_module_meta_genes_hgnc, recursive = T, use.names=F) 
    df_meta_module_genes <- data.frame(data, run, cell_cluster, module, ensembl, hgnc, pkMs, row.names = NULL)
    colnames(df_meta_module_genes) <- c("data", "run", "cell_cluster", "module", "ensembl", "hgnc", paste0("p", fuzzyModMembership))
  } else {
    df_meta_module_genes <- data.frame(data, run, cell_cluster, module, hgnc, pkMs, row.names = NULL)
    colnames(df_meta_module_genes) <- c("data", "run", "cell_cluster", "module", "hgnc", paste0("p", fuzzyModMembership))
  }
  if (fuzzyModMembership=="kME") df_meta_module_genes[["gene_loadings"]] <- gene_loadings
  
  ##########################################################################
  ############################# OUTPUT TABLES ##############################
  ##########################################################################
  
  message("Writing outputs to disk")
  
  ##################### WRITE OUT TABLES OF MODULE GENES ####################
  ###########################################################################
  
  write.csv(df_meta_module_genes, file = gzfile(sprintf("%s%s_%s_cell_cluster_module_genes.csv.gz",tables_dir, data_prefix, run_prefix)), quote = F, row.names = F)
  
  ####################### SAVE STRINGDB PPI OUTPUT #########################
  ##########################################################################
  
  
  # output 
  # fun = function(module_PPI_uniq,cell_cluster) write.csv(as.matrix(module_PPI_uniq, ncol=2), file=gzfile(sprintf("%s%s_%s_%s_STRINGdb_output_all.csv.gz", tables_dir, data_prefix, run_prefix, cell_cluster)), row.names = F, quote = F)
  # args = list(module_PPI_uniq=list_module_PPI_uniq, 
  #             cell_cluster= sNames_PPI[match(names(list_module_PPI_uniq),sNames_PPI)])
  # invisible(safeParallel(fun=fun, args=args, outfile=outfile, tables_dir=tables_dir, data_prefix=data_prefix, run_prefix))
  
  list_module_PPI_uniq <- mapply(FUN = function(module_PPI_uniq, cell_cluster) {
    if (!is.null(nrow(module_PPI_uniq))) module_PPI_uniq[["cell_cluster"]] <- rep(cell_cluster, nrow(module_PPI_uniq))
    return(module_PPI_uniq)
  },module_PPI_uniq=list_module_PPI_uniq, cell_cluster=names(list_module_PPI_uniq), SIMPLIFY=F)
  module_PPI_uniq_df <- Reduce(x=list_module_PPI_uniq, f = rbind)
  colnames(module_PPI_uniq_df) <- c("colors","q.value","expected.interactions", "cell_cluster")
  write.csv(module_PPI_uniq_df, file=gzfile(sprintf("%s%s_%s_STRINGdb_output_all.csv.gz", tables_dir, data_prefix, run_prefix)), row.names = F, quote = F)
  
  # Now significantly PPI enriched modules
  # convert the module PPI dataframe columns from list to numeric 
  list_module_PPI_signif_uniq %>% Filter(f=nrow) -> list_module_PPI_signif_uniq_f
  
  if (length(list_module_PPI_signif_uniq_f) > 0) {
    list_module_PPI_signif_uniq_f <- mapply(FUN = function(module_PPI_signif_uniq_f, cell_cluster) {
      if (!is.null(nrow(module_PPI_signif_uniq_f))) module_PPI_signif_uniq_f[["cell_cluster"]] <- rep(cell_cluster, nrow(module_PPI_signif_uniq_f))
      return(module_PPI_signif_uniq_f)
    },module_PPI_signif_uniq_f=list_module_PPI_signif_uniq_f, cell_cluster=names(list_module_PPI_signif_uniq_f), SIMPLIFY=F)
    module_PPI_signif_uniq_df_f <- Reduce(x=list_module_PPI_signif_uniq_f, f = rbind)
    colnames(module_PPI_signif_uniq_df_f) <- c("colors","q.value","expected.interactions","cell_cluster")
    
    write.csv(module_PPI_signif_uniq_df_f, file=gzfile(sprintf("%s%s_%s_STRINGdb_output_signif.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=F, quote=F)
    
    # fun = function(x,y) write.csv(as.matrix(x, ncol=2), file=gzfile(sprintf("%s%s_%s_%s_STRINGdb_output_signif.csv.gz", tables_dir, data_prefix, run_prefix, y)), 
    #                               row.names=F, 
    #                               quote = F)
    # args=list(x=list_module_PPI_signif_uniq_f, 
    #           y= sNames_PPI[match(names(list_module_PPI_signif_uniq_f),sNames_PPI)])
    # 
    # invisible(mapply(function(x,y) write.csv(as.matrix(x, ncol=2), file=gzfile(sprintf("%s%s_%s_%s_STRINGdb_output_signif.csv.gz", tables_dir, data_prefix, run_prefix, y)), 
    #                                                  row.names=F, 
    #                                                  quote = F), 
    #                      x=list_module_PPI_signif_uniq_f, 
    #                      y= sNames_PPI[match(names(list_module_PPI_signif_uniq_f),sNames_PPI)], 
    #                      SIMPLIFY = F))
    
  }
  
  ################## SAVE KMs FOR ENRICHED MODULES #########################
  ##########################################################################
  
  # Prepare dfs with a gene column followed by kMEs / kIMs 
  list_kMs_meta_out <- lapply(list_kMs_meta, function(kMs) {
    kMs<-cbind(genes=rownames(kMs), kMs)
    rownames(kMs) <- NULL
    kMs
  })
  
  
  # invisible(mapply(function(x,y) write.csv(x, file=gzfile(sprintf("%s%s_%s_%s_kMs.csv.gz", tables_dir, data_prefix, run_prefix, y)), 
  #                                                  row.names=F, quote = F), 
  #                  list_kMs_meta_out, sNames_meta, SIMPLIFY = F))
  
  # save kMs as a single outer join 
  kMs_meta_out <- Reduce(f = function(df1, df2) {
    dplyr::full_join(df1, df2, by="genes")}, 
    x=list_kMs_meta_out)
  
  write.csv(kMs_meta_out, file=gzfile(sprintf("%s%s_%s_kMs_full_join.csv.gz", tables_dir, data_prefix, run_prefix)), 
            row.names=F, quote = F)
  
  ######################## OUTPUT GSEA RESULTS ##############################
  ###########################################################################
  
  if (!is.null(list_genesets_path)) {
    try({GSEA_df <- arrange(GSEA_df, desc(scaled.score))  
    invisible(write.csv(GSEA_df, file=gzfile(sprintf("%s%s_%s_GSEA_df.csv.gz", tables_dir, data_prefix, run_prefix)), 
                        row.names=F, quote = F))
    })
  }
  ######################## OUTPUT MAGMA RESULTS #############################
  ###########################################################################
  
  if (!is.null(magma_gwas_dir)) {
    invisible(write.csv(magma.p.all, file=gzfile(sprintf("%s%s_%s_magma.p.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=F, quote = F))
    invisible(write.csv(magma.emp.p.all, file=gzfile(sprintf("%s%s_%s_magma.emp.p.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=F, quote = F))
    invisible(write.csv(magma.p.fdr.log, file=gzfile(sprintf("%s%s_%s_magma.fdr.log.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=F, quote = F))
    invisible(write.csv(magma.r.all, file=gzfile(sprintf("%s%s_%s_magma.r.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=F, quote = F))
    
    if(!is.null(magma.p.fdr.log.sig)) invisible(write.csv(magma.p.fdr.log.sig, file=gzfile(sprintf("%s%s_%s_magma.fdr.log.sig.csv.gz", tables_dir, data_prefix, run_prefix)), 
                                                          row.names=F, quote = F))
    if(!is.null(magma.r.sig))invisible(write.csv(magma.r.sig, file=gzfile(sprintf("%s%s_%s_magma.r.sig.csv.gz", tables_dir, data_prefix, run_prefix)), 
                                                 row.names=F, quote = F))
  }
  
  ####################### OUTPUT METADATA CORRELATION RESULTS ###############
  ###########################################################################
  
  if (!is.null(metadata_corr_col) & !is.null(metadata)) {
    invisible(write.csv(corr_fdr.log, file=gzfile(sprintf("%s%s_%s_all_metadata_corr_logfdr.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=T, quote = F))
    invisible(write.csv(corr_rho, file=gzfile(sprintf("%s%s_%s_all_metadata_rho.csv.gz", tables_dir, data_prefix, run_prefix)), row.names=T, quote = F))
    
    # invisible(mapply(function(x,y) write.csv(x, file=gzfile(sprintf("%s%s_%s_%s_metadata_corr_rho.csv.gz", tables_dir, data_prefix, run_prefix, y)), row.names=T, quote = F), 
    #                      list_mod_metadata_corr_rho, 
    #                      sNames_gwas, SIMPLIFY = F))
    # invisible(mapply(function(x,y) write.csv(x, file=gzfile(sprintf("%s%s_%s_%s_metadata_corr_logfdr.csv.gz", tables_dir, data_prefix, run_prefix, y)), row.names=T, quote = F), 
    #                      list_mod_metadata_corr_fdr.log, 
    #                      sNames_gwas, SIMPLIFY = F))
  }
  
  ################## OUTPUT CELL MODULE EMBEDDINGS MATRIX ###################
  ###########################################################################
  
  invisible(write.csv(x=as.data.frame(cellModEmbed_mat_annot), file=gzfile(sprintf("%s%s_%s_%s_cellModEmbed.csv.gz", tables_dir, data_prefix, run_prefix, fuzzyModMembership)), row.names=F, quote = F))
  
  ######## OUTPUT MODULE LEFT SINGULAR COMPONENTS (U) FOR EACH MODULE (U) ###
  ########################################################################### 
  
  saveRDS(object = list_u_meta, file=sprintf("%s%s_%s_list_list_module_u.RDS.gz", RObjects_dir, data_prefix, run_prefix), compress = "gzip")
  
  ##################### WRITE PARAMETERS AND STATS TO FILES ################
  ##########################################################################
  
  t_finish <- as.character(Sys.time())
  
  sumstats_celltype_df <- matrix(data = NA_real_, nrow=length(sNames_0), ncol=19) %>% data.frame 
  colnames(sumstats_celltype_df) = c("run_prefix", 
                                     "subset", 
                                     "n_cells", 
                                     "n_genes", 
                                     "softpower", 
                                     "SFT.R.sq", 
                                     "median.k.", 
                                     "plot_label_final", 
                                     "n_genes_reassign", 
                                     "prop_genes_t.test_fail", 
                                     "prop_genes_assign",
                                     "prop_genes_assign_PPI",
                                     "n_modules",
                                     "n_modules_PPI_enriched",
                                     "prop_genes_mapped_to_ortholog",
                                     "n_modules_GSEA_enriched",
                                     "n_modules_gwas_enriched",
                                     "n_modules_meta_enriched")
  
  sumstats_celltype_df[["run_prefix"]] = rep(run_prefix, times=length(sNames_0))
  sumstats_celltype_df[["subset"]] = sNames_0
  sumstats_celltype_df[["n_cells"]] = n_cells_subsets
  sumstats_celltype_df[["n_genes"]] = n_genes_subsets
  sumstats_celltype_df[["softPower"]][sNames_0 %in% sNames_2] = sapply(list_sft, function(x) x$Power)
  sumstats_celltype_df[["SFT.R.sq"]][sNames_0 %in% sNames_2] = sapply(list_sft, function(x) x$SFT.R.sq)
  sumstats_celltype_df[["median.k."]][sNames_0 %in% sNames_2] = sapply(list_sft, function(x) x$median.k.)
  sumstats_celltype_df[["plot_label_final"]][sNames_0 %in% sNames_5 ] = unlist(x=list_plot_label_final, use.names = F)
  if (kM_reassign) sumstats_celltype_df[["n_genes_reassign"]][sNames_0 %in% sNames_5] = sapply(list_reassign_log, function(rlog) if(!is.null(rlog)) nrow(rlog) else 0, simplify=T)
  if (kM_signif_filter) sumstats_celltype_df[["prop_genes_t.test_fail"]][sNames_0 %in% sNames_5] = sapply(list_geneMod_t.test, function(x) {if (!is.null(x)) {sum(as.numeric(!x$signif))/length(x$signif)} else {NA_real_}}, simplify =T)
  sumstats_celltype_df[["prop_genes_assign"]][sNames_0 %in% sNames_5] = sapply(list_colors_all, function(x) round(sum(x!="grey")/length(x),2), simplify=T)
  sumstats_celltype_df[["prop_genes_assign_PPI"]][sNames_0 %in% sNames_PPI] = sapply(list_colors_PPI_uniq, function(x) round(sum(x!="grey")/length(x),2), simplify=T)
  sumstats_celltype_df[["n_modules"]][sNames_0 %in% sNames_5] = sapply(list_colors_all, function(x) length(unique(as.character(x)))-1, simplify=T)
  sumstats_celltype_df[["n_modules_PPI_enriched"]][sNames_0 %in% sNames_PPI] = sapply(list_colors_PPI_uniq, function(x) length(unique(as.character(x)))-1, simplify=T)
  if (!is.null(magma_gwas_dir)) sumstats_celltype_df[["prop_genes_mapped_to_ortholog"]][sNames_0 %in% sNames_PPI] = vec_MMtoHsmapping_prop.mapped
  if (!is.null(list_genesets_path)) sumstats_celltype_df[["n_modules_GSEA_enriched"]] = n_modules_GSEA_enriched
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) sumstats_celltype_df[["n_modules_gwas_enriched"]] = n_modules_gwas_enriched
  if (!is.null(metadata) & !is.null(metadata_corr_col))  sumstats_celltype_df[["n_modules_meta_enriched"]] =n_modules_meta_enriched
  
  
  # sumstats_celltype_df <- data.frame(run_prefix = rep(run_prefix, times=length(sNames_0)),
  #                                    subset = sNames_0, 
  #                                    n_cells = n_cells_subsets,
  #                                    n_genes = n_genes_subsets,
  # softPower = sapply(list_sft, function(x) x$Power),
  # SFT.R.sq = sapply(list_sft, function(x) x$SFT.R.sq),
  # median.k. = sapply(list_sft, function(x) x$median.k.),
  # plot_label_final = ifelse(test=sNames_0 %in% sNames_5, yes= unlist(x=list_plot_label_final, use.names = F), no=NA),
  # n_genes_reassign = if (kM_reassign) ifelse(test=sNames_0 %in% sNames_5, yes = sapply(list_reassign_log, function(rlog) if(!is.null(rlog)) nrow(rlog) else 0, simplify=T), no = 0) else numeric(length=sNames_0),
  # prop_genes_t.test_fail = if (kM_signif_filter) ifelse(test= sNames_0 %in% sNames_5, yes=sapply(list_geneMod_t.test, function(x) {if (!is.null(x)) {sum(as.numeric(!x$signif))/length(x$signif)} else {NA_real_}}, simplify =T), no=NA) else rep(NA, times=length(sNames_0)),
  # prop_genes_assign =  ifelse(test = sNames_0 %in% sNames_5, yes = sapply(list_colors_all, function(x) round(sum(x!="grey")/length(x),2), simplify=T), no = 0),
  # prop_genes_assign_PPI = ifelse(test = sNames_0 %in% sNames_PPI, yes = sapply(list_colors_PPI_uniq, function(x) round(sum(x!="grey")/length(x),2), simplify=T), no = 0),
  # n_modules = ifelse(test=sNames_0 %in% sNames_5, yes=sapply(list_colors_all, function(x) length(unique(as.character(x)))-1, simplify=T), no = 0),
  # n_modules_PPI_enriched = ifelse(test=sNames_0 %in% sNames_PPI, yes=sapply(list_colors_PPI_uniq, function(x) length(unique(as.character(x)))-1, simplify=T), no = 0),
  # prop_genes_mapped_to_ortholog = if (!is.null(magma_gwas_dir)) if (data_organism=="mmusculus") ifelse(test=sNames_0 %in% sNames_PPI, yes=vec_MMtoHsmapping_prop.mapped, no=NA) else rep(NA, length(sNames_0)) else rep(NA, length(sNames_0)), 
  # n_modules_GSEA_enriched = n_modules_GSEA_enriched,
  # n_modules_gwas_enriched = if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) n_modules_gwas_enriched else rep(NA, times=length(sNames_0)),
  # n_modules_meta_enriched = if (!is.null(metadata) & !is.null(metadata_corr_col)) n_modules_meta_enriched else rep(NA, times=length(sNames_0)),
  # row.names = NULL, stringsAsFactors = F)
  
  sumstats_run <- c(run_prefix = run_prefix,
                    t_start = t_start,
                    t_finish = t_finish,
                    mean_percent_mito = mean(percent.mito, na.rm=T),
                    mean_percent_ribo = mean(percent.mito, na.rm=T),
                    n_celltypes = length(sNames_0),
                    n_celltypes_used = length(sNames_5),
                    prop_genes_mapped_to_ensembl = if (map_genes_to_ensembl) round(sum(!is.na(mapping$ensembl))/nrow(mapping),2) else NA_real_,
                    prop_genes_assign_mean = round(mean(sumstats_celltype_df$prop_genes_assign, na.rm=T),2),
                    prop_genes_assign_PPI_mean = round(mean(sumstats_celltype_df$prop_genes_assign_PPI, na.rm=T),2),
                    prop_genes_t.test_fail_mean = if (kM_signif_filter) tryCatch({sum((sapply(list_geneMod_t.test, function(x)  {if (!is.null(x)) {sum(as.numeric(!x$signif))/length(x$signif)} else {NA_real_}}, simplify=T)))/ sum(sapply(list_geneMod_t.test, function(x) length(x$signif), simplify =T))}, error = function(err) {NA_real_}) else NA,
                    n_modules_total = sum(sumstats_celltype_df$n_modules),
                    n_modules_PPI_enriched_total = sum(sumstats_celltype_df$n_modules_PPI_enriched), 
                    prop_genes_mapped_to_ortholog = if (!is.null(magma_gwas_dir)) if (data_organism=="mmusculus") round(mean(sumstats_celltype_df$prop_genes_mapped_to_ortholog, na.rm=T),2) else NA else NA,
                    n_modules_GSEA_enriched = if(!is.null(list_genesets_path)) sum(sumstats_celltype_df$n_modules_GSEA_enriched, na.rm=T) else NA_real_,
                    n_modules_magma_enriched = if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) sum(n_modules_gwas_enriched, na.rm = T) else NA,
                    n_modules_meta_enriched = if (!is.null(metadata) & !is.null(metadata_corr_col)) sum(sumstats_celltype_df$n_modules_meta_enriched, na.rm = T) else NA) 
  
  params_run <- opt[-length(opt)]# all the user options/defaults
  
  # Convert sumstats and params from list to one-row data.frame
  matrix(unlist(sumstats_run), nrow=1, dimnames=list(NULL, c(names(sumstats_run)))) %>% as.data.frame(row.names=NULL, stringsAsFactors=F) -> sumstats_run_df
  matrix(unlist(params_run), nrow=1, dimnames=list(NULL, c(names(params_run)))) %>% as.data.frame(row.names=NULL, stringsAsFactors=F) -> params_run_df
  
  # make path strings
  sumstats_celltype_path = sprintf("%s%s_%s_sumstats_celltype.tab", log_dir, data_prefix, run_prefix)
  sumstats_run_path = sprintf("%s%s_%s_sumstats_run.tab", log_dir, data_prefix, run_prefix)
  params_run_path = sprintf("%s%s_%s_params_run.tab", log_dir, data_prefix, run_prefix)
  
  # set append param values
  append_sumstats_celltype = F
  append_sumstats_run = F
  append_params_run = F
  
  # write to file
  write.table(sumstats_celltype_df, file=sumstats_celltype_path, quote = F, sep = "\t", row.names=F, append = append_sumstats_celltype, col.names = !append_sumstats_celltype)
  write.table(sumstats_run_df, file=sumstats_run_path, quote = F, sep = "\t", row.names=F, append = append_sumstats_run, col.names = !append_sumstats_run)
  write.table(params_run_df, file=params_run_path, quote = F, sep = "\t", row.names=F, append = append_params_run, col.names = !append_params_run)
  
  ######################################################################
  ################# SAVE SESSION IMAGE AND FINISH ######################
  ######################################################################
  
  message("Saving final session image")
  
  save.image(file=sprintf("%s%s_%s_final_session_image.RData.gz", RObjects_dir, data_prefix, run_prefix), compress = "gzip")
  
  checkpoint_paths <- dir(path = scratch_dir, pattern = sprintf("%s_%s_checkpoint_?_image", data_prefix, run_prefix), full.names = T)
  if (length(checkpoint_paths)>0) for (ch in checkpoint_paths) if (file.exists(ch)) file.remove(ch)
  
  message("Script DONE!")
}

