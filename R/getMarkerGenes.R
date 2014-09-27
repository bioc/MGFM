## Function to get marker genes.  Input parameters: data.mat: The microarray data matrix with probe sets corresponding
## to rows and samples corresponding to columns.  chip: chip name samples2compare: A character vector with the sample
## names to be compared (e.g. c('liver', 'lung', 'brain')). By default all samples are used.  annotate: A boolean value.
## If TRUE the gene symbol and the entrez gene id are shown.  score.cutoff: An integer value to filter the marker genes.
## Default is 1 ('no filtering').  output: A list with marker genes associated with each of the given sample types.
getMarkerGenes <- function(data.mat, samples2compare = "all", annotate = TRUE, chip = NULL, score.cutoff = 1) {
    if (length(samples2compare) > 1) {
        if (any(!samples2compare %in% colnames(data.mat))) {
            stop("The samples to compare should be included in the data matrix!!")
        } else {
            ii <- which(colnames(data.mat) %in% samples2compare)
            data.mat <- data.mat[, ii]
        }
    }
    markers.list <- list()
    rep.vec <- table(colnames(data.mat))
    res.list <- apply(data.mat, 1, .isMarker, rep.vec = rep.vec)
    res.list[sapply(res.list, is.null)] <- NULL
    mar.len <- length(unlist(res.list))
    samples.vec <- unname(unlist(res.list))[seq(1, mar.len, 2)]
    scores.vec <- round(as.numeric(unname(unlist(res.list))[seq(2, mar.len, 2)]), digits = 2)
    if (length(which(scores.vec <= score.cutoff)) == 0) {
        message("No markers found using the given score cut-off!!!")
        return(list())
    }
    
    mar.ps.vec <- names(res.list)
    names(scores.vec) <- mar.ps.vec
    u.snames <- unique(colnames(data.mat))
    if (annotate) {
        if (is.null(chip)) {
            stop("The chip name is required to map probe sets to genes")
        }
        if (length(grep("\\.db$", chip))) 
            chip <- substr(chip, 1, nchar(chip) - 3)
        annot.df <- .get.annotation(rownames(data.mat), chip = chip)
        if (dim(data.mat)[1] != dim(annot.df)[1]) {
            warning("Number of probe sets on the chip is not equal to the number of probe sets in         the data matrix!!!")
        }
        for (i in seq_along(u.snames)) {
            inds <- which(samples.vec == u.snames[i])
            sort.scores <- sort(scores.vec[inds])
            sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
            ps.inds <- match(names(sort.scores), annot.df[, "PROBEID"])
            markers.list[[i]] <- paste(names(sort.scores), annot.df[ps.inds, "SYMBOL"], annot.df[ps.inds, "ENTREZID"], 
                unname(sort.scores), sep = " : ")
        }
    } else {
        for (i in seq_along(u.snames)) {
            inds <- which(samples.vec == u.snames[i])
            sort.scores <- sort(scores.vec[inds])
            sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
            markers.list[[i]] <- paste(names(sort.scores), unname(sort.scores), sep = " : ")
        }
    }
    names(markers.list) <- paste(u.snames, "markers", sep = "_")
    return(markers.list)
} 
