.checkIfPkgInstalled <- function(pkg) {
    if (!suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))) {
        stop("'", pkg, "' package is not currently installed.\n", "  You first need to install it, which you can do with:\n", 
            "      library(BiocInstaller)\n", "      biocLite(\"", pkg, "\")")
    }
}

## Function to get the mapping of genes to a given vector of probe sets
.get.annotation <- function(IDs = "", chip = NULL, cols.vec = c("SYMBOL", "ENTREZID")) {
    if (is.null(chip)) 
        stop("Invalid chip name 'NULL', please provide a valid chip name!")
    
    .checkIfPkgInstalled(paste(chip, "db", sep = "."))
    
    ann.df <- suppressWarnings(select(get(paste(chip, "db", sep = ".")), keys = IDs, columns = cols.vec))
    dupl.ps <- unique(ann.df[, 1][duplicated(ann.df[, 1])])
    dupl.inds <- which(ann.df[, 1] %in% dupl.ps)
    dupl.df <- ann.df[dupl.inds, ]
    if (length(cols.vec) == 2) {
        temp.list <- lapply(split(seq(nrow(dupl.df)), dupl.df$PROBEID), function(.d) {
            c(PROBEID = dupl.df$PROBEID[.d[1]], SYMBOL = paste(dupl.df$SYMBOL[.d], collapse = ","), ENTREZID = paste(dupl.df$ENTREZID[.d], 
                collapse = ","))
        })
    } else {
        temp.list <- lapply(split(seq(nrow(dupl.df)), dupl.df$PROBEID), function(.d) {
            c(PROBEID = dupl.df$PROBEID[.d[1]], SYMBOL = paste(dupl.df$SYMBOL[.d], collapse = ","), ENTREZID = paste(dupl.df$ENTREZID[.d], 
                collapse = ","), GENENAME = paste(dupl.df$GENENAME[.d], collapse = ","))
        })
    }
    temp.df <- data.frame(t(sapply(temp.list, c)))
    result.df <- rbind(ann.df[-dupl.inds, ], temp.df)
    m.inds <- match(IDs, result.df[, 1])
    result.df <- result.df[m.inds, ]
    rownames(result.df) <- IDs
    return(result.df)
}


### Function to test if a given probe set is a potential marker This function is designed for internal use
.isMarker <- function(named.vec, rep.vec) {
    sort.vec <- sort(named.vec, decreasing = TRUE)
    sv.len <- length(sort.vec)
    names.vec <- names(sort.vec)
    rep.fe <- unname(rep.vec[names.vec[1]])
    poss.marker <- length(unique(names.vec[1:rep.fe])) == 1
    if (poss.marker) {
        mean.vec <- c()
        sort.num <- as.numeric(sort(named.vec, decreasing = TRUE))
        mean.vec <- c(mean.vec, mean(sort.num[1:rep.fe]))
        start.p <- rep.fe + 1
        rep.se <- unname(rep.vec[names.vec[start.p]])
        end.pos <- start.p + rep.se - 1
        cp.found <- FALSE
        while (end.pos <= sv.len && !cp.found) {
            len.sa <- length(unique(names.vec[start.p:end.pos]))
            if (sum(unique(names.vec[start.p:end.pos]) %in% names.vec[(end.pos + 1):sv.len]) == 0 | end.pos == sv.len) {
                mean.vec <- c(mean.vec, mean(sort.num[start.p:end.pos]))
                cp.found <- TRUE
            } else {
                
                end.pos <- start.p + sum(rep.vec[unique(names.vec[start.p:end.pos])]) - 1
                
            }
            
        }
        return(c(names.vec[1], (mean.vec[2]/mean.vec[1])))
    }
    
} 
