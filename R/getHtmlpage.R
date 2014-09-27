getHtmlpage <- function(markers.list, chip, directory = getwd()) {
    
    annot.df <- suppressWarnings(select(get(paste(chip, "db", sep = ".")), keys = keys(get(paste(chip, "db", sep = "."))), 
        columns = c("SYMBOL", "ENTREZID", "GENENAME")))
    
    for (i in seq_along(markers.list)) {
        u.list <- unlist(strsplit(markers.list[[i]], split = " : "))
        ps.vec <- u.list[seq(1, length(u.list), 4)]
        scores1 <- u.list[seq(4, length(u.list), 4)]
        names(scores1) <- ps.vec
        
        ps.inds <- which(annot.df[, 1] %in% ps.vec)
        temp.ps <- unname(annot.df[ps.inds, 1])
        new.inds <- match(temp.ps, ps.vec)
        scores2 <- scores1[new.inds]
        sort.scores2 <- sort(scores2, decreasing = FALSE)
        ii <- which(annot.df[, 1] %in% names(scores2))
        temp.df <- annot.df[ii, ]
        temp.df <- cbind(annot.df[ii, ], SCORE = unname(sort.scores2[temp.df[, 1]]))
        res.df <- temp.df[order(temp.df$SCORE), ]
        
        suppressWarnings(htmlpage(res.df[, c(1, 2, 3)], paste(directory, paste(names(markers.list[i]), "html", sep = "."), 
            sep = "/"), " ", table.head = c("Probe Set", "Gene Symbol", "Entrez Id", "Gene Name", "Score"), othernames = res.df[, 
            c(4, 5)], repository = list("affy", "gb", "en")))
    }
}
 
