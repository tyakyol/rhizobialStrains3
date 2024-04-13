rarefyFilter = function (x, min = 0) 
{
    keep = c()
    if (min < 0) {
        stop("Min should be either 0 or positive.")
    }
    if (min == 0) {
        min = min(colsums = apply(x, 2, sum))
        print(paste("Rarefy to minimum count", min))
        keep = c(1:ncol(x))
    }
    else {
        colsums = apply(x, 2, sum)
        if (min(colsums) < min) {
            for (j in 1:ncol(x)) {
                if (colsums[j] >= min) {
                  keep = c(keep, j)
                }
            }
            print(paste("Number of columns", ncol(x)))
            print(paste("Keeping ", length(keep), " columns with column sums equal or above", 
                min))
            x = x[, keep]
        }
    }
    rar = t(vegan::rrarefy(t(x), min))
    res = list(rar, keep)
    names(res) = c("rar", "colindices")
    return(res)
}
