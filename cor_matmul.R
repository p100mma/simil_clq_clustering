corfast<- function (D) 
{
    stopifnot(any(class(D) %in% c("matrix")))
    D_means <- colMeans(D)
    D_sds <- colSds(D)
    tZ <- (t(D) - D_means)/(D_sds)
    (1/(nrow(D) - 1)) * (tZ %*% t(tZ))
}

corpval<- function (r, n) 
{
    t_stat = sqrt(n - 2) * r/sqrt(1 - r^2)
    2 * pmin(pt(t_stat, n - 2), pt(t_stat, n - 2, lower.tail = FALSE))
}

fastPearsonData<- function (D, NA.warn = TRUE, p.adjust.method = "holm") 
{
    if (NA.warn) 
        if (any(is.na(D))) 
            warn("warning: D contains NAs. Consider using rcorrData")
    r = corfast(D)
    ltr <- lower.tri(r)
    corData <- list()
    corData$N = ncol(D)
    corData$r = r[ltr]
    corData$P = corpval(corData$r, nrow(D))
    corData$Pa = p.adjust(corData$P, method = p.adjust.method)
    corData
}

similarity_matrix<- function (corData, addLoops = TRUE, power = 2, level = 0.05, 
    return.lower.tri = FALSE) 
{
    if (!return.lower.tri) {
        S <- matrix(nrow = corData$N, ncol = corData$N)
        if (addLoops) 
            diag(S) = 1
        else diag(S) = 0
        ltr = lower.tri(S)
    }
    S_ltr = corData$r
    S_ltr[corData$Pa >= level] = 0
    S_ltr = abs(S_ltr)^power
    if (!return.lower.tri) {
        S[ltr] = S_ltr
        S <- t(S)
        S[ltr] = S_ltr
        return(S)
    }
    else {
        return(S_ltr)
    }
}

