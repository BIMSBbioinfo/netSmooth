random.walk <- function(f0, A, lambda, max.iter=1000, tol=1e-3) {
    e <- Inf
    f <- f0
    for(k in 1:max.iter){
        f_next <- lambda * A %*% f + (1-lambda) * f0
        e <- Matrix::norm(f_next-f)
        f <- f_next
        if(e < tol){
            return(f)
        }
    }
    warning("Did not converge. Make sure `A` is row-normalized. Try increasing `max.iter`.")
    return(f)
}

l1.normalize.rows <- function(A) {
    rs = Matrix::rowSums(A)
    factors = 1 / replace(rs, rs==0, 1)
    return(Matrix::Diagonal(x=factors) %*% A)
}

