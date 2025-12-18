#' Use clustered coefficients models in area-level models
#' @import sae
#' @import cluster
#' @import scclust
#' @import dplyr
#' @export
#' @param dom 1:m, m is the number of domains
#' @param y direct estimator
#' @param x covaraites including intercept
#' @param index_d the index with the cluster structure
#' @param vardir the direct variance
#' @param lamvec the sequence of tuning parameters lambda
#' @examples
#' library(CLAreaModel)
#' data(simdata)
#'lamvec <- 10^seq(-1.05, -0.5, length.out=50)
#'res_clust <- Est_clust_area2(dom = 1:nrow(df_area),
#'                              y = df_area$y, x = cbind(1, df_area$x),
#'                              index_d = 1:2, vardir = df_area$vardir,
#'                              lamvec = lamvec)
#'estCL <- res_clust$area_estm ### area estimate
#'clusterm <- res_clust$cluster_estm ## estimated number
#'
Est_clust_area2 <- function(dom, y, x, index_d, vardir, lamvec,
                           nu = 1, gam = 3,
                           lam0 = 0.001, maxiter= 500, tol = 1e-4,
                           seed = 1234, init_control = NULL, beta00 = NULL)
{
  n0 <- length(y)
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  res_bic <- Clust_area_pw_bic2(indexy=dom, y, x, index_d, vardir,
                               lamvec, nu = nu, gam = gam,
                               lam0 = lam0, maxiter= maxiter, tol = tol,
                               seed = seed, init_control = init_control,beta00 = beta00)

  npp <- length(unique(dom))

  bic_df <- res_bic$res_df %>%
    mutate(mbic = loss0 + log(npp)*log(n0)*(Khat*len_d + nz +1))

  res_fm <- res_bic$out[[which.min(bic_df$mbic)]]
  area_estm <- res_fm$area_est

  fit_fm <- NULL
  tryCatch({
    fit_fm <- fit_clust_area(dom = dom, y = y, x = x, index_d = index_d,
                             vardir = vardir, cluster = res_fm$cluster)
    area_estm <- fit_fm$area_est_ada
  }, error = function(e) {
    area_estm <- res_fm$area_est
  })


  out <- list(bic_df = bic_df,
              resfm = res_fm,
              fitm = fit_fm,
              cluster_estm = res_fm$cluster,
              area_estm = area_estm,
              outlist = res_bic$out)

  return(out)

}



