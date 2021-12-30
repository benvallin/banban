#' Adjust choose_on_aic() output for multiple comparisons
#'
#' @description get_sign_mod() adjusts the p value related to the F statistic of each "best model" provided by choose_on_aic().
#'
#' The function first filters out every best model containing only the intercept (1) as independent variable (dependent variable ~ mean(dependent variable)).
#'
#' Adjustment for multiple comparisons then relies either on false discovery rate (FDR) using the Benjamini-Hochberg procedure, or on one of the following adjustment methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'
#' Best models receive a label "best_is_sign" = TRUE if the p value related to their F statistic reaches statistical significance after adjustment (passed FDR correction, or adjusted value <= 0.05) AND the p value of the t statistic related to the independent variable of interest <= 0.05.
#'
#' @param tibble a tibble produced by choose_on_aic().
#' @param method a string indicating the adjustment method to be used. Must be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none".
#' @param FDR a numeric indicating the q-value threshold to use for FDR.
#'
#' @return A tibble containing, for each combination of dependent variable, independent variable of interest and covariates, i) the corresponding best model, ii) its associated statistics (F p value, t p value) and iii) the label "best_is_sign" = TRUE/FALSE.
#'
#' Combinations of dependent variable, independent variable of interest and covariates for which the best model contains only the intercept as independent variable are discarded, so the output may be shorter than the input tibble.
#'
#' @export
#'
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom rlang !!
#' @importFrom rlang enexpr
#' @importFrom stats p.adjust
#' @importFrom stringr str_c
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#'
#' @examples
#'
get_sign_mod <- function(tibble, method, FDR = 0.2) {

  ### START mcp.project ###

  reject=function(sorted,criticals){
    m=length(sorted)
    stopifnot( length(criticals) == m )
    indicators= sorted<criticals # Marking p-values below the critical values

    if(!any(indicators)) return(list(cutoff=0,cut.index=0))

    cut.index=max((1:m)[indicators])

    cutoff=sorted[cut.index] #The highest rejected p-value

    return( list(cutoff=cutoff,cut.index=cut.index) )
  }

  bh.adjust=function(sorted,m,m0,constant=1){
    adjusted=rep(NA,m)
    temp.min=sorted[m]
    min.ind=rep(0,m)
    for (i in m:1) {
      temp= min(m0*sorted[i]*constant / i,1)
      if  ( temp <= temp.min  ) {
        temp.min = temp
        min.ind[i]=1
      }
      adjusted[i]=temp.min
    }
    return(adjusted)
  }

  bh=function(sorted, q, m, adjust = F, m0 = m, pi0, constant = 1){
    # cat("Calling bh \n")
    flush.console
    if (missing(m0) & !missing(pi0))
      m0 = pi0 * m
    else {
      criticals = (1:m) * q/(m0 * constant)
      cutoff = reject(sorted, criticals)
      rejected = sorted <= cutoff$cutoff
      adjusted = rep(NA, m)
      if (adjust)
        adjusted = bh.adjust(sorted, m = m, m0 = m0, constant = constant)
      multiple.pvals = data.frame(original.pvals = sorted,
                                  criticals = criticals, rejected = rejected, adjusted.pvals = adjusted)
      output = list(Cutoff = cutoff, Pvals = multiple.pvals)
      return(output)
    }
  }

  fdr=function(x,q,method='BH',pi0,lambda=0.5, m0){
    #Input checking
    stopifnot( is.numeric(x) , is.numeric(q) , any(x<1) , any(x>0) , q<1 , q>0 )
    m=length(x)
    if(!missing(m0)) stopifnot( is.numeric(m0) , m0<m , m0>=0 )
    else if(!missing(m0) & !missing(pi0)) stop("Only m0 or pi0 can be specified.")
    else if (method!='Oracle' && !missing(m0) ) stop("m0 can only be specified for the Oracle method.")
    else if(missing(m0) & !missing(pi0)) m0=pi0*m

    ranks=rank(x)
    sorted=sort(x)

    # Selecting the propoer procedure
    output=switch(method,
                  'BH'=bh(sorted=sorted,q=q,adjust=TRUE,m=m),
                  'General Dependency'=bh(sorted=sorted,q=q,m=m,constant=log(m),adjust=TRUE),
                  'Oracle'=bh(sorted=sorted,q=q,m=m,m0=m0,adjust=TRUE),
                  'BH Adaptive'=adaptive.bh(sorted=sorted,q=q,m=m),
                  'ST Adaptive'=adaptive.st(sorted=sorted,q=q,m=m,lambda),
                  'Two Stage'=two.stage(sorted=sorted,q=q,m=m),
                  'Multiple Step Down'=multiple.down(sorted=sorted,q=q,m=m),
                  'Multiple Step Up'=stop("Not implemented yet.") ,#look in file "multiple up (unused).r",
                  'Smoother'=stop("Not implemented yet.") ,
                  'Bootstrap'=stop("Not implemented yet."),
                  'Median ST Adaptive'=stop("Not implemented yet.")
    )

    #Arranging output
    output$Pvals=output$Pvals[ranks,]
    output=list(method=method,q=q,Cutoff=output[['Cutoff']],Pvals=output[['Pvals']])
    class(output)=c(class(output),'fdr')
    return(output)
  }

  ### END mcp.project ###

  method <- enexpr(method)

  output <- tibble %>% filter(!is.na(best_mod_f_p_value))

  if (method == "BH") {
    output <- output %>%
      nest(data = everything()) %>%
      mutate(fdr = map(data, ~ fdr(.x$best_mod_f_p_value, q = FDR, method = "BH")$Pvals)) %>%
      unnest(cols = c(data, fdr)) %>%
      rename(pass_fdr = rejected) %>%
      mutate(fdr_result = map_chr(pass_fdr, ~ if_else(isTRUE(.x),
                                                      str_c("sign (FDR = ", FDR, ")"),
                                                      str_c("ns (FDR = ", FDR, ")"))),
             best_is_sign = case_when(best_mod_t_p_value <= 0.05 & pass_fdr == TRUE ~ TRUE,
                                      TRUE ~ FALSE))
  } else {
    output <- output %>%
      mutate(adj_best_mod_f_p_value = p.adjust(best_mod_f_p_value, method = !!method),
             best_is_sign = case_when(best_mod_t_p_value <= 0.05 & adj_best_mod_f_p_value <= 0.05 ~ TRUE,
                                      TRUE ~ FALSE))
  }
}
