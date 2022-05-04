#' Title
#'
#' @param series 
#' @param indicator 
#' @param differencing Not implemented yet
#' @param conversion 
#' @param conversion.obsposition 
#' @param outliers 
#'
#' @return
#' @export
#'
#' @examples
denton.modelbased<-function(series, indicator, differencing=1, conversion=c("Sum", "Average", "Last", "First", "UserDefined"), conversion.obsposition=1, 
                            outliers=NULL){
  conversion=match.arg(conversion)
  jseries=rjd3toolkit::ts_r2jd(series)
  if (is.null(outliers)){
    odates=.jcast(.jnull(), "[Ljava/lang/String;")
    ovars=.jnull("[D")
  }else{
    odates=.jarray(names(outliers))
    ovars=.jarray(as.numeric(outliers))
  }
  jindicator<-rjd3toolkit::ts_r2jd(indicator)
  jrslt<-.jcall("demetra/benchmarking/r/TemporalDisaggregation", "Ljdplus/tempdisagg/univariate/ModelBasedDentonResults;",
                "processModelBasedDenton", jseries, jindicator, as.integer(1), conversion, as.integer(conversion.obsposition), odates, ovars)
  # Build the S3 result
  estimation<-list(
    disagg=rjd3toolkit::proc_ts(jrslt, "disagg"),
    edisagg=rjd3toolkit::proc_ts(jrslt, "edisagg"),
    biratio=rjd3toolkit::proc_ts(jrslt, "biratio"),
    ebiratio=rjd3toolkit::proc_ts(jrslt, "ebiratio")
  )
  likelihood<-rjd3toolkit::proc_likelihood (jrslt, "ll")
  
  return(structure(list(
    estimation=estimation,
    likelihood=likelihood),
    class="JD3MBDenton"))
}
