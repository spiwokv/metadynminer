#' Calculate probability of free energy surface
#'
#' `prob.fes` calculates probability from free energy in a fes object.
#'
#' @param inputfes fes object.
#' @param temp temperature in Kelvins (default 300).
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default).
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' print(prob(tfes))
prob.fes<-function(inputfes, temp=300, eunit="kJ/mol") {
  if(class(inputfes)=="fes") {
    if(eunit=="kJ/mol") {
      if(fes1$dimension==1) {
        cfes<-list(fes=exp(-1000*inputfes$fes/8.314/temp), hills=fes1$hills, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
      }
      if(fes1$dimension==2) {
        cfes<-list(fes=exp(-1000*inputfes$fes/8.314/temp), hills=fes1$hills, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
      }
    } else if (eunit=="kJ/mol") {
      if(fes1$dimension==1) {
        cfes<-list(fes=exp(-1000*4.184*inputfes$fes/8.314/temp), hills=fes1$hills, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
      }
      if(fes1$dimension==2) {
        cfes<-list(fes=exp(-1000*4.184*inputfes$fes/8.314/temp), hills=fes1$hills, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
      }
    } else {
      stop("Error: Wrong eunit")
    }
    class(cfes) <- "fes"
    return(cfes)
  } else {
    stop("Error: Input must be fes object")
  }
}
