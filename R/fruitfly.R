#' Fruit fly data 
#'
#'This data set consist of the eggs laid from 667 female fruit flies. This is a subset of the original dataset collected in Dr. Carey's laboratory at University of California, Davis. 
#'A description of the experiment can be found in Carey et al.(1998). We processed the original data set  avaiable at https://anson.ucdavis.edu/~mueller/data/medfly1000.html to include the ones that lived past 30 days. 
#'The objective is to explore the relationship of egg laying pattern and longevity of fruit flies. 
#'
#
#'
#' @format A list with 3 fields:
#' \describe{
#'   \item{eggs}{a matrix of  667 by 30 representing the number of eggs laid at each day.}
#'   \item{lifetime}{an integer giving the lifetime (in days).}
#'   \item{eggs30}{an integer of the total number of eggs laid past 30 days.}
#'   ...
#' }
#' @source Carey, J.R., Liedo, P., Müller, H.G., Wang, J.L., Chiou, J.M. (1998). Relationship of age patterns of fecundity to mortality, longevity, and lifetime reproduction in a large cohort of Mediterranean fruit fly females. J. of Gerontology --Biological Sciences 53, 245-251.
"fruitfly"