#' Medicial decision data
#'
#' A dataset containing the binary responses and associated response times for
#' an experiment where experts and novices were presented with images of cells
#' and each participant had to decide whether the pictured cell was "blast" or
#' "non-blast"
#' 
#' @docType data
#' @keywords dataset
#' @name med_dec
#' @usage med_dec
#'
#' @format A data frame with 11000 rows and 9 variables:
#' \describe{
#'   \item{id}{identification number of the participant}
#'   \item{group}{expertise of participant; either "expert" or "novice"}
#'   \item{block}{block number}
#'   \item{trial}{index of trial for each participant}
#'   \item{classification}{true classification of the pictured cell; i.e. the correct response}
#'   \item{difficulty}{adjudged difficulty of the task for the particular image}
#'   \item{response}{response given by the participant; either "blast" or "non-blast"}
#'   \item{rt}{the response time associated with the response, in seconds}
#'   \item{stimulus}{the image file used for the specific trial}
#' }
#' 
#' @source Henrik
"med_dec"