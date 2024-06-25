#' @title runApp
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description Shiny App
#' 
#' @import DT
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import shinythemes
#' 
#' @param ui UI of the Shiny app
#' @param server Server of the Shiny app
#'
#' @return none
#' 
#' @example ./inst/examples/example.R
#' 
#' @export

runPDACMOC <- function(){
  shinyApp(ui, server, options = list(host = "0.0.0.0", port = 1995, launch.browser = FALSE))
}
