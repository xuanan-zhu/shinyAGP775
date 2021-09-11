#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @noRd
app_ui <- function(request) {
  tagList(
    dashboardPage(skin = "blue",
                  dashboardHeader(title = "AGP775"),
                  dashboardSidebar(
                    sidebarMenu(menuItem("global alignments",tabName = "GA"),
                                menuItem("Affine Gap Penalties",tabName = "AGP")),
                    textInput("data1","seq1:"),
                    textInput("data2","seq2:")
                  ),
                  dashboardBody(
                    tabItems(
                      tabItem(tabName = "GA",
                              numericInput("m1","march :score:",value=0),
                              numericInput("rm1","dismarch :score:",value=0),
                              numericInput("gap1","gap score:",value=0),
                              actionButton("submit1","submit now"),
                              tableOutput("r1"),
                              tableOutput("r1seq")
                              
                      ),
                      tabItem(tabName = "AGP",
                              numericInput("m2","march :score:",value=0),
                              numericInput("rm2","dismarch :score:",value=0),
                              numericInput("gap2","gap score:",value=0),
                              numericInput("gapp2","gapgap score:",value=0),
                              actionButton("submit2","submit now"),
                              tableOutput("r2"),
                              tableOutput("r2seq")
                      )
                    )
                  )
    )
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'shinyAGP775'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

