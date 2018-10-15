library(shiny)
library(shinydashboard)
rm(ui)
rm(server)

ui <- dashboardPage(
   dashboardHeader(title = "dynamic sidebar menus"),
   dashboardSidebar(
    sidebarMenu(id = "tab",
       menuItem("1", tabName = "1"),
       menuItem("2", tabName = "2")),
    uiOutput("out1")),

    dashboardBody(
      tabItems(
         tabItem(tabName="igvAndTable",
            fluidRow(
               column(width=9, id="igvColumn",
                  igvShinyOutput('igvShiny', height="800px")
                  ),
               column(width=3, id="dataTableColumn", title="mtcars",
                  DTOutput("table")
                  #actionButton("button2", "mainPanel button")
                  )
               ) # fluidRow
            ), # tabItem 1
         tabItem(tabName="buildModels",
             h4("Build trena regulatory model from DNase footprints in the genomic region currently displayed in the IGV view, with these constraints:"),
             br(),
             fluidRow(
                column(2, offset=1, checkboxGroupInput("footprintDatabases", "Footprint Databases",
                                            footprintDatabases, selected=footprintDatabases[1])),
                column(2, checkboxGroupInput("intersectWithRegions", "Intersect footprints with:",
                                             c("GeneHancer" = "genehancer",
                                               "Encode DHS" = "encodeDHS"))),
                column(3, selectInput("expressionSet", "Choose gene expression dataset:",
                                      expressionMatrixNames))),
             br(),
             fluidRow(column(width=1, offset=3, actionButton("buildModelButton", "Build model")))
            ) # tabItem 2
         ) # tabItems
      ) # dashboardBody
) # dashboardPage

server <- function(input, output) {

  output$out1 <- renderUI({

    if (input$tab == "1") {

      dyn_ui <- tabsetPanel(id = "tabset_id", selected = "t1",
                            tabPanel("t1", value = "t1"),
                            tabPanel("t2", value = "t2"))

    }
    if (input$tab == "2") {

      dyn_ui <- list(selectInput("s1", label = "Select", choices = letters[1:3]),
                     selectInput("s2", label = "Select2", choices = letters[4:6]))
    }
    return(dyn_ui)
  })
}

app <- shinyApp(ui, server)
