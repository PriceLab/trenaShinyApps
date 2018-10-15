## app.R ##
library(shiny)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title="dashboard 1"),
   dashboardSidebar(
     sidebarMenu(id="sidebarMenu",
        menuItem("One", tabName = "tabOne", icon = icon("dashboard")),
        menuItem("Two", tabName = "tabTwo", icon = icon("th"))
        ),
     conditionalPanel(
        condition = "input.sidebarMenu == 'tabTwo'",
        selectInput("period", "Period:",
                    choices = list("Years" = 1, "Months" = 2))
        ) # conditionalPanel
     ), # dashboarSidebar

  dashboardBody(
    tabItems(
       tabItem(tabName = "tabOne",
          h2("Tab one")),
       tabItem(tabName = "tabTwo",
          h2("Tab two"))
       ) # tabItems
     ) # dashboardBody
  ) # dashboarPage

server <- function(input, output) { }

app <- shinyApp(ui, server)
