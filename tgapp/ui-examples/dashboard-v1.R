## app.R ##
library(shiny)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title="dashboard 1"),
   dashboardSidebar(
    sidebarMenu(
      menuItem("One", tabName = "tabOne", icon = icon("dashboard")),
      menuItem("Two", tabName = "tabTwo", icon = icon("th"))
      )
  ),
  dashboardBody(
    tabItems(

      tabItem(tabName = "tabOne",
        h2("Tab one")
        ),
      tabItem(tabName = "tabTwo",
        h2("Tab two")
      )
    )
  )
)

server <- function(input, output) { }

app <- shinyApp(ui, server)
