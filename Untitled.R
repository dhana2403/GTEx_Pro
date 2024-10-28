# Load the shiny package
library(shiny)

# Define the UI
ui <- fluidPage(
  titlePanel("Simple Shiny App"),
  sidebarLayout(
    sidebarPanel(),
    mainPanel(
      h1("Hi")  # Display "Hi" as a heading
    )
  )
)

# Define the server logic (not needed for this simple app)
server <- function(input, output) {}

# Run the application 
shinyApp(ui = ui, server = server)
