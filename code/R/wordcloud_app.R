library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("k-mer word clouds"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      checkboxInput(inputId = "test1",
                    label = "Show word cloud for all samples",
                    value = FALSE),
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "bins",
                  label = "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30),
      numericInput(inputId = "test",
                   label = "test",
                   value = 1,
                   step = 10)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  ),
  hr(),
  fluidRow(column(3, verbatimTextOutput("value"))),
)

# Define server logic required to compute and render wordcloud ----
server <- function(input, output) {

  output$distPlot <- renderPlot({
    
    
    
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#007bc2", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
    
  })
  
  output$value <- renderPrint({ input$test1 })
  
}

shinyApp(ui = ui, server = server)