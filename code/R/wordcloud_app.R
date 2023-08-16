library(shiny)

# NOTE: This app currently only works when launched from pgc-perf-opt.Rproject.
# TODO: Publish to shinyapps.io with data loaded from somewhere else.

# Define UI for app that draws a word cloud
ui <- fluidPage(
  # App title
  titlePanel("k-mer word clouds"),
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      # Select strat_size, hardcoded available strat_sizes only
      selectInput(inputId = "strat_size",
                  label = "Select strat_size",
                  choices = list("100 = 500 samples"=100,
                                 "250 = 1239 samples"=250,
                                 "500 = 2458 samples"=500,
                                 "750 = 3687 samples"=750,
                                 "1000 = 4877 samples"=1000,
                                 "1500 = 7224 samples"=1500,
                                 "2000 = 9408 samples"=2000)),
      radioButtons(inputId = "k_value",
                   label = "Select k-mer size",
                   choices = list(3,5,7),
                   inline = TRUE),
      checkboxInput(inputId = "show_all",
                    label = "Show word cloud for full dataset",
                    value = FALSE),
      checkboxInput(inputId = "strain",
                    label = "Select")
      
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
  
  df <- kmer_df %>% dplyr::select(!(strain:ncol(kmer_df)))

  output$distPlot <- renderPlot({
    
    # remove metadata columns
    test <- t(df[1,])
    test2 <- test[order(test[,1],decreasing=TRUE),]
    wordcloud(words = names(test2), freq = test2, min.freq = 1,
              max.words=200, random.order=FALSE, rot.per=0.35,
              colors=brewer.pal(8, "Dark2"))
  })
  
  output$value <- renderPrint({ input$test1 })
  
}

shinyApp(ui = ui, server = server)