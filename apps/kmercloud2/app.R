library(tidyverse)
library(RColorBrewer)
library(wordcloud)
library(wordcloud2)

# Define UI for app that draws a word cloud
ui <- fluidPage(
  shinyjs::useShinyjs(),
  # App title
  titlePanel("k-mer word clouds"),
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      # Select strat_size, hardcoded available strat_sizes only
      p("Note: ", strong("data/kmers/"), "must be populated with the prescribed, pre-generated kmer files."),
      p("Only the top 200 most frequent k-mers will be displayed."),
      selectInput(inputId = "strat_size",
                  label = "Select stratum size",
                  choices = list(100,250,500,750,1000,1500,2000)),
      textOutput(outputId = "numsamples"),
      hr(),
      radioButtons(inputId = "k",
                   label = "Select k-mer size",
                   choices = list(3,5,7),
                   inline = TRUE),
      tags$div(title="This corresponds to the strain metadata attribute. You can type here.",
               selectizeInput("sample_name", label = "Select sample name", choices = NULL)
      ),
      tags$div(title="Whether to compute the wordcloud using the mean of the entire k-mer matrix",
               checkboxInput(inputId = "show_all",
                             label = "Show wordcloud for full dataset",
                             value = FALSE)
      ),
      hr(),
      tableOutput(outputId = "summary")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      wordcloud2::wordcloud2Output(outputId = "wordcloud")

    )
  )
)

# Define server logic required to compute and render wordcloud ----
server <- function(input, output, session) {
  kmer_df <- reactiveVal()
  
  data <- reactive({
    path <- sprintf("data/kmer_%s_%s.csv", input$k, input$strat_size)
    kmer_df(readr::read_csv(path))
    strains <- lapply(as.list(dplyr::select(kmer_df(), strain)), sort)
    updateSelectizeInput(session, "sample_name", choices = strains, server = TRUE)
    output$numsamples <- renderText(sprintf("This stratum size yields %s samples.", nrow(kmer_df())))
  })
  
  figure <- reactive({
    data()
    df <- kmer_df()
    if(!is.null(input$sample_name)) {
      if(!input$show_all) {
        sample <- df %>%
          dplyr::filter(strain == input$sample_name) %>%
          dplyr::select(!(strain:length(df)))
        sample <- t(sample)
        sample <- sample[order(sample,decreasing=TRUE),]
        summary_tbl <- as.data.frame(names(sample))
        summary_tbl["frequency"] <- sample
        if (length(colnames(summary_tbl)) == 2)
          colnames(summary_tbl) <- c("kmer", "frequency")
        output$summary <- renderTable({summary_tbl})
        if(length(sample)>0) {
          wordcloud2::wordcloud2(summary_tbl,
                                 color=RColorBrewer::brewer.pal(8, "Dark2"))
        }
      } else {
        sample <- df %>%
          dplyr::select(!(strain:length(df))) %>%
          dplyr::summarise(dplyr::across(dplyr::everything(), mean))
        sample <- t(sample[1,])
        sample <- sample[order(sample,decreasing=TRUE),]
        summary_tbl <- as.data.frame(names(sample))
        summary_tbl["frequency"] <- sample
        colnames(summary_tbl) <- c("kmer", "frequency")
        output$summary <- renderTable({summary_tbl})
        if(length(sample)>0) {
          wordcloud2::wordcloud2(summary_tbl,
                                 color=RColorBrewer::brewer.pal(8, "Dark2"),
                                 shuffle=FALSE)
        }
      }
    }
  })
  
  observe({
    if((input$show_all == TRUE)) {
      shinyjs::disable("sample_name")
    } else {
      shinyjs::enable("sample_name")
    }
  })
  
  output$wordcloud <- wordcloud2::renderWordcloud2({figure()})
  
}

shinyApp(ui = ui, server = server)