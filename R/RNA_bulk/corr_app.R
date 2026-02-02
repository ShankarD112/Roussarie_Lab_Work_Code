library(shiny)
library(shinyWidgets)
library(ggplot2)
library(reshape2)
library(DT)
library(pheatmap)
library(readr)

options(shiny.maxRequestSize = 1000 * 1024^2)

# Shiny app for exploring gene ↔ proteoform correlations.
# Expected input columns include:
#   - Gene, Proteoform, Pvalue, FDR_BH, Correlation, CorrelationMethod


# UI
ui <- fluidPage(
  titlePanel("Gene-Proteoform Correlations Viewer"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "file", 
        "Upload Correlation File:", 
        accept = c(".txt", ".csv")
      ),
      sliderInput(
        "pval_thresh", 
        "P-value Threshold:", 
        min = 0, 
        max = 0.05, 
        value = 0.05, 
        step = 0.01
      ),
      sliderInput(
        "fdr_thresh", 
        "FDR Threshold:", 
        min = 0, 
        max = 0.05, 
        value = 0.05, 
        step = 0.01
      ),
      sliderInput(
        "cor_thresh", 
        "Correlation Threshold:", 
        min = -1, 
        max = 1, 
        value = c(-1, 1),
        step = 0.1
      ),
      pickerInput(
        "method", 
        "Select Correlation Method:", 
        choices = c("pearson", "spearman", "kendall"),
        selected = "pearson",
        multiple = FALSE
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Heatmap",
          plotOutput("heatmap")
        ),
        tabPanel(
          "Filtered Table",
          dataTableOutput("table")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    file <- input$file$datapath
    read_delim(file, delim = "\t")
  })
  
  filtered_data <- reactive({
    req(data())
    df <- data()
    df <- df[df$Pvalue < input$pval_thresh &
               df[[paste0("FDR_BH")]] < input$fdr_thresh &
               df$CorrelationMethod == input$method &
               df$Correlation >= input$cor_thresh[1] &
               df$Correlation <= input$cor_thresh[2], ]
    df
  })
  
  heatmap_data <- reactive({
    req(filtered_data())
    df <- filtered_data()
    mat <- acast(df, Gene ~ Proteoform, value.var = "Correlation", fun.aggregate = mean, fill = NA)
    mat
  })
  
  output$heatmap <- renderPlot({
    req(heatmap_data())
    pheatmap(heatmap_data(),
             color = colorRampPalette(c("blue", "white", "red"))(50),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             main = paste("Heatmap for", input$method, "Correlations"))
  })
  
  output$table <- renderDataTable({
    req(filtered_data())
    datatable(
      filtered_data(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
}

# App
shinyApp(ui = ui, server = server)
