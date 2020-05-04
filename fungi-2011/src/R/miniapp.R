#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)

# Define UI for application that draws a histogram
ui <- fluidPage(
    plotlyOutput("plot"),
    verbatimTextOutput("click"),
    verbatimTextOutput("brush")
)

server <- function(input, output, session) {
    frame1 <- data.frame()
    frame2 <- data.frame()
    
    
    df3D <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
    
    nms <- row.names(pca$x)
    
    output$plot <- renderPlotly({
        p <- ggplot(df3D, aes(x = PC1, y = PC2)) + geom_point()
        ggplotly(p) %>% layout(dragmode = "lasso")
    })
    
    output$click <- renderPrint({
        d <- event_data("plotly_click")
        if (!is.null(d)) {
            frame1 <<- frame1[is.null(frame1$pointNumber), ] # Optional line to remove the previous selections
            frame1 <<- rbind(frame1, d) 
        }
        frame1
    })
    
    output$brush <- renderPrint({
        d <- event_data("plotly_selected")
        if (!is.null(d)) {
            frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
            frame2 <<- rbind(frame2, d)
        }
        frame2
        
    })
    
}

shinyApp(ui, server)