# CUSTOM1 Tab
tabPanel(span("Fungi explorer", title="CUSTOM1"), value = "tab_custom1",
   h2("PCA"),
   p("Interactive Principal Component Analysis"),
   sidebarLayout(
     sidebarPanel(
       selectInput("custom1_conditionSelect", "Select condition", c("ND+NE", 
                                                                    "ND", 
                                                                    "NE")),
       # 
       selectInput("custom1_propertySelect", "Select property", c()),
       selectInput("custom1_subpropertySelect", "Select subproperty", c()),
       selectizeInput("custom1_selectGroup", "Select group", choices=c("week", "treatment"), multiple=TRUE, selected = c("week", "treatment"),options = list(maxItems = 2)),
       # tags$hr(style="border-color: black;"),
       actionButton('custom1_loadBtn', "Load"),
       width=4
     ),
     
     mainPanel(
       fluidRow(
         h4("2D PCA"),
         plotlyOutput("custom1_pcaPlot"),
           h4("3D PCA"), 
          plotlyOutput("custom1_pca3dPlot")
       )
     )
   )
)
         