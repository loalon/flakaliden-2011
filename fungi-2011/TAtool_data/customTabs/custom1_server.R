print("tab 1 server code")


observeEvent(input$custom1_conditionSelect, {

  properties <- c(customData1$taxa, "Trophic.Mode", "Guild")
  updateSelectInput(session, 'custom1_propertySelect', choices = properties )
  if("species" %in% properties)
    updateSelectInput(session, 'custom1_propertySelect', selected = "species" )

  
})

observeEvent(input$custom1_propertySelect, {
  
  thisValues <- sort(unique(customData1$megaData[[input$custom1_propertySelect]]))
  updateSelectInput(session, 'custom1_subpropertySelect', choices = thisValues)
  
  if(input$custom1_propertySelect == "species")
    updateSelectInput(session, 'custom1_subpropertySelect', selected = "Cortinarius glaucopus" )
})



observeEvent(input$custom1_loadBtn, {

  print(paste("condition:", input$custom1_conditionSelect))
  print(paste("property:", input$custom1_propertySelect))
  print(paste("subproperty:", input$custom1_subpropertySelect))
  print(paste("group:", input$custom1_selectGroup))

  mini <- isolate(customData1$megaData[customData1$megaData[[input$custom1_propertySelect]] == input$custom1_subpropertySelect,])

  if(isolate(input$custom1_conditionSelect) == 'control'){

    mini <- mini[grep("W\\d+C",colnames(mini))] 
  } else if (isolate(input$custom1_conditionSelect) == 'fertilised'){
    print("aqui fertilised")
    mini <- mini[grep("W\\d+F",colnames(mini))] 
  } else {
    mini <- mini[grep("W\\d+[CF]",colnames(mini))] 
  }

  miniMeta <- data.frame(treatment=substr(colnames(mini), 4, 4),
                         week = substr(colnames(mini), 2, 3))
  miniMeta$treatment <- ifelse(miniMeta$treatment == "C", "ND", "NE")
  

  withProgress(message = 'Processing', detail = "Might take a while",value = 0, {
    
    vsd <- tryCatch({
      print("dds")
      dds <- DESeqDataSetFromMatrix(mini, miniMeta, design = ~1)
      dds <- dds[(rowSums(counts(dds)) > 0),]
      dds <- dds[,(colSums(counts(dds)) > 0)]
      dds$treatment
      dds <- estimateSizeFactors(dds, type="poscounts")
      print("transform")
      vsd <- varianceStabilizingTransformation(dds)
      miniMeta <- data.frame(treatment=vsd$treatment,
                             week = vsd$week)
      
      
      vsd <- t(assay(vsd))
    },error=function(e){
      getDialog("This data can't be processed")
      return(NULL)
    })

  })
  

  
  #dds <- estimateSizeFactors(dds)
  # sizes <- sizeFactors(ddsMatrix)
  # pander(sizes)
  # boxplot(sizes, main="Sequencing libraries size factor")
  
  # vsd.QA <- varianceStabilizingTransformation(dds, blind=TRUE)
  # assay(vsd.QA) <- assay(vsd.QA) - min(assay(vsd.QA))
    output$custom1_pcaPlot <- renderPlotly({
      print("plot")
      withProgress(message = 'Plotting 2D PCA', detail = "Might take a while",value = 0, {
        if(!is.null(vsd)) {
         ggplotly (
           plotInteractivePCA(vsd, mutate_all(miniMeta, as.character), 
                     group=isolate(input$custom1_selectGroup)) 

           )
        } else{ 
          plot.new()
        }
      })
    })
    
    output$custom1_pca3dPlot <- renderPlotly({
      withProgress(message = 'Plotting 3D PCA', detail = "Might take a while", value = 0, {
      if(!is.null(vsd)) {
        plotInteractivePCA(vsd, mutate_all(miniMeta, as.character), group=isolate(input$custom1_selectGroup), model="3D")
      } else{ 
        plot.new()
      }
      })
    })
  
  
})

#getDialog("after")

