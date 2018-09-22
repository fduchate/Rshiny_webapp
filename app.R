# Load packages ----
#install.packages("shiny")
#install.packages("epinet")
#install.packages("ggtree")
#install.packages("shinythemes")
setwd("D:/PhD/GitHub/Rshiny_webapp/Rshiny_webapp/")
library(shiny)
library(epinet)
library(ggtree)
library(ape)
library(igraph)
library(visNetwork) 
library(shinycssloaders)
library(shinythemes)
library(shinyjs)
library(DT)
library(shinydashboard)
library(gsubfn)
library(Biostrings)
library(NbClust)
library(cluster)
library(dbscan)
library(factoextra)
library(data.table)
source("small_network.R")
########In dashboard######
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Epidemia simulation",icon = icon("laptop", lib = "font-awesome"), tabName = "episim",
             menuSubItem("Sim. 2 pop", tabName = "2pop"),
             menuSubItem("Sim. 4 pop. circ.", tabName = "4popcirc"),
             startExpanded = TRUE),
    menuItem("Phylogenetic analysis", tabName = "widgets",icon = icon("exchange alt",lib = "font-awesome"))
  )
)

body <- dashboardBody(shinyjs::useShinyjs(),
                      tabItems(
                        tabItem(tabName = "2pop",
                            fluidRow(column(width = 4,
                                box(
                                  title = "Inputs", status = "primary", solidHeader = TRUE,width = NULL,
                                  #        "Box content here", br(), "More box content",
                                  numericInput("num", label = h3("Number of hosts in population one"), value = 10,min = 0),
                                  numericInput("num2", label = h3("Number of hosts in population two"), value = 10,min = 0),
                                  actionButton("goButton", "Create the transmission network tree"),
                                  actionButton("checkbox", "Plot the transmission network"),
                                  actionButton("makesequence","Simulate the sequences")),
                                p("Timestamp: ",
                                  span(id = "time", date())
                                )),
                                column(width = 7,
                                tabBox(title = "Output",width = NULL,
                                       side = "right", height = "600px",selected = "Summary",
                                       tabPanel("Epidemia graph",textOutput("check_network2"),visNetworkOutput("network_graph")),
                                       tabPanel("Epidemia table", DT::dataTableOutput("epi_network"),downloadButton('downloadData','Download the epidemia network')),
                                       tabPanel("Population structure", DT::dataTableOutput("tpop_struct_red")),
                                       tabPanel("Summary", textOutput("selected_var"),textOutput("check_network"),textOutput("status"),textOutput("make_sequence"),
                                                imageOutput("name_tree")%>% withSpinner(color="#0dc5c1"),
                                                uiOutput("download"),uiOutput("download_dna"))
                                )#))
                        ))),
                        tabItem(tabName = "4popcirc",
                                fluidRow(column(width = 4,
                                box(title = "Inputs",width = NULL, status = "primary", solidHeader = TRUE,
                                    numericInput("num14", label = h3("Number of hosts in population one"), value = 10,min = 0),
                                    numericInput("num24", label = h3("Number of hosts in population two"), value = 10,min = 0),
                                    numericInput("num34", label = h3("Number of hosts in population three"), value = 10,min = 0),
                                    numericInput("num44", label = h3("Number of hosts in population four"), value = 10,min = 0),
                                    actionButton("goButton4", "Create the transmission network tree"),
                                    actionButton("checkbox4", "Plot the transmission network"),
                                    actionButton("makesequence4","Simulate the sequences")),
                                p("Timestamp: ",
                                  span(id = "time", date())
                                )),
                                column(width = 7,
                                tabBox(title = "Output", width = NULL,
                                       side = "right", height = "600px",selected = "Summary",
                                       tabPanel("Epidemia graph",textOutput("check_network24"),visNetworkOutput("network_graph4")),
                                       tabPanel("Epidemia table", DT::dataTableOutput("epi_network_circular"),downloadButton('downloadData4','Download the epidemia network')),
                                       tabPanel("Population structure", DT::dataTableOutput("pop_struct_red_circular")),
                                       tabPanel("Summary", textOutput("selected_var4"),textOutput("check_network4"),textOutput("status4"),textOutput("make_sequence4"),
                                                imageOutput("name_tree4")%>% withSpinner(color="#0dc5c1"),
                                                uiOutput("download4"),uiOutput("download_dna4")))))),
                        
                        tabItem(tabName = "widgets",
                                fluidRow(column(width = 4,
                                box(
                                  title = "Inputs",width = NULL, status = "primary", solidHeader = TRUE,
                                  actionButton("treeforestsequences","Subsampling the sequences"),                                  
                                  a(id = "toggleAdvanced", "Show/hide advanced info", href = "#"),
                                  shinyjs::hidden(
                                    div(id = "advanced",
                                        numericInput("numPerCat","Max number of sample per category after sampling",30),
                                        numericInput("repsPerCat","Number of repetitions per joint sampling",10),
                                        numericInput("jn", "Joint number per category", 30),
                                        numericInput("jr", "Joint repetitions", 30)
                                    )),
                                  actionButton("performbssvs","Perform a BSSVS analysis")),
                                p("Timestamp: ",
                                  span(id = "time", date())
                                )),
                                column(width = 7,
                                tabBox(
                                  title = "Output",width = NULL,side = "right",selected = "DNA sequences sampled",
                                  tabPanel("BSSVS analysis output",textOutput("status_dna3"),uiOutput("download_BSSVS"),DT::dataTableOutput("BSSVS_result_table"),visNetworkOutput("BSSVS_graph")),
                                  tabPanel("Trait sequences sampled",textOutput("status_dna2"),uiOutput("download_traits_table_samp"),DT::dataTableOutput("traits_table_samp")),
                                  tabPanel("DNA sequences sampled",textOutput("sampled_txt_output"),uiOutput("download_dna_sampled"),uiOutput("download_trees_sampled"),textOutput("status_dna"),DT::dataTableOutput("table_continent_samp"))
                                ))))))
 
# Put them together into a dashboardPage

ui <-   dashboardPage(skin = "black",
                      dashboardHeader(title = "Epi.Tree Friends"),
                      sidebar,
                      body
)

########In Server#################
server <- function(input, output) {
  stampval <- reactiveVal()
  whichtab <- reactiveVal()
  max_Nhosttree <- reactiveVal()
  dna_seqtree <- reactiveVal()
#####First Tab panel ####### 
    #####List of reactive values#######
  value_1 <- reactiveVal()
  value_2 <- reactiveVal()
  max_Nhost <- reactiveVal()
  status <- reactiveVal()
  status_dna <- status_dna2 <-  reactiveVal()
  status_dna3 <- reactiveVal()
  dna_seq <- reactiveVal()
  traits_table_samp <- reactiveVal()
  table_continent_samp <- reactiveVal()
  list_trees_samp <- reactiveVal()
  BSSVS_result <- reactiveVal()
  
    #####Initialize values and entry data check ###########
  
  output$status <- renderText({status()})
  output$status_dna <- renderText({status_dna()})
  output$status_dna2 <- renderText({status_dna2()})
  output$status_dna3 <- renderText({status_dna3()})
  observeEvent({input$num & input$num2 },
               {output$name_tree <- renderPlot(0)
               output$make_sequence <- NULL
               output$BSSVS_result_table <- NULL
               output$download <- NULL
               output$download_dna <- NULL
               output$download_dna_sampled <- NULL
               output$download_traits_table_samp <- NULL
               output$download_bssvs_result <- NULL
               output$traits_table_samp <- NULL
               output$sampled_txt_output <- NULL
               output$table_continent_samp <- NULL
               output$download_trees_sampled <- NULL
               output$downloadBSSVS_result <- NULL
               output$download_BSSVS <- NULL
               stampval(as.numeric(format(Sys.Date(), "%m%j%H")))
               print(as.numeric(format(Sys.Date(), "%m%j%H")))
               print(typeof(as.numeric(format(Sys.Date(), "%m%j%H"))))
               #print( as.integer(format(Sys.Date(), "%Y%m%d%j%H")))
               })
  
   
  observeEvent({input$num & input$num2 },{              
    output$network_graph <- NULL
    value_1(input$num)
    value_2(input$num2)
    if(is.na(input$num))
    {
      value_1 <- 0
    }
    if(is.na(input$num2))
    {
      value_2 <- 0
    }
    
    max_Nhost(max(value_1(), value_2()))
    #print("before_desabled")
    shinyjs::disable("makesequence")
    shinyjs::disable("treeforestsequences")
    shinyjs::disable("performbssvs")
    status("Please create the transmission network tree in order to create the corresponding sequences.")
    status_dna("Please create or upload sequences in order to start the subsampling.")
    status_dna2("Please create or upload sequences in order to start the subsampling.")
    status_dna3("Please create and subsample the sequences to perform a BSSVS analysis.")
    #print("before req")
    
    # print((value_1, value_2))
    value <- value_1()+ value_2()
    #print(paste("the value is",value,sep=" "))
    #              req(input$num & input$num2)
    #print(value < 2)
    if ((value) < 2)
    { 
      #print("dans le warning")
      showModal(modalDialog(
        title = "Warning!",
        "You need at least two host to calculate a network."))
    }
    output$check_network <- output$check_network2 <- renderText({
      # print("before warning")
      #print(length(epi_network()[,1]))
      #print(length(epi_network()[complete.cases(epi_network()),1]))
      validate(
        need(!is.null(epi_network()),'There is no epidemia, change a parameter1'))
      validate(
        need(length(epi_network()[complete.cases(epi_network()),1]) > 5 , 'There is no epidemia, change a parameter2'))
      paste("There is an epidemia, we can produce a transmission network")
    })
  })
  
  observeEvent(input$goButton ,{
    # print("before_enable_sequences")
    shinyjs::enable("makesequence")
    status("You can create the sequences.")} )
  
  observeEvent(input$makesequence ,{
    shinyjs::enable("treeforestsequences")
    status_dna("There are sequences you can subsample.")
    status_dna2("There are sequences you can subsample.")
    status_dna3("Please subsample the sequences to be able to perform a BSSVS analysis.")
  })
  observeEvent(input$treeforestsequences,{
    shinyjs::enable("performbssvs")
  })
  
  
    #####original tables and values to remove ####### 
  
  pop_struct <- reactive({
    dual_population_structure(max_Nhost())
  })
  
  dyadcov <- reactive({
    validate(
      need(max_Nhost() > 1, "Please select a number of hosts greater than one."))
    dual_population_dyadcov(max_Nhost())
  })
  
  rmvalues <- reactive({
    # print("dans rmvalues")
    element_to_rmv(max_Nhost(),value_1(),value_2())
  })
  
    #####Tables reduced #####
  
  dyadcov_red <- reactive({
    #print("dams reactive dyadcov reduce")
    epidemia_dyadcov_reduced(dyadcov(),rmvalues())
  })
  
  pop_struct_red <- reactive({
    #print("avant pop reduced")
    population_structure_reduced(pop_struct(),rmvalues())
  })
    #####Reactive on network ########################
  
  
  epi_network <- reactive({
    # print("dans epi net")
    #print(paste("values to rm:",rmvalues(),sep=""))
    network_creation(dyadcov(),max_Nhost(),value_1(),value_2(),rmvalues())})
  
  output$tpop_struct_red <- DT::renderDataTable({
    pop_struct_red()})
  
  output$epi_network <- DT::renderDataTable({
    epi_network()})
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Epidemia",input$num, "hosts.csv", sep = "_")
    },
    content = function(file) {
      write.csv(epi_network(), file, row.names = FALSE)
    })
  
  output$check_network <- output$check_network2 <- renderText({
    # print("before warning")
    #print(length(epi_network()[,1]))
    #print(length(epi_network()[complete.cases(epi_network()),1]))
    validate(
      need(!is.null(epi_network()),'There is no epidemia, change a parameter1'))
    validate(
      need(length(epi_network()[complete.cases(epi_network()),1]) > 5 , 'There is no epidemia, change a parameter2'))
    paste("There is an epidemia, we can produce a transmission network")
  })
  
  
    #####Plot the transmission network########## 
  observeEvent(input$checkbox,{
    if (length(epi_network()[complete.cases(epi_network()),1]) <=5)
    { 
      showModal(modalDialog(
        title = "Warning!",
        "The epidemia is too small to produce a tranmission network."))
    }
    output$network_graph <- renderVisNetwork({
      validate(need(length(epi_network()[complete.cases(epi_network()),1]) > 5,
                    "Please change a parameter to simulate an epidemia"))  
      isolate(network_visual(epi_network(),pop_struct_red(),value_1(),value_2()))})
  })
    #####Press buttons ###################
  observeEvent(input$goButton,{#print("ici8")
    # TIf the network is too small, sends a warning pop-up
    if (length(epi_network()[complete.cases(epi_network()),1]) <=5||max_Nhost() <=2)
    { 
      showModal(modalDialog(
        title = "Warning!",
        "The epidemia is too small to produce a tranmission tree."))
    }
    else{
      output$check_network <- NULL
      plot1 <-make_tree(epi_network(),pop_struct_red(),max_Nhost(),1,stampval()) 
      output$name_tree <- renderPlot(plot1)
      output$downloadPlot <- downloadHandler(
        filename = 
          function() { 
            paste("Epidemia_network",max_Nhost(),'_hosts.png', sep='') 
          },
        content = function(file) {
          ggsave(file, plot = plot1, device = "png",dpi = 1080, width = 10, units = "in")
        },
        contentType = 'image/png'
      )
      
      output$download <- renderUI({downloadButton('downloadPlot', 'Download the tranmission tree.')})
    }
  })
  
  observeEvent(input$makesequence,{
    whichtab(1) 
    #output$make_sequence <-renderText(
    dna_sequences <-(withProgress(message = 'Creating the sequences', 
                                  value = 0, { make_sequences(epi_network(),pop_struct_red(),max_Nhost(),stampval())}))
    
    #print(dna_sequences)
    dna_seq(dna_sequences)
    output$downloadsequences <- downloadHandler(
      filename = function() {
        paste("Epidemia",max_Nhost(), "hosts_sequences.fas", sep = "_")
      },
      content = function(file) {
        write.dna(dna_seq(),file,format="fasta", nbcol=-1, colsep="")
      })
    output$download_dna <- renderUI({downloadButton('downloadsequences', 'Download the DNA sequences.')})
  })
 
  
 
#####Second tab panel#####
    ####List of reactive values####
  value_14 <- reactiveVal()
  value_24 <- reactiveVal()
  value_34 <- reactiveVal()
  value_44 <- reactiveVal()
  max_Nhost4 <- reactiveVal()
  dna_seq4 <- reactiveVal()
  status4 <- reactiveVal()
  
    ####Initialize values and entry data check####
  output$status4 <- renderText({status4()}) 
  observeEvent({input$num14 & input$num24 &input$num34 & input$num44 },{   
  #output$name_tree4 <- renderPlot(0)
  output$make_sequence4 <- NULL
  output$BSSVS_result_table4 <- NULL
  output$download4 <- NULL
  output$download_dna4 <- NULL
  output$download_dna_sampled4 <- NULL
  output$download_traits_table_samp4 <- NULL})
  
  observeEvent({input$num14 & input$num24 &input$num34 & input$num44 },{
    shinyjs::disable("treeforestsequences")
    shinyjs::disable("performbssvs")
    output$network_graph4 <- NULL
    output$name_tree4 <- NULL
    value_14(input$num14)
    #  print(value_1())
    value_24(input$num24)
    #   print(value_2())
    value_34(input$num34)
    #  print(value_1())
    value_44(input$num44)
    #   print(value_2())
    if(is.na(input$num14))
    {
      #   print("num is na")
      value_14 <- 0
      #   print(value_1)
    }
    if(is.na(input$num24))
    {
      # print("num is na")
      value_24 <- 0
      # print(value_2)
    }
    if(is.na(input$num34))
    {
      #   print("num is na")
      value_34 <- 0
      #   print(value_1)
    }
    if(is.na(input$num44))
    {
      # print("num is na")
      value_44 <- 0
      # print(value_2)
    }
    max_Nhost4(max(value_14(),value_24(),value_34(),value_44()))
    status4("Please create the transmission network tree in order to create the corresponding sequences.")
    value <- value_14()+ value_24()+value_34()+ value_44()
    #print(paste("the value is",value,sep=" "))
    #              req(input$num & input$num2)
    #print(value < 2)
    if ((value) < 2)
    { 
      #print("dans le warning")
      showModal(modalDialog(
        title = "Warning!",
        "You need at least two host to calculate a network."))
    }
    shinyjs::disable("makesequence4")
    output$check_network4 <- output$check_network24 <- renderText({
      # print("before warning")
      #print(length(epi_network()[,1]))
      #print(length(epi_network()[complete.cases(epi_network()),1]))
      validate(
        need(!is.null(epi_network_circular()),'There is no epidemia, change a parameter1'))
      validate(
        need(length(epi_network_circular()[complete.cases(epi_network_circular()),1]) > 5 , 'There is no epidemia, change a parameter2'))
      paste("There is an epidemia, we can produce a transmission network")
    })})

  
  observeEvent(input$goButton4 ,{
    # print("before_enable_sequences")
    shinyjs::enable("makesequence4")
    status4("You can create the sequences.")} )
  
  observeEvent(input$makesequence4 ,{
    shinyjs::enable("treeforestsequences")
    status_dna("There are sequences you can subsample.")
    status_dna2("There are sequences you can subsample.")
    status_dna3("Please subsample the sequences to be able to perform a BSSVS analysis.")
  })  
    #####original tables and values to remove ####### 
  pop_struct_circular4 <- reactive({
    circular_population_structure(max_Nhost4())
  })
  
  dyadcov_circular4 <- reactive({
    validate(
      need(max_Nhost4() > 1, "Please select a number of hosts greater than one."))
    circular_dyadcov(max_Nhost4())
  })
  
  rmvalues_circular4 <- reactive({
    # print("dans rmvalues")
    element_to_rmv_circular(max_Nhost4(),value_14(),value_24(),value_34(),value_44())
  }) 
  
  pop_struct_circular4_reduced <-reactive({
    #print("in pop struct red")
    #print(rmvalues_circular4())
    population_structure_reduced(pop_struct_circular4(),rmvalues_circular4())
  })
    #####Reactive on network ########################
  epi_network_circular <- reactive({
    network_creation_circular(dyadcov_circular4(),max_Nhost4(),value_14(),value_24(),value_34(),value_44(),rmvalues_circular4())})
  
  output$pop_struct_red_circular <- DT::renderDataTable({
    pop_struct_circular4_reduced()})
  
  output$epi_network_circular <- DT::renderDataTable({
    epi_network_circular()})  
  
  output$downloadData4 <- downloadHandler(
    filename = function() {
      paste("Epidemia",input$num, "hosts.csv", sep = "_")
    },
    content = function(file) {
      write.csv(epi_network_circular(), file, row.names = FALSE)
    }) 
  
  output$check_network4 <- output$check_network24 <- renderText({
    # print("before warning")
    #print(length(epi_network()[,1]))
    #print(length(epi_network()[complete.cases(epi_network()),1]))
    validate(
      need(!is.null(epi_network_circular()),'There is no epidemia, change a parameter1'))
    validate(
      need(length(epi_network_circular()[complete.cases(epi_network_circular()),1]) > 5 , 'There is no epidemia, change a parameter2'))
    paste("There is an epidemia, we can produce a transmission network")
  })
    #####Plot the transmission network########## 
  observeEvent(input$checkbox4,{
    print("in check box")
    if (length(epi_network_circular()[complete.cases(epi_network_circular()),1]) <=5)
    { 
      showModal(modalDialog(
        title = "Warning!",
        "The epidemia is too small to produce a tranmission network."))
    }
    output$check_network24 <- NULL
    output$network_graph4 <- renderVisNetwork({
      validate(need(length(epi_network_circular()[complete.cases(epi_network()),1]) > 5,
                    "Please change a parameter to simulate an epidemia"))  
      isolate(network_visual_circular(epi_network_circular(),pop_struct_circular4_reduced(),value_14(),value_24(),value_34(),value_44()))})
  })
    #####Press buttons ###################
  observeEvent(input$goButton4,{#print("ici8")
    # TIf the network is too small, sends a warning pop-up
    if (length(epi_network_circular()[complete.cases(epi_network_circular()),1]) <=5||max_Nhost4() <=2)
    { 
      showModal(modalDialog(
        title = "Warning!",
        "The epidemia is too small to produce a tranmission tree."))
    }
    else{
      output$check_network4 <- NULL
      plot1 <-make_tree(epi_network_circular(),pop_struct_circular4_reduced(),max_Nhost4(),1,stampval()) 
      output$name_tree4 <- renderPlot(plot1)
      output$downloadPlot4 <- downloadHandler(
        filename = 
          function() { 
            paste("Epidemia_network",max_Nhost4(),'_hosts.png', sep='') 
          },
        content = function(file) {
          ggsave(file, plot = plot1, device = "png",dpi = 1080, width = 10, units = "in")
        },
        contentType = 'image/png'
      )
      
      output$download <- renderUI({downloadButton('downloadPlot4', 'Download the tranmission tree.')})
    }
  })

  observeEvent(input$makesequence4,{
    whichtab(2)
    print(paste("tab:",whichtab(),sep=""))
    #output$make_sequence <-renderText(
    dna_sequences4 <-(withProgress(message = 'Creating the sequences', 
                                  value = 0, { make_sequences(epi_network_circular(),pop_struct_circular4_reduced(),max_Nhost4(),stampval())}))
    status4("The sequences are created")
    #print(dna_sequences)
    dna_seq4(dna_sequences4)
    output$downloadsequences4 <- downloadHandler(
      filename = function() {
        paste("Epidemia",max_Nhost4(), "hosts_sequences.fas", sep = "_")
      },
      content = function(file) {
        write.dna(dna_seq4(),file,format="fasta", nbcol=-1, colsep="")
      })
    output$download_dna4 <- renderUI({downloadButton('downloadsequences4', 'Download the DNA sequences.')})
  })  
 
#####Third tab panel #####
  shinyjs::onclick("toggleAdvanced",
                   shinyjs::toggle(id = "advanced", anim = TRUE))
    ####Sample sequences####
  observeEvent(input$treeforestsequences,{
    #output$make_sequence <-renderText(
   # print("avant_sampling")
    if(whichtab()==1)
    {
      print("in which tab")
      print(max_Nhost())
      max_Nhosttree(max_Nhost()) 
      dna_seqtree(dna_seq())
      #max_Nhosttree() <- max_Nhost()
      #dna_seqtree() <- dna_seq()
    }else{
      print(max_Nhost4())
      max_Nhosttree(max_Nhost4()) 
      dna_seqtree(dna_seq4())
      #max_Nhosttree() <- max_Nhost4()
      #dna_seqtree() <- dna_seq4()
    }
    list[dna_sequences_sampled, trait_table_samp,list_tree_samp] <-(withProgress(message = 'Subsample the sequences', 
                                  value = 0, { make_tree_of_forest_sequences(max_Nhosttree(),dna_seqtree(),input$jn,input$jr,input$numPerCat,input$repsPerCat,stampval())}))
    #print("apres sampling avant bssvs")
    #print(trait_table_samp)
    
    status_dna("The sequences are sampled")
    status_dna2(NULL)
   # print(dna_sequences_sampled)
    #print(trait_table_samp)
    dna_seq(dna_sequences_sampled)
    traits_table_samp(trait_table_samp)
    list_trees_samp(list_tree_samp)
   # print(list_trees_samp())
    num_seq <- length(traits_table_samp()[,1])
    #print(num_seq)
    type_traits <- unique(traits_table_samp()[,2])
    ntraits <- length(type_traits)
    
    output$sampled_txt_output <- renderText({
      paste("There is ",num_seq," sequences in ",ntraits," groups")
    })
    
    table_continent_sampled <-aggregate(data.frame(count = traits_table_samp()[,2]), list(value = traits_table_samp()[,2]), length)
    colnames(table_continent_sampled) <- c("Trait","Number sequences")
    #print(table_continent)
    table_continent_samp(table_continent_sampled)
    
    output$traits_table_samp <- DT::renderDataTable({traits_table_samp()})
    output$table_continent_samp <- DT::renderDataTable({table_continent_samp()})
    
    output$downloadtraits_table_samp <- downloadHandler(
      filename = function() {
        paste("Epidemia",max_Nhost(), "hosts_sequences_sampled_trait_table.csv", sep = "_")
      },
      content = function(file) {
        write.csv(traits_table_samp(), file, row.names = FALSE)
      })
    output$download_traits_table_samp <- renderUI({downloadButton('downloadtraits_table_samp', 'Download the DNA sequences trait table.')})
    
    output$downloadsequences_sampled <- downloadHandler(
      filename = function() {
        paste("Epidemia",max_Nhost(), "hosts_sequences_sampled.fas", sep = "_")
      },
      content = function(file) {
        write.dna(dna_seq(),file,format="fasta", nbcol=-1, colsep="")
      })
    output$downloadstrees_sampled <- downloadHandler(
      filename = function() {
        paste("Epidemia",max_Nhost(), "hosts_sequences_sampled_tree_set.trees", sep = "_")
      },
      content = function(file) {
        writeLines(list_trees_samp(),file)
      })
    output$download_dna_sampled <- renderUI({downloadButton('downloadsequences_sampled', 'Download the DNA sequences.')})
    output$download_trees_sampled <- renderUI({downloadButton('downloadstrees_sampled', 'Download the corresponding tree set.')})
  })
  
    ####Perform BBSVS####
  observeEvent(input$performbssvs,{
    withProgress(message = 'Performing the BSSVS analysis', 
                  value = 0,make_bssvs_anal(dna_seq(),traits_table_samp(),max_Nhost(),stampval()))
    bf_table <- Calculate_bf(max_Nhost())
    #print(bf_table)
    BSSVS_result(bf_table)
    output$BSSVS_result_table <- DT::renderDataTable({BSSVS_result()})

    output$downloadBSSVS_result <- downloadHandler(
      filename = function() {
        paste("Epidemia",max_Nhost(), "hosts_sequences_sampled_BSSVS_result.csv", sep = "_")
      },
      content = function(file) {
        write.csv(BSSVS_result(), file, row.names = FALSE)
      })
    output$download_BSSVS <- renderUI({downloadButton('downloadBSSVS_result', 'Download the BSSVS table.')})
#    output$BSSVS_graph <- renderVisNetwork({
#      isolate(plot_bssvs(BSSVS_result()))})
  })
  
}

shinyApp(ui, server)