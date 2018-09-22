
#library(epinet)
#library(ggtree)
#library(ape)
#library(igraph)
#library(visNetwork) 
setwd("D:/PhD/GitHub/Rshiny_webapp/Rshiny_webapp/")
source("modif_transmission_tree_2.R")
Rpath <- getwd()
source(paste(Rpath,"modif_transmission_tree_2.R",sep="/"))
source(paste(Rpath,"traits_of_the_forest_tree_set.R",sep = "/"))



getEl	<- function( line, sep=",", ind=-1, final=FALSE, reconstruct=FALSE, ex=-1, fromEnd=FALSE ) {
  els	<- strsplit(line, sep)[[1]]
  
  if (ind[1] != -1) {
    if (fromEnd) {
      ind <- length(els)-(ind-1)
    }
  }
  
  if (final) {
    return( els[length(els)] )
  } else {
    
    if (reconstruct) {
      if (ex[1] > 0) {
        if (fromEnd) {
          ex <- length(els)-(ex-1)
        }
        ind <- setdiff((1:length(els)),ex)
      }
      
      newLine <- els[ind[1]]
      if (length(ind) > 1) {
        for (i in 2:length(ind)) {
          newLine <- paste(newLine, els[ind[i]], sep=sep)
        }
      }
      return ( newLine )
    } else {
      if ( ind[1] == -1 ) {
        return( els )
      } else {
        return( els[ind] )
      }
    }
  }
}


#### Dual population ####
dual_population_structure <- function(number_elements_node)
{
  N = number_elements_node*2
  population_structure <- data.frame(cbind(seq(1, N, by=1),rep(c("A","B"), each=N/2)))
  colnames(population_structure) <- c("node","population_structure")
  return(population_structure)
}

dual_population_dyadcov <- function(number_elements_node)
{
  #print("make dyadcov")
  N = number_elements_node*2
  mycov1 = data.frame(id = 1:N,inpopulation_structure1 = rep(1:0, each=N/2),inpopulation_structure = rev(rep(1:0, each=N/2)))
  
  # make matrix 
  dyadCov1 = BuildX(mycov1, binaryFunc = "manhattan",binaryCol = list(c(2,3)))
  return(dyadCov1)
}

element_to_rmv <- function(number_elements_node,Na_ori,Nb_ori)
{
  #print("in elem to remove")
  N_ori <- Na_ori+Nb_ori
  Na <-Nb <- number_elements_node
  #Nd = number_elements_node/2
  
  N_final= (Na+Nb)
  Na_start <- 1
  Na_end <- Na
  Nb_start <-Na_end+1
  Nb_end <- Nb_start+Nb-1
  
  #Na_ori <- 15
  ####start remove A elements ##
  Na_difference <- number_elements_node-Na_ori
  Na_values <- Na_start:Na_end
  Na_values_to_remove <- sample(Na_values,Na_difference)
  ####start remove B elements ###
  Nb_difference <- number_elements_node-Nb_ori
  Nb_values <- Nb_start:Nb_end
  Nb_values_to_remove <- sample(Nb_values,Nb_difference)
  ####start remove C elements ###
  
  ###
  values_to_remove <- c(Na_values_to_remove,Nb_values_to_remove)
  #print(paste("in fct value to rm:",values_to_remove,sep=""))
  #print(typeof(values_to_remove))
  return(values_to_remove)
}

epidemia_dyadcov_reduced <- function(dual_population_dyadcov,values_to_remove)
{
  #print("dans dyadcov reduce")
  if(length(values_to_remove) == 0)
  {
    #print("No values to remove")
    #name_dyadcov <- paste(file_name,"reduced_dyadcov.txt",sep="")
    #write.table(circular_epidemia_dyadcov,name_dyadcov,quote = FALSE,row.names = FALSE,sep="\t")
    return(dual_population_dyadcov)
  }
  else{
    #print("there are values to remove")
    #print(values_to_remove)
    #print(typeof(values_to_remove))
    #print(paste("lenght of dyadic cov:",length(dual_population_dyadcov),sep=""))
    in_col_1 <- which(dual_population_dyadcov[,1] %in% values_to_remove )
    in_col_2 <- which(dual_population_dyadcov[,2] %in% values_to_remove )
    to_remove <- as.numeric(c(in_col_1,in_col_2))
    #print(paste("lenght of to remove:",length(to_remove),sep=""))
    
    dual_population_dyadcov_reduced<- dual_population_dyadcov[-to_remove,]
    num_host <- length(unique(dual_population_dyadcov_reduced[,2])) +1
    host_num <- seq(1:num_host)
    
    new_host_num <- list()
    new_host_conn <- list()
    for (i in 1:(num_host-1))
    {
      new_host_num <- c(new_host_num,rep(host_num[i],each =(num_host - i )))
      new_host_conn <- c(new_host_conn,(i+1):num_host)
    }
    
    new_host_num <- data.table(new_host_num)
    new_host_conn <- data.table(new_host_conn)
    
    new_dyncov <- cbind(new_host_num,new_host_conn,dual_population_dyadcov_reduced[,3:4])
    #name_dyadcov <- paste(file_name,"reduced_dyadcov.txt",sep="")
    df_new_dyncov <- data.frame(matrix(unlist(new_dyncov)))
    #write.table(df_new_dyncov,name_dyadcov,quote = FALSE,row.names = FALSE,sep="\t")
    return(df_new_dyncov)
  }
} 

population_structure_reduced <-function(population_structure,values_to_remove)
{
  #print("in pop reduced")
  if(length(values_to_remove) == 0)
  {
    #name_pop_struct <- paste(file_name,"reduced_pop_struct.txt",sep="")
    #write.table(population_structure,name_pop_struct,quote = FALSE,row.names = FALSE,sep="\t")
    return(population_structure)
  }
  else{
    result_population_structure <- population_structure[-values_to_remove,]
    len_result_population_structure <- length(result_population_structure[,1])
    host_num <- seq(1:len_result_population_structure)
    result_population_structure[,1] <- host_num
    #name_pop_struct <- paste(file_name,"reduced_pop_struct.txt",sep="")
    #write.table(result_population_structure,name_pop_struct,quote = FALSE,row.names = FALSE,sep="\t")
    return(result_population_structure)
  }
}

network_creation <- function(circular_epidemia_dyadcov,number_elements_node,Na_ori,Nb_ori,values_to_remove)
{
  #print(values_to_remove)
  N_ori <- Na_ori+Nb_ori
  Na <-Nb  <- number_elements_node
  N_final= (Na+Nb)
  
  old_value <- c(1:N_final)
  old_value <- old_value[-values_to_remove]
  new_value <- c(1:length(old_value))
  value_traduction <- cbind(old_value,new_value)  
  
  population_structure_final <- circular_epidemia_dyadcov
  #explain the relation between distance and the resistance
  eta <- c(0,-0.3)
  net <- SimulateDyadicLinearERGM(N = N_final, dyadiccovmat = population_structure_final, eta = eta)
  if(length(values_to_remove) == 0)
  {
    epi<- SEIR.simulator(M = net, N = N_ori, beta = 1,ki = 1, thetai = 2,thetae = 2, ke = 1,latencydist="gamma")
    return(epi)
  }
  else{
    
    #remove unnecessary values
    in_col_1_net <- which(net[,1] %in% values_to_remove )
    in_col_2_net <- which(net[,2] %in% values_to_remove )
    to_remove_net <- as.numeric(c(in_col_1_net,in_col_2_net))
    net_reduced_net <- net[-to_remove_net,]
    
    net_reduced_traducted_1 <- as.character(1:length(value_traduction[,2] ))[ match(net_reduced_net[,1], value_traduction[,1] ) ] 
    net_reduced_traducted_2 <- as.character(1:length(value_traduction[,2] ))[ match(net_reduced_net[,2], value_traduction[,1] ) ] 
    net_reduced_traducted <-cbind(net_reduced_traducted_1,net_reduced_traducted_2)
    #print(net_reduced_traducted)
    #simulate the epidemia
    epi_reduced <- SEIR.simulator(M = net_reduced_traducted, N = N_ori, beta = 1,ki = 1, thetai = 2,thetae = 2, ke = 1,latencydist="gamma")
    epi_table <- as.data.frame(epi_reduced)
    
    return(epi_table)
  }
}

network_visual <- function(network_dual,population_structure,Na_ori,Nb_ori)
{
  modif_epi_table <- network_dual
  min_epi_table <- min(modif_epi_table[,3])
  modif_epi_table[,3] <- modif_epi_table[,3]- min_epi_table
  modif_epi_table[,4] <- modif_epi_table[,4]- min_epi_table
  modif_epi_table[,5] <- modif_epi_table[,5]- min_epi_table
  max_epi_table <- max(modif_epi_table[,3])
  
  epi_network <- cbind(modif_epi_table[,2],modif_epi_table[,1])
  
  epi_network<-   epi_network[complete.cases(epi_network), ]
  colnames(epi_network) <- c("from","to")
  
  net1 <- graph_from_data_frame(epi_network)
  
  #colr_a <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =Na_ori )
  #colr_b <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =Nb_ori )
  
  
  nodes <- data.frame(population_structure)
  #print(nodes)
  colnames(nodes) <-c("id","population_structure")
  links <- data.frame(epi_network)
  
  vis.nodes <- nodes
  vis.links <- links
  
  vis.nodes$shape  <- "dot"  
  vis.nodes$group <- vis.nodes$population_structure # add groups on nodes 
  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  vis.nodes$title  <- vis.nodes$population_structure # Text on click
  vis.nodes$label  <- vis.nodes$id # Node label
  vis.nodes$borderWidth <- 2 # Node border width
  
  vis.nodes$color.background <- c("#8dd3c7", "#ffffb3", "#bebada","#fb8072")[nodes$population_structure]
  vis.nodes$color.border <- "black"
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"
  
  vis.links$color <- "gray"    # line color  
  vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- TRUE    # should the edges be curved?
  vis.links$shadow <- TRUE
  network_visual_output <-   visNetwork(vis.nodes, vis.links,main="Tranmission between the infected hosts", submain="and the different populations",height = "1080px")%>%
    visGroups(groupname = "A", color = "#8dd3c7") %>%
    visGroups(groupname = "B", color = "#ffffb3")%>%
    visGroups(groupname = "C", color = "#bebada")%>%
    visGroups(groupname = "D", color = "#fb8072")%>%
    visPhysics(enabled = FALSE) %>% visEdges(smooth = FALSE) %>%
    visOptions(selectedBy = "group",highlightNearest = TRUE) %>%
    visLegend( position = "right", main = "Population")%>% 
    visExport(type = "png", name = "export-network", 
              float = "left", label = "Download network picture", style= "") 
  
  return(network_visual_output)
}

#########################
#### four population circular####
##create a table which give the population structure for a circulare pop 
circular_population_structure <- function(number_elements_node)
{
  N = number_elements_node*4
  population_structure <- cbind(c(1:N),rep(c("A","B","C","D"),each =number_elements_node ))
  colnames(population_structure) <- c("node","population_structure")
  return(population_structure)
}

#create the dyad cov matrix for a circular population
circular_dyadcov	<- function( number_elements_node)
{
  N = number_elements_node*4
  ####create the population structure####
  #circulare population structure between 4 populations#
  num = 1:N
  population_structure_temp <- data.table()
  population_structure_temp2 <- data.table()
  #create the links between the different "hosts"
  for(i in 1:N)
  {
    population_structure_temp <- c(population_structure_temp,rep(num[i],(N-i)))
    if (i+1 <= N)
    {
      population_structure_temp2 <- c(population_structure_temp2,(i+1):N)
    }
  }
  
  population_structure_temp3 <- data.frame(cbind(population_structure_temp,population_structure_temp2))
  
  population_structure_temp4<- data.frame(id = numeric(length(population_structure_temp2)))
  interept = rep(1,length(population_structure_temp3[,1]))
  #create the distance between te different "hosts"
  for (i in 1:length(population_structure_temp2))
  {
    leftside <-as.numeric(as.character(unlist(population_structure_temp3[i,1])))
    rightside <- as.numeric(as.character(unlist(population_structure_temp3[i,2])))
    multiple_N2_leftside <- (leftside %% number_elements_node ==0)  
    divide_N2_leftside <- leftside %/% number_elements_node
    multiple_N2_rightside <- (rightside %% number_elements_node ==0)
    divide_N2_rightside  <- rightside %/% number_elements_node
    est_multiple <- multiple_N2_rightside == multiple_N2_leftside
    value <- rightside- leftside
    #  print(value)
    if(value <= (number_elements_node-1))
    {
      if (multiple_N2_leftside ==  FALSE & multiple_N2_rightside ==  FALSE & divide_N2_leftside == divide_N2_rightside)
      {
        population_structure_temp4[i,1]<- 0
      }else if(multiple_N2_leftside == FALSE & multiple_N2_rightside ==  TRUE & divide_N2_leftside != divide_N2_rightside )# & leftside %% number_elements_node == 0 )
      {
        population_structure_temp4[i,1]<- 0
      }else if(multiple_N2_leftside == FALSE & multiple_N2_rightside ==  FALSE & divide_N2_leftside != divide_N2_rightside )# & leftside %% number_elements_node == 0 )
      {
        population_structure_temp4[i,1]<- 1
      }
      if (multiple_N2_leftside ==  TRUE)
      {
        population_structure_temp4[i,1]<- 1
      }
    }else {
      if (multiple_N2_rightside)
      {
        divide_N2_rightside  <- (rightside -1) %/% number_elements_node
      }
      if (multiple_N2_leftside)
      {
        divide_N2_leftside  <- (leftside -1) %/% number_elements_node
      }
      if (divide_N2_rightside == divide_N2_leftside +2)
      {
        #    if(multiple_N2_rightside)
        #    {
        #      population_structure4[i,1]<- 1
        #    }else{
        population_structure_temp4[i,1]<- 10
        #    }
      }else{
        #    if(multiple_N2_rightside)
        #   {
        #      population_structure4[i,1]<- 2
        #    }else{
        population_structure_temp4[i,1]<- 1
      }
    }
  }
  population_structure_final <- cbind(population_structure_temp3,interept,population_structure_temp4) 
  # this is the dyad cov matrix where we will simulate the epidemia
  population_structure_final <- data.matrix(population_structure_final)
  return(population_structure_final)
  
}

element_to_rmv_circular <- function(number_elements_node,Na_ori,Nb_ori,Nc_ori,Nd_ori)
{
  N_ori <- Na_ori+Nb_ori+Nc_ori+Nd_ori
  Na <-Nb <- Nc <- Nd <- number_elements_node
  #Nd = number_elements_node/2
  
  N_final= (Na+Nb+Nc+Nd)
  Na_start <- 1
  Na_end <- Na
  Nb_start <-Na_end+1
  Nb_end <- Nb_start+Nb-1
  Nc_start <-Nb_end+1
  Nc_end <- Nc_start+Nc-1
  Nd_start <-Nc_end+1
  Nd_end <- Nd_start+Nd-1
  
  #Na_ori <- 15
  ####start remove A elements ##
  Na_difference <- number_elements_node-Na_ori
  Na_values <- Na_start:Na_end
  Na_values_to_remove <- sample(Na_values,Na_difference)
  ####start remove B elements ###
  Nb_difference <- number_elements_node-Nb_ori
  Nb_values <- Nb_start:Nb_end
  Nb_values_to_remove <- sample(Nb_values,Nb_difference)
  ####start remove C elements ###
  Nc_difference <- number_elements_node-Nc_ori
  Nc_values <- Nc_start:Nc_end
  Nc_values_to_remove <- sample(Nc_values,Nc_difference)
  ####start remove D elements ###
  Nd_difference <- number_elements_node-Nd_ori
  Nd_values <- Nd_start:Nd_end
  Nd_values_to_remove <- sample(Nd_values,Nd_difference)
  ###
  values_to_remove <- c(Na_values_to_remove,Nb_values_to_remove,Nc_values_to_remove,Nd_values_to_remove)
  #print(paste("in fct value to rm:",values_to_remove,sep=""))
  #print(typeof(values_to_remove))
  return(values_to_remove)
}

network_creation_circular <- function(circular_epidemia_dyadcov,number_elements_node,Na_ori,Nb_ori,Nc_ori,Nd_ori,values_to_remove)
{
  N_ori <- Na_ori+Nb_ori+Nc_ori+Nd_ori
  Na <-Nb <- Nc <- Nd <- number_elements_node
  N_final= (Na+Nb+Nc+Nd)
  
  old_value <- c(1:N_final)
  old_value <- old_value[-values_to_remove]
  new_value <- c(1:length(old_value))
  value_traduction <- cbind(old_value,new_value)  
  
  population_structure_final <- circular_epidemia_dyadcov
  #explain the relation between distance and the resistance
  eta <- c(0,-0.3)
  net <- SimulateDyadicLinearERGM(N = N_final, dyadiccovmat = population_structure_final, eta = eta)
  if(length(values_to_remove) == 0)
  {
    epi<- SEIR.simulator(M = net, N = N_ori, beta = 1,ki = 1, thetai = 2,thetae = 2, ke = 1,latencydist="gamma")
#    name_epi <- paste(file_name,"reduced_network.txt",sep="")
#    write.table(epi,name_epi,quote = FALSE,row.names = FALSE,sep="\t")
    return(epi)
  }
  else{
    
    #remove unnecessary values
    in_col_1_net <- which(net[,1] %in% values_to_remove )
    in_col_2_net <- which(net[,2] %in% values_to_remove )
    to_remove_net <- as.numeric(c(in_col_1_net,in_col_2_net))
    net_reduced_net <- net[-to_remove_net,]
    
    net_reduced_traducted_1 <- as.character(1:length(value_traduction[,2] ))[ match(net_reduced_net[,1], value_traduction[,1] ) ] 
    net_reduced_traducted_2 <- as.character(1:length(value_traduction[,2] ))[ match(net_reduced_net[,2], value_traduction[,1] ) ] 
    net_reduced_traducted <-cbind(net_reduced_traducted_1,net_reduced_traducted_2)
    
    #simulate the epidemia
    epi_reduced <- SEIR.simulator(M = net_reduced_traducted, N = N_ori, beta = 1,ki = 1, thetai = 2,thetae = 2, ke = 1,latencydist="gamma")
    epi_table <- as.data.frame(epi_reduced)
#    name_epi <- paste(file_name,"reduced_network.txt",sep="")
#    write.table(epi_table,name_epi,quote = FALSE,row.names = FALSE,sep="\t")
    return(epi_table)
  }
}

population_structure_reduced <-function(population_structure,values_to_remove)
{
  if(length(values_to_remove) == 0)
  {
 #   name_pop_struct <- paste(file_name,"reduced_pop_struct.txt",sep="")
   # write.table(population_structure,name_pop_struct,quote = FALSE,row.names = FALSE,sep="\t")
    return(population_structure)
  }
  else{
    result_population_structure <- population_structure[-values_to_remove,]
    len_result_population_structure <- length(result_population_structure[,1])
    host_num <- seq(1:len_result_population_structure)
    result_population_structure[,1] <- host_num
    #name_pop_struct <- paste(file_name,"reduced_pop_struct.txt",sep="")
    #write.table(result_population_structure,name_pop_struct,quote = FALSE,row.names = FALSE,sep="\t")
    return(result_population_structure)
  }
}

network_visual_circular <- function(epidemia_network,result_population_reduced,Na_ori,Nb_ori,Nc_ori,Nd_ori)
{
 #print("in creation network")
  circular_epidemia_network <- epidemia_network
  
  #starting some network vizualisation 
  test <- circular_epidemia_network[,-is.na(circular_epidemia_network[,3])]
  
  modif_epi_table <- circular_epidemia_network
  min_epi_table <- min(modif_epi_table[,3])
  modif_epi_table[,3] <- modif_epi_table[,3]- min_epi_table
  modif_epi_table[,4] <- modif_epi_table[,4]- min_epi_table
  modif_epi_table[,5] <- modif_epi_table[,5]- min_epi_table
  #max_epi_table <- max(modif_epi_table[,3])
  
  epi_network <- cbind(modif_epi_table[,2],modif_epi_table[,1])
 # print("ici2")
  epi_network<-   epi_network[complete.cases(epi_network), ]
  colnames(epi_network) <- c("from","to")
  
  net1 <- graph_from_data_frame(epi_network)
  
  #colr_a <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =Na_ori )
  #colr_b <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =Nb_ori )
  #colr_c <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =Nc_ori )
  #colr_d <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =Nd_ori )
  #print("ici3")
  nodes <- data.frame(result_population_reduced)
  print("in network visuall circular")
  print(nodes)
  colnames(nodes) <-c("id","population_structure")
  links <- data.frame(epi_network)
  #print(nodes)
  vis.nodes <- nodes
  vis.links <- links
  
  vis.nodes$shape  <- "dot"  
  vis.nodes$group <- vis.nodes$population_structure # add groups on nodes 
  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  vis.nodes$title  <- vis.nodes$population_structure # Text on click
  vis.nodes$label  <- vis.nodes$id # Node label
  vis.nodes$borderWidth <- 2 # Node border width
  #print("ici4")
  #print(c("#8dd3c7", "#ffffb3", "#bebada","#fb8072")[nodes$population_structure])
  #print(unique(nodes$population_structure))
  
  rep1 <- length(which(nodes$population_structure == "A"))
  rep2 <- length(which(nodes$population_structure == "B"))
  rep3 <- length(which(nodes$population_structure == "C"))
  rep4 <- length(which(nodes$population_structure == "D"))

  col1 <- rep("#8dd3c7",rep1)
  col2 <- rep("#ffffb3",rep2)
  col3 <- rep("#bebada",rep3)
  col4 <- rep("#fb8072",rep4)
  print(c(col1,col2,col3,col4))
    
 #print(rep(c("A","B","C","D"),c(rep1,rep2,rep3,rep4)))
  vis.nodes$color.background <- c(col1,col2,col3,col4)
  #vis.nodes$color.background <- c("#8dd3c7", "#ffffb3", "#bebada","#fb8072")[nodes$population_structure]
  vis.nodes$color.border <- "black"
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"
  
  vis.links$color <- "gray"    # line color  
  vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- TRUE    # should the edges be curved?
  vis.links$shadow <- TRUE
  #print(vis.nodes)
  #print("ici5")
  network_visual_output <-   visNetwork(vis.nodes, vis.links,main="Tranmission between the infected hosts", submain="and the different populations",height = "1080px")%>%
    visGroups(groupname = "A", color = "#8dd3c7") %>%
    visGroups(groupname = "B", color = "#ffffb3")%>%
    visGroups(groupname = "C", color = "#bebada")%>%
    visGroups(groupname = "D", color = "#fb8072")%>%
    visPhysics(enabled = FALSE) %>% visEdges(smooth = FALSE) %>%
    visOptions(selectedBy = "group",highlightNearest = TRUE) %>%
    visLegend( position = "right", main = "Population")%>% 
    visExport(type = "png", name = "export-network", 
              float = "left", label = "Download network picture", style= "") 
  
  return(network_visual_output)
}

##########################
####Function for all####


make_tree <- function(network_dual,population_structure,number_elements_node,value,datenum)
{
  non_NA_network <-   data.frame(network_dual[complete.cases(network_dual[,3:length(network_dual[1,])]),])
  ok_name <- which(population_structure[,1] %in% matrix(non_NA_network[,1]) ) 
  result_population_reduced_non_NA <- population_structure[ok_name,]
  
  file_name <- paste("simulation",number_elements_node,datenum,sep = "_")
  #print(number_elements_node)
  epi <- network_dual
  #print("A")
  tree_population_structure <- modif_transmission_tree_2(epi,result_population_reduced_non_NA)
  tree_result_population2 <- tree_population_structure
  
  for(i in 1:length(ok_name))
  {
    query <- paste("\\(",result_population_reduced_non_NA[i,1],"\\[",sep="")
    replacement <- paste("\\(",i,"[",sep="")
    tree_result_population2 <- gsub(query, replacement,tree_result_population2)
    query2 <- paste("\\,",result_population_reduced_non_NA[i,1],"\\[",sep="")
    replacement2 <- paste("\\,",i,"[",sep="")
    tree_result_population2 <- gsub(query2, replacement2,tree_result_population2)
  }
  name_file_simple <- paste(file_name,"simple_tree.trees",sep="_")
  write(tree_result_population2,name_file_simple)
  
  name <- paste(result_population_reduced_non_NA[,1],"/",result_population_reduced_non_NA[,2],sep="")
  
  position <- 1:length(ok_name)
  traduction <- paste(position ,name, sep=" ")
  traduction_last <- traduction[length(traduction)]
  traduction <- traduction[-length(traduction)]
  name_tree <- paste(file_name,"tree.trees",sep="_")
  fileConn<- name_tree
  line1<-("#NEXUS")
  line2<-("Begin taxa;")
  line3<-(paste("Dimensions ntax=",length(name),";"))
  line4<-("Taxlabels")
  line5<-(name)
  line6<-(";")
  line7<-("End;")
  line8<-("Begin trees;")
  line9<-("Translate")
  line10<-(paste(traduction,",",sep=""))
  line11<-(traduction_last)
  line12<-(";")
  line13<-(paste("tree TREE1 =",tree_result_population2,";",sep=" ")) 
  line14<-("End;")
  #print("E")
  value <- value +1
  writeLines(c(line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12,line13,line14),fileConn)
  
  beast <- read.beast(name_tree)
  #print("F")
  p <- ggtree(beast,size=1, ladderize=FALSE,aes(color=location)) + theme_tree2(legend.position='left')
  p <- p + geom_tippoint(aes(color=location), size=3, alpha=.75) +
    scale_color_brewer("Legend", palette="Set1")+
    theme(panel.grid.major.x = element_line(color="black", size=0.2) ,
          panel.grid.minor.x = element_line(color="grey", size=0.1),
          axis.text.x= element_text(color = "steelblue", size = 10 ),
          legend.background = element_rect(colour = 'steelblue',fill="gray90", size=1),
          legend.title = element_text(colour = 'steelblue', size = 16),
          legend.text = element_text(colour="steelblue"),
          legend.key = element_rect(colour = 'steelblue',fill = "gray90", size =1))
  return(p)  
}

make_sequences <- function(network_dual,population_structure,number_elements_node,datenum) 
{  incProgress(0.1, detail = paste("Doing part", 1))
  Rpath <- getwd()
  epi_table  <- as.data.frame(network_dual)
  #print(network_dual)
  file_name <- paste("simulation",number_elements_node,datenum,sep = "_")
  print("ici")
  name_tree <- paste(file_name,"tree.trees",sep="_")
  tree_location <- paste(Rpath,name_tree,sep="/")
  name_sequence <- paste(file_name,".fasta",sep="")
  sequence_location <- paste(Rpath,name_sequence,sep="/")
  query <- paste("java -jar ",Rpath, "/pibuss.jar -treeFile ",tree_location," -from 1 -to 2000 -every 1 -branchSubstitutionModel HKY -HKYsubstitutionParameterValues 1.0 -clockRateModel strictClock -strictClockParameterValues 0.003 -siteRateModel gammaSiteRateModel -gammaSiteRateModelParameterValues 4.0 0.5 0.0 : ",sequence_location,sep="" )
  incProgress(0.2, detail = "Doing part 2 : calling piBuss")
  system(query)

  #return(query)
  incProgress(0.7, detail = paste("Doing part", 7))
  #rootname 	<- "D:/PhD/data/Simulation_analysis/"
  #fext		<- "name_sequence"
  fname		<- paste(Rpath,name_sequence,sep="/")
  seqs 		<- read.dna( fname, format="fasta", as.matrix=FALSE)
  taxa	<- as.matrix(attributes(seqs)$names )
  sequence_id <- apply(taxa, 1, getEl, ind=1, sep="/")
  length(sequence_id)
  sequence_pop <- apply(taxa, 1, getEl, ind=2, sep="/")
  taxa_in_epi <- match(as.numeric(sequence_id),epi_table$`Node ID`)
  newname_epi <- paste(epi_table$`Node ID`,population_structure[match(epi_table$`Node ID`,population_structure[,1]),2],epi_table$Rtime + abs(epi_table$Etime),sep="/" )
  new_name <- newname_epi[taxa_in_epi]
  pop_struct_ok <- population_structure[taxa_in_epi,]
  #print(paste("leng_new_names",length(new_name)))
  #print(paste("leng_pop_struct",length(population_structure[,2])))
  trait_table <- cbind(new_name,as.character(pop_struct_ok[,2]))
  colnames(trait_table) <- c("taxa","trait")
  write.table(trait_table,paste(file_name,"_trait_table.txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  newSeqs <- seqs
  #print(new_name)
  attributes(newSeqs)$names <- new_name
  #sname			<- paste("newnames_test_",new_name,sep="")
  #sname <- paste("test",sname,sep="_")
  incProgress(0.9, detail = paste("Doing part", 3))
  #write.dna(newSeqs, file=sname, format="fasta", nbcol=-1, colsep="")
  return(newSeqs)
}

###### VIZUAL CREATION FUNCTIONS #######

network_animation <- function(epidemia_network,number_elements_node,population_structure)#,Na_ori,Nb_ori,Nc_ori,Nd_ori)
{
  modif_epi_table <- epidemia_network
  first_line_modif_epi_table <- modif_epi_table[1,]
  modif_epi_table<-   modif_epi_table[complete.cases(modif_epi_table), ]
  if (length(modif_epi_table[,1]) <=3 )
  {
    stop("Network to small to produce an animation ")
    # break
  }
  modif_epi_table <- rbind(first_line_modif_epi_table,modif_epi_table)
  min_epi_table <- min(modif_epi_table[,3])
  modif_epi_table[,3] <- modif_epi_table[,3]- min_epi_table
  modif_epi_table[,4] <- modif_epi_table[,4]- min_epi_table
  modif_epi_table[,5] <- modif_epi_table[,5]- min_epi_table
  
  max_epi_table <- max(modif_epi_table[,3])
  epi_table_new <- modif_epi_table[-1,]
  min_epi_table <- min(epi_table_new[,3])
  #max_epi_table <- max(modif_epi_table[,3])
  colr_a <- rep("#8dd3c7",each =number_elements_node )
  colr_b <- rep("#ffffb3",each =number_elements_node )
#  colr_c <- rep("#bebada",each =Nc_ori )
 # colr_d <- rep("#fb8072",each =Nd_ori )
  
  colr <- c(colr_a,colr_b)#,colr_c,colr_d)
  #  colr <- rep(c("#8dd3c7","#ffffb3","#bebada","#fb8072"),each =number_elements_node )
  epi_network <- cbind(modif_epi_table[,2],modif_epi_table[,1])
  colnames(epi_network) <- c("from","to")
  epi_network <- epi_network[-1,]
  net1 <- as.network(epi_network)
  net1 %v% "population_structure" <- colr[net1 %v% "vertex.names"]
  activate.edges(net1,at=0)
  
  step_time <- round((max_epi_table - floor(min_epi_table))/10,digits = 2)
  
  vs <- data.frame(onset=0, terminus=max_epi_table, vertex.id=1:length(epi_table_new[,3]))
  es <- data.frame(onset=epi_table_new[,3], terminus=max_epi_table, 
                   head=epi_table_new[,2],
                   tail=epi_table_new[,1])
  as.data.frame(net1)
  net1.dyn <- networkDynamic(base.net=net1, edge.spells=es, vertex.spells=vs)
  population_structure <- as.data.frame(net1.dyn)
  
  
  start_value <- (floor(min_epi_table)-step_time)
  if(start_value < 0)
  {
    start_value = 0
  }
  
  #  filmstrip(net1.dyn, displaylabels=F, mfrow=c(1, 5),
  #          slice.par=list(start=start_value, end=max_epi_table , interval=step_time, 
  #                         aggregate.dur=step_time, rule='any'))
  
  compute.animation(net1.dyn, animation.mode = "kamadakawai",
                    slice.par=list(start=start_value, end=max_epi_table, interval=step_time, 
                                   aggregate.dur=step_time, rule='any'))
  
  output <- render.d3movie(net1.dyn, usearrows = T, 
                           displaylabels = T, label=net1 %v% "vertex.names",
                           bg="#ffffff", vertex.border="#333333",
                           vertex.col = net1.dyn %v% "population_structure",
                           edge.col = '#55555599',
                           vertex.tooltip = paste("<b>Name:</b>", (net1 %v% 'vertex.names') , "<br>",
                                                  "<b>Type:</b>", (net1 %v% 'population_structure')),
                           launchBrowser=T, filename="Media-Network-Dynamic.html",
                           render.par=list(tween.frames = 30, show.time = T),
                           plot.par=list(mar=c(0,0,0,0)),output.mode='htmlWidget')
  
  return(output)
}

create_small_network <- function(number_elements_node)
{
#number_elements_node <- 5
value <- 1
  file_name <- paste("simulation",number_elements_node,sep = "_")
  
  population_structure <- dual_population_structure(number_elements_node)
  
  dyadcov2 <- dual_population_dyadcov(number_elements_node)
  
  network_dual <- network_creation_dual(dyadcov2,number_elements_node)
  
  make_tree(network_dual,population_structure,number_elements_node,value) 
}

make_tree_of_forest_sequences <- function(number_elements_node,seqs,jn,jr,numPerCat,repsPerCat,datenum)
{ 
  incProgress(0.1, detail = paste("Doing part", 1))
Rpath <- getwd()  
jointTrait1		<- "trait"
sampleType 		<- "trait"
#numPerCat 	<- 10
#repsPerCat <- 10
#jn <- 25
#jr <- 10
#path 			<- "D:/PhD/data/Simulation_analysis/"
#name 			<- paste("newnames_",file_name,sep="")
file_name_ori <- paste("simulation",number_elements_node,datenum,sep = "_")
traitFileName <- paste(file_name_ori,"_trait_table",sep="")
#print("avant_le_call_fct_R")
list[seqs, trait_table,list_tree] <- tree_of_the_forest(Rpath,seqs,file_name_ori,traitFileName,jointTrait1,sampleType,numPerCat,repsPerCat,jn,jr)
return(list(first=seqs, second=trait_table,third = list_tree))
#file_name <- paste(file_name_ori,"_nj_tree_set",sep="")
}

make_sequences_df <- function(seqs)
{
  #print("make sequence df")
  #name_sequence <- paste(file_name,".fas",sep="")
  #  sequences_name <- paste("newnames_",name_sequence,sep="")
  #sequences_name <- paste(name_sequence,sep="")
 # print(seqs)
  test <-  matrix(0, nrow = length(seqs))
  #test[1,1]
  for( i in 1:length(seqs))
  {
    test[i,1] <- paste(unlist(as.character(seqs[i])),collapse='')
  }
  
  taxa	<- as.matrix( attributes(seqs)$names )
  
  df <- cbind(taxa,test)
  df <- data.frame(df)
  return(df)
}

change_taxa <- function(sequence_df,trait_table_samp)
{
  #print("change taxa")
  #  Trait_table_name <- paste(file_name,"_Trait_table.txt",sep="")
  #rait_table_name <- paste(file_name,"_traits.txt",sep="")
 # Trait_table <- read.table(Trait_table_name,sep="\t", header = TRUE)
  Trait_table <- trait_table_samp
  line1 <- paste('\t\t<taxon id="',sequence_df[1,1],'">',sep="")
  line2 <- paste("\t \t \t",'<attr name="Trait">',sep="")
  line3 <- paste("\t \t \t","\t",Trait_table[1,2],sep="")
  line4 <- paste("\t \t \t</attr>","\t \t</taxon>",sep="\n")
  
  final_txt <- paste(line1,line2,line3,line4,sep="\n")
  
  for(i in 2:length(sequence_df[,1]))
  {
    line1 <- paste('\t \t <taxon id="',sequence_df[i,1],'">',sep="")
    line2 <- paste("\t \t \t",'<attr name="Trait">',sep="")
    line3 <- paste("\t \t \t","\t",Trait_table[i,2],sep="")
    line4 <- paste("\t \t \t</attr>","\t \t</taxon>",sep="\n")
    
    text_new <- paste(line1,line2,line3,line4,sep="\n")
    final_txt <- paste(final_txt,text_new,sep="\n")
  }
  return(final_txt)
}

change_Trait <- function(trait_table_samp)
{#print("change trait")
  #  Trait_table_name <- paste(file_name,"_Trait_table.txt")
  #Trait_table_name <- paste(file_name,"_traits.txt",sep="")
  #Trait_table <- read.table(Trait_table_name,sep="\t", header = TRUE)
  Trait_table <- data.frame(trait_table_samp)
  states <- sort(unique(Trait_table[,2]))
  states <- gsub(" ", "", states, fixed = TRUE)
  #print(states)
  all_lines <- paste('\t\t<state code="',states[1],'"/>',sep="")
  
  for(i in 2:length(states))
  {
    #  i=1
    line1 <- paste('\t\t<state code="',states[i],'"/>',sep="")
    all_lines <- paste(all_lines,line1,sep="\n")
  }
  #all_lines <- gsub(" ", "", all_lines, fixed = TRUE)
 # print(all_lines)
  return(all_lines)
}

change_BSSVS_xml <- function(taxa_lines,trait_table_samp,Trait_lines,num,datenum)
{print("change bssvs file")
  #  Trait_table_name <- paste(file_name,"_Trait_table.txt")
  #Trait_table_name <- paste(file_name,"_traits.txt",sep="")
  #Trait_table <- read.table(Trait_table_name,sep="\t", header = TRUE)
  Trait_table <- trait_table_samp
  print(Trait_table)
  tree_set_name <- paste('"',Rpath,"/simulation_",num,"_",datenum,'_nj_tree_set.trees.txt"',sep="")
  #print(tree_set_name)
  bssvs_name <- paste(Rpath,"/simulation_",num,"_",datenum,"_BSSVS_analysis",sep="")
  print(bssvs_name)
  states <- sort(unique(Trait_table[,2]))
  n_states <- length(unique(states))
  
  xml_lines <- readLines("bssvs_template_file.xml")
  #print(xml_lines)
  xml_lines <- paste(xml_lines,"\n",sep="")
  xml_lines <- paste(unlist(xml_lines), collapse ="\n")
  
  xml_lines_ok_taxa <- gsub("taxa to change",taxa_lines,xml_lines)
  xml_lines_ok_Trait<- gsub("Trait to change",Trait_lines,xml_lines_ok_taxa)
  
  new_dim1 <- paste('<parameter id="Trait.frequencies" dimension="',n_states,'"/>',sep="")
  xml_new_dim1 <- gsub('<parameter id="Trait.frequencies" dimension="change_dim"/>',new_dim1,xml_lines_ok_Trait)
  
  new_dim2 <- paste('<parameter id="Trait.rates" dimension="',n_states*(n_states-1),'" value="1.0" lower="0.0"/>',sep="")
  xml_new_dim2 <- gsub('<parameter id="Trait.rates" dimension="change_dim" value="1.0" lower="0.0"/>',new_dim2,xml_new_dim1)
  
  new_dim3 <- paste('<parameter id="Trait.indicators" dimension="',n_states*(n_states-1),'" value="1.0"/>',sep="")
  xml_new_dim3 <- gsub('<parameter id="Trait.indicators" dimension="change_dim" value="1.0"/>',new_dim3,xml_new_dim2)
  
  new_dim4 <- paste('<parameter id="Trait.root.frequencies" dimension="',n_states,'"/>',sep="")
  xml_new_dim4 <- gsub('<parameter id="Trait.root.frequencies" dimension="change_dim"/>',new_dim4,xml_new_dim3)
  
  xml_new_dim4 <- gsub('"\t',"",xml_new_dim4)
  xml_new_dim4 <- gsub('>"',">",xml_new_dim4)
  
  xml_new_dim4 <- gsub('new_tree_set',tree_set_name,xml_new_dim4)
  xml_new_dim4 <- gsub("bssvs_template_file",bssvs_name,xml_new_dim4)
  
  xml_new_dim4 <- gsub("bssvs_template_tree",tree_set_name,xml_new_dim4)
  
  xml_new_dim4 <- gsub("number_states",length(states),xml_new_dim4)
  
  file_bssvs_name <- paste(bssvs_name,".xml",sep="")
  cat(xml_new_dim4,file = file_bssvs_name)
  return(file_bssvs_name)
  #print(bssvs_name)
  # setwd("C:/PhD software/BEAST v1.8.3/bin")
  #  bssvs_file <- paste("D:/PhD/data/Simulation_analysis/",file_bssvs_name,sep="")
  #  call <- paste("beast -beagle_gpu -beagle_order 1", bssvs_file)
  #  system(call)
  # setwd("D:/PhD/data/Simulation_analysis/")
  #return(file_bssvs_name)
 
  
}

make_bssvs_anal <- function(seqs,trait_table_samp,num,datenum)
{  incProgress(0.1, detail = "Perform prim.analysis")
  sequence_df <- make_sequences_df(seqs)
  taxa_changed <- change_taxa(sequence_df,trait_table_samp)
  trait_changed <- change_Trait(trait_table_samp)
  bssvs_name_run <- change_BSSVS_xml(taxa_changed,trait_table_samp,trait_changed,num,datenum)
  #print("avant_print_bssvs_xml")
  incProgress(0.1, detail = "Call to BEAST.")
  call <- paste("java -jar ",Rpath,"/beast.jar -overwrite -beagle_CPU -beagle_order 1 ", bssvs_name_run,sep="")
  incProgress(0.1, detail = "Running BEAST ... It might take a while.")
  #print(call)
  system(call)
  #print(bssvs_name_run)
  #Calculate_bf(Rpath,name,trait,ext,burn_in)
  #plot_bssvs(Rpath,name,trait)  
}

getStatesFromXML <- function( xmlName ) {
  xmlLines <- readLines(xmlName)
  is	   <- grep("<!\\-\\- Number Of States =", xmlLines)
  ie	   <- grep("</generalDataType", xmlLines)
  states   <- xmlLines[(is+1):(ie-1)]
  states   <- apply(as.matrix(states), 1, getEl, ind=2, sep="=")
  states   <- gsub("\"", "", states)
  states   <- gsub("/>", "", states)
  return( states )
}

calcBF	<- function( p, q ) {
  upper <- p/(1-p)
  lower <- q/(1-q)
  return( upper/lower )
}

calcBF_from_indicators	<- function( probAccept, numberOfStates=numberOfStates ) {
  p 	<- probAccept
  q 	<- calc_qk( numberOfStates )
  BF	<- calcBF(p, q)
  return( BF )
}

calc_qk	<- function( K, eta=log(2) ) {
  upper <- eta*(K-1)
  lower <- K*(K-2)/2
  return( upper/lower )
}

Calculate_bf <- function(num)
{
  trait <-"Trait"
  burn_in <- 10
  # load states from beauti xml file
  xmlName	<- paste(Rpath,"/simulation_",num,"_BSSVS_analysis.xml",sep="")
  states	<- getStatesFromXML( xmlName )
  
  Original_table <- read.table( paste(Rpath,"/simulation_",num,"_BSSVS_analysis.",trait,".rates.log",sep=""), header=TRUE, sep="\t")
  
  lines_to_remove <- round(length(Original_table[,1])/burn_in)
  Original_table <- Original_table[-(1:lines_to_remove),]
  
  n_col <- length(Original_table[1,])
  indic =  (n_col-2)/2
  numRows	<- length(Original_table[,1])
  numStates	<- length(states)
  
  #get qk
  qk = (log(2)+(indic -1))/(indic*(indic-2)/2)
  
  #selection of indicator columns
  location_rate_original <- Original_table[,1 :indic +1]
  names_location_rate_original <- colnames(location_rate_original)
  Indicator_original <- Original_table[,(indic + 2) :((indic * 2)+1 )]
  names_Indicator_original <- colnames(Indicator_original)
  
  fromState	<- apply(as.matrix(names_location_rate_original), 1, getEl, ind=3, sep="\\.")
  toState	<- apply(as.matrix(names_Indicator_original), 1, getEl, ind=4, sep="\\.")
  
  RelativeRates		<- array(0, indic)
  meanIndic			<- array(0, indic)
  
  for (j in 1:indic )
  {
    location_rate_original_col <- location_rate_original[,j]
    location_rate_original_col <- as.numeric(as.character(location_rate_original_col))
    inds				<- which( Indicator_original[j] == 1)
    RelativeRates[j] <- mean(location_rate_original_col)
    meanIndic[j]		<- length(inds)/numRows
    
  }
  
  bayesFactor	<- calcBF_from_indicators(meanIndic, numberOfStates=numStates)
  tbl 		<- cbind(fromState,toState,RelativeRates,meanIndic,bayesFactor)
 # print(tbl)
  
  outname	<- paste(path,name,".",trait,".mean_rates_with_bf.txt",sep="")
  #write.table( tbl, outname,
   #            sep="\t", col.names=TRUE, row.names=FALSE)
  
  #print( paste("Written results to",outname) )
  
  bfThres	<- 3
  sigRates	<- which(bayesFactor >= bfThres)
 # write.table( tbl[sigRates,], outname,
  #             sep="\t", col.names=TRUE, row.names=FALSE)
  #print(tbl[sigRates,])
  return(tbl[sigRates,])
}

plot_bssvs <- function(BSSVS_result)
{
  bssvs_table <- read.table("D:/PhD/R_file/plot_tranmission_network/Epidemia_10_hosts_sequences_sampled_BSSVS_result.csv",header=TRUE,sep=",")
 # bssvs_table_name <- paste(path,name,".",trait,".mean_rates_with_bf.txt",sep="")
  #bssvs_table <- read.table(bssvs_table_name,header=TRUE,sep="\t")
 # bssvs_table <- BSSVS_result
  df <- data.frame(bssvs_table[,5])
  df <- data.frame(as.numeric(gsub("Inf", "1000", as.matrix(df))))

  
  gap_stat <- clusGap(df, FUN = hcut, K.max = 10, B = 10)
  number_K <- with(gap_stat,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
  cluster_df <- kmeans(df,number_K)
  cluster_df <- data.frame(df,cluster_df$cluster)
  
  ncolor <- data.frame(palette(rainbow(number_K)))
  if(length(ncolor[,1])!=number_K)
  {
    ncolor <- data.frame(palette(rainbow(number_K)))
  }
  palette_color <- data.frame(cbind(ncolor,1: number_K))
  colnames(palette_color) <- c("color","id")
  
  #plot(test)
  epi_network <- bssvs_table[,1:2]
  
  links_net <- data.frame(epi_network)
  
  nodes_net <- data.frame(1:length(unique(links_net[,1])),cbind(as.character(unlist(unique(links_net[,1])))))
  colnames(nodes_net) <-c("id","population_structure")
  
  links_net$from<- nodes_net$id[match(links_net$fromState ,nodes_net$population_structure)]
  links_net$to <- nodes_net$id[match(links_net$toState,nodes_net$population_structure)]
  links_net$color<- palette_color$color[match(cluster_df$cluster_df.cluster,palette_color$id)]
  
  legend_edge<- data.frame()
  for(i in 1:number_K)
  {
    ii <- which(cluster_df$cluster_df.cluster == i)
    min_value <- round(min(df$bssvs_table...5.[ii]))
    max_value <- round(max(df$bssvs_table...5.[ii]))
    if(min_value == max_value)
    {
      interval_color <- max_value
    }else
    {
      interval_color <- paste(min_value,max_value,sep="_")
    }
    legend_edge[i,1] <- interval_color
  }
  
  legend_edge <- cbind(ncolor,legend_edge)
  
  colnames(legend_edge) <- c("color","label")
  links_net$interval<- legend_edge$id[match(links_net$color,legend_edge$color)]
  
  vis.nodes <- data.frame(nodes_net$id)
  colnames(vis.nodes) <-"id"
  vis.links <- links_net[,3:4] 
  
  vis.nodes$shape  <- "dot"  
  vis.nodes$group <- nodes_net$population_structure # add groups on nodes 
  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  vis.nodes$title  <- nodes_net$population_structure # Text on click
  #vis.nodes$label  <- nodes_net$population_structure # Node label
  vis.nodes$borderWidth <- 2 # Node border width
  
  vis.nodes$color.background <- c("#8dd3c7", "#ffffb3", "#bebada","#fb8072")[nodes_net$population_structure]
  vis.nodes$color.border <- "black"
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"
  
  vis.links$color <- links_net$color# line color  
  vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- TRUE    # should the edges be curved?
  vis.links$shadow <- TRUE
  visNetwork(vis.nodes, vis.links,main="BSSVS output",submain="Well supported transition routes")%>%
    visGroups(groupname = "A", color = "#8dd3c7") %>%
    visGroups(groupname = "B", color = "#ffffb3")%>%
    visGroups(groupname = "C", color = "#bebada")%>%
    visGroups(groupname = "D", color = "#fb8072")%>%
    #  visOptions(selectedBy = "group",highlightNearest = TRUE) %>%
    #  visEdges(shadow = TRUE,
    #           arrows =list(to = list(enabled = TRUE, scaleFactor = 2)),
    #           color = list(color = "lightblue", highlight = "red"))%>%
    visLegend( position = "right", main = "Legend",addEdges = legend_edge,ncol = 2)
}
