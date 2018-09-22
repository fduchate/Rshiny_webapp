modif_transmission_tree_2 <-function (epi,house,node = epi[1, 1], lastrow = 1, infectedatrow = 1, 
          rightchild = FALSE) 
{
  nex = ""
  atleaf = FALSE
  if (lastrow == nrow(epi)) {
    atleaf = TRUE
  }
  else {
    nextrow = match(node, epi[(lastrow + 1):nrow(epi), "Parent"])
    atleaf = is.na(nextrow)
  }
  if (atleaf) {
    if (rightchild) {
      test2 <- paste("[&location=",house[which( house[,1] == node),2],"]",sep="")
      return(paste( node,test2, ":", epi[infectedatrow, 
                                                      "Rtime"] - epi[infectedatrow, "Itime"]+ 
                   epi[infectedatrow, "Itime"] - epi[lastrow, "Etime"], 
                   sep = ""))
    }
    else {
#      print(node)
      test2 <- paste("[&location=",house[which( house[,1] == node),2],"]",sep="")
#      print(test2)
      return(paste(node,test2, ":", epi[infectedatrow, 
                                                 "Rtime"] - epi[lastrow, "Etime"], sep = ""))
    }
  }
  else {
    nextrow = nextrow + lastrow

    if (rightchild == TRUE) {
      leftside <- modif_transmission_tree_2(epi,house, node = node,lastrow = nextrow, infectedatrow = infectedatrow)
      left_tips_name <- gsub("\\(","",(sub(":.*$", "",leftside)))
      left_tips_name <- sub("\\[.*$", "",left_tips_name)
                   
      left_annotation <-paste("[&location=",house[which( house[,1] == left_tips_name),2],"]",sep="")
      
      rightside <- modif_transmission_tree_2(epi,house, node = epi[nextrow,"Node ID"], lastrow = nextrow, infectedatrow = nextrow,rightchild = TRUE)

   
      return(paste("(", leftside,
                   ",", rightside,
                   ")",
                   left_annotation,
                   ":", 
                   epi[lastrow,"Itime"] - epi[lastrow, "Etime"]+epi[nextrow, "Etime"] - epi[lastrow, "Itime"], 
                   sep = ""))
    }
    else {
      leftside <- modif_transmission_tree_2(epi,house, node = node,lastrow = nextrow, infectedatrow = infectedatrow)
      left_tips_name <- gsub("\\(","",(sub(":.*$", "",leftside)))
      left_tips_name <- sub("\\[.*$", "",left_tips_name)
      left_annotation <-paste("[&location=",house[which( house[,1] == left_tips_name),2],"]",sep="")

      rightside <- modif_transmission_tree_2(epi,house, node = epi[nextrow,"Node ID"], lastrow = nextrow, infectedatrow = nextrow,rightchild = TRUE)
            return(paste("(",  leftside,
                         ",",  rightside,
                         ")",
                   left_annotation,
                   ":", 
                   epi[nextrow,"Etime"] - epi[lastrow, "Etime"], sep = ""))
    }
  }
}
#test <- modif_transmission_tree_2(epi,house)

#gsub( ":.*$", "", modif_transmission_tree_2(epi, node = node,lastrow = nextrow, infectedatrow = infectedatrow))