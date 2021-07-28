


#'@Relative  Which method should be used to estimate contribution, can be "PandomPoints","Volume" or "Percentage". Note that values of contribution differ between methods 
#'@ if TRUE returns a vector of length = nrow(hv@RandomPoints) with the number of the species which the random point was assigned.


kernel.spp.contrib<-function(comm = NULL,trait = NULL,relative= "Percentage", method="gaussian", level="alpha", abund=FALSE, verbose=FALSE, parallel=getOption("mc.cores"), return.hv=TRUE){
  


  #Auxiliary functions are loaded now ----
  
  
  #Contribution to alpha 
  spp.contrib<-function(hv,rel=relative){
    
    require(tidyverse)
    spp.density<-lapply(1:nrow(hv@Data),
                        function(x) sapply(1:nrow(hv@RandomPoints),
                                           function(i) sqrt(sum((hv@Data[x,] - hv@RandomPoints[i,])^2))))
    #merge species column wise
    spp.columns<-do.call(cbind,spp.density)
    colnames(spp.columns)<-rownames(hv@Data)
    #assign the random points to the closest species (minimum euclidean distance)
    neighbours.dist<-colnames(spp.columns)[apply(spp.columns,1,which.min)]
    
    #estimate the contribution to functional beta diversity using:
    #PointDensity
    #Volume 
    #NumberofRandomPoints
    
    spp.contrib<-c(table(neighbours.dist))
    
    if (rel=="RandomPoints") result<-data.frame(spp.contrib) else 
      if(rel=="Volume") result<-data.frame(spp.contrib/hv@PointDensity) else 
        if(rel=="Percentage") {result<- data.frame(spp.contrib)/hv@PointDensity/hv@Volume}  #testing with volume of the union}
    
    result<-result %>% rownames_to_column()
    colnames(result)<-c("Spp",hv@Name)
    
    output<-result
    return(output)
  }
  
  #Contribution to beta
  spp.contrib.beta<-function(hvlist,rel.=relative){

    #Assign trait data to elements
    hvlist@HVList$Union@Data<- unique(rbind(hvlist@HVList$HV1@Data,hvlist@HVList$HV2@Data))
    hvlist@HVList$Unique_1@Data<- hvlist@HVList$HV1@Data
    hvlist@HVList$Unique_2@Data<- hvlist@HVList$HV2@Data
    hvlist@HVList$Intersection@Data<- unique(rbind(hvlist@HVList$HV1@Data,hvlist@HVList$HV2@Data))
    
    #estimate contribution for each element
    #union.contrib<-spp.contrib(hvlist@HVList$Union,relative) #Real Union from hypervolume_set
    unique1.contrib<-spp.contrib(hvlist@HVList$Unique_1,rel=rel.) #contribution to unique comm A
    unique2.contrib<-spp.contrib(hvlist@HVList$Unique_2,rel=rel.) #contribution to unique comm B
    shared.contrib<-spp.contrib(hvlist@HVList$Intersection,rel=rel.) #shared contribution of comm A and B
    
    #Combine the results into a single data frame
    comm.contrib<-full_join(unique1.contrib,
                            full_join(unique2.contrib,shared.contrib,by="Spp"),
                            by="Spp")
    names(comm.contrib)<-c("Spp","Unique_1","Unique_2","Shared") #assign names to the columns
    comm.contrib[is.na(comm.contrib)]<-0 #Convert NAs to zero
    comm.contrib$Union<- comm.contrib$Shared+comm.contrib$Unique_1+comm.contrib$Unique_2 #Union has to be build from the sum of components
    
    Total<-c()
    Brepl<-c()
    Brich<-c()
    for (i in 1:nrow(comm.contrib)){
      Total[i]<-comm.contrib$Shared[i] / comm.contrib$Union[i]
      Brepl[i]<-2*min(comm.contrib$Unique_1[i],comm.contrib$Unique_2[i]) / comm.contrib$Union[i]
      Brich[i]<-abs(comm.contrib$Unique_1[i]-comm.contrib$Unique_2[i]) / comm.contrib$Union[i]
    }
    #Divide by the total value so we can get percentage for each beta component
    spp.contrib.components<-data.frame((cbind(Total,Brepl,Brich)*comm.contrib$Union)/sum(comm.contrib$Union))
    rownames(spp.contrib.components)<- comm.contrib$Spp
    
    return(spp.contrib.components)
  }
  
  
  
  # Commands to test suitability of data begin now ---
  
  if (!(class(comm)[1] %in% c("Hypervolume", "data.frame", 
                              "matrix","HypervolumeList","List"))) stop("A Hypervolume, or a sites x species matrix or data.frame is needed as input data.")
  if (class(comm)[1]  %in% c("data.frame", "matrix")) {
    
    if (!isTRUE(ncol(comm) == nrow(trait))) 
      stop("Number of species in comm and trait matrices are different")
    if (any(is.na(c(comm,trait)))) 
      stop("The function cannot be computed with missing values. Please remove observations with missing values.")
    if (class(comm)[1] == "data.frame") 
      comm <- as.matrix(comm)
    if (class(trait)[1] == "data.frame") 
      trait = as.matrix(trait)
    if (is.null(colnames(comm))) 
      colnames(comm) = paste(rep("Sp", ncol(comm)), 1:ncol(comm), 
                             sep = "")
    if (is.null(rownames(trait))) 
      rownames(trait) = paste(rep("Sp", nrow(trait)), 1:nrow(trait), 
                              sep = "")
    comm2 = comm[rowSums(comm) > 1, ]
    if (nrow(comm2) != nrow(comm)) 
      warning(paste("In the site x species matrix (comm), one or more rows contain 0 or 1 species.\n  These rows have been removed prior to hypervolume estimation."))
    comm <- comm2
    nComm <- nrow(comm)
    subComm <- comm[1, ]
    subTrait <- trait[comm[1, ] > 0, ]
    subComm <- subTrait[rep(1:nrow(subTrait), times = subComm[comm[1, ] > 0]), ]
    if (abund) {
      subComm <- lapply(1:nrow(comm), function(s) comm[s, ])
      subTrait <- lapply(1:nrow(comm), function(s) trait[comm[s, ] > 0, ])
      subComm <- lapply(1:length(subComm), function(s) subTrait[[s]][rep(1:nrow(subTrait[[s]]), times = subComm[[s]][comm[s,] > 0]), ])}
    else {
      subComm <- lapply(1:nrow(comm), function(s) trait[comm[s, ] > 0, ])
    }

    #Checking availability of parallel
 parallel<-5   
    if (is.null(parallel)) 
      parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
      if (.Platform$OS.type == "unix" && !hasClus) { #both OSX and Linux are "Unix" and run Parallel
        hv <- do.call(hypervolume_join, pbmclapply(1:length(subComm), 
                                                                function(x,...) do.call(paste0("hypervolume_",
                                                                                                match.arg(method,c("svm","box","gaussian"))),
                                                                                         list(subComm[[x]],
                                                                                              name=rownames(comm)[x],
                                                                                              verbose=FALSE,...)),
                                                                mc.cores = parallel,mc.style = "ETA"))
      }
      else {
        if (!hasClus) {
          parallel <- makeCluster(parallel)
        }
        clusterExport(cl=parallel,c("subComm","method","verbose","comm","..."))
        clusterEvalQ(cl=parallel,c(library("hypervolume")))
        list.results <-pblapply(cl=parallel,
                                      X=1:length(subComm),
                                      fun = function(x,...) do.call(paste0("hypervolume_",
                                                                           match.arg(method,c("svm","box","gaussian"))),
                                                                    list(subComm[[x]],
                                                                         name=rownames(comm)[x],
                                                                         verbose=verbose,...)))
        if (!hasClus) 
          stopCluster(parallel)
        hv<-do.call(hypervolume_join,list.results)
      }
    } else{  #NO PARALLEL AT ALL
      
      hv<-do.call(hypervolume_join,lapply(1:length(subComm),
                           function(x) do.call(paste0("hypervolume_",match.arg(tolower(method),c("svm","box","gaussian"))),
                                               list(subComm[[x]],
                                                    name=rownames(comm)[x],
                                                    verbose=verbose))))
    }
  }  else {
    
    hv <- comm
  }
  
  
  #Contribution to hypervolume start now
 results<-list() 
  
  if(class(hv)[1] == "Hypervolume") { 
    
    if (is.null(rownames(hv@Data))) 
      rownames(hv@Data) <- paste(rep("Sp", nrow(hv@Data)), 
                                 1:nrow(hv@Data), sep = "")
    if(is.null(hv@Data)) hv@Data <- trait
    if(is.null(hv@Data)) stop("Trait data is missing from hypervolume and/or not given as data.frame")
    results$Contribution<-  spp.contrib(hv, rel=relative) }
  
  if(class(hv)[1] == "HypervolumeList" && level == "alpha")  {
    output<- lapply(1:length(hv@HVList),function(x) spp.contrib(hv@HVList[[x]],relative))
    output<- output %>% reduce(full_join,by="Spp")
    output[is.na(output)]<-0
    results$Contribution<-output
    if (isTRUE(return.hv))  results$hypervolumes<-hv
  } 
  
  if(class(hv)[1] == "HypervolumeList" && level == "beta"){
    
    if(is.null(hv@HVList$Intersection)){
      
      nComm<-length(hv@HVList)
      pairwise.contrib<-list()
      z=0
      for (i in 1:nComm){
        hyper =  hv@HVList[[i]]
        for (j in i:nComm){
          
          if (i == j){
            next
          } else {
            z=z+1
            {hyper2 <- hv@HVList[[j]]
            hyperSet <- hypervolume_set(hyper, hyper2, check.memory = FALSE, 
                                        verbose = FALSE, num.points.max = 10000)
            pairwise.contrib[[z]]<-spp.contrib.beta(hyperSet,rel.=relative)
            names(pairwise.contrib)[z]<-paste("Pair:",hyper@Name,"and",hyper2@Name)}
            print(z)
            print(paste("Pair:",hyper@Name,"and",hyper2@Name))
          }        
        }
      }
      results$Contribution<-pairwise.contrib
    }  else { results$Contribution<- spp.contrib.beta(hvlist,relative) }
  }  
  

  return(results)}

#
#
# #Species contribution to functional beta diversity (based on hypervolumes)
# # #
# library("BAT")
 library("hypervolume")
 library("pbmcapply")
comm <- rbind(c(1,1,0,0), c(0,0,1,1), c(1,1,1,0))
rownames(comm) <- c("Community_1","Community_2","Community_3")
colnames(comm) <- c("Sp_1","Sp_2","Sp_3","Sp_4")

#community 1 and 2 = Complete Turnover (taxonomic)
#community 2 and 3 = Richness differences (Taxonomic)

trait <- scale(cbind(c(2.2,4.4,6.1,8.3),c(0.5,1,0.5,0.4),c(0.7,1.2,0.5,0.4)))
rownames(trait) <- c("Sp_1","Sp_2","Sp_3","Sp_4")
colnames(trait) <- c("Trait_1","Trait_2","Trait_3")

test<-kernel.spp.contrib(comm=comm,trait=trait,relative="RandomPoints",level="alpha",parallel=3)
