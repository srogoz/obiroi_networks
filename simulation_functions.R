spread_SI_anttype_ab_subsets <-function(edge_list, spreading_prob, seed_ant, duration, present_ants){
  #SI-spreading process on timeprdered edgelist with a seed ant
  
  
  
  
  infection_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  a_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(a_hist)<-present_ants
  b_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(b_hist)<-present_ants
  colnames(infection_hist)<-present_ants
  #SI, no recovery infect seed ant
  
  infection_hist[1:duration, seed_ant]<-1
  
  #assign genotye of seedant
  seed_idx <- match(seed_ant, present_ants)
 
  gen <- genotype_match[[seed_ant]]
  
  if (gen == "a") {
    a_hist[1:duration, seed_idx] <- 1
  } else if (gen == "b") {
    b_hist[1:duration, seed_idx] <- 1
  }
  
  #order interactions by time
  timeordered <- edge_list[order(edge_list$onset),]
  number_interactions<-nrow(timeordered)
  
  
  for (i in seq_len(number_interactions)){
    time_step = timeordered$onset[i]
    
    #respect simulation duration
    if (time_step > duration) next
    
    #interaction-partners
    partner1<-timeordered$tail[i]
    partner2<-timeordered$head[i]
    
    #collect infected ants
    infected_ants<-present_ants[infection_hist[time_step, ] == 1]
    
    
    #if both infected, skip
    
    if (partner2 %in% infected_ants && partner1 %in% infected_ants){
      next
    }
    
    # transmission from 1 to 2 with prob
    if (partner1 %in% infected_ants && !(partner2 %in% infected_ants)){
      if (runif(1) < spreading_prob) {
        # once infected, stays infected for all future times
        seed_idx <- match(partner2, present_ants)
        
        infection_hist[time_step:duration, seed_idx] <-1
        if(genotype_match[[partner2]] == "a"){
          a_hist[time_step:duration, seed_idx] <-1
        }else{
          b_hist[time_step:duration, seed_idx] <-1
        }
      }
    }
    #transmission from 2 to 1
    if (partner2 %in% infected_ants && !(partner1 %in% infected_ants)){
      if(runif(1) < spreading_prob){
        seed_idx <- match(partner1, present_ants)
        infection_hist[time_step:duration, seed_idx ] <-1
        if(genotype_match[[partner1]] == "a"){
          a_hist[time_step:duration, seed_idx] <-1
        }else{
          b_hist[time_step:duration, seed_idx] <-1
        }
      }
    }
    #end simulation if everyone is infected  
    if (length(infected_ants) == length(present_ants)) break
    saturated<-timeordered$terminus[i]
    
  }
  
  prevalence<-rowSums(infection_hist)
  prevalence_a<-rowSums(a_hist)
  prevalence_b<-rowSums(b_hist)
  return(list(infection_hist = infection_hist, 
              prevalence_a= prevalence_a,
              prevalence_b= prevalence_b,
              prevalence = prevalence,
              saturated = saturated)
  )
  
}

spread_SI_anttype_all_subsets <-function(edge_list, spreading_prob, seed_ant, duration, present_ants, anttype, treatment){
  #SI-spreading process on timeprdered edgelist with a seed ant
   infection_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  other_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(other_hist)<-present_ants
  b_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(b_hist)<-present_ants
  colnames(infection_hist)<-present_ants
  #SI, no recovery infect seed ant
  
  infection_hist[1:duration, seed_ant]<-1
  
  #assign genotye of seedant
  seed_idx <- match(seed_ant, present_ants)
  
  gen <- anttype_match[[treatment]][[seed_ant]]
  
  if (gen == anttype) {
    other_hist[1:duration, seed_idx] <- 1
  } else if (gen == "b") {
    b_hist[1:duration, seed_idx] <- 1
  }
  
  #order interactions by time
  timeordered <- edge_list[order(edge_list$onset),]
  number_interactions<-nrow(timeordered)
  
  
  for (i in seq_len(number_interactions)){
    time_step = timeordered$onset[i]
    
    #respect simulation duration
    if (time_step > duration) next
    
    #interaction-partners
    partner1<-timeordered$tail[i]
    partner2<-timeordered$head[i]
    
    #collect infected ants
    infected_ants<-present_ants[infection_hist[time_step, ] == 1]
    
    
    #if both infected, skip
    
    if (partner2 %in% infected_ants && partner1 %in% infected_ants){
      next
    }
    
    # transmission from 1 to 2 with prob
    if (partner1 %in% infected_ants && !(partner2 %in% infected_ants)){
      if (runif(1) < spreading_prob) {
        # once infected, stays infected for all future times
        seed_idx <- match(partner2, present_ants)
        
        infection_hist[time_step:duration, seed_idx] <-1
        if(anttype_match[[treatment]][[partner2]] == anttype){
          other_hist[time_step:duration, seed_idx] <-1
        }else{
          b_hist[time_step:duration, seed_idx] <-1
        }
      }
    }
    #transmission from 2 to 1
    if (partner2 %in% infected_ants && !(partner1 %in% infected_ants)){
      if(runif(1) < spreading_prob){
        seed_idx <- match(partner1, present_ants)
        infection_hist[time_step:duration, seed_idx ] <-1
        if(anttype_match[[treatment]][[partner1]] == anttype){
          other_hist[time_step:duration, seed_idx] <-1
        }else{
          b_hist[time_step:duration, seed_idx] <-1
        }
      }
    }
    #end simulation if everyone is infected  
    if (length(infected_ants) == length(present_ants)) break
    saturated<-timeordered$terminus[i]
    
  }
  
  prevalence<-rowSums(infection_hist)
  prevalence_other<-rowSums(other_hist)
  prevalence_b<-rowSums(b_hist)
  return(list(infection_hist = infection_hist, 
             prevalence_other= prevalence_other,
              prevalence_b= prevalence_b,
              prevalence = prevalence,
              saturated = saturated)
  )
  
}
