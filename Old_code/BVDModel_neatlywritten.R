###############################Code Content################################
#1a. Creating function
 #Code inside of function#
  #1. Setting parameters
  #2. Movement generators
  #3. BVD model with SimInf
  #4. Result

#2. Entering variable

#~#~#~#~#~#~#~#~#~#~#~#
#########Libraries#######
library(SimInf)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(systemfonts)
library(hrbrthemes)
library(viridisLite)
library(viridis)
library(purrr)
library(readr)
library(plyr)
library(reshape)
library(reshape2)
library(scales)


#################################Build function#################################
BVD_model <- function(allowed_import_A, animal_vaccine_A,animal_vaccine_B,
                      farm_vaccine_A, farm_vaccine_B){
  ###############################1.SETTING PARAMETERS############################
  
  #TIME
  #Moves over the whole period of time
  #total movement 128928
  number_of_moves <- 128928
  #number of time steps in month
  end <- 9*12
  tspan <- seq(from = 1, to= end)
  
  #NUMBER OF FARMS IN EACH GROUP
  # how many farms you want of each type? The first part initialises a vector
  number_of_farms <- c(0,0,0)
  # give the name to the group
  names(number_of_farms) <- c('A','B','C')
  # choose how many farms you want of each type
  #A 3891; B 1360; C 743
  number_of_farms['A'] <- 3891
  number_of_farms['B'] <- 1360
  number_of_farms['C'] <- 743
  # count the total
  total_farms <- sum(number_of_farms)
  #number of farms within Scotland
  Scotland_farms <- as.numeric(number_of_farms['A']+number_of_farms['B'])
  
  #NUMBER OF ANIMALS IN EACH FARM
  # n of animal in each type of farm
  animal_A_farms <- 100
  animal_B_farms <- 100
  animal_C_farms <- 100
  #setting the prevalence for farms in group C
  Type_C_prevalence <- 0.3
  Animal_C_prevalence <- 0.02
  
  #VACCINATION
  #Enter efficacy
  w <- 0.8
  #Coverage
  #coverage of vaccinated animals in a farm
  u <- c(0,0)
  #giving names to the farms' group
  names(u) <- c('A','B')
  #Choose animals vaccination coverage in a year
  u['A'] <- animal_vaccine_A
  u['B'] <- animal_vaccine_B
  #vaccination with efficacy
  uw <- u*w
  
  #coverage of vaccinated farms in each group
  #group A
  u_group_A <- farm_vaccine_A
  #group B
  u_group_B <- farm_vaccine_B
  
  #TRANSMISSION RATE PER MONTH
  # Beta  : transmission rate per month (0.03)
  # Mu    : removal rate (0.08/0.04/0.0132)
  # Delta : TI animal give birth to PI  
  # gamma : rate from recovered to susceptible 
  # Reference of rate from paper :
  # "Farm productive contexts and the dynamics of bovine viral diarrhea (BVD) transmission"
  # W     : Vaccine efficacy
  # u     : percentage individual vaccinated/unit time 
  # b     : birth rate
  # m     : mortality rate not related to BVD
  # s     : number of susceptible animal introduce of S in a node less than 10
  gdata <- c(beta = 0.03,
             delta = 0.00400756,
             Mu = 0.0132,
             gamma = 0.0833,
             u = 0,
             w = w,
             b = 0.02563349,
             m = 0.02563349,
             s = 100)
  
  #MOVEMENT PROPORTION MATRIX
  # create matrix 3 by 3
  move_prob <- c(rep(0,9))
  dim(move_prob) <- c(3,3)
  # then we add names so move_prob['A','B'], for example, can be called 
  rownames(move_prob) <- c('A','B','C')
  colnames(move_prob) <- c('A','B','C')
  # set the proportion of movement between the different types
  # These need to add to 1
  move_prob['C','A'] <- 0.04*(allowed_import_A)
  
  move_prob['A','A'] <- 0.22
  move_prob['A','B'] <- 0.06
  move_prob['A','C'] <- 0
  move_prob['B','A'] <- 0.04+((0.04-move_prob['C','A'])/2)
  move_prob['B','B'] <- 0.61
  move_prob['B','C'] <- 0
  
  move_prob['C','B'] <- 0.03+((0.04-move_prob['C','A'])/2)
  move_prob['C','C'] <- 0
  
  #########################2. MOVEMENT GENERATOR##################################
  #Creating stochastic movement from node origin to destination
  
  # make the farms. farms[['A']] will return a vector of farm numbers
  farms <- c(0,0,0)
  names(farms) <- c('A','B','C')
  #This is for events in SimInf
  # this part could be coded better. If you had more types then you might want to 
  # write a for loop instead of writing a new line for each type
  farms['A'] <- list(1:number_of_farms['A'])
  farms['B'] <- list((number_of_farms['A']+1):(number_of_farms['A']+number_of_farms['B']))
  farms['C'] <- list((number_of_farms['A']+number_of_farms['B']+1):total_farms)
  # start with empty lists
  node <- NULL
  dest <- NULL
  
  # compute the number of moves from each group
  for( node_type in c('A','B','C')){
    for( dest_type in c('A','B','C')){
      # calculate how many moves there will be between each group 
      moves <- move_prob[node_type,dest_type]*number_of_moves
      
      # choose random node farms of type 'node_type'
      node_nodes <- sample(farms[[node_type]],moves,replace=TRUE)
      # and random dest farms of type 'dest_type'
      dest_nodes <- sample(farms[[dest_type]],moves,replace=TRUE)
      
      # add these to the list to make the dataframe 
      node <-c(node,node_nodes)
      dest <-c(dest,dest_nodes)
      
    }
  }
  
  #print(node)
  #print(dest)
  
  # if 'moves' is not an integer in the above code, then 'node_nodes' will be 
  # than 'moves' and, in the end, 'node' and 'dest' will be smaller than
  # 'number_of_nodes'. To make them the same length we do
  number_of_moves <- length(node)
  
  # choose random times for these to happen
  time <- sample(1:end,number_of_moves,replace=TRUE)
  
  # now create the dataframe needed for SimInf
  events <- data.frame(
    event      = rep("extTrans", number_of_moves), ## Event "extTrans" is a movement between nodes
    time       = time,      ## The time that the event happens
    node       = node,      ## In which node does the event occur
    dest       = dest,      ## Which node is the destination node
    n          = rep(1, number_of_moves),      ## How many animals are moved
    proportion = rep(0, number_of_moves),      ## This is not used when n > 0
    select     = rep(1, number_of_moves),      ## 1st column in the select matrix (see below)
    shift      = rep(0, number_of_moves))      ## Not used in this example
  
  # remove movements that don't go anywhere i.e. node and dest are the same
  events<-events[events$node!=events$dest,]
  
  #Introduce import to group A and B
  #In here we change the SELECT column of nodes from group C as Enter event
  events$event <- ifelse(events$node>Scotland_farms, "enter", "extTrans")
  
  #Change the origin node of all Enter events as destination node
  #As an enter event directly transfer animal to the node 
  events$node <- ifelse(events$event=="extTrans", events$node,events$dest)
  
  #Change the destination node of all Enter events as 0
  #As an enter event does not need destination node
  events$dest <- ifelse(events$event=="extTrans", events$dest,0)
  
  #Randomly introduce susceptible and PI animals to the nodes with enter event
  #Number of PI animals introduced based on animal prevalence
  #The number 1/2/3 can be referred to E matrix 
  events$select[events$event == "enter"] <- sample(c(2,3),length(events$event[which(events$event=="enter")]),
                                                   prob = c((1-Animal_C_prevalence), Animal_C_prevalence),
                                                   replace = TRUE)
  
  # put the dataframe in order of time
  events <- events[order(events$time),]
  
  #VACCCINE COVERAGE
  #Vaccinating only once a year
  # if the animal coverage is set then this will run, 
  # otherwise there is no vaccination in a group 
  
  # group A
  if(animal_vaccine_A > 0){
    vaccination_A <- data.frame(event="intTrans",
                                time = rep(seq(from=12, to=end, 12),ceiling(u_group_A*number_of_farms['A'])),
                                node = sample(1:number_of_farms['A'],ceiling(u_group_A*number_of_farms['A']),replace=FALSE),
                                dest=0, n=0,
                                proportion= as.numeric(uw['A']), select=4,shift=1)
  }
  
  #group B
  if(animal_vaccine_B > 0 ){
    vaccination_B <- data.frame(event="intTrans", 
                                time=rep(seq(from=12, to=end, 12),ceiling(u_group_B*number_of_farms['B'])),
                                node=sample(number_of_farms['A']+1:number_of_farms['B'],ceiling(u_group_B*number_of_farms['B']),replace=FALSE),
                                dest=0,n=0,
                                proportion=as.numeric(uw['B']),select=4,shift=1)
  }
  
  #########################3. BVD MODEL WITH SIMINF###############################
  
  #COMPARTMENT
  compartment <- c("S", "PI", "PIL", "R", "V")
  #Define matrix E to handle which compartment to sample during a scheduled event
  E <- matrix(c(1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,1,0),
              nrow = 5, ncol=4,
              dimnames = list(c("S", "PI", "PIL", "R", "V"),
                              c("1", "2", "3","4")))
  
  #Define the shift matrix N to process internal transfer events: vaccination
  N <- matrix(c(4,3,2,1,0), nrow = 5, ncol = 1,
              dimnames = list(c("S","PI","PIL","R","V"),"1"))
  
  #ARRAY: NUMBER OF ANIMALS IN EACH COMPARTMENT
  #creating nodes with number of animals in each compartment
  #susceptible 
  S <- c(rep(animal_A_farms,number_of_farms['A']),
         rep(animal_B_farms,number_of_farms['B']))
  
  #final data.frame of animals in each node
  u0 <- data.frame(S = S,
                   PI = rep(0,Scotland_farms),
                   PIL = rep(0, Scotland_farms),
                   R = rep(0, Scotland_farms),
                   V = rep(0, Scotland_farms))
  
  #NETWORK FLOW
  transition <- c("S -> PI > 0 ? beta*S*PI/(S+PIL+PI+R+V) : 0 -> PIL",
                  "@ -> S < 10 ? s : 0 -> S",
                  "@ -> b*S -> S",
                  "S -> m*S -> @",
                  "@ -> b*R -> R",
                  "R -> m*R -> @",
                  "PIL -> beta*PIL -> R",
                  "@ -> delta*PIL -> PI",
                  "PI -> Mu*PI -> @",
                  "R -> gamma*R -> S")
  
  #RUNNING THE MODEL
  # type of model run according to the condition given
  # if there is vaccination in a group therefore the model will add table for vaccination
  if(animal_vaccine_B > 0){
    if(animal_vaccine_A >0){
      Model_vac <- mparse(transitions = transition,
                          compartments = compartment,
                          gdata = gdata, u0 = u0, tspan = tspan,
                          E = E, N = N,
                          events = rbind(events,vaccination_A,vaccination_B))
    } else {
      Model_vac <- mparse(transitions = transition,
                          compartments = compartment,
                          gdata = gdata, u0 = u0, tspan = tspan,
                          E = E, N = N,
                          events = rbind(events,vaccination_A))
    }
  } else{
    Model_vac <- mparse(transitions = transition,
                        compartments = compartment,
                        gdata = gdata, u0 = u0, tspan = tspan,
                        E = E, N = N,
                        events = events)
  }
  
  set_num_threads(1)
  result <- run(Model_vac)
  
  #################################4. RESULT######################################
  
  #Run trajectory
  infected <- trajectory(result)
  
  #NODES GROUPING
  Scotland_farms <- as.numeric(Scotland_farms)
  #rename nodes with group[A,B,C]
  infected$group[infected$node <= number_of_farms['A']] <- "GroupA"
  infected$group[infected$node > number_of_farms['A'] & infected$node <= Scotland_farms] <- "GroupB"
  #infected$group[infected$node > Scotland_farms] <- "GroupC"
  
  # GRAPH AVERAGE FOR ALL GROUPS
  # mean PI
  # this is for taking the number of animal in each compartement at each time step
  PI_mean <- ggplot(infected, aes(x = time, y = PI, color = group, group = group))+
    geom_line(stat = "summary", fun = "mean")+
    scale_x_continuous("Time[month]")+
    scale_y_continuous("Number of PI")+
    scale_color_viridis(discrete = TRUE)+
    scale_color_manual(values=c("#FFC602", "#3CAEA3", "#ED553B"))+
    theme_bw()
  
  #PREVALENCE
  #PI Prevalence based on group
  #create label for months and group
  month <- unique(infected$time)
  group <- unique(infected$group)
  
  #Calculating total number of PI at each time step in each group
  # creating empty variable
  PI_A<- 0
  PI_B<- 0
  PI_C<- 0
  
  # summing PI based on time and group
  for(i in group) {
    for(j in month){
      if (i == "GroupA"){
        PI_A[j] <- sum(infected$PI[which(infected$time == j &
                                           infected$group == "GroupA")])
      } else if ( i == "GroupB"){
        PI_B[j] <-  sum(infected$PI[which(infected$time == j &
                                            infected$group == "GroupB")])
      } else {
        PI_C[j] <- sum(infected$PI[which(infected$time == j &
                                           infected$group == "GroupC")])
      }
    }
  }
  
  # removing empty value
  total_PI <- na.omit(data.frame(PI_A,PI_B,PI_C))
  
  #counting total number of animals at each time step in each group
  infected <- infected %>% mutate(n_animal = S+PI+PIL+R+V)
  
  # creating empty variable for total animals
  n_A <- 0
  n_B <- 0
  n_C <- 0
  
  # sum the total number of animal based on their group 
  for(i in group) {
    for(j in month){
      if (i == "GroupA"){
        n_A[j] <- sum(infected$n_animal[which(infected$time == j &
                                                infected$group == "GroupA")])
      } else if ( i == "GroupB"){
        n_B[j] <-  sum(infected$n_animal[which(infected$time == j &
                                                 infected$group == "GroupB")])
      } else {
        n_C[j] <- sum(infected$n_animal[which(infected$time == j &
                                                infected$group == "GroupC")])
      }
    }
  }
  
  # remove empty value
  total_animal <- na.omit(data.frame(n_A,n_B,n_C))
  
  #merge the number of PI and total number of animals
  PI_prev <- cbind(total_PI, total_animal)
  #calculate the prevalence
  PI_prev <- transform(PI_prev, prev_A = (PI_prev$PI_A/PI_prev$n_A)*100,
                       prev_B=(PI_prev$PI_B/PI_prev$n_B)*100,
                       prev_C=(PI_prev$PI_C/PI_prev$n_C)*100)
  
  #SELECTING WHAT KIND OF RESULT YOU WANT
  
  #Result: mean for every compartement in each time step
  #fin_result <- data.frame(S= with(PI_mean$data, mean(S[group=="GroupA"])),
  # PI = with(PI_mean$data, mean(PI[group=="GroupA"])),
  #  PIL = with(PI_mean$data, mean(PIL[group=="GroupA"])),
  #  R = with(PI_mean$data, mean(R[group=="GroupA"])),
  #  V= with(PI_mean$data, mean(V[group=="GroupA"])),
  #  Prev_PI = mean(PI_prev$prev_A),
  #  import = allowed_import)
  
  #Result: all animal in EACH NODE each time point in each group (HEAVY)
  #all_node <- PI_mean[["data"]]
  #all_node$total <- rowSums(all_node[,3:7])
  #all_node$import <- ifelse(all_node$group =="GroupA", allowed_import, 1+(1-allowed_import))
  #all_node$vaccine <- 0
  #final_result <- all_node
  
  
  #Result: MEAN of all animals in each time point in each group
  #S_dat <- ggplot_build(S_mean)
  #PIL_dat <- ggplot_build(PIL_mean)
  #V_dat <- ggplot_build(V_mean)
  #R_dat <- ggplot_build(R_mean)
  #PI_data <- ggplot_build(PI_mean)
  
  #fin_result <- data.frame(S= do.call("rbind", S_dat$data) %>% select(x,group,y),
  #PI = do.call("rbind", PI_data$data) %>% select(y),
  #PIL = do.call("rbind", PIL_dat$data) %>% select(y),
  #R = do.call("rbind", R_dat$data) %>% select(y),
  #V= do.call("rbind", V_dat$data) %>% select(y),
  #Prev_PI_A = PI_prev$prev_A, 
  #Prev_PI_B = PI_prev$prev_B,
  #import = allowed_import)
  
  #names(fin_result)[names(fin_result)=="y"] <- "PI"
  #names(fin_result)[names(fin_result)=="y.1"] <- "PIL"
  #names(fin_result)[names(fin_result)=="y.2"] <- "R"
  #names(fin_result)[names(fin_result)=="y.3"] <- "V"
  
  #final_result <- fin_result
  
  #Result: ALL NUMBER of animals in each time points in each group
  #take data from ggplot
  total_animal <- PI_mean[["data"]]
  #sum based on group and time
  total_animal <- total_animal %>% 
    group_by(group,time) %>% 
    dplyr :: summarise(S=sum(S, na.rm = TRUE),PI=sum(PI, na.rm = TRUE),
                       TI=sum(PIL, na.rm = TRUE),R=sum(R, na.rm = TRUE),
                       V=sum(V, na.rm=TRUE))
  #add the sum of all animals in each time
  total_animal$total <- rowSums(total_animal[,3:7])
  
  #adding vaccine information
  total_animal$import<- ifelse(total_animal$group=="GroupA", allowed_import_A, 1+((1-allowed_import_A)/2))
  total_animal$vaccine <- ifelse(total_animal$group =="GroupA", "0", animal_vaccine_B)
  fin_result <- total_animal
  
}
################################### Entering variables########################
#set the proportion of import that you allowed to enter group A
allowed_import_A <- 1

#Enter animal vaccination coverage
animal_vaccine_A <- 1
animal_vaccine_B <- 1

#Enter farm vaccination coverage
farm_vaccine_A <- 0.8
farm_vaccine_B <- 0.8

#number of iteration
iteration <- 100

end_result<- vector("list", iteration)
#run the function
for (i in 1:iteration) {
  end_result[[i]] <- BVD_model(allowed_import_A, animal_vaccine_A, animal_vaccine_B,
                               farm_vaccine_A, farm_vaccine_B)
}
end_result <- do.call("rbind", end_result)

# save to working directory

setwd("D:/Roslin computer/Scotland_BVD/R model/Output")

# setwd("D:/WORK/IDOH - EMJMD/4 Sem - Internship/Data/Result/Scenarios/2/New Parameter")
write.csv(end_result, paste0("farm_vaccine_A",farm_vaccine_A,"farm_vaccine_B",farm_vaccine_B,".csv"),
          row.names = FALSE)