library(SimInf)
library(ggplot2)

##########################BVD model function ################################

BVD_model <- function(allowed_import_A, animal_vaccine_A,animal_vaccine_B,
                      farm_vaccine_A, farm_vaccine_B){
  ##########################1.SETTING PARAMETERS############################
  
  ## ##TIME
  
  ## Total movements over the whole period of time
  number_of_moves <- 1458258
  
  ## number of time steps in month
  end <- 10*12
  tspan <- seq(from = 1, to= end)
  
  ## ## NUMBER OF FARMS IN EACH GROUP
  
  ## Create empty value
  number_of_farms <- c(0,0,0)
  
  ## Assign name to each group
  names(number_of_farms) <- c('A','B','C')
  # group C represent non-Scottish farms
  
  ## How many farms for each group
  ## A 3891; B 1360; C 743
  number_of_farms['A'] <- 32129
  number_of_farms['B'] <- 5918
  number_of_farms['C'] <- 10
  
  ## Count the total farms 
  total_farms <- sum(number_of_farms)
  
  ## Define number of farms within Scotland
  Scotland_farms <- as.numeric(number_of_farms['A']+number_of_farms['B'])
  
  ## ## NUMBER OF ANIMALS IN EACH FARM
  
  ## Number of animal in each group
  animal_A_farms <- 100
  animal_B_farms <- 100
  animal_C_farms <- 100
  
  ## Setting the PI prevalence for farms in group C
  Type_C_prevalence <- 0.3
  Animal_C_prevalence <- 0.02
  
  ## ## VACCINATION
  
  ## Enter vaccine efficacy
  w <- 0.8
  
  ## Coverage
  # Coverage of vaccinated animals in a farm
  u <- c(0,0)
  
  ## Assign names to the farms' group
  names(u) <- c('A','B')
  
  ## Choose animals vaccination coverage in a year
  u['A'] <- animal_vaccine_A
  u['B'] <- animal_vaccine_B
  
  ## Vaccination efficacy
  uw <- u*w
  
  ## Coverage of vaccinated farms in each group
  #group A
  u_group_A <- farm_vaccine_A
  #group B
  u_group_B <- farm_vaccine_B
  
  
  ## ## PARAMETERS
  # Beta  	: transmission rate per month
  # Mu    	: PI removal rate
  # Delta 	: TI become PI rate. 
  # gamma 	: rate from recovered to susceptible 
  # w     	: Vaccine efficacy
  # u     		: percentage individual vaccinated/unit time 
  # b     		: birth rate
  # m		: mortality rate not related to BVD
  # a		: rate of TI animals giving birth to normal calves
  # epsilon 	: TI animals to recovered
  
  gdata <- c(beta = 0.5,
             delta = 0.0094,
             Mu = 0.083,
             gamma = 0.0833,
             u = 0,
             w = w,
             b = 0.0246,
             m = 0.0246,
             s = 100,
             a = 0.0098,
             epsilon = 0.189)
  
  
  ## ## MOVEMENT PROPORTION MATRIX
  ## Create matrix 3 by 3
  move_prob <- c(rep(0,9))
  dim(move_prob) <- c(3,3)
  
  ## Then we add names so move_prob['A','B']
  rownames(move_prob) <- c('A','B','C')
  colnames(move_prob) <- c('A','B','C')
  
  ## Set the proportion of movement between the different types
  # these need to add to 1
  
  move_prob['C','A'] <- 0.04*(allowed_import_A)
  
  move_prob['A','A'] <- 0.66
  move_prob['A','B'] <- 0.03-(1-allowed_import_A)*0.02
  move_prob['A','C'] <- 0.03+(1-allowed_import_A)*0.02
  move_prob['B','A'] <- 0.11+(1-allowed_import_A)*0.02
  move_prob['B','B'] <- 0.08
  move_prob['B','C'] <- 0.04-(1-allowed_import_A)*0.02
  
  move_prob['C','B'] <- 0.01+(1-allowed_import_A)*0.02
  move_prob['C','C'] <- 0
  
  #####################2. MOVEMENT GENERATOR###########################
  
  ## ##Creating stochastic movement from node origin to destination
  
  ## Make the farms. farms[['A']] will return a vector of farm numbers
  farms <- c(0,0,0)
  names(farms) <- c('A','B','C')
  
  ## This is for events in SimInf
  farms['A'] <- list(1:number_of_farms['A'])
  farms['B'] <- list((number_of_farms['A']+1):(number_of_farms['A']+number_of_farms['B']))
  farms['C'] <- list((number_of_farms['A']+number_of_farms['B']+1):total_farms)
  
  
  ## Create an empty matrix to accomodate number of moves from each group
  # generated from the probability
  moves <- matrix(nrow = 3, ncol = 3)
  
  # then we add names so move_prob['A','B'],
  rownames(moves) <- c('A','B','C')
  colnames(moves) <- c('A','B','C')
  
  ## Calculate the number of moves from each group
  # and put into matrix
  for( node_type in c('A','B','C')){
    for( dest_type in c('A','B','C')){
      # calculate how many moves there will be between each group 
      moves[node_type,dest_type] <- ceiling(move_prob[node_type,dest_type]*number_of_moves)
    }
  }
  
  # create empty list
  node <- NULL
  dest <- NULL
  
  # Loop to randomly select nodes origin and destination
  for (node_name in c('A','B','C')) {
    for (dest_name in c('A','B','C')) {
      
      # choose random node of origin
      # select the origin based on farm list
      # and the number of moves based on the row matrix
      node_nodes <- sample(farms[[node_name]],moves[node_name,dest_name],replace = TRUE)
      
      # choose random destination
      # select the destination based on farm list
      # and the number of moves based on the row matrix
      dest_nodes <- sample(farms[[dest_name]],moves[node_name,dest_name],replace = TRUE)
      
      # add to the empty list to make data frame
      node <- c(node, node_nodes)
      dest <- c(dest, dest_nodes)
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
  
  ##  Now create the dataframe needed for SimInf
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
  
  # Introduce import to group A and B
  #In here we change the SELECT column of nodes from group C as Enter event
  events$event <- ifelse(events$node>Scotland_farms, "enter", "extTrans")
  
  # Change the origin node of all Enter events as destination node
  # As an enter event directly transfer animal to the node 
  events$node <- ifelse(events$event=="extTrans", events$node,events$dest)
  
  # Change the destination node of all Enter events as 0
  # As an enter event does not need destination node
  events$dest <- ifelse(events$event=="extTrans", events$dest,0)
  
  # Randomly introduce susceptible and PI animals to the nodes with enter event
  # Number of PI animals introduced based on animal prevalence
  # The number 1/2/3 can be referred to E matrix 
  events$select[events$event == "enter"] <- sample(c(2,3),length(events$event[which(events$event=="enter")]),
                                                   prob = c((1-Animal_C_prevalence), Animal_C_prevalence),
                                                   replace = TRUE)
  
  # put the dataframe in order of time
  events <- events[order(events$time),]
  
  ## ## VACCCINE COVERAGE
  ## Vaccinating only once a year
  #group B
  if(animal_vaccine_B > 0 ){
    # if there is vaccination in place: 0 mean there is no vaccination and 
    # 1 mean all animals get vaccinated
    
    vaccination_B <- data.frame(event="intTrans", 
                                time=rep(seq(from=12, to=end, 12),ceiling(u_group_B*number_of_farms['B'])),
                                node=sample(number_of_farms['A']+1:number_of_farms['B'],ceiling(u_group_B*number_of_farms['B']),replace=FALSE), 
                                # we randomly choose susceptible animals in all farms in Group B
                                dest=0,n=0,
                                proportion=as.numeric(uw['B']),select=4,shift=1)
    # select and shift is assigning these vaccinated animals in SimInf matrix
    # see SimInf CRAN Project: specification of scheduled events
  }
  
  # Similar like command above but for 
  # group A
  if(animal_vaccine_A > 0){
    vaccination_A <- data.frame(event="intTrans",
                                time = rep(seq(from=12, to=end, 12),ceiling(u_group_A*number_of_farms['A'])),
                                node = sample(1:number_of_farms['A'],ceiling(u_group_A*number_of_farms['A']),replace=FALSE),
                                dest=0, n=0,
                                proportion= as.numeric(uw['A']), select=4,shift=1)
  }
  
  ######################3. BVD MODEL WITH SIMINF ############################
  
  ## ## COMPARTMENT
  ## Creating model's compartment 
  compartment <- c("S", "PI", "PIL", "R", "V")
  
  ## Define matrix E to handle which compartment to sample during a scheduled event
  # Matrix E is used to specify scheduled events between disease compartments 
  E <- matrix(c(1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,1,0), 
              # the binary value is to allow animals from a compartment,
              # to go to certain compartments
              
              nrow = 5, ncol=4, 
              dimnames = list(c("S", "PI", "PIL", "R", "V"),
                              c("1", "2", "3","4")))
  
  # Matrix N is to process internal transfer events, which is vaccination
  N <- matrix(c(4,3,2,1,0), nrow = 5, ncol = 1,
              dimnames = list(c("S","PI","PIL","R","V"),"1"))
  
  ## ## NUMBER OF ANIMALS IN EACH COMPARTMENT
  ## Creating nodes with number of animals in each compartment
  # at the beginning of the simulation, all animals in all farms are susceptible 
  S <- c(rep(animal_A_farms,number_of_farms['A']),
         rep(animal_B_farms,number_of_farms['B']),
         rep(animal_C_farms,number_of_farms['C']))
  
  ## Create data frame for the start of the simulation
  u0 <- data.frame(S = S,
                   PI = rep(0,total_farms),
                   PIL = rep(0, total_farms),
                   R = rep(0, total_farms),
                   V = rep(0, total_farms))
  
  ## ## TRANSITION BETWEEN COMPARTMENTS
  transition <- c("S -> PI > 0 ? beta*S*PI/(S+PIL+PI+R+V) : 0 -> PIL",
                  "@ -> b*S -> S",
                  "S -> m*S -> @",
                  "PIL -> epsilon*PIL -> R",
                  "@ -> delta*PIL -> PI",
                  "PI -> Mu*PI -> @",
                  "R -> gamma*R -> S",
                  "@ -> a*PIL -> S",
                  "@ -> S < 10 ? S+s : 0 -> S")
  
  ## ## RUNNING THE MODEL
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
                          events = rbind(events,vaccination_B))
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
  
  #############################4. RESULT################################
  
  ## ## RUN TRAJECTORY
  infected <- trajectory(result)
  
  ## ##NODES GROUPING
  ## we then labels the nodes' ID into assigned group
  Scotland_farms <- as.numeric(Scotland_farms)
  #rename nodes with group[A,B,C]
  infected$group[infected$node <= number_of_farms['A']] <- "GroupA"
  infected$group[infected$node > number_of_farms['A'] & infected$node <= Scotland_farms] <- "GroupB"
  # since we are only interested in Group A and B, we ignore imports (group C)
  
  ## ## GRAPH AVERAGE FOR ALL GROUPS
  #mean PI
  PI_mean <- ggplot(infected, aes(x = time, y = PI, color = group, group = group))+
    geom_line(stat = "summary", fun = "mean")+
    scale_x_continuous("Time[month]")+
    scale_y_continuous("Number of PI")+
    scale_color_viridis(discrete = TRUE)+
    scale_color_manual(values=c("#FFC602", "#3CAEA3", "#ED553B"))+
    theme_bw()
  
  
  ## ## RESULT
  ## All number of animals in each time points in each group
  # take data from ggplot
  total_animal <- PI_mean[["data"]]
  #sum based on group and time
  total_animal <- total_animal %>% 
    group_by(group,time) %>% 
    dplyr :: summarise(S=sum(S, na.rm = TRUE),PI=sum(PI, na.rm = TRUE),
                       TI=sum(PIL, na.rm = TRUE),R=sum(R, na.rm = TRUE),
                       V=sum(V, na.rm=TRUE))
  
  #add the sum of all animals in each time
  total_animal$total <- rowSums(total_animal[,3:7])
  
  ## Adding vaccine and import information
  total_animal$import<- ifelse(total_animal$group=="GroupA", allowed_import_A, 1+((1-allowed_import_A)/2))
  total_animal$vaccine <- ifelse(total_animal$group =="GroupA", "0", animal_vaccine_B)
  fin_result <- total_animal
  
}


###################################Entering variables########################
sc_name <- "sc0"

## NOTE
## For baseline scenario: sc0
# allowed_import is 1
# animal_vaccine_A and animal_vaccine_B are 0
# farm_vaccine_A and farm_vaccine_B are 0

## For import restriction: sc1
# allowed_import is 0.5, meaning 50% reduction
# animal_vaccine_A and animal_vaccine_B are 0
# farm_vaccine_A and farm_vaccine_B are 0

## For targeted vaccination: sc2
# allowed_import is 1
# animal_vaccine_A is 0 and animal_vaccine_B is 0.8, meaning 80% animal get vaccinated
# farm_vaccine_A is 0 and farm_vaccine_B is 1, meaning all farms get vaccinated. 

## For combined: sc3
# allowed_import is 0.5
# animal_vaccine_A is 0 and animal_vaccine_B is 0.8
# farm_vaccine_A is 0 and farm_vaccine_B is 1

## ##set the proportion of import that you allowed to enter group A
## 0 is when there is total import restrictions
## 1 is when there is no import restriction to group A
allowed_import_A <- 0.5

## ## Enter animal vaccination coverage
## 0 is when there is no vaccination at all
## 1 is when all animals get vaccinated
animal_vaccine_A <- 0
animal_vaccine_B <- 0.8

## ## Enter farm vaccination coverage
## 0 is when there is no vaccination at all
## 1 is when all farms get vaccinated
farm_vaccine_A <- 0
farm_vaccine_B <- 1

#number of iteration
iteration <- 100

end_result<- vector("list", iteration)
#run the function
for (i in 1:iteration) {
  end_result[[i]] <- BVD_model(allowed_import_A, animal_vaccine_A,animal_vaccine_B,
                               farm_vaccine_A, farm_vaccine_B)
}

end_result <- do.call("rbind", end_result)


write.csv(end_result, paste0(sc_name,".csv"),
          row.names = FALSE)

#################################ANALYSIS##################################
## ## READ ALL RESULTS -----

## baseline scenario
sc0 <- read.table("sc0.csv", header = TRUE,
                  sep = ",", dec = ".")

## sc 1
sc1 <- read.table("sc1.csv", header = TRUE,
                  sep = ",", dec = ".")

## sc 2
sc2 <- read.table("sc2.csv", header = TRUE,
                  sep = ",", dec = ".")

## sc 3
sc3 <- read.table("sc3.csv", header = TRUE,
                  sep = ",", dec = ".")

## add scenario information 
sc0$scenario <- 0
sc0 <- na.omit(sc0)

sc1$scenario <- 1
sc1 <- na.omit(sc1)

sc2$scenario <- 2
sc2<- na.omit(sc2)

sc3$scenario <- 3
sc3 <- na.omit(sc3)

## bind all data
data <- rbind(sc0,sc1,sc2,sc3)

## level the scenario
data$scenario <- as.factor(data$scenario)

## ## RESULT OVER TIME-----
# we want to create graph to see the dynamic of disease in PI and TI population 

## GROUP A
# Analysis in group A with average of all nodes
# taking only group A in the data
data_A <- data %>% filter(group == "GroupA")

# over time analysis in group A with max and min value
# based on group

## baseline scenario: sc0
data_A_0 <- data_A %>% filter(scenario==0) %>% select(2:7)

# mean, max, min of rows on each iteration 
data_A_0_mean <- data_A_0 %>% group_by(time) %>% 
  dplyr::summarize(S=mean(S), PI=mean(PI),
                   TI=mean(TI), R=mean(R),
                   V= mean(V))
data_A_0_max <- data_A_0 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_A_0_min <- data_A_0 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, min
data_A_0_mean <- melt(data_A_0_mean, id.vars = 'time', variable_name = 'series')
data_A_0_max <- melt(data_A_0_max, id.vars = 'time', variable_name = 'series')
data_A_0_min <- melt(data_A_0_min, id.vars = 'time', variable_name = 'series')

# bind into one data frame and give unique name
data_A_0_p <- cbind(data_A_0_mean,data_A_0_max,data_A_0_min)
colnames(data_A_0_p)<- make.unique(colnames(data_A_0_p),sep="_")

#taking only PI 
data_A_0_p <- data_A_0_p %>% filter(variable == "PI")

# Ceate the graph
Sc0_A_fig <-ggplot(data_A_0_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Baseline scenario")+
  coord_cartesian(ylim = c(0,1000))

Sc0_A_fig

## Import restriction: sc1
# taking only sc1
data_A_1 <- data_A %>% filter(scenario==1) %>% select(2:7)
# creating df with value of: average, max, min
data_A_1_mean <- data_A_1 %>% group_by(time) %>%
  dplyr::summarize(S=median(S), PI=median(PI),
                   TI=median(TI), R=median(R),
                   V= median(V))
data_A_1_max <- data_A_1 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_A_1_min <- data_A_1 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, and min 
data_A_1_mean <- melt(data_A_1_mean, id.vars = 'time', variable_name = 'series')
data_A_1_max <- melt(data_A_1_max, id.vars = 'time', variable_name = 'series')
data_A_1_min <- melt(data_A_1_min, id.vars = 'time', variable_name = 'series')

#bind into one df and give unique name
data_A_1_p <- cbind(data_A_1_mean, data_A_1_max, data_A_1_min)
colnames(data_A_1_p) <- make.unique(colnames(data_A_1_p), sep = "_") 

#taking only PI
data_A_1_p <- data_A_1_p %>% filter(variable=="PI")

# creating graph
sc1_A_fig <- ggplot(data_A_1_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Scenario 1")+
  coord_cartesian(ylim = c(0,1000))

sc1_A_fig

## Targeted vaccination: sc2 
# taking only sc2
data_A_2 <- data_A %>% filter(scenario==2) %>% select(2:7)

# creating df with value of: average, max, min
data_A_2_mean <- data_A_2 %>% group_by(time) %>%
  dplyr::summarize(S=median(S), PI=median(PI),
                   TI=median(TI), R=median(R),
                   V= median(V))
data_A_2_max <- data_A_2 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_A_2_min <- data_A_2 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, and min 
data_A_2_mean <- melt(data_A_2_mean, id.vars = 'time', variable_name = 'series')
data_A_2_max <- melt(data_A_2_max, id.vars = 'time', variable_name = 'series')
data_A_2_min <- melt(data_A_2_min, id.vars = 'time', variable_name = 'series')

#bind into one df and give unique name
data_A_2_p <- cbind(data_A_2_mean, data_A_2_max, data_A_2_min)
colnames(data_A_2_p) <- make.unique(colnames(data_A_2_p), sep = "_") 

# taking only PI 
data_A_2_p <- data_A_2_p %>% filter(variable=="PI")

# creating graph
sc2_A_fig <- ggplot(data_A_2_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Scenario 2")+
  coord_cartesian(ylim = c(0,1000))

sc2_A_fig

## combined: sc3
# taking only sc3
data_A_3 <- data_A %>% filter(scenario==3) %>% select(2:7)

# creating df with value of: average, max, min
data_A_3_mean <- data_A_3 %>% group_by(time) %>%
  dplyr::summarize(S=median(S), PI=median(PI),
                   TI=median(TI), R=median(R),
                   V= median(V))
data_A_3_max <- data_A_3 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_A_3_min <- data_A_3 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, and min 
data_A_3_mean <- melt(data_A_3_mean, id.vars = 'time', variable_name = 'series')
data_A_3_max <- melt(data_A_3_max, id.vars = 'time', variable_name = 'series')
data_A_3_min <- melt(data_A_3_min, id.vars = 'time', variable_name = 'series')

#bind into one df and give unique name
data_A_3_p <- cbind(data_A_3_mean, data_A_3_max, data_A_3_min)
colnames(data_A_3_p) <- make.unique(colnames(data_A_3_p), sep = "_") 

# taking only PI 
data_A_3_p <- data_A_3_p %>% filter(variable=="PI")

# creating graph
sc3_A_fig <- ggplot(data_A_3_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Scenario 3")+
  coord_cartesian(ylim = c(0,1000))

sc3_A_fig

## ## Merge all the graph
A_figures <- ggarrange(Sc0_A_fig, sc1_A_fig, sc2_A_fig, sc3_A_fig,
                       nrow = 2, ncol = 2, common.legend = TRUE,
                       labels = c("a","b","c","d"))

A_figures


## ## creating one graph for all scenario
## group A

#taking only PI sc 0
data_A_0_p <- data_A_0_p %>% filter(variable == "PI")
data_A_0_p$scenario <- "0"

#taking only PI sc 1
data_A_1_p <- data_A_1_p %>% filter(variable=="PI")
data_A_1_p$scenario <- "1"

# taking only PI sc 2
data_A_2_p <- data_A_2_p %>% filter(variable=="PI")
data_A_2_p$scenario <- "2"

# taking only PI sc 3
data_A_3_p <- data_A_3_p %>% filter(variable=="PI")
data_A_3_p$scenario <- "3"

data_A_scenario <- rbind(data_A_0_p, data_A_1_p,
                         data_A_2_p, data_A_3_p)


data_A_scenario_p <- ggplot(data_A_scenario, aes(time, value, group=scenario, fill=scenario))+
  geom_line(aes(linetype=scenario,color=scenario))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  geom_point(aes(shape=scenario))+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A")+
  coord_cartesian(ylim = c(0,250))
data_A_scenario_p


## GROUP B
# Analysis in group A with average of all nodes
# taking only group A in the data
data_B <- data %>% filter(group == "GroupB")

## Baseline scenario:sc0 
data_B_0 <- data_B %>% filter(scenario==0) %>% select(2:7)

# mean, max, min of rows on each iteration 
data_B_0_mean <- data_B_0 %>% group_by(time) %>% 
  dplyr::summarize(S=mean(S), PI=mean(PI),
                   TI=mean(TI), R=mean(R),
                   V= mean(V))
data_B_0_max <- data_B_0 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_B_0_min <- data_B_0 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, min, max
data_B_0_mean <- melt(data_B_0_mean, id.vars = 'time', variable_name = 'series')
data_B_0_max <- melt(data_B_0_max, id.vars = 'time', variable_name = 'series')
data_B_0_min <- melt(data_B_0_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_B_0_p <- cbind(data_B_0_mean,data_B_0_max,data_B_0_min)
colnames(data_B_0_p)<- make.unique(colnames(data_B_0_p),sep="_")

# taking only PI population
data_B_0_p <- data_B_0_p %>% filter(variable=="PI")

# creating graph
Sc0_B_fig <-ggplot(data_B_0_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Baseline scenario B")+
  coord_cartesian(ylim = c(0,1000))


## import restriction: sc1
# taking only sc1
data_B_1 <- data_B %>% filter(scenario==1) %>% select(2:7)
# creating df with value of: average, max, min
data_B_1_mean <- data_B_1 %>% group_by(time) %>%
  dplyr::summarize(S=median(S), PI=median(PI),
                   TI=median(TI), R=median(R),
                   V= median(V))
data_B_1_max <- data_B_1 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_B_1_min <- data_B_1 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, and min 
data_B_1_mean <- melt(data_B_1_mean, id.vars = 'time', variable_name = 'series')
data_B_1_max <- melt(data_B_1_max, id.vars = 'time', variable_name = 'series')
data_B_1_min <- melt(data_B_1_min, id.vars = 'time', variable_name = 'series')

#bind into one df and give unique name
data_B_1_p <- cbind(data_B_1_mean, data_B_1_max, data_B_1_min)
colnames(data_B_1_p) <- make.unique(colnames(data_B_1_p), sep = "_") 

# taking only PI population
data_B_1_p <- data_B_1_p %>% filter(variable=="PI")

# creating graph
sc1_B_fig <- ggplot(data_B_1_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Scenario 1 B")+
  coord_cartesian(ylim = c(0,1000))


## targeted vaccination: sc2
# taking only sc2
data_B_2 <- data_B %>% filter(scenario==2) %>% select(2:7)

# creating df with value of: average, max, min
data_B_2_mean <- data_B_2 %>% group_by(time) %>%
  dplyr::summarize(S=median(S), PI=median(PI),
                   TI=median(TI), R=median(R),
                   V= median(V))
data_B_2_max <- data_B_2 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_B_2_min <- data_B_2 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, and min 
data_B_2_mean <- melt(data_B_2_mean, id.vars = 'time', variable_name = 'series')
data_B_2_max <- melt(data_B_2_max, id.vars = 'time', variable_name = 'series')
data_B_2_min <- melt(data_B_2_min, id.vars = 'time', variable_name = 'series')

#bind into one df and give unique name
data_B_2_p <- cbind(data_B_2_mean, data_B_2_max, data_B_2_min)
colnames(data_B_2_p) <- make.unique(colnames(data_B_2_p), sep = "_") 

# taking only PI population
data_B_2_p <- data_B_2_p %>% filter(variable=="PI")

# creating graph
sc2_B_fig <- ggplot(data_B_2_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Scenario 2 B")+
  coord_cartesian(ylim = c(0,1000))


## combined approach: sc3
# taking only sc3
data_B_3 <- data_B %>% filter(scenario==3) %>% select(2:7)

# creating df with value of: average, max, min
data_B_3_mean <- data_B_3 %>% group_by(time) %>%
  dplyr::summarize(S=median(S), PI=median(PI),
                   TI=median(TI), R=median(R),
                   V= median(V))
data_B_3_max <- data_B_3 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_B_3_min <- data_B_3 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))

# melt according to: mean, max, and min 
data_B_3_mean <- melt(data_B_3_mean, id.vars = 'time', variable_name = 'series')
data_B_3_max <- melt(data_B_3_max, id.vars = 'time', variable_name = 'series')
data_B_3_min <- melt(data_B_3_min, id.vars = 'time', variable_name = 'series')

#bind into one df and give unique name
data_B_3_p <- cbind(data_B_3_mean, data_B_3_max, data_B_3_min)
colnames(data_B_3_p) <- make.unique(colnames(data_B_3_p), sep = "_") 

# Taking only PI population
data_B_3_p <- data_B_3_p %>% filter(variable=="PI")

# creating graph
sc3_B_fig <- ggplot(data_B_3_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Scenario 3 B")+
  coord_cartesian(ylim = c(0,1000))



## ## merge the graph
B_figures <- ggarrange(Sc0_B_fig, sc1_B_fig, sc2_B_fig, sc3_B_fig,
                       nrow = 2, ncol = 2, common.legend = TRUE,
                       labels = c("a","b","c","d"))
B_figures


data_B_scenario <- rbind(data_B_0_p, data_B_1_p,
                         data_B_2_p, data_B_3_p)


## ## combine two graph (A and B)  in one window
PI_all_sc <- ggarrange(data_A_scenario_p, data_B_scenario_p,
                       nrow = 2, ncol = 1, common.legend = TRUE,
                       legend = "right",
                       labels = c("A", "B")) 

PI_all_sc

## ## PREVALENCE ----
## we want to see the prevalence at the end of simulation

# take only the last time step or the last month
data_end <- data %>% filter(time==120)

## TI prevalence
data_end$TI_prev <- data_end$TI/data_end$total
# plot the table
ggplot(data_end, aes(scenario,TI_prev, fill=group))+
  geom_boxplot()+
  ggtitle("at t = 120: \nTI prevalence")+
  theme_minimal()+
  ylab("%")+
  xlab("Scenario")


## PI prevalence
data_end$PI_prev <- (data_end$PI/data_end$total)*100
# plot the table
ggplot(data_end, aes(scenario,PI_prev, fill=group))+
  geom_boxplot()+
  ggtitle("at t = 120: \nPI prevalence")+
  theme_minimal()+ 
  ylab("%")
xlab("Scenario")

