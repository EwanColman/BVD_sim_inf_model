# MODEL DESCRIPTION
# this is the model with adjusted parameter based on
# Anonymized_data
# Changed parameters: movement probabilities, birth and mortality rate, 
# and rate of newborn PI

# In this model vaccination is scheduled once a year, at the end of the year
# or at month 12

# Number of movements and farms is extracted from Anonymized_data
# from 2008-2017.
# The time step in model is monthly time step and run for 10 years or
# 120 months, mimicing the movement condition from original data

# Prevalence of PI animals is taken from literature, which is 2% of animals
# are PIs animals. This might create over estimation of introduced PIs to 
# Scotland. 

########LIBRARIES#########
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

#############################VACCINATION ONCE A YEAR############################
#################################Build function#################################
BVD_model <- function(allowed_import_A, animal_vaccine_A,animal_vaccine_B,
                      farm_vaccine_A, farm_vaccine_B){
  ###############################1.SETTING PARAMETERS############################
  
  #TIME
  #Moves over the whole period of time
  #total movement 128928
  number_of_moves <- 1458258
  #number of time steps in month
  end <- 10*12
  tspan <- seq(from = 1, to= end)
  
  #NUMBER OF FARMS IN EACH GROUP
  # how many farms you want of each type? The first part initialises a vector
  number_of_farms <- c(0,0,0)
  # give the name to the group
  names(number_of_farms) <- c('A','B','C')
  # choose how many farms you want of each type
  #A 3891; B 1360; C 743
  number_of_farms['A'] <- 32129
  number_of_farms['B'] <- 5918
  number_of_farms['C'] <- 10
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
  # Mu    : PI removal rate
  # Delta : TI become PI rate. Calculated from reduced fertility.
  # gamma : rate from recovered to susceptible 
  # Reference of rate from paper :
  # "Farm productive contexts and the dynamics of bovine viral diarrhea (BVD) transmission"
  # W     : Vaccine efficacy
  # u     : percentage individual vaccinated/unit time 
  # b     : birth rate
  # m     : mortality rate not related to BVD
  # a    : rate of TI animal giving birth to healthy cow
  gdata <- c(beta = 0.03,
             delta = 0.00787921,
             Mu = 0.0132,
             gamma = 0.0833,
             u = 0,
             w = w,
             k = 0.03*(1-0.8),
             b = 0.02138261,
             m = 0.02138261,
             s = 100,
             a = 0.003527)
  
  #MOVEMENT PROPORTION MATRIX
  # create matrix 3 by 3
  move_prob <- c(rep(0,9))
  dim(move_prob) <- c(3,3)
  # then we add names so move_prob['A','B'], for example, can be called 
  rownames(move_prob) <- c('A','B','C')
  colnames(move_prob) <- c('A','B','C')
  # set the proportion of movement between the different types
  # These need to add to 1
  move_prob['C','A'] <- 0.03*(allowed_import_A)
  
  move_prob['A','A'] <- 0.46
  move_prob['A','B'] <- 0.04-(1-allowed_import_A)*0.03
  move_prob['A','C'] <- 0.02+(1-allowed_import_A)*0.03
  move_prob['B','A'] <- 0.07+(1-allowed_import_A)*0.03
  move_prob['B','B'] <- 0.27
  move_prob['B','C'] <- 0.07-(1-allowed_import_A)*0.03
  
  move_prob['C','B'] <- 0.04+(1-allowed_import_A)*0.03
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
  
  
  # create empty matrix to accomodate number of moves from each group
  # generated from the probability
  moves <- matrix(nrow = 3, ncol = 3)
  # then we add names so move_prob['A','B'], for example, can be called 
  rownames(moves) <- c('A','B','C')
  colnames(moves) <- c('A','B','C')
  
  # calculate the number of moves from each group
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
  
  # loop to randomly select nodes origin and destination
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
  #group B
  if(animal_vaccine_B > 0 ){
    vaccination_B <- data.frame(event="intTrans", 
                                time=rep(seq(from=12, to=end, 12),ceiling(u_group_B*number_of_farms['B'])),
                                node=sample(number_of_farms['A']+1:number_of_farms['B'],ceiling(u_group_B*number_of_farms['B']),replace=FALSE),
                                dest=0,n=0,
                                proportion=as.numeric(uw['B']),select=4,shift=1)
  }
  # group A
  if(animal_vaccine_A > 0){
    vaccination_A <- data.frame(event="intTrans",
                                time = rep(seq(from=12, to=end, 12),ceiling(u_group_A*number_of_farms['A'])),
                                node = sample(1:number_of_farms['A'],ceiling(u_group_A*number_of_farms['A']),replace=FALSE),
                                dest=0, n=0,
                                proportion= as.numeric(uw['A']), select=4,shift=1)
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
         rep(animal_B_farms,number_of_farms['B']),
         rep(animal_C_farms,number_of_farms['C']))
  
  #final data.frame of animals in each node
  u0 <- data.frame(S = S,
                   PI = rep(0,total_farms),
                   PIL = rep(0, total_farms),
                   R = rep(0, total_farms),
                   V = rep(0, total_farms))
  
  #NETWORK FLOW
  transition <- c("S -> PI > 0 ? beta*S*PI/(S+PIL+PI+R+V) : 0 -> PIL",
                  "@ -> b*S -> S",
                  "S -> m*S -> @",
                  "PIL -> beta*PIL -> R",
                  "@ -> delta*PIL -> PI",
                  "PI -> Mu*PI -> @",
                  "R -> gamma*R -> S",
                  "@ -> a*PIL -> S",
                  "@ -> S < 10 ? S+s : 0 -> S")
  
  #RUNNING THE MODEL
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
  
  #################################4. RESULT######################################
  
  #Run trajectory
  infected <- trajectory(result)
  
  #NODES GROUPING
  Scotland_farms <- as.numeric(Scotland_farms)
  #rename nodes with group[A,B,C]
  infected$group[infected$node <= number_of_farms['A']] <- "GroupA"
  infected$group[infected$node > number_of_farms['A'] & infected$node <= Scotland_farms] <- "GroupB"
  #infected$group[infected$node > Scotland_farms] <- "GroupC"
  
  #GRAPH AVERAGE FOR ALL GROUPS
  #mean PI
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
###################################Entering variables########################
#set the proportion of import that you allowed to enter group A
allowed_import_A <- 1

#Enter animal vaccination coverage
animal_vaccine_A <- 0
animal_vaccine_B <- 0.8

#Enter farm vaccination coverage
farm_vaccine_A <- 0
farm_vaccine_B <- 1

#number of iteration
iteration <- 100

end_result<- vector("list", iteration)
#run the function
for (i in 1:iteration) {
  end_result[[i]] <- BVD_model(allowed_import_A, animal_vaccine_A, animal_vaccine_B,
                               farm_vaccine_A, farm_vaccine_B)
}
end_result <- do.call("rbind", end_result)
setwd("D:/WORK/IDOH - EMJMD/4 Sem - Internship/Data/Result/Scenarios/5/5_anon_data")
write.csv(end_result, paste0("farm_vaccine_A",farm_vaccine_A,"farm_vaccine_B",farm_vaccine_B,"2.csv"),
          row.names = FALSE)

#################################ANALYSIS#######################################
# read all result
setwd("D:/WORK/IDOH - EMJMD/4 Sem - Internship/Data/Result/Scenarios/5/5_anon_data")
sc0 <- read.table("3_farm_vaccine_A0farm_vaccine_B0.csv", header = TRUE,
                  sep = ",", dec = ".")
sc1 <- read.table("farm_vaccine_A0farm_vaccine_B02.csv", header = TRUE,
                  sep = ",", dec = ".")
sc2 <- read.table("farm_vaccine_A0farm_vaccine_B12.csv", header = TRUE,
                  sep = ",", dec = ".")
sc3 <- read.table("farm_vaccine_A0farm_vaccine_B12.csv", header = TRUE,
                  sep = ",", dec = ".")
# add scenario information 
sc0$scenario <- 0
sc0 <- na.omit(sc0)

sc1$scenario <- 1
sc1 <- na.omit(sc1)

sc2$scenario <- 2
sc2<- na.omit(sc2)

sc3$scenario <- 3
sc3 <- na.omit(sc3)

# bind all data
data <- rbind(sc0,sc1,sc2,sc3)

# level the scenario
data$scenario <- as.factor(data$scenario)
################################Checking the distribution#######################
# taking several times (t)
# group A
data_end_A <- data %>% filter(time==120) %>% filter(group=="GroupA")
ggplot(data_end_A, aes(PI))+
  geom_histogram(aes(color=scenario, fill=scenario),
                 position = "identity", bins=30, alpha= 0.4)
# result evenly distributed but not perfectly bellshape

# t = x
data %>% 
  filter(time==120) %>% filter(group=="GroupA") %>%
  ggplot(aes(PI))+
  geom_histogram(aes(color=scenario, fill=scenario),
                 position = "identity", bins=30, alpha=0.4)
# qq plot 
data %>% 
  filter(time==120) %>% filter(group=="GroupA") %>%
  ggplot(aes(sample=PI))+
  geom_qq()+geom_qq_line()+
  facet_wrap(~scenario, scales = "free_y")

#group B
data_end_B <- data %>% filter(time==120) %>% filter(group=="GroupB")
ggplot(data_end_B, aes(PI))+
  geom_histogram(aes(color=scenario, fill=scenario),
                 position = "identity", bins=30, alpha= 0.4)
# t = x
data %>% 
  filter(time==60) %>% filter(group=="GroupB") %>%
  ggplot(aes(PI))+
  geom_histogram(aes(color=scenario, fill=scenario),
                 position = "identity", bins=30, alpha=0.4)
# qq plot 
data %>% 
  filter(time==120) %>% filter(group=="GroupA") %>%
  ggplot(aes(sample=PI))+
  geom_qq()+geom_qq_line()+
  facet_wrap(~scenario, scales = "free_y")

# conclusion:
# over all the distribution is normal (bell shaped)
# however, there are certain time point that showed skewed distribution
#########################At TIME 120############################################
# take only the last time step or the last month
data_end <- data %>% filter(time==120)
data_end <- ifelse(data_end$scenario==2, 3, ifelse())
# plot the table
ggplot(data_end, aes(scenario,PI_prev, fill=group))+
  geom_boxplot()+
  ggtitle("at t = 120: \nnumber of PI")+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Scenario")+
  scale_x_discrete(limits=c("0", "1", "3", "2"),
                   labels=c("3"="2", "2"="3"))


########################OVER TIME POPULATION#################################
# analysis in group A with average of all nodes
# taking only group A in the data
data_A <- data %>% filter(group == "GroupA")
# plot group A over time
ggplot(data_A_mean,
       aes(time, TI, color=scenario))+
  geom_line()

data_A_mean <- data_A %>% group_by(time, scenario) %>%
  dplyr :: summarise(S=mean(S, na.rm = TRUE),PI=mean(PI, na.rm = TRUE),
                     TI=mean(TI, na.rm = TRUE),R=mean(R, na.rm = TRUE),
                     V=mean(V, na.rm=TRUE), total=mean(total,na.rm=TRUE))
data_A_mean_1 <- data_A_mean %>% filter(scenario=="1")
data_A_mean_1<- subset(data_A_mean_1, select = -c(scenario, total, S))
data_A_mean_1_melt <- melt(data_A_mean_1, id.vars = 'time', variable.name = 'series')
ggplot(data_A_mean_1_melt, aes(time,value))+
  geom_line(aes(color = series))+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Over time Group A \nScenario 1")+
  coord_cartesian(ylim = c(0,850))+
  xlim(0,120)

data_A_mean_2 <- data_A_mean %>% filter(scenario=="2")
data_A_mean_2 <- subset(data_A_mean_2, select = -c(scenario,total,S))
data_A_mean_2_melt <- melt(data_A_mean_2, id.vars = 'time', variable.name='series')
ggplot(data_A_mean_2_melt, aes(time,value))+
  geom_line(aes(color = series))+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Over time Group A \nScenario 2")+
  coord_cartesian(ylim = c(0,850))

data_A_mean_0 <- data_A_mean %>% filter(scenario=="0")
data_A_mean_0 <- subset(data_A_mean_0, select = -c(scenario,total,S))
data_A_mean_0_melt <- melt(data_A_mean_0, id.vars = 'time', variable.name='series')
ggplot(data_A_mean_0_melt, aes(time,value))+
  geom_line(aes(color = series))+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Overtime Group A \nControl")+
  coord_cartesian(ylim = c(0,850))


# overtime analysis in group A with max and min value
# based on group
# control
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
data_A_0_mean <- melt(data_A_0_mean, id.vars = 'time', variable_name = 'series')
data_A_0_max <- melt(data_A_0_max, id.vars = 'time', variable_name = 'series')
data_A_0_min <- melt(data_A_0_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_A_0_p <- cbind(data_A_0_mean,data_A_0_max,data_A_0_min)
colnames(data_A_0_p)<- make.unique(colnames(data_A_0_p),sep="_")

# creating graph
ggplot(data_A_0_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A \nBaseline scenario")+
  coord_cartesian(ylim = c(0,1000))
 data_A_0_p <- data_A_0_p %>% filter(!variable=="S")
#taking only PI 
data_A_0_p <- data_A_0_p %>% filter(variable == "PI")

# scenario 1
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

# creating graph
ggplot(data_A_1_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A \nScenario 1")+
  coord_cartesian(ylim = c(0,1000))
data_A_1_p <- data_A_1_p %>% filter(!variable=="S")
#taking only PI
data_A_1_p <- data_A_1_p %>% filter(variable=="PI")

# scenario 2
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

# creating graph
ggplot(data_A_2_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A \nScenario 3")+
  coord_cartesian(ylim = c(0,1000))
data_A_2_p <- data_A_2_p %>% filter(!variable=="S")
# taking only PI 
data_A_2_p <- data_A_2_p %>% filter(variable=="PI")

# scenario 3
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

# creating graph
ggplot(data_A_3_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A \nScenario 2")+
  coord_cartesian(ylim = c(0,1000))
data_A_3_p <- data_A_3_p %>% filter(!variable=="S")
# taking only PI 
data_A_3_p <- data_A_3_p %>% filter(variable=="PI")


#GROUP B
# control 
# taking only group A in the data
data_B <- data %>% filter(group == "GroupB")

# overtime analysis in group B with max and min value
# based on group
# control
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
data_B_0_mean <- melt(data_B_0_mean, id.vars = 'time', variable_name = 'series')
data_B_0_max <- melt(data_B_0_max, id.vars = 'time', variable_name = 'series')
data_B_0_min <- melt(data_B_0_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_B_0_p <- cbind(data_B_0_mean,data_B_0_max,data_B_0_min)
colnames(data_B_0_p)<- make.unique(colnames(data_B_0_p),sep="_")

# creating graph
ggplot(data_B_0_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group B \nBaseline scenario")+
  coord_cartesian(ylim = c(0,1000))
data_B_0_p <- data_B_0_p %>% filter(!variable=="S")
# taking only PI
data_B_0_p <- data_B_0_p %>% filter(variable =="PI")

# scenario 1
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

# creating graph
ggplot(data_B_1_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group B \nScenario 1")+
  coord_cartesian(ylim = c(0,1000))
data_B_1_p <- data_B_1_p %>% filter(!variable=="S")
#taking only PI
data_B_1_p <- data_B_1_p %>% filter(variable=="PI")

# scenario 2
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

# creating graph
ggplot(data_B_2_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group B \nScenario 3")+
  coord_cartesian(ylim = c(0,1000))
data_B_2_p <- data_B_2_p %>% filter(!variable=="S")
#taking only PI
data_B_2_p <- data_B_2_p %>% filter(variable=="PI")

# scenario 3
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

# creating graph
ggplot(data_B_3_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group B \nScenario 2")+
  coord_cartesian(ylim = c(0,1000))
data_B_3_p <- data_B_3_p %>% filter(!variable=="S")
#taking only PI
data_B_2_p <- data_B_2_p %>% filter(variable=="PI")

################################Statistic over time ############################
# this statistic aim to compare the mean of PI animal
# in each grout at each time step
# check distribution
hist(data_B_0$PI[data_B_0$time==120])
#well seems okay

# t test for sc0
t.test(data_A_0$PI[data_A_0$time==120], 
       data_B_0$PI[data_B_0$time==120],
       alternative = "two.sided",
       var.equal = TRUE)

# relative difference of sc1 A and B
# at time t=12
# group A
Rel_dif_A_1 <-(data_A_0_mean$value[data_A_0_mean$variable=="PI"]-
                 data_A_1_mean$value[data_A_1_mean$variable=="PI"])/data_A_0_mean$value[data_A_0_mean$variable=="PI"]
Rel_dif_A_1[120]*100

# group B
Rel_dif_B_1 <-(data_B_0_mean$value[data_B_0_mean$variable=="PI"]-
                 data_B_1_mean$value[data_B_1_mean$variable=="PI"])/data_B_0_mean$value[data_B_0_mean$variable=="PI"]
Rel_dif_B_1[120]*100


# t test for group A sc1 vs sc0
t.test(data_A_0$PI[data_A_0$time==120], 
       data_A_3$PI[data_A_3$time==120],
       alternative = "two.sided",
       var.equal = TRUE)

# t test for group A sc3 vs sc1
t.test(data_A_1$PI[data_A_1$time==120], 
       data_A_2$PI[data_A_2$time==120],
       alternative = "two.sided",
       var.equal = TRUE)
#################PI Prevalence######################################
data$PI_prev <- data$PI/data$total
setwd("D:/WORK/IDOH - EMJMD/4 Sem - Internship/Data/Result/Scenarios/Result/Result_anon_data")
write.csv(data,"sc0-3_with_PI_Prev.csv", row.names = FALSE, quote = FALSE)

#looking at the highest PI prev in group A in sc0
data %>% filter(group=="GroupA") %>% 
  filter(!scenario==0) %>% filter(time==120) %>% summarize(max(PI_prev))
#looking at the highest PI prev in group B in sc0
data %>% filter(group=="GroupB") %>% 
  filter(!scenario==0) %>% filter(time==120) %>% summarize(min(PI_prev))

#statistic PI analyasis (ANOVA) between scenarios in group A
data_A_PI <- data %>% filter(group=="GroupA")
write.csv(data_A_PI, "sc0-2_with_PI_A_Prev.csv", row.names = FALSE, quote = FALSE)
one_way_PI_A <- aov(PI_prev~scenario, data=data_A_PI)
summary(one_way_PI_A)
kruskal.test(PI_prev~scenario, data=data_A_PI)
TukeyHSD(one_way_PI_A)
pairwise.wilcox.test(data_A_PI$PI_prev, data_A_PI$scenario,
                     p.adjust.method = "BH")
PI_A_lm<-lm(PI_prev~scenario,data = data_A_PI)
summary(PI_A_lm)

data_12 <- rbind(sc0,sc1,sc2)
data_12$PI_prev <- data_12$PI/data_12$total
data_12$scenario<- as.factor(data_12$scenario)
data_12_A_t <- data_12 %>% filter(group=="GroupA") %>% 
  filter(time==120)
PI_A_lm <- aov(PI~scenario, data = data_12_A_t)
summary(PI_A_lm)
tukey_A_test<-TukeyHSD(PI_A_lm)
kruskal.test(PI~scenario, data=data_12_A_t)
plot(tukey_A_test)

# creating df with value of: average, max, min
data_A_PI_mean <- data_A_PI %>% group_by(time,scenario) %>%
  dplyr::summarize(PI_prev=median(PI_prev))
data_A_PI_max <- data_A_PI %>% group_by(time,scenario) %>% 
  dplyr::summarize(PI_prev_max=max(PI_prev))
data_A_PI_min <- data_A_PI %>% group_by(time,scenario) %>% 
  dplyr::summarize(PI_prev_min=min(PI_prev))

#bind into one df and give unique name
data_A_PI_p <- cbind(data_A_PI_mean, data_A_PI_max, data_A_PI_min)
colnames(data_A_PI_p) <- make.unique(colnames(data_A_PI_p), sep = "_") 
data_A_PI_p <- data_A_PI_p %>% dplyr::rename(scenario=scenario...2)

# creating graph for prevalence group A
ggplot(data_A_PI_p, aes(time...1, PI_prev, group=scenario, fill=scenario))+
  geom_line(aes(color=scenario))+
  geom_ribbon(aes(ymax=PI_prev_max, ymin=PI_prev_min), alpha=0.3)+
  theme_minimal()+
  ylab("%")+
  xlab("Time (month)")+
  ggtitle("PI prevalence Group A")+
  coord_cartesian(ylim = c(0,0.0016))+
  scale_color_discrete(limits=c("0", "1", "3", "2"),
                   labels=c("3"="2", "2"="3"))+
  scale_fill_discrete(limits=c("0", "1", "3", "2"),
                       labels=c("3"="2", "2"="3"))

data_A_PI_p %>% filter(time...1==120)

#statistic PI analysis between scenarios in group B
data_B_PI <- data %>% filter(group == "GroupB")
write.csv(data_B_PI, "sc0-2_with_PI_B_Prev.csv", row.names = FALSE, quote = FALSE)
one_way_PI_B <- aov(PI_prev~scenario, data=data_B_PI)
summary(one_way_PI_B)
kruskal.test(PI_prev~scenario, data=data_B_PI)
TukeyHSD(one_way_PI_B)

data_B_PI %>% filter(scenario==2) %>% summarize(max(PI))

# creating df with value of: average, max, min
data_B_PI_mean <- data_B_PI %>% group_by(time,scenario) %>%
  dplyr::summarize(PI_prev=median(PI_prev))
data_B_PI_max <- data_B_PI %>% group_by(time,scenario) %>% 
  dplyr::summarize(PI_prev_max=max(PI_prev))
data_B_PI_min <- data_B_PI %>% group_by(time,scenario) %>% 
  dplyr::summarize(PI_prev_min=min(PI_prev))

#bind into one df and give unique name
data_B_PI_p <- cbind(data_B_PI_mean, data_B_PI_max, data_B_PI_min)
data_B_PI_p <- data_B_PI_p %>% dplyr::rename(scenario=scenario...2)

# creating graph for prevalence group B
ggplot(data_B_PI_p, aes(time...1, PI_prev, group=scenario, fill=scenario))+
  geom_line(aes(color=scenario))+
  geom_ribbon(aes(ymax=PI_prev_max, ymin=PI_prev_min), alpha=0.3)+
  theme_minimal()+
  ylab("%")+
  xlab("Time (month)")+
  ggtitle("PI prevalence Group B")+
  coord_cartesian(ylim = c(0,0.0016))+
  scale_color_discrete(limits=c("0", "1", "3", "2"),
                       labels=c("3"="2", "2"="3"))+
  scale_fill_discrete(limits=c("0", "1", "3", "2"),
                      labels=c("3"="2", "2"="3"))

data_B_PI_p %>% filter(time...1==120)

############################VACCINATION ANALYSIS################################
# read vaccination result
setwd("D:/WORK/IDOH - EMJMD/4 Sem - Internship/Data/Result/Scenarios/2/2_anon_data")
sc0 <- read.table("3_farm_vaccine_A0farm_vaccine_B0.csv", header = TRUE,
                  sep = ",", dec = ".")
sc4 <- read.table("farm_vaccine_A1farm_vaccine_B1.csv", header = TRUE,
                  sep = ",", dec = ".")
sc5 <- read.table("farm_vaccine_A0.8farm_vaccine_B0.8.csv", header = TRUE,
                  sep = ",", dec = ".")

# add scenario information 
sc4$scenario <- 1
sc4 <- na.omit(sc4)

sc5$scenario <- 2
sc5<- na.omit(sc5)

# bind all data
data_vac <- rbind(sc0,sc4,sc5)

# level the scenario
data_vac$scenario <- as.factor(data_vac$scenario)

# overtime analysis in group A with max and min value
# based on group
# scenario 3
data_vac_A_1 <- data_vac %>% filter(scenario==1) %>% 
  filter(group=="GroupA") %>% select(2:7)
# mean, max, min of rows on each iteration 
data_vac_A_1_mean <- data_vac_A_1 %>% group_by(time) %>% 
  dplyr::summarize(S=mean(S), PI=mean(PI),
                   TI=mean(TI), R=mean(R),
                   V= mean(V))
data_vac_A_1_max <- data_vac_A_1 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_vac_A_1_min <- data_vac_A_1 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))
data_vac_A_1_mean <- melt(data_vac_A_1_mean, id.vars = 'time', variable_name = 'series')
data_vac_A_1_max <- melt(data_vac_A_1_max, id.vars = 'time', variable_name = 'series')
data_vac_A_1_min <- melt(data_vac_A_1_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_vac_A_1_p <- cbind(data_vac_A_1_mean,data_vac_A_1_max,data_vac_A_1_min)
colnames(data_vac_A_1_p)<- make.unique(colnames(data_vac_A_1_p),sep="_")
# creating graph
ggplot(data_vac_A_1_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A \nvaccine analysis 1")+ 
  coord_cartesian(ylim = c(0,5e+06))
  coord_cartesian(ylim = c(0,1000))
 
data_vac_A_1_p <- data_vac_A_1_p %>% filter(!variable=="S")

# scenario 4
data_vac_A_2 <- data_vac %>% filter(scenario==2) %>% 
  filter(group=="GroupA") %>% select(2:7)
# mean, max, min of rows on each iteration 
data_vac_A_2_mean <- data_vac_A_2 %>% group_by(time) %>% 
  dplyr::summarize(S=mean(S), PI=mean(PI),
                   TI=mean(TI), R=mean(R),
                   V= mean(V))
data_vac_A_2_max <- data_vac_A_2 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_vac_A_2_min <- data_vac_A_2 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))
data_vac_A_2_mean <- melt(data_vac_A_2_mean, id.vars = 'time', variable_name = 'series')
data_vac_A_2_max <- melt(data_vac_A_2_max, id.vars = 'time', variable_name = 'series')
data_vac_A_2_min <- melt(data_vac_A_2_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_vac_A_2_p <- cbind(data_vac_A_2_mean,data_vac_A_2_max,data_vac_A_2_min)
colnames(data_vac_A_2_p)<- make.unique(colnames(data_vac_A_2_p),sep="_")
# creating graph
ggplot(data_vac_A_2_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Group A \nvaccine analysis 2")+
  coord_cartesian(ylim = c(0,5e+06))
  coord_cartesian(ylim = c(0,1000))
  
data_vac_A_2_p <- data_vac_A_2_p %>% filter(!variable=="S")

# Group B
# scenario 3
data_vac_B_3 <- data_vac %>% filter(scenario==3) %>% 
  filter(group=="GroupB") %>% select(2:7)
# mean, max, min of rows on each iteration 
data_vac_B_3_mean <- data_vac_B_3 %>% group_by(time) %>% 
  dplyr::summarize(S=mean(S), PI=mean(PI),
                   TI=mean(TI), R=mean(R),
                   V= mean(V))
data_vac_B_3_max <- data_vac_B_3 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_vac_B_3_min <- data_vac_B_3 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))
data_vac_B_3_mean <- melt(data_vac_B_3_mean, id.vars = 'time', variable_name = 'series')
data_vac_B_3_max <- melt(data_vac_B_3_max, id.vars = 'time', variable_name = 'series')
data_vac_B_3_min <- melt(data_vac_B_3_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_vac_B_3_p <- cbind(data_vac_B_3_mean,data_vac_B_3_max,data_vac_B_3_min)
colnames(data_vac_B_3_p)<- make.unique(colnames(data_vac_B_3_p),sep="_")
# creating graph
ggplot(data_vac_B_3_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Over time Group B \nvaccine analysis 1")+
  coord_cartesian(ylim = c(0,5e+06))
  coord_cartesian(ylim = c(0,1000))
data_vac_B_3_p <- data_vac_B_3_p %>% filter(!variable=="S")

# scenario 4
data_vac_B_4 <- data_vac %>% filter(scenario==4) %>% 
  filter(group=="GroupB") %>% select(2:7)
# mean, max, min of rows on each iteration 
data_vac_B_4_mean <- data_vac_B_4 %>% group_by(time) %>% 
  dplyr::summarize(S=mean(S), PI=mean(PI),
                   TI=mean(TI), R=mean(R),
                   V= mean(V))
data_vac_B_4_max <- data_vac_B_4 %>% group_by(time) %>% 
  dplyr::summarize(Smax=max(S), PImax=max(PI), TImax=max(TI),
                   Rmax=max(R), Vmax=max(V))
data_vac_B_4_min <- data_vac_B_4 %>% group_by(time) %>% 
  dplyr::summarize(Smin=min(S), PImin=min(PI), TImin=min(TI),
                   Rmin=min(R), Vmin=min(V))
data_vac_B_4_mean <- melt(data_vac_B_4_mean, id.vars = 'time', variable_name = 'series')
data_vac_B_4_max <- melt(data_vac_B_4_max, id.vars = 'time', variable_name = 'series')
data_vac_B_4_min <- melt(data_vac_B_4_min, id.vars = 'time', variable_name = 'series')

# bind into one df and give unique name
data_vac_B_4_p <- cbind(data_vac_B_4_mean,data_vac_B_4_max,data_vac_B_4_min)
colnames(data_vac_B_4_p)<- make.unique(colnames(data_vac_B_4_p),sep="_")
# creating graph
ggplot(data_vac_B_4_p, aes(time, value, group=variable, fill=variable))+
  geom_line(aes(color=variable))+
  geom_ribbon(aes(ymax=value_1, ymin=value_2), alpha=0.3)+
  theme_minimal()+
  ylab("n of cattle")+
  xlab("Time (month)")+
  ggtitle("Over time Group B \nvaccine analysis 2")+
  coord_cartesian(ylim = c(0,5e+06))
  coord_cartesian(ylim = c(0,1000))
data_vac_B_4_p <- data_vac_B_4_p %>% filter(!variable=="S")

####box plot at t=120
# take only the last time step or the last month
data_vacc_end <- data_vac %>% filter(time==120)
data_vacc_end$PI_prev <- data_vacc_end$PI/data_vacc_end$total

#rename the scenario
data_vacc_end$scenario[data_vacc_end$scenario=="3"] <- "1"
data_vacc_end$scenario[data_vacc_end$scenario=="4"] <- "2"
data_vacc_end$scenario <- as.factor(data_vacc_end$scenario)
# plot the table
ggplot(data_vacc_end, aes(scenario,PI_prev, fill=group))+
  geom_boxplot()+
  ggtitle("at t = 120: \nPI prevalence")+
  theme_minimal()+
  ylab("%")+
  xlab("Analysis")+
  scale_y_continuous(breaks = seq(0,0.0016, by=0.0004))+
  coord_cartesian(ylim = c(0, 0.0017))

################### vaccination statistic #################################

# group A
# take only group A
data_vacc_end_A <- data_vacc_end %>% filter(group=="GroupA")

# check histogram
hist(data_vacc_end_A$PI_prev[data_vacc_end_A$scenario==2])
# ok normal

# to see if prevalence of VA1 and VA2 in group A is different
# t test for comparing 

t.test(data_vacc_end_A$PI_prev[data_vacc_end_A$scenario==1],
       data_vacc_end_A$PI_prev[data_vacc_end_A$scenario==2],
       alternative = "two.sided",
       var.equal = TRUE)


##########################CHI sq for movement ###################################

# testing movement between group if they are independent 
origin <- c(rep("import",96651), rep("domestic", 1213539))
dest <- c(rep("A",32421), rep("B",64230), rep("A",775500), rep("B",438039))
data <- data.frame(origin,dest)
table(data$dest, data$origin)
chisq.test(data$dest, data$origin, correct=FALSE)
