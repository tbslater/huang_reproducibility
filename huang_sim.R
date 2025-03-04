###########################################################################
############### HUANG ET AL. - REPRODUCIBILITY ANALYSIS ###################
###########################################################################

## Load relevant packages ##
library(simmer)
library(simmer.plot)
library(parallel)
library(dplyr)
library(plotly)
library("gridExtra")

## Set seed for simulation ##
set.seed(14112024)

## Parameter names ##
paramNames <- c("ct", "angio_inr", "angio_ir", "stroke_staff", "ed_staff",
                "angio_staff", "ir", "inr", "angio_staff_night", "ir_night",
                "inr_night", "ed_pt", "st_pt", "ais_pt", "ecr_pt", "inr_pt", 
                "eir_pt", "ir_pt", "shifts", "nsim", "run_t")

###########################################################################
########################## SIMULATION FUNCTION ############################
###########################################################################

simulate_nav <- function(ct = 2, angio_inr = 1, angio_ir = 1, 
                         stroke_staff = 1, ed_staff = 10, angio_staff = 6,
                         ir = 2, inr = 1, angio_staff_night = 3, 
                         ir_night = 1, inr_night = 1, ed_pt = 107700, 
                         st_pt = 750, ais_pt = 450, ecr_pt = 58, 
                         inr_pt = 104, eir_pt= 468, ir_pt = 3805,  
                         shifts = c(8,17), nsim = 10, run_t = 365) {
  
  ########################################################################
  ######### A function to perform DES applied to ECR service ############# 
  ########################################################################
  ## Inputs: 
  ## ct:                [int] no. of CT scan machines (default = 2)
  ## angio_inr:         [int] no. of INR+IR angio. machines (default = 1)
  ## angio_ir:          [int] no. of IR angio. machines (default = 1)
  ## stroke_staff:      [int] no. of stroke  staff (default = 1)
  ## ed_staff:          [int] no. of ED staff (default = 10)
  ## angio_staff:       [int] no. of angio. staff (default = 10)
  ## ir:                [int] no. of IRs (default = 1)
  ## inr:               [int] no. of INRs (default = 1)
  ## angio_staff_night: [int] no. of angio. night staff (default = 3)
  ## ir_night:          [int] no. of night IRs (default = 1)
  ## inr_night:         [int] no. of night INRs (default = 1)
  ## ed_pt:             [int] no. of ED patients (default = 107,000)
  ## st_pt:             [int] no. of stroke patients (default = 750)
  ## ais_pt:            [int] no. of AIS patients (default = 450)
  ## ecr_pt:            [int] no. of ECR patients (default = 58)
  ## inr_pt:            [int] no. of elective INR patients (default = 300)
  ## eir_pt:            [int] no. of emergency IR patients (default = 1,000)
  ## ir_pt:             [int] no. of elective IR patients (default = 4,000)
  ## shifts:            [vec] shift start/end hours (default = c(8,17))
  ## nsim:              [int] no. of simulations to run (default = 1)
  ## run_t:             [int] run time of simulation (default = 10,000)
  #######################################################################
  ## Outputs:
  ##                    [lst] list of arrivals/resources monitoring stats
  #######################################################################
  
  ## MODEL PARAMETERS ##
  
  # Physical resources #
  CT = ct
  ANGIO_INR = angio_inr
  ANGIO_IR = angio_ir 
  
  # Human resources #
  ED_STAFF = ed_staff
  ST_STAFF = stroke_staff
  ANGIO_STAFF = angio_staff
  INR = inr
  IR = ir
  ANGIO_STAFF_NIGHT = angio_staff_night
  INR_NIGHT = inr_night
  IR_NIGHT = ir_night
  
  # Proportions of stroke pt #
  PROB_STROKE = st_pt / ed_pt 
  PROB_AIS = ais_pt / st_pt   
  PROB_ECR = ecr_pt / ais_pt 
  
  # Inter-arrival t #
  year2min = 525600 
  I_ED  = round(year2min/ed_pt)
  I_ST  = round(year2min/st_pt)
  I_AIS = round(year2min/ais_pt)
  I_ECR = round(year2min/ecr_pt)
  I_INR = round(year2min/inr_pt)
  I_EIR = round(year2min/eir_pt)
  I_IR  = round(year2min/ir_pt)
  
  # Day shift start and end # 
  T_START = shifts[1] * 60
  T_END = shifts[2] * 60
  
  # Simulation set-up #
  RUN_T = run_t * 60 * 24
  N_SIM = nsim
  
  #######################################################################
  
  ## MODEL ##
  
  ecr_traj <- trajectory("ecr trajectory") %>% # ECR patient trajectory
    seize("angio_inr", 1) %>% # angio_INR machine in use
    seize("inr", 1) %>% # INR in use
    seize("angio_staff", 3) %>% # 3 x angio_staff in use
    timeout(function() rnorm(1,120,60)) %>% # delay ~N(120,3600)
    release("angio_inr", 1) %>% # angio_INR machine available
    release("inr", 1) %>% # INR available
    release("angio_staff", 3) # 3 x angio_staff available
  
  ais_traj <- trajectory("LVO trajectory") %>% # LVO patient trajectory
    set_attribute("priority", 3) %>%
    set_prioritization(values = c(3, 3, FALSE)) %>%
    branch(option = function() sample(1:2, 1, prob = # sample from 2 paths
                                        c(PROB_ECR, (1-PROB_ECR))), 
           continue = c(TRUE,FALSE), # 1: cont. to main trajectory
           ecr_traj, # 1: ECR trajectory
           
           trajectory("tpa only") %>% # 2: TPA trajectory
             timeout(1) # delay 1 min
    )
  
  nonstroke_traj <- trajectory(name = "non stroke traj") %>% # non-stroke
    branch(option = function() sample(1:2, 1, # sample from 2 paths
                                      prob = c(.9, .1)), 
           continue = c(FALSE,FALSE), # neither cont. to main trajectory
           trajectory("discharge path") %>% # 1: discharge trajectory
             timeout(1) %>% # delay 1 min
             leave(prob=1), # leave system
           
           trajectory("ct review") %>% # CT trajectory
             seize("ct",1) %>% # CT scan machine in use
             timeout(20) %>% # 20 minute delay 
             release("ct",1) %>% # CT scan machine available
             leave(prob=1) # leave system
    )
  
  stroke_traj <- trajectory("stroke trajectory") %>% # stroke trajectory
    set_attribute("priority", 3) %>% # set priority attribute
    set_prioritization(values = c(3, 3, FALSE)) %>% # 3x priority, 3x not
    
    seize("stroke_doctor",1) %>% # stroke doc. in use
    timeout(function() rnorm(1, 30, 10)) %>% # delay ~N(30,100) mins
    release("stroke_doctor",1) %>% # stroke doc. available
    
    seize("ct",1) %>% # CT scan machine in use
    timeout(function() rnorm(1, 20,10)) %>% # delay ~N(20,100) mins
    release("ct",1) %>% # CT scan machine available
    
    branch(option = function() sample(1:2, 1, # sample from 2 paths
                                      prob = c(PROB_AIS, (1-PROB_AIS))), 
           continue = c(FALSE,TRUE), # 2: cont. to main trajectory
           ais_traj, # 1: AIS trajectory
           
           trajectory("not ais") %>% # 2: not AIS trajectory
             timeout(1) %>% # delay 1 min
             leave(prob=1) # leave system
    )
  
  new_patient_traj <- trajectory(name = "new patient's path") %>% # NP traj.
    seize("ed_staff", 1) %>% # 1 x ED staff in use
    timeout(function() rnorm(1, 20,10)) %>% # delay ~N(20,100) mins
    release("ed_staff", 1) %>% # 1 x ED staff available
    
    branch(option = function() sample(1:2, 1, # sample from 2 paths
                                      prob = c(PROB_STROKE, (1-PROB_STROKE))),
           continue = c(F,T), # 2:  cont. to main trajectory
           stroke_traj, # 1: stroke trajectory
           
           nonstroke_traj # 2: non-stroke trajectory
    )
  
  ir_traj <- trajectory("ir traj") %>% # IR trajectory
    seize("door", 1)  %>% # door in use
    release("door", 1) %>% # door available
    
    seize("angio_staff", 1) %>% # 1x angio_staff in use
    timeout(function() rnorm(1, 20,10)) %>% # delay ~N(20,100) mins
    release("angio_staff", 1) %>% # 1x angio_staff available
    # comment in/out the below 3 lines for exclusive scenario
    simmer::select(resources = c("angio_ir", "angio_inr"),
                   policy = "shortest-queue") %>% # choose quickest option
    seize_selected(amount = 1) %>% # IR/INR machine in use
    #seize("angio_ir", 1) %>% # IR machine in use
    seize("ir", 1) %>% # IR in use
    seize("angio_staff", 3) %>% # 3x angio_staff in use
    timeout(function() rnorm(1, 60,30)) %>% # delay ~N(60,900) mins
    # comment in/out the below 2 lines for exclusive scenario
    release_selected(amount=1) %>% # IR/INR machine available
    #release("angio_ir", 1) %>% # IR machine available
    release("ir", 1) %>% # IR available
    release("angio_staff",3) # 3x angio_staff available
  
  inr_traj <- trajectory("inr traj") %>% # INR trajectory
    seize("door", 1)  %>% # door in use
    release("door", 1) %>% # door available
    
    seize("angio_staff", 1) %>% # 1x angio_staff in use
    timeout(function() rnorm(1, 20,10)) %>% # delay ~N(20,100) mins
    release("angio_staff", 1) %>% # 1x angio_staff available
    
    seize("angio_inr", 1) %>% # angio_INR machine in use
    seize("inr", 1) %>% # INR in use
    seize("angio_staff", 3) %>% # 3x angio_staff in use
    timeout(function() rnorm(1, 60,30)) %>% # delay ~N(60,900) mins
    release("inr",1) %>% # INR available
    release("angio_staff",3) %>% # 3x angio_staff available
    release("angio_inr", 1) # angio_INR machine available
  
  eir_traj <- trajectory("eir traj") %>% # EIR trajectory
    seize("angio_staff", 1) %>% # 1x angio_staff in use
    timeout(function() rnorm(1, 20,10)) %>% # delay ~N(20,100) mins
    release("angio_staff", 1) %>% # 1x angio_staff available
    
    simmer::select(resources = c("angio_ir", "angio_inr"), 
                   policy = "shortest-queue") %>% # choose quickest option
    seize_selected(amount = 1) %>% # IR/INR machine in use
    seize("angio_staff", 3) %>% # 3x angio_staff in use
    seize("ir", 1) %>% # IR in use
    timeout(function() rnorm(1, 60,30)) %>% # delay ~N(60,900) mins
    release_selected(amount=1) %>% # IR/INR machine available
    release("ir",1) %>% # IR available
    release("angio_staff",3) # 3x angio_staff available
  
  #######################################################################
  
  ## MODEL SIMULATION ##
  
  # Progress bar #
  # progress <- shiny::Progress$new()
  # on.exit(progress$close())
  # progress$set(message = "Running simulation", value = 0)
  
  # Capacity schedule #
  DOOR_SCHEDULE <- schedule(c(T_START, T_END), 
                            c(Inf, 0), # never ends
                            period = 1440) # repeat every day 
  STAFF_SCHEDULE <- schedule(c(T_START, T_END), 
                             c(ANGIO_STAFF, ANGIO_STAFF_NIGHT), # shifts
                             period = 1440) # repeat every day
  IR_SCHEDULE <- schedule(c(T_START, T_END), 
                          c(IR, IR_NIGHT), # shift rotations
                          period = 1440) # repeat every day
  INR_SCHEDULE <- schedule(c(T_START, T_END), 
                           c(INR, INR_NIGHT), # shift rotations
                           period = 1440) # repeat every day
  
  env <- lapply(1:N_SIM, function(i) { # apply func. for each sim.
    # progress$inc(1/N_SIM)     
    simmer() %>% # initialise a simulation environment
      add_resource("door", DOOR_SCHEDULE) %>% # add door
      
      add_resource("ct", capacity = CT) %>% # add CT
      add_resource("angio_inr", capacity = ANGIO_INR) %>% # add INR machine
      add_resource("angio_ir", capacity = ANGIO_IR) %>% # add IR machine
      
      add_resource("ed_staff", capacity = ED_STAFF) %>% # add ED staff
      add_resource("stroke_doctor", capacity = ST_STAFF) %>% # add stroke staff
      add_resource("angio_staff", STAFF_SCHEDULE) %>% # add angio. staff
      add_resource("inr", INR_SCHEDULE) %>% # add INR 
      add_resource("ir", IR_SCHEDULE) %>% # add IR
      
      add_generator("pt_ed", new_patient_traj, # ED patient arrival
                    function() rpois(1, I_ED) ) %>%
      add_generator("pt_inr", inr_traj, # INR patient arrival
                    function() rpois(1, I_INR) ) %>%
      add_generator("pt_eir", eir_traj, priority = 1, # EIR patient arrival
                    function() rpois(1, I_EIR) ) %>%
      add_generator("pt_ir", ir_traj, # IR patient arrival
                    function() rpois(1, I_IR) ) %>%
      
      run(RUN_T) #%>% # run simulation until RUN_T
    #wrap()
  })
  
  #######################################################################
  
  ## RETURN SIMULATION OUTPUT ##
  
  resources <- get_mon_resources(env) # resources monitoring stats
  
  monArr <- get_mon_arrivals(env, per_resource = TRUE) # arriv monitoring stats
  
  arrivals <- data.frame(get_mon_arrivals(env, per_resource = TRUE)) %>%
    transmute(name,
              resource, 
              replication, 
              start_time, 
              wait_time = end_time - start_time - activity_time) %>% 
    mutate(category = gsub("pt_([a-z]+)\\d+","\\1", name))
  
  arrivalsCounts <- arrivals %>% # no. of each arrivals
    select(name, category, replication) %>%
    group_by(replication) %>%
    unique(.) %>%
    summarize(total = n(),
              totalEd = sum(category == "ed"),
              totalElectInr = sum(category == "inr"),
              totalElectIr = sum(category == "ir"),
              totalEmergIr = sum(category == "eir")
    )
  
  arrivalsAngioInr <- subset(arrivals, resource == "angio_inr") %>%
    group_by(replication) %>%  # execute subsequent commands for each group
    summarize(total = n(),
              totalStroke = sum(category == "ed"),
              totalElectInr = sum(category == "inr"),
              totalElectIr = sum(category == "ir"),
              totalEmergIr = sum(category == "eir")
    )
  
  arrivalsAngioIr <- subset(arrivals, resource == "angio_ir") %>%
    group_by(replication) %>%  #execute subsequent commands for each group
    summarize(total = n(),
              totalStroke = sum(category == "ed"),
              totalElectInr = sum(category == "inr"),
              totalElectIr = sum(category == "ir"),
              totalEmergIr = sum(category == "eir")
    )
  
  print(arrivalsCounts)
  print(arrivalsAngioInr)
  print(arrivalsAngioIr)
  print("")
  
  list_containing_output <- list(arrivals, resources)
  
  return(list_containing_output)

}

#############################################################################
########################### PLOTTING FUNCTIONS ##############################
#############################################################################

plot_nav <- function(list_containing_output) {
 
  ###########################################################################
  ############## Function to plot the simulation results ####################
  ###########################################################################
  ## INPUTS: ##
  ## list_containing_output: [lst] list of results 
  ###########################################################################
  ## OUTPUTS: ##
  ## grid of plots
  ###########################################################################
  arrivals <- data.frame(list_containing_output[[1]]) 
  resources <- list_containing_output[[2]] 
  
  # s1 <- subset(arrivals, resource == "angio_inr") %>%
  #   group_by(replication, category) %>% 
  #   summarize("fast <20 min" = sum(wait_time < 20), 
  #             "medium 20-40 min" = sum(wait_time >= 20 & wait_time <= 40), 
  #             "slow >40 min" = sum(wait_time > 40)) %>%
  #   ungroup %>%
  #   dplyr::select(replication, category, 
  #                 "fast <20 min", "slow >40 min", "medium 20-40 min") %>%
  #   tidyr::gather(key=statistic, value = count, -replication, -category) %>% 
  #   filter(category != "inr" & category != "ir" & category != "eir") %>%
  #   ggplot(aes(factor(0),count)) + 
  #   geom_boxplot() + 
  #   facet_wrap(category~statistic) + 
  #   ylab("Patient count") + 
  #   ggtitle("Stroke patient wait time at angio INR") + 
  #   theme(plot.title = element_text(size=18, face = "bold"), 
  #         axis.text = element_text(size=12), 
  #         axis.title = element_text(size=14), 
  #         strip.text = element_text(size=14, face = "bold"),  
  #         axis.title.x=element_blank(), 
  #         axis.text.x=element_blank(), 
  #         axis.ticks.x=element_blank())
  # 
  # s2 <- subset(arrivals, resource == "angio_inr") %>%
  #   group_by(replication, category) %>%
  #   summarize("fast <20 min" = sum(wait_time < 20), 
  #             "medium 20-40 min" = sum(wait_time >= 20 & wait_time <= 40), 
  #             "slow >40 min" = sum(wait_time > 40)) %>%
  #   ungroup %>%
  #   dplyr::select(replication, category, 
  #                 "fast <20 min", "slow >40 min", "medium 20-40 min") %>%
  #   tidyr::gather(key=statistic, value = count, -replication, -category) %>%
  #   filter(category != "ed" & category != "ir" & category != "eir") %>%
  #   ggplot(aes(factor(0),count)) + 
  #   geom_boxplot() + 
  #   facet_wrap(category~statistic) + 
  #   ylab("Patient count") + 
  #   ggtitle("INR patient wait time at angio INR") + 
  #   theme(plot.title = element_text(size=18, face = "bold"), 
  #         axis.text = element_text(size=12), 
  #         axis.title = element_text(size=14), 
  #         strip.text = element_text(size=14, face = "bold"),  
  #         axis.title.x = element_blank(), 
  #         axis.text.x=element_blank(), 
  #         axis.ticks.x=element_blank())
  # 
  # if (sum(arrivals$resource == "angio_ir") == 0)
  #   s3 <- NULL
  # else {
  #   s3 <- filter(arrivals, resource == "angio_ir") %>%
  #     group_by(replication, category) %>% 
  #     summarize("fast <20 min" = sum(wait_time < 20), 
  #               "medium 20-40 min" = sum(wait_time >= 20 & wait_time <= 40), 
  #               "slow >40 min" = sum(wait_time > 40)) %>% 
  #     ungroup %>%
  #     dplyr::select(replication, category, 
  #                   "fast <20 min", "slow >40 min", "medium 20-40 min") %>%
  #     tidyr::gather(key=statistic, value = count, -replication, -category) %>%
  #     filter(category != "ed" & category != "inr" & category != "ir") %>%
  #     ggplot(aes(factor(0),count)) + 
  #     geom_boxplot() + 
  #     facet_wrap(category~statistic) +
  #     ylab("Patient count") +
  #     ggtitle("Emergency IR patient wait time at angio IR") +
  #     theme(plot.title = element_text(size=18, face = "bold"),
  #           axis.text = element_text(size=12), 
  #           axis.title = element_text(size=14), 
  #           strip.text = element_text(size=14, face = "bold"),
  #           axis.title.x=element_blank(), 
  #           axis.text.x=element_blank(), 
  #           axis.ticks.x=element_blank())
  # }
  # 
  # if (sum(arrivals$resource == "angio_ir") == 0)
  #   s4 <- NULL
  # else {
  #   s4<- subset(arrivals, resource == "angio_ir") %>%
  #     group_by(replication, category) %>%
  #     summarize("fast <20 min" = sum(wait_time < 20), 
  #               "medium 20-40 min" = sum(wait_time >= 20 & wait_time <= 40), 
  #               "slow >40 min" = sum(wait_time > 40)) %>%
  #     ungroup %>%
  #     dplyr::select(replication, category, 
  #                   "fast <20 min", "slow >40 min", "medium 20-40 min") %>%
  #     tidyr::gather(key=statistic, value = count, -replication, -category) %>%
  #     filter(category != "ed" & category != "inr" & category != "eir") %>%
  #     ggplot(aes(factor(0),count)) +
  #     geom_boxplot() +
  #     facet_wrap(category~statistic) +
  #     ylab("Patient count") +
  #     ggtitle("Non-emergency IR patient wait time at angio IR") +
  #     theme(plot.title = element_text(size=18, face = "bold"),
  #           axis.text = element_text(size=12),
  #           axis.title = element_text(size=14),
  #           strip.text = element_text(size=14, face = "bold"),  
  #           axis.title.x=element_blank(), 
  #           axis.text.x=element_blank(), 
  #           axis.ticks.x=element_blank())
  # }
  # 
  n1 <- subset(arrivals, category != "ir" & category != "eir" &
                 category != "inr" & resource != "door") %>%
    ggplot() +
    geom_density(aes(x=wait_time) ) +
    # facet_wrap(~resource) +
    ggtitle("Overview of Stroke patient wait time at various resources") +
    theme(plot.title = element_text(size=18, face = "bold"),
          axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          strip.text = element_text(size=14, face = "bold")) +
    ylab("Patient count") + scale_y_continuous(limits=c(0,80))
  # 
  # note that this not data.frame ob
  s2 <- plot(resources, metric ="usage", c("angio_inr", "inr"),
             items=c("server", "queue"), steps = TRUE)
  s3 <- plot(resources, metric ="usage", c("angio_ir", "ir"),
             items=c("server", "queue"), steps = TRUE)

  scat1 <- subset(arrivals, resource == "angio_inr") %>%
    ggplot(aes(x=start_time, y=wait_time, color=category)) +
    geom_point() +
    ggtitle("angio inr wait t vs run time")
  scat2 <- subset(arrivals, resource == "angio_ir") %>%
    ggplot(aes(x=start_time, y=wait_time, color=category)) +
    geom_point() +
    ggtitle("angio ir wait t vs run time")
  
  a <- subset(resources, resource != "door" )
  u1 <- plot(a, metric ="utilization") + 
    theme(plot.title = element_text(size=18, face = "bold"), 
          axis.text = element_text(size=12), 
          axis.title = element_text(size=14)) + 
    ylab("Median occupancy ratio") + ggtitle("Resource utilisation")
  
  myPlotList <- plyr::compact(list(n1))
  do.call(grid.arrange,  myPlotList)
}
  
server <- function(input, output) {
  
  getParams <- function(prefix) {
    input[[paste0(prefix, "_recalc")]]
    
    params <- lapply(paramNames, function(p) {
      input[[paste0(prefix, "_", p)]]
    })
    
    names(params) <- paramNames
    return(params)
  }
  
  output$a_distPlot <- renderPlot({
    input$goButton
    plot_nav(isolate(do.call(simulate_nav, getParams("a"))))
  })

}





