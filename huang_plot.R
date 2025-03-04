################################
## REPRODUCE FIG.2 - BASELINE ##
################################

## LOAD RELEVANT PACKAGES ##
library(tidyverse)
library(ggthemes)
library(extrafont)

rel_probs <- function(sim_output, resource){
  
  ####################################################################
  ######## Function to obtain binned relative probabilities ##########
  ####################################################################
  ## Parameters ##                                                  ##
  ## sim_output:  output from simulation [df]                       ##
  ## resource:    resource to investigate [character]               ##
  ####################################################################
  ## Returns ##                                                     ##
  ## [df] ##                                                        ##
  ####################################################################
  
  ## PREPARE DATA ##
  data <- sim_output[[1]] # extract arrivals data
  data$wait_time[which(data$wait_time < 0)] <- 0 # -ve set to 0
  data <- data %>% mutate(wait = wait_time > 0) # wait ? (bool)
  
  ## ASSIGN TO BINS ##
  
  # Wait time > 0 for resource #
  x1 <- data$wait_time[which(data$resource == resource & 
                               data$wait == TRUE)]
  # Wait time = 0 for resource # 
  x2 <- data$wait_time[which(data$resource == resource &
                             data$wait == FALSE)]
  
  # Seq. inputs #
  from <- 0
  to <- max(x1)
  by <- to/500
  
  # Bins and midpoints #
  x <- cut(x1, seq(from,to,by), include.lowest=TRUE) # cut x1 into bins
  vec <- seq(from+by/2, to-by/2, by) # midpoints
  names(vec) <- levels(x) # matching names to levels
  cut_vec <- unname(vec[x]) # map names of bins to midpoints
  
  # Transform to relative probabilities #
  y <- table(cut_vec)/length(x2) # take ratios
  df2 <- as.data.frame( (2/(1+exp(-y)) - 1) ) # sigmoid transform
  df2[,1] <- as.numeric(as.character(df2[,1])) # change to num.
  df2 <- rbind(c(0,1), df2, c(200,1e-10)) # add value for wait_time = 0 & 200
  df2 <- cbind(df2, rep(resource, times = nrow(df2))) # add resource 
  colnames(df2)[3] <- 'resource'
  
  ## RETURN DATA ##
  df2
}


plot_func <- function(sim_lst, scenarios){
  
  ##############################################################
  ############## Function to reproduce figure 2 ################
  ##############################################################
  
  ## PARAMETERS ##
  ## df_lst:    list of simulated data frames [list] ##
  ## scenarios: vector of scenarios simulated [vec]  ##
  
  ## OUTPUT ##
  ## plot ## 
  
  ## PREPARE DATA ##
  
  n <- length(sim_lst) # length of list
  plot_df <- data.frame() # empty df for storage
  
  resources <- c('angio_inr', 'angio_staff', 'ct', 
                    'ed_staff', 'inr', 'stroke_doctor')
  
  for (i in 1:n) { # for each df in df_lst
    
    sim_output <- sim_lst[[i]] # reassign name
     
    rel_prob_df <- data.frame() # empty df for storage
    
    ## Obtain rel. probs for required resources ##
    
    for (j in resources) { # for each resource
      rel_prob_df <- rbind(rel_prob_df, rel_probs(sim_output,j))
    }
    
    ## Assign scenario names to data ##
    rel_prob_df <- cbind(rel_prob_df, 
                         rep(scenarios[i], nrow(rel_prob_df)))
    colnames(rel_prob_df)[ncol(rel_prob_df)] <- 'scenario'
    plot_df <- rbind(plot_df, rel_prob_df)
  }
  
  ## PRODUCE PLOT ##
  
  plt <- ggplot(plot_df) + 
    geom_line(aes(x = cut_vec, y = Freq, col = resource),
              linewidth = 0.75, alpha = 0.9) +
    facet_wrap(vars(factor(scenario, c(scenarios)))) +
    xlim(0,200) +
    coord_trans(x = 'sqrt', y = 'sqrt') + 
    xlab('Patient wait time (min)') + 
    ylab('Standardized density of patients in queue') + 
    guides(colour = guide_legend(nrow = 1)) +
    theme_minimal() +
    theme(text = element_text(family = 'Centaur'),
          axis.title.x = element_text(colour = 'black'), 
          axis.title.y = element_text(colour = 'black'), 
          strip.text = element_text(colour = 'black', 
                                    face = 'bold'),
          legend.position = 'bottom',
          legend.direction = 'horizontal', 
          legend.title = element_blank())
  
  plt
  
}


plot_func_hours <- function(sim_lst1, sim_lst2, sim_lst3, scenarios){
  
  res <- 'angio_inr'
  n <- length(scenarios) # no. of scenarios
  plot_df <- data.frame() # empty df for storage 
  
  for (i in 1:n) { # for each scenario
    
    rel_df1 <- rel_probs(sim_lst1[[i]], res)
    rel_df1 <- cbind(rel_df1, rep('5pm', nrow(rel_df1)))
    colnames(rel_df1)[ncol(rel_df1)] <- 'shift_end'
    
    rel_df2 <- rel_probs(sim_lst2[[i]], res)
    rel_df2 <- cbind(rel_df2, rep('6pm', nrow(rel_df2)))
    colnames(rel_df2)[ncol(rel_df2)] <- 'shift_end'
    
    rel_df3 <- rel_probs(sim_lst3[[i]], res)
    rel_df3 <- cbind(rel_df3, rep('7pm', nrow(rel_df3)))
    colnames(rel_df3)[ncol(rel_df3)] <- 'shift_end'
    
    df <- rbind(rel_df1, rel_df2, rel_df3)
    df <- cbind(df, rep(scenarios[i], nrow(df)))
    colnames(df)[ncol(df)] <- 'scenario'
    plot_df <- rbind(plot_df, df)
    
  }

  
  ## Produce plot ##
  plt <- ggplot(plot_df) + 
    geom_line(aes(x = cut_vec, y = Freq, col = shift_end), 
              linewidth = 0.75, alpha = 0.9) +
    facet_wrap(vars(factor(scenario, scenarios))) +
    xlim(0,200) + 
    coord_trans(x = 'sqrt', y = 'sqrt') + 
    xlab('Patient wait time (min)') + 
    ylab('Standardized density of patients in queue') + 
    guides(colour = guide_legend(nrow = 1)) +
    theme_minimal() +
    theme(text = element_text(family = 'Centaur'),
          strip.text = element_text(colour = 'black', 
                                    face = 'bold'),
          legend.position = 'bottom',
          legend.direction = 'horizontal', 
          legend.title = element_blank())
  
  plt
  
}


plot_wait <- function(df_base, df_lst, scenarios){
  
  ##########################################################
  
  ## Function to plot the mean disability-free life added 
  ## compared to baseline                                 
  
  ## PARAMETERS ##
  ## df_base:   baseline scenario df [df]
  ## df_lst:    list of dataframes for each scenario [lst]
  ## scenarios: vector of scenarios [vec]
  
  ## OUTPUT ##
  ## plot
  
  ##########################################################
  
  # Average wait time for baseline #
  x0 <- mean(
    df_base$wait_time[which(df_base$category == 'ed' &
                            df_base$resource == 'angio_inr')]
    ) 
  
  # Average wait time for all other scenarios #
  n <- length(df_lst)
  x <- numeric(n) # empty vector for storing
  
  for (i in 1:n){
    
    x[i] <- mean(
      df_lst[[i]]$wait_time[which(df_lst[[i]]$category == 'ed' &
                              df_lst[[i]]$resource == 'angio_inr')]
    ) 
  }
  
  # Form df of results #
  df <- as.data.frame(cbind((x0-x) * 4.2, scenarios))
  colnames(df) <- c('added_time', 'scenario')
  df$added_time <- as.numeric(df$added_time)
  df
  
  # Create plot #
  ggplot(df) +
    geom_col(aes(x = factor(scenario, scenarios),
                 y = added_time)) +
    xlab('Scenarios') +
    ylab('Mean disability-free life added (days)')
  
}

df_base <- baseline[[1]]
df_lst <- list(exclusive[[1]], 
               double[[1]], 
               exclusive_hr[[1]])
scenarios <- c('Exclusive-use', 
               'Two angio INRs', 
               'Exclusive-use and +1hr work')
plot_wait(df_base, df_lst, scenarios)

##############################################################################

## Fig. 5 Plot ##

plot_util <- function(df_lst, scenarios, res){
  
  ##########################################################
  
  ## Function to plot the utilization % for a particular
  ## resource under each scenario
  
  ## PARAMETERS ##
  ## df_lst:    list of dataframes for each scenario [lst]
  ## scenarios: vector of scenarios [vec]
  ## res:       resources to isolate [str]
  
  ## OUTPUT ##
  ## plot
  
  ##########################################################
  
  ###############################
  ## OBTAIN UTILIZATION VALUES ##
  ###############################
  
  n <- length(df_lst) # no. of scenarios
  df <- data.frame() # empty data frame for storing
  
  for (i in 1:n) { # for each scenario
    
    ## Extract only the required resource ## 
    df_lst[[i]] <- df_lst[[i]][
      which(df_lst[[i]]$resource == res),]
    
    ## Code from simmer to obtain utilization values ## 
    df_lst[[i]] <- df_lst[[i]] %>%
      group_by(.data$resource, .data$replication) %>%
      mutate(dt = lead(.data$time) - .data$time) %>%
      mutate(capacity = ifelse(.data$capacity < .data$server, 
                               .data$server, .data$capacity)) %>%
      mutate(dt = ifelse(.data$capacity > 0, .data$dt, 0)) %>%
      mutate(in_use = .data$dt * .data$server / .data$capacity) %>%
      summarise(utilization = sum(.data$in_use, na.rm = TRUE) / 
                  sum(.data$dt, na.rm=TRUE)) %>%
      summarise(Q25 = stats::quantile(.data$utilization, .25),
                Q50 = stats::quantile(.data$utilization, .5),
                Q75 = stats::quantile(.data$utilization, .75))
    
    df <- rbind(df, c(df_lst[[i]]$Q25,
                      df_lst[[i]]$Q50,
                      df_lst[[i]]$Q75, 
                      scenarios[i])
                )
  }
  
  #################
  ## CREATE PLOT ##
  #################
  
  colnames(df) <- c('Q25', 'Q50', 'Q75', 'scenario')
  df[1:3] <- sapply(df[,1:3], as.numeric)
  q50_vals <- paste(round(100*df$Q50), "%", sep='')
  
  ggplot(df) +
    aes(x = scenario, y = Q50, ymin = Q25, ymax = Q75) +
    geom_col() +
    geom_errorbar(width = .25, color = "black") +
    geom_text(aes(label = q50_vals), vjust = -1) +
    scale_y_continuous(labels = scales::percent, 
                       limits = c(0, 0.4), 
                       breaks = seq(0, 2, .2)) +
    xlab('Scenario') +
    ylab('Utilization')
  
}

df_lst <- list(baseline[[2]],
               exclusive[[2]], 
               double[[2]])
scenarios <- c('Baseline',
               'Exclusive-use', 
               'Two angio INRs')
res <- 'angio_inr'
plot_util(df_lst, scenarios, res)





