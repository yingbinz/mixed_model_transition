# load data and library


library(dplyr)
library(purrr)
library(tidyr)
library(lme4)
library(texreg)



#### preprocess data####

## compute intervals between observations

#studentid: student id.  
#msoffsetfromstart: timestamp of the observation.  
#affect: affect observation.


raw_data <- read.csv("example_data.csv")

raw_data$sid2 <- c(raw_data[-1,"studentid"], NA)
raw_data$t2 <- c(raw_data[-1,"msoffsetfromstart"], NA)

raw_data$interval <- with(raw_data, ifelse(studentid == sid2, (t2-msoffsetfromstart)/1000, NA))
summary(raw_data$interval)


## convert affect sequence to frequencies of transitions
## Note in the dataset, engagement is named CONCENTRATING.

three_list <- list()

studentids <- unique(raw_data$studentid)
affect_pairs <- expand.grid(c('CONCENTRATING','CONFUSED','FRUSTRATED','BORED'),
                            c('CONCENTRATING','CONFUSED','FRUSTRATED','BORED'),
                            c('CONCENTRATING','CONFUSED','FRUSTRATED','BORED'))
colnames(affect_pairs) <- c("Conditional", "Given", "Target")
for (i in 1:length(studentids)) {
  single_sequence <- filter(raw_data, studentid == studentids[i])
  three_seq <- cbind(single_sequence$affect, c(single_sequence[-1, "affect"], NA), c(single_sequence[c(-1,-2), "affect"], NA, NA), single_sequence$interval, unlist(c(single_sequence[-1, "interval"], NA)))
  
  three_seq <- as.data.frame(three_seq)
  colnames(three_seq) <- c("Conditional", "Given", "Target", "first_interval", "second_interval")
  
  three_seq$studentid <- studentids[i]
  three_list[[i]] <- three_seq
}

base_df <- bind_rows(three_list) %>% 
  # drop invalid transitions due to data preprocessing and transitions involving unknown states
  na.omit()  %>% 
  mutate(
    # dummy coding the conditional affective states
    c_engagement = ifelse(Conditional == "CONCENTRATING" , 1, 0),
    c_confus = ifelse(Conditional == "CONFUSED" , 1, 0),
    c_frustration = ifelse(Conditional == "FRUSTRATED" , 1, 0),
    c_bored = ifelse(Conditional == "BORED" , 1, 0),
    
    # dummy coding the given affective states
    g_engagement = ifelse(Given == "CONCENTRATING" , 1, 0),
    g_confus = ifelse(Given == "CONFUSED" , 1, 0),
    g_frustration = ifelse(Given == "FRUSTRATED" , 1, 0),
    g_bored = ifelse(Given == "BORED" , 1, 0),
    
    # dummy coding the target affective states
    t_engagement = ifelse(Target == "CONCENTRATING", 1, 0),
    t_confus = ifelse(Target == "CONFUSED", 1, 0),
    t_frustration = ifelse(Target == "FRUSTRATED", 1, 0),
    t_bored = ifelse(Target == "BORED" , 1, 0))

# convert interval columns into numeric formats
base_df$first_interval <- as.numeric(base_df$first_interval)
base_df$second_interval <- as.numeric(base_df$second_interval)

base_df <- base_df %>% 
  ## drop transitions with long intervals. Here the threshold is 180
  filter(first_interval <= 180 & second_interval <= 180) %>% 
  ## drop columns no longer used
  select(-first_interval, -second_interval)

## the code below is for filtering students with too few affect observations
observation_threshold <- 10
nfreq <- base_df %>%
  group_by(studentid) %>%
  summarise(n = n()) %>%
  filter(n > observation_threshold)

base_df <- base_df %>%
  filter(studentid %in% unlist(nfreq$studentid))


#### transition analysis - compare marginal, conditional, and interaction models####

## to CONCENTRATING
three_df <- base_df %>% mutate(
  # interaction
  engagement_confus = ifelse(Conditional == "CONCENTRATING" & Given == "CONFUSED", 1, 0),
  confusion_confus = ifelse(Conditional == "CONFUSED" & Given == "CONFUSED", 1, 0),
  frustration_confus = ifelse(Conditional == "FRUSTRATED" & Given == "CONFUSED", 1, 0),
) %>% dplyr::select(-Conditional, -Given, -Target)

## estimating the marginal transition strengths, i.e., the model with only given states as predictors
engagement_marginal <- glmer(
  ## The left side is the dummy coding variable for the target state, while the left side contains the given states. For '(1 | studentid)', the '1' represents the random effect of the intercept, while 'studentid' indicates the id column.
  t_engagement ~ g_engagement + g_confus + (1 | studentid), 
  data = three_df, ## The data used
  family = binomial, ## Set the distribution as binomial
  control = glmerControl(optimizer ="Nelder_Mead"))

## adding conditional states to the model to estimate the conditional transition strengths
engagement_marginal_conditional <- update(engagement_marginal, 
                                          .~.+c_engagement + c_confus + c_frustration)

## adding interactions between given and conditional states
engagement_interaction <- update(engagement_marginal_conditional,
                                 .~.+ engagement_confus + confusion_confus + frustration_confus)

screenreg(list(engagement_marginal, engagement_marginal_conditional, engagement_interaction))

anova(engagement_marginal_conditional, engagement_interaction)
anova(engagement_marginal, engagement_marginal_conditional)



## to CONFUSED
three_df <- base_df %>% mutate(
  # interaction
  engagement_engagement = ifelse(Conditional == "CONCENTRATING" & Given == "CONCENTRATING", 1, 0),
  confusion_engagement = ifelse(Conditional == "CONFUSED" & Given == "CONCENTRATING", 1, 0),
  frustration_engagement = ifelse(Conditional == "FRUSTRATED" & Given == "CONCENTRATING", 1, 0),
  
  engagement_frustration = ifelse(Conditional == "CONCENTRATING" & Given == "FRUSTRATED", 1, 0),
  confusion_frustration = ifelse(Conditional == "CONFUSED" & Given == "FRUSTRATED", 1, 0),
  frustration_frustration = ifelse(Conditional == "FRUSTRATED" & Given == "FRUSTRATED", 1, 0)
) %>% dplyr::select(-Conditional, -Given, -Target)



confusion_marginal <- glmer(t_confus ~ g_engagement + g_confus + g_frustration + (1 | studentid), 
                 data = three_df, 
                 family = binomial,
                 control = glmerControl(optimizer ="Nelder_Mead"))

confusion_marginal_conditional <- update(confusion_marginal, 
                                         .~.+c_engagement + c_confus + c_frustration)

confusion_interaction <- update(confusion_marginal_conditional, 
                                         .~.+ engagement_engagement + confusion_engagement + 
                                  frustration_engagement + engagement_frustration + 
                                  confusion_frustration + frustration_frustration)


screenreg(list(confusion_marginal, confusion_marginal_conditional, confusion_interaction))

anova(confusion_marginal_conditional, confusion_interaction)
anova(confusion_marginal_conditional, confusion_marginal)



## to FRUSTRATED
three_df <- base_df %>% mutate(
  # interaction
  engagement_confus = ifelse(Conditional == "CONCENTRATING" & Given == "CONFUSED", 1, 0),
  confusion_confus = ifelse(Conditional == "CONFUSED" & Given == "CONFUSED", 1, 0),
  frustration_confus = ifelse(Conditional == "FRUSTRATED" & Given == "CONFUSED", 1, 0),

  engagement_bored = ifelse(Conditional == "CONCENTRATING" & Given == "BORED", 1, 0),
  confusion_bored = ifelse(Conditional == "CONFUSED" & Given == "BORED", 1, 0),
  frustration_bored = ifelse(Conditional == "FRUSTRATED" & Given == "BORED", 1, 0)
) %>% dplyr::select(-Conditional, -Given, -Target)

fru_marginal <- glmer(t_frustration ~ g_confus + g_frustration + g_bored + (1 | studentid), 
                 data = three_df, 
                 family = binomial,
                 control = glmerControl(optimizer ="Nelder_Mead"))

fru_marginal_conditional <- update(fru_marginal, .~.+c_engagement + c_confus + c_frustration)

fru_interaction <- update(fru_marginal_conditional, 
                          .~.+ engagement_confus + confusion_confus +
                            frustration_confus + engagement_bored + 
                            confusion_bored +frustration_bored)

screenreg(list(fru_marginal, fru_marginal_conditional, fru_interaction))

anova(fru_marginal_conditional, fru_interaction)
anova(fru_marginal_conditional, fru_marginal)


## to BORED
three_df <- base_df %>% mutate(
  # interaction
  engagement_frustration = ifelse(Conditional == "CONCENTRATING" & Given == "FRUSTRATED", 1, 0),
  confusion_frustration = ifelse(Conditional == "CONFUSED" & Given == "FRUSTRATED", 1, 0),
  
  frustration_frustration = ifelse(Conditional == "FRUSTRATED" & Given == "FRUSTRATED", 0, 1)
) %>% dplyr::select(-Conditional, -Given, -Target)

bor_marginal <- glmer(t_bored ~ g_frustration + g_bored + (1 | studentid), 
                 data = three_df, 
                 family = binomial,
                 control = glmerControl(optimizer ="Nelder_Mead"))

bor_marginal_conditional <- update(bor_marginal, .~.+c_engagement + c_confus + c_frustration)

bor_interaction <- update(bor_marginal_conditional, .~.+engagement_frustration + confusion_frustration + 
    frustration_frustration)

screenreg(list(bor_marginal, bor_marginal_conditional, bor_interaction))

anova(bor_marginal_conditional, bor_interaction)
anova(bor_marginal_conditional, bor_marginal)


#### prediction####

library(caret)
library(pROC)

# create a function to run 10 folds cross-validation and compute the performance metric

get_cv_metrics <- function(input_model, your_data, num_folds=10){
  formula <-  input_model@call$formula
  outcome <- as.character(formula)[2]
  
  # Create the group-wise folds using groupKFold()
  set.seed(230806)  # For reproducibility
  folds <- groupKFold(your_data$studentid, k = num_folds)
  
  # Initialize an empty vector to store the cross-validated predictions
  cv_predictions <- rep(NA, nrow(your_data))
  
  metrics <- matrix(NA, nrow = num_folds, ncol = 4)
  colnames(metrics) <- c("accuracy", "kappa", "F1","AUC")
  
  # Implement group-wise cross-validation
  for (i in 1:num_folds) {
    # Get the indices for the current fold
    train_indices <- folds[[i]]
    
    # Split data into training and testing based on the fold indices
    train_data <- your_data[train_indices, ]
    test_data <- your_data[-train_indices, ]
    
    # Fit the mixed-effects model using glmer on the training data
    model <- glmer(formula, data = train_data, family = binomial)  # Modify family as needed
     
    # Predict on the test data for this fold
    predictions <- predict(model, test_data, type = "response", allow.new.levels = TRUE)
    
    pred_class <- as.factor(predictions > 0.5)
    true_labels <- as.factor(as.logical(test_data[,outcome]))    
    
    c_matrix <- confusionMatrix(pred_class, true_labels, mode = "everything",
                  positive="TRUE")
    
    metrics[i,] <- c(c_matrix$overall[c( 'Accuracy','Kappa')], c_matrix$byClass['F1'], auc(roc(true_labels, predictions, direction = '<')))
  }
  return(metrics)
}

metric_list <- list()


## to engagement
three_df <- base_df %>% mutate(
  
  engagement_confus = ifelse(Conditional == "CONCENTRATING" & Given == "CONFUSED", 1, 0),
  confusion_confus = ifelse(Conditional == "CONFUSED" & Given == "CONFUSED", 1, 0),
  frustration_confus = ifelse(Conditional == "FRUSTRATED" & Given == "CONFUSED", 1, 0),
) %>% dplyr::select(-Conditional, -Given, -Target)


engagement_marginal_conditional_prediction <- get_cv_metrics(engagement_marginal_conditional, three_df)
engagement_marginal_prediction <- get_cv_metrics(engagement_marginal, three_df)

colMeans(engagement_marginal_conditional_prediction)
colMeans(engagement_marginal_prediction)


#### adding external factor####

## create a synthetic external factor

# Create a synthetic factor of prior knowledge, ranging from 1-20.
set.seed(2024)
simulated_prior_knowledge <- sample(1:20, length(unique(base_df$studentid)), replace = TRUE)
# centering the external factor to ease interpretation
simulated_prior_knowledge <- simulated_prior_knowledge - mean(simulated_prior_knowledge)

prior_knowledge <- data.frame(studentid = unique(base_df$studentid),
                              prior_knowledge = simulated_prior_knowledge)


## examine the influence of the external factor on transitions
three_df <- base_df %>% 
  dplyr::select(-Conditional, -Given, -Target) %>% 
  left_join(prior_knowledge, by = 'studentid')

## estimating the transition strength without pripr knowledge
engagement_marginal_conditional <- glmer(
  t_engagement ~ +c_engagement + c_confus + c_frustration+g_engagement + g_confus + (1 | studentid), 
  data = three_df, 
  family = binomial, 
  control = glmerControl(optimizer ="Nelder_Mead"))

## adding pripr knowledge to the model 
engagement_marginal_conditional_with_prior_knowledge <- update(engagement_marginal_conditional, 
                                          .~.+prior_knowledge + prior_knowledge*g_confus)

screenreg(list(engagement_marginal_conditional,
          engagement_marginal_conditional_with_prior_knowledge))



#### reference states ####

## to CONCENTRATING


three_df <- base_df %>% select(-Conditional, -Given, -Target)

# all other transitions as reference
engagement_marginal_all_ref <- glmer(t_engagement ~ g_confus + (1 | studentid), 
                 data = three_df, 
                 family = binomial,
                 control = glmerControl(optimizer ="Nelder_Mead"))

# control self-transitions
engagement_marginal <- update(engagement_marginal_all_ref, .~.+g_engagement )

# only one transition as reference: boredom -> engagement
engagement_marginal_one_ref <- update(engagement_marginal, .~.+g_engagement+g_frustration )

screenreg(list(engagement_marginal_all_ref, engagement_marginal, engagement_marginal_one_ref))


