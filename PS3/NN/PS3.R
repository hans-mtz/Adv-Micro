setwd(paste0(here::here(),"/PS3/NN"))
# install.packages("ggm", dependencies = TRUE)
library(ggm)
# ps <- powerset(1:3)

## Read data ----

# data <- read.csv("data/ABKK_full_experiment.csv")
data_a <- read.csv("data/ABKK_attributes.csv", stringsAsFactors = TRUE)

# table(data$Latitud)
## Prepare data ----
library(leaps)
library(car)
library(recipes)

df <- recipe(choice ~., data = data_a) %>%
  # step_num2factor(c(Order,Latitud,Longitud,Cost,Gender,
  #                   Education,Ethnicity,Marital_status,
  #                   Income, Employment,firs_question))# %>%
  step_scale(has_type(match = "double")) %>%
  step_dummy(has_type(match = "integer"), one_hot = TRUE) %>%
  step_dummy(all_nominal(), one_hot = TRUE) %>%
  prep(.,data=data_a) %>%
  bake(.,new_data=data_a)

# X <- df %>% select(-starts_with("binary"),-choice) %>% as.matrix()
# y <- df %>% select(choice) %>% pull

## ----
dnn <- df %>% select(-contains("binary"))
write.csv(dnn, "dnn.csv")

## regssubset ----
# sub_m <- regsubsets(x=X,y=y, method = "exhaustive",
#                      all.best=TRUE, nbest = 10, really.big = TRUE)

## My first Neural Net ----

# # Helper packages
# library(dplyr)         # for basic data wrangling
#
# # Modeling packages
# library(keras)         # for fitting DNNs
# library(tfruns)        # for additional grid search & model training functions
#
# # Modeling helper package - not necessary for reproducibility
# library(tfestimators)  # provides grid search & model training interface
