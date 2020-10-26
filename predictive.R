# A script to create predictive model for Bcl-xL inhibition
# Created by: svavil
# Created on: 2020-10-24
# Edited on: 2020-10-26

library(tidyverse)
library(readxl)
library(randomForest)
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.minor = element_blank())
data_block_1 <- read_xlsx("bcl-xl_1_descriptors.xlsx") %>% 
  rename(IC = `Inhibition(10.9uM)`)
inconsistent_rows <- data_block_1[xor(M_ACTIVITY_OUTCOME == "Inactive", data_block_1$IC <= 10.9),]
stopifnot(nrow(inconsistent_rows) == 0)data_block_1$PUBCHE

ggplot(data = data_block_1) + 
  geom_point(aes(x = IC, y = Norm_Act, color = PUBCHEM_ACTIVITY_OUTCOME)) + 
  labs(x = "Inhibitory concentration (uM)", 
       y = "Normalized activity", 
       color = "Active?")

ggplot(data = data_block_1) + 
  geom_point(aes(x = doxrub_similarity, y = Norm_Act, color = PUBCHEM_ACTIVITY_OUTCOME)) + 
  labs(x = "Fingerprint similarity to doxorubicin", 
       y = "Normalized activity", 
       color = "Active?")

ggplot(data = data_block_1) + 
  geom_point(aes(x = exact_mol_weight, y = Norm_Act, color = PUBCHEM_ACTIVITY_OUTCOME)) + 
  labs(x = "Molecular weight (g/mol)", 
       y = "Normalized activity", 
       color = "Active?")

ggplot(data = data_block_1) + 
  geom_point(aes(x = ring_count, y = Norm_Act, color = PUBCHEM_ACTIVITY_OUTCOME)) + 
  labs(x = "Ring count", 
       y = "Normalized activity", 
       color = "Active?")

ggplot(data = data_block_1) + 
  geom_point(aes(x = reference_similarity, y = doxrub_similarity))

predictors <- data_block_1 %>%
  select(-PUBCHEM_ACTIVITY_SCORE, -PUBCHEM_SID, -PUBCHEM_ACTIVITY_OUTCOME, -IC, -SMILES, -Norm_Act, -PUBCHEM_CID)

model.regression <- lm(data_block_1$IC~., data = predictors)
predicted.IC <- predict.lm(model.regression, newdata = predictors)
ggplot() + 
  geom_point(aes(x = data_block_1$IC, y = predicted.IC)) + 
  geom_line(aes(x = data_block_1$IC, y = data_block_1$IC), color = "red", size = 1)

model.class <- glm(data_block_1$PUBCHEM_ACTIVITY_OUTCOME %>% recode("Inactive" = 0, "Active" = 1) ~ ., data = predictors, family = "binomial")
predicted.outcome <- predict.glm(model.class, newdata = predictors, family = "binomial") > 0
table(predicted.outcome, data_block_1$PUBCHEM_ACTIVITY_OUTCOME)

rf.model <- randomForest(data_block_1$IC ~ ., data = predictors, ntree = 1000)
rf.predictions <- predict(rf.model, newdata = predictors)
ggplot() + 
  geom_point(aes(x = data_block_1$IC, y = rf.predictions)) + 
  geom_line(aes(x = data_block_1$IC, y = data_block_1$IC), color = "red", size = 1)
