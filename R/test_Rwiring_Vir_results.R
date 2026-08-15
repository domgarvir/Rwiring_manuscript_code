library(Rwiring)
library(dplyr)

interactions <- readr::read_csv("DATA/Ilist_byyear.csv")

names(interactions) <- c(
  "Order",
  "Site",
  "Time",
  "Origin_node",
  "Destiny_node",
  "Strength"
)

interactions_2 <- interactions %>% dplyr::select(-Order) %>% 
  mutate(Origin_guild = "Pollinator",
         Destiny_guild = "Plant")

# Prepare the interaction table
prepared <- prepare_interactions(
  interactions_2,
  directed = FALSE,
  multipartite = TRUE
)

summary(prepared)

# Compare two annual networks
changes <- quantify_network_changes(
  prepared,
  site_1 = "aznalcazar",
  time_1 = 2015,
  site_2 = "aznalcazar",
  time_2 = 2016
)

# Display a summary
summary(changes)

# Link-level classifications
head(changes$links)

# Node-level classifications
head(changes$nodes)

# Number of gained, lost and permanent nodes by guild
changes$nodes %>% dplyr::select(Guild,Node_status) %>% 
  group_by(Guild,Node_status) %>% count()
