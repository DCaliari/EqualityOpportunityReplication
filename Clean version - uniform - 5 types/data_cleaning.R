#---------------------------------------
# Load libraries
#---------------------------------------
library("tidyverse")
library(numDeriv)
library(parallel)
#---------------------------------------



#---------------------------------------
# Download and load the data
#---------------------------------------
url <- "https://www.dropbox.com/scl/fi/dywc9u23bufx1b3cuxmhx/October2025.csv?rlkey=pmc8d4vajfr9eith98kuq9pku&st=yjlf3fmm&dl=1"

dest <- "data_mobility.csv"
download.file(url, destfile = dest, mode = "wb")

df_raw <- read_csv("data_mobility.csv")
#---------------------------------------



#---------------------------------------
# Pre-processing of data
#---------------------------------------

# Rename variables to avoid "."
df_raw <- df_raw %>%
  rename_with(~ str_replace_all(.x, "\\.", "_"))

df <- df_raw %>%
  select(participant_code,
         socmob_October2025_36_player_stars, #PlayerStars
         socmob_October2025_36_player_parent, #PlayerParentStatus
         socmob_October2025_36_player_child, #PlayerChildStatus
         questions_October2025_1_player_top_third, #BeliefTop
         questions_October2025_1_player_middle_third, #BeliefMiddle
         questions_October2025_1_player_bottom_third, #BeliefBottom
         #socmob_October2025_1_player_q_order, #OrderParts
         socmob_October2025_1_player_total_order, #OrderQuestions
         socmob_October2025_1_player_position_order, #LeftRight
         #starts_with("questions_October2025_1_player_error_q_"), #ControlQuestions
         #starts_with("socmob_October2025_1_player_error_q_"),  #ControlQuestions
         matches("^socmob_October2025_.*_player_button_1$"), #Choices
         matches("^socmob_October2025_.*_player_button_2$") #Choices
         )

# Expands the question type and the left and right position of the alternatives into multiple columns
df <- df %>%
 separate(
    socmob_October2025_1_player_total_order,
    into = paste0("question_", 1:36),
    sep = " "
  ) %>%
  separate(
    socmob_October2025_1_player_position_order,
    into = paste0("right_left_", 1:36),
  )

# Change the name of the variables capturing the button pressing and puts the question number at the end
# Rename variables to get a shorter name
# Change status labels from Red, yellow, Green to High, Medium, Low
df <- df %>%
  rename_with(
    ~ str_replace(
        .x,
        "^socmob_October2025_(\\d+)_player_button_1$",
        "button_L_\\1"
      ),
    matches("^socmob_October2025_\\d+_player_button_1$")
  ) %>%
  rename_with(
    ~ str_replace(
        .x,
        "^socmob_October2025_(\\d+)_player_button_2$",
        "button_R_\\1"
      ),
    matches("^socmob_October2025_\\d+_player_button_2$")
  ) %>%
  rename(
    player_stars = socmob_October2025_36_player_stars,
    parent_status = socmob_October2025_36_player_parent,
    child_status = socmob_October2025_36_player_child,
    belief_top = questions_October2025_1_player_top_third,
    belief_middle = questions_October2025_1_player_middle_third,
    belief_bottom = questions_October2025_1_player_bottom_third
  ) %>%
  mutate(parent_status = replace_na(parent_status, "Spectators"),
    parent_status = recode(parent_status,
      "red" = "Low",
      "yellow" = "Medium",
      "green" = "High"
    ),
    parent_status = factor(
      parent_status,
      levels = c("Spectators", "High", "Medium", "Low")
    )
  )

# Convert the data table to long
df_long <- df %>%
  pivot_longer(
    cols = matches("^(right_left_|button_L_|button_R_|question_)\\d+$"),
    names_to = c(".value", "q"),
    names_pattern = "(right_left|button_L|button_R|question)_(\\d+)"
  ) %>%
  mutate(q = as.integer(q))

# Filter pout the questions that have not been seen
df_long <- df_long %>%
  filter(button_L == 1 | button_R == 1)

# Add the meritocracy factor
df_long <- df_long %>%
  mutate(
    meritocracy = case_when(
      str_detect(question, "_po_") ~ "positive",
      str_detect(question, "_ne_") ~ "negative",
      str_detect(question, "_no_") ~ "no",
      TRUE ~ NA_character_
    )
  )

# Add the ordered altenatives
df_long <- df_long %>%
  mutate(
    alternatives = case_when(
      str_detect(question, "q1_") ~ "A1_A4",
      str_detect(question, "q2_") ~ "A4_A7",
      str_detect(question, "q3_") ~ "A1_A7",
      str_detect(question, "q4_") ~ "A2_A3",
      str_detect(question, "q5_") ~ "A2_A3",
      str_detect(question, "q6_") ~ "A2_A7",
      str_detect(question, "q7_") ~ "A3_A7",
      str_detect(question, "q8_") ~ "A1_A7",
      str_detect(question, "q9_") ~ "A1_A4",
      str_detect(question, "q10_") ~ "A4_A7",
      str_detect(question, "q11_") ~ "A1_A7",
      str_detect(question, "q12_") ~ "A2_A3",
      TRUE ~ NA_character_
    )
  )

# Include the choice variable
df_long <- df_long %>%
  mutate(
    choice = if_else(
      !xor(button_R == 0, right_left == 0),
      word(alternatives, 1, sep = "_"),  # first alternative "_"
      word(alternatives, 2, sep = "_")   # second alternative "_"
    )
  )
  
rm(list = setdiff(ls(), c("df", "df_long")))
#---------------------------------------