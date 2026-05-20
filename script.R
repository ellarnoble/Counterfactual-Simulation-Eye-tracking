library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(arm)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# (1) Load data
###############################################################

# Relative path: RData must be in same directory as script
path <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "RData")

files <- list.files(path, full.names = TRUE)

# Create df using all files
df_long <- files |>
  map_dfr(function(f) {
    e <- new.env()
    load(f, envir = e)
    obj <- get(ls(e)[1], envir = e)
    obj
  })


df_long |> count(condition, participant)

# (2) Select 10 participants per condition with least data loss
###############################################################
# data loss = proportion of frames where v.max is NA
data_loss <- df_long |>
  group_by(condition, participant) |>
  summarise(
    pct_loss = mean(is.na(v.max)) * 100,
    .groups = "drop"
  )

# keep 10 with lowest data loss per condition
keep_participants <- data_loss |>
  group_by(condition) |>
  slice_min(order_by = pct_loss, n = 10) |>
  ungroup() |>
  select(condition, participant)

# filter main data to only those participants
df_long <- df_long |>
  semi_join(keep_participants, by = c("condition", "participant"))

# Filter data samples for only the first instance of the clip
df_long <- df_long |>
  filter(instance == 1)

# check
df_long |> count(condition, participant)

# (3) Reproduce main findings from study 
#################################################################
# Filter data to pre-collision
df_pre <- df_long |>
  filter(frame < collision.time)

# summarise per participant per condition
df_summary <- df_pre |>
  group_by(condition, participant) |>
  summarise(
    pct_B_counterfactual = mean(v.max == "B.counterfactual", na.rm = TRUE) * 100
  )

# box plot comparing conditions
ggplot(df_summary, aes(x = condition, y = pct_B_counterfactual)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(y = "% B counterfactual looks", x = "Condition")

# Check condition means and sd
df_summary |> 
  group_by(condition) |> 
  summarise(mean = mean(pct_B_counterfactual),
            sd = sd(pct_B_counterfactual))

# ANOVA model: Is condition a significant 
# predictor of counterfactual gazes?

aov_result <- aov(pct_B_counterfactual ~ condition, data = df_summary)
summary(aov_result)

# post hoc
TukeyHSD(aov_result)

# (4) Time series analysis
###############################################################

# calculate mean B-counterfactual probability per frame per condition
df_timecourse <- df_long |>
  group_by(condition, frame) |>
  summarise(
    mean_B_counterfactual = mean(B.counterfactual, na.rm = TRUE),
    se = sd(B.counterfactual, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# find average collision time
mean_collision <- df_long |> 
  distinct(clip, collision.time) |> 
  summarise(mean(collision.time)) |> 
  pull()

# min and max collision time across clips for instance 1
collision_range <- df_long |>
  filter(instance == 1) |>
  distinct(clip, collision.time) |>
  summarise(
    min_collision = min(collision.time),
    max_collision = max(collision.time),
    mean_collision = mean(collision.time))

collision_range

# Plot conditions time series data 
ggplot(df_timecourse, aes(x = frame, y = mean_B_counterfactual, 
                          colour = condition)) +
  geom_vline(xintercept = collision_range$min_collision, 
             linetype = "dotted", colour = "black") +
  geom_vline(xintercept = collision_range$max_collision, 
             linetype = "dotted", colour = "black") +
  geom_line() +
  geom_ribbon(aes(ymin = mean_B_counterfactual - se,
                  ymax = mean_B_counterfactual + se,
                  fill = condition), alpha = 0.25, colour = NA) +
  labs(x = "Frame", y = "Mean B-counterfactual look probability",
       fill = "Condition") +
  scale_colour_manual(values = c("causality" = "cornflowerblue", 
                                 "counterfactuals" = "rosybrown",
                                 "observations" = "seagreen")) +
  scale_fill_manual(values = c("causality" = "cornflowerblue", 
                               "counterfactuals" = "rosybrown",
                               "observations" = "seagreen"),
                    labels = c("causality" = "Causality",
                               "counterfactuals" = "Counterfactual",
                               "observations" = "Outcome")) +
  guides(colour = "none") +
  theme_classic()

# (5) AR(1) Model for Time Series Data (prior to collision)
##############################################################

# create binary outcome from l.max
df_long <- df_long |>
  mutate(B_cf = as.integer(l.max == "B.counterfactual"))

# create lagged binary outcome
df_long <- df_long |>
  arrange(participant, clip, instance, frame) |>
  group_by(participant, clip, instance) |>
  mutate(B_cf_lag = dplyr::lag(B_cf)) |>
  ungroup()

# filter to pre-collision frames
df_pre <- df_long |>
  filter(frame < collision.time)

# AR(1) logistic GLMM on pre-collision frames only
timemodel <- glmer(B_cf ~ condition + B_cf_lag +
                     (1 | participant) + (1 | clip),
                   data = df_pre |> filter(!is.na(B_cf_lag)),
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"))

summary(timemodel)

# (6) Model Diagnostics
#############################################################
icc(timemodel)

r2(timemodel)

# (7) Onset Frame Analysis 
###############################################################
# Determine when probabilities of B-counterfactual looks start
# to significantly differ across conditions
# Causality versus Observations
frame_effects_ca.o <- df_long |>
  filter(condition %in% c("causality", "observations")) |>
  group_by(frame) |>
  summarise(
    t_stat = t.test(B.counterfactual ~ condition)$statistic,
    p_val = t.test(B.counterfactual ~ condition)$p.value
  )

ca.o_onset_frame <- df_long |>
  filter(condition %in% c("causality", "observations")) |>
  group_by(frame) |>
  summarise(p_val = t.test(B.counterfactual ~ condition)$p.value,
            .groups = "drop") |>
  filter(!is.na(p_val)) |>
  mutate(sig = p_val < .05,
         run_id = cumsum(!sig)) |>
  group_by(run_id) |>
  mutate(run_length = cumsum(sig)) |>
  ungroup() |>
  filter(run_length >= 5) |>
  slice(1) |>
  mutate(ca.o_onset_frame = frame - 5 + 1) |>
  pull(ca.o_onset_frame)


# get p-value at onset frame
frame_effects_ca.o |>
  filter(frame == ca.o_onset_frame) |>
  select(frame, t_stat, p_val)


# Counterfactual versus Observations
frame_effects_co.o <- df_long |>
  filter(condition %in% c("counterfactuals", "observations")) |>
  group_by(frame) |>
  summarise(
    t_stat = t.test(B.counterfactual ~ condition)$statistic,
    p_val = t.test(B.counterfactual ~ condition)$p.value
  )

co.o_onset_frame <- df_long |>
  filter(condition %in% c("counterfactuals", "observations")) |>
  group_by(frame) |>
  summarise(p_val = t.test(B.counterfactual ~ condition)$p.value,
            .groups = "drop") |>
  filter(!is.na(p_val)) |>
  mutate(sig = p_val < .05,
         run_id = cumsum(!sig)) |>
  group_by(run_id) |>
  mutate(run_length = cumsum(sig)) |>
  ungroup() |>
  filter(run_length >= 5) |>
  slice(1) |>
  mutate(co.o_onset_frame = frame - 5 + 1) |>
  pull(co.o_onset_frame)

# get p-value at onset frame
frame_effects_co.o |>
  filter(frame == co.o_onset_frame) |>
  select(frame, t_stat, p_val)


# Causality versus Counterfactual
frame_effects_cc <- df_long |>
  filter(condition %in% c("causality", "counterfactuals")) |>
  group_by(frame) |>
  summarise(
    t_stat = t.test(B.counterfactual ~ condition)$statistic,
    p_val = t.test(B.counterfactual ~ condition)$p.value
  )

cc_onset_frame <- df_long |>
  filter(condition %in% c("causality", "counterfactuals")) |>
  group_by(frame) |>
  summarise(p_val = t.test(B.counterfactual ~ condition)$p.value,
            .groups = "drop") |>
  filter(!is.na(p_val)) |>
  mutate(sig = p_val < .05,
         run_id = cumsum(!sig)) |>
  group_by(run_id) |>
  mutate(run_length = cumsum(sig)) |>
  ungroup() |>
  filter(run_length >= 5) |>
  slice(1) |>
  mutate(cc_onset_frame = frame - 5 + 1) |>
  pull(cc_onset_frame)

# get p-value at onset frame
frame_effects_cc |>
  filter(frame == cc_onset_frame) |>
  select(frame, t_stat, p_val)

# (8) Session Info for Reproducibility
###############################################################
sink("sessionInfo.txt")
sessionInfo()
sink()
