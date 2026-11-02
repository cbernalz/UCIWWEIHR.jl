library(tidyverse)
library(here)
library(patchwork)
library(ggh4x)
library(viridis)
library(scales)

pp_hosp_inc_prev = read_csv(
  here(
    "UCIMaterials/LabWork/VolodLab/UCIWWEIHR.jl/model_output_1.csv"
  )
)

ww_data = read_csv(
  here(
    "UCIMaterials/LabWork/VolodLab/UCIWWEIHR.jl/simulated_data.csv"
  )
)

# Posterior Predictives - Wastewater
ww_pp_plot = pp_hosp_inc_prev %>%
  select(-iteration, -chain) %>%
  pivot_longer(
    everything()
  ) %>%
  filter(
    str_detect(name, "data_wastewater")
  ) %>%
  mutate(
    obstime = str_extract(
      name,
      "(?<=\\[)\\d+(?=\\])"
    ),
    obstime = as.numeric(obstime)
  ) %>%
  select(
    -name
  ) %>%
  group_by(
    obstime
  ) %>%
  summarize(
    lower_50th = quantile(value, 0.25),
    lower_80th = quantile(value, 0.1),
    lower_95th = quantile(value, 0.025),
    upper_50th = quantile(value, 0.75),
    upper_80th = quantile(value, 0.9),
    upper_95th = quantile(value, 0.975),
    median = median(value)
  ) %>%
  arrange(
    as.numeric(obstime)
  ) %>%
  mutate(
    date = c(
      ww_data$dates[1:90],
      seq.Date(
        from = as.Date("202-02-27"), 
        to = as.Date("2025-03-12"), 
        by = "day"
      )
    ),
    model = "EIHR-h-ww"
  ) %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = date,
      ymin = lower_95th,
      ymax = upper_95th,
      fill = model
    ),
    alpha = 0.3
  ) +
  geom_line(
    aes(
      x = date,
      y = median,
      color = model
    ),
  ) +
  geom_line(
    data = ww_data,
    aes(
      x = dates,
      y = log_ww_conc
    ),
    color = "black"
  ) +
  geom_point(
    data = ww_data,
    aes(
      x = dates,
      y = log_ww_conc,
      shape = "Fit"
    ),
    color = "black",
    size = 2.5
  ) +
  geom_point(
    data = ww_data,
    aes(
      x = dates,
      y = log_ww_conc,
      shape = "Forecast"
    ),
    color = "red",
    size = 2.5
  ) +
  geom_vline(
    xintercept = as.Date("2025-11-23"),
    linetype = "dashed"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels =  "%b %Y"
  ) +
  scale_fill_manual(
    breaks = c("EIHR-h-ww"),
    values = c("dodgerblue")
  ) +
  scale_color_manual(
    breaks = c("EIHR-h-ww"),
    values = c("dodgerblue")
  ) +
  scale_shape_manual(
    breaks = c("Fit", "Forecast"),
    values = c(16, 1)
  ) +
  labs(
    title = "Wastewater Posterior Predictive",
    y = "Log Conc",
    x = "Date (Fake Currently)",
    fill = "Model",
    color = "Model",
    shape = "Period"
  ) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust=1),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    #plot.background = element_rect(fill = "whitesmoke"),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    text = element_text(size = 18)
  )
ww_pp_plot


# Posterior Predictives - Prevalence
# hprev_pp_plot = pp_hosp_inc_prev %>%
#   select(-iteration, -chain) %>%
#   pivot_longer(
#     everything()
#   ) %>%
#   filter(
#     str_detect(name, "data_hosp_prev")
#   ) %>%
#   mutate(
#     obstime = str_extract(
#       name,
#       "(?<=\\[)\\d+(?=\\])"
#     ),
#     obstime = as.numeric(obstime)
#   ) %>%
#   select(
#     -name
#   ) %>%
#   group_by(
#     obstime
#   ) %>%
#   summarize(
#     lower_50th = quantile(value, 0.25),
#     lower_80th = quantile(value, 0.1),
#     lower_95th = quantile(value, 0.025),
#     upper_50th = quantile(value, 0.75),
#     upper_80th = quantile(value, 0.9),
#     upper_95th = quantile(value, 0.975),
#     median = median(value)
#   ) %>%
#   arrange(
#     as.numeric(obstime)
#   ) %>%
#   mutate(
#     date = seq(as.Date("2025-09-01"), by = "day", length.out = 84),
#     model = "EIHR-h-ww"
#   ) %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = date,
#       ymin = lower_95th,
#       ymax = upper_95th,
#       fill = model
#     ),
#     alpha = 0.3
#   ) +
#   geom_line(
#     aes(
#       x = date,
#       y = median,
#       color = model
#     ),
#   ) +
#   geom_line(
#     data = ww_data %>% filter(period=="Fit"),
#     aes(
#       x = date,
#       y = log_conc
#     ),
#     color = "black"
#   ) +
#   geom_point(
#     data = ww_data %>% filter(period=="Fit"),
#     aes(
#       x = date,
#       y = log_conc,
#       shape = "Fit"
#     ),
#     color = "black",
#     size = 2.5
#   ) +
#   geom_point(
#     data = ww_data %>% filter(period=="Forecast"),
#     aes(
#       x = date,
#       y = log_conc,
#       shape = "Forecast"
#     ),
#     color = "red",
#     size = 2.5
#   ) +
#   geom_vline(
#     xintercept = as.Date("2025-11-23"),
#     linetype = "dashed"
#   ) +
#   scale_x_date(
#     date_breaks = "1 month",
#     date_labels =  "%b %Y"
#   ) +
#   scale_fill_manual(
#     breaks = c("EIHR-h-ww"),
#     values = c("dodgerblue")
#   ) +
#   scale_color_manual(
#     breaks = c("EIHR-h-ww"),
#     values = c("dodgerblue")
#   ) +
#   scale_shape_manual(
#     breaks = c("Fit", "Forecast"),
#     values = c(16, 1)
#   ) +
#   labs(
#     title = "Prevalence Posterior Predictive",
#     y = "",
#     x = "Date (Fake Currently)",
#     fill = "Model",
#     color = "Model",
#     shape = "Period"
#   ) +
#   theme_bw() +
#   theme(
#     #axis.text.x = element_text(angle = 45, hjust=1),
#     #axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     #plot.background = element_rect(fill = "whitesmoke"),
#     strip.background = element_blank(),
#     panel.border = element_rect(colour = "black", fill = NA),
#     plot.title = element_text(face = "bold"),
#     plot.subtitle = element_text(face = "bold"),
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(face = "bold"),
#     strip.text = element_text(face = "bold"),
#     legend.text = element_text(face = "bold"),
#     legend.title = element_text(face = "bold"),
#     legend.position = "bottom",
#     text = element_text(size = 18)
#   )
# hprev_pp_plot
# 
# 
# ## Posterior Predictives - Incidence
# hinc_pp_plot = pp_hosp_inc_prev %>%
#   select(-iteration, -chain) %>%
#   pivot_longer(
#     everything()
#   ) %>%
#   filter(
#     str_detect(name, "data_hosp_inc")
#   ) %>%
#   mutate(
#     obstime = str_extract(
#       name,
#       "(?<=\\[)\\d+(?=\\])"
#     ),
#     obstime = as.numeric(obstime)
#   ) %>%
#   select(
#     -name
#   ) %>%
#   group_by(
#     obstime
#   ) %>%
#   summarize(
#     lower_50th = quantile(value, 0.25),
#     lower_80th = quantile(value, 0.1),
#     lower_95th = quantile(value, 0.025),
#     upper_50th = quantile(value, 0.75),
#     upper_80th = quantile(value, 0.9),
#     upper_95th = quantile(value, 0.975),
#     median = median(value)
#   ) %>%
#   arrange(
#     as.numeric(obstime)
#   ) %>%
#   mutate(
#     date = seq(as.Date("2025-09-01"), by = "day", length.out = 84),
#     model = "EIHR-h-ww"
#   ) %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = date,
#       ymin = lower_95th,
#       ymax = upper_95th,
#       fill = model
#     ),
#     alpha = 0.3
#   ) +
#   geom_line(
#     aes(
#       x = date,
#       y = median,
#       color = model
#     ),
#   ) +
#   geom_line(
#     data = ww_data %>% filter(period=="Fit"),
#     aes(
#       x = date,
#       y = log_conc
#     ),
#     color = "black"
#   ) +
#   geom_point(
#     data = ww_data %>% filter(period=="Fit"),
#     aes(
#       x = date,
#       y = log_conc,
#       shape = "Fit"
#     ),
#     color = "black",
#     size = 2.5
#   ) +
#   geom_point(
#     data = ww_data %>% filter(period=="Forecast"),
#     aes(
#       x = date,
#       y = log_conc,
#       shape = "Forecast"
#     ),
#     color = "red",
#     size = 2.5
#   ) +
#   geom_vline(
#     xintercept = as.Date("2025-11-23"),
#     linetype = "dashed"
#   ) +
#   scale_x_date(
#     date_breaks = "1 month",
#     date_labels =  "%b %Y"
#   ) +
#   scale_fill_manual(
#     breaks = c("EIHR-h-ww"),
#     values = c("dodgerblue")
#   ) +
#   scale_color_manual(
#     breaks = c("EIHR-h-ww"),
#     values = c("dodgerblue")
#   ) +
#   scale_shape_manual(
#     breaks = c("Fit", "Forecast"),
#     values = c(16, 1)
#   ) +
#   labs(
#     title = "Incidence Posterior Predictive",
#     y = "Log Conc",
#     x = "Date (Fake Currently)",
#     fill = "Model",
#     color = "Model",
#     shape = "Period"
#   ) +
#   theme_bw() +
#   theme(
#     #axis.text.x = element_text(angle = 45, hjust=1),
#     #axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     #plot.background = element_rect(fill = "whitesmoke"),
#     strip.background = element_blank(),
#     panel.border = element_rect(colour = "black", fill = NA),
#     plot.title = element_text(face = "bold"),
#     plot.subtitle = element_text(face = "bold"),
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(face = "bold"),
#     strip.text = element_text(face = "bold"),
#     legend.text = element_text(face = "bold"),
#     legend.title = element_text(face = "bold"),
#     legend.position = "bottom",
#     text = element_text(size = 18)
#   )
# hinc_pp_plot
