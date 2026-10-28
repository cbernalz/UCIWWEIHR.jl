library(tidyverse)
library(here)
library(patchwork)
library(formattable)
library(paletteer)
library(sns) 
library(scoringutils)

custom_theme = theme_bw() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.background = element_rect(fill = "whitesmoke"),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

pp_hosp_inc_prev = read_csv(
  here(
    "UCIMaterials/LabWork/VolodLab/UCIWWEIHR.jl/model_output_1.csv"
  )
)

params = read_csv(
  here(
    "UCIMaterials/LabWork/VolodLab/UCIWWEIHR.jl/model_output_2.csv"
  )
)
## Priors

## Reading in data from path in preliminaries, ww is assumed to be included.
df = read_csv(
  here(
    "UCIMaterials/LabWork/VolodLab/UCIWWEIHR.jl/simulated_data.csv"
  )
)

params_of_inter = c("E_init", "I_init", "H_init", 
                    "gamma", "nu", "epsilon", 
                    "rho_gene", "sigma_ww", "sigma_hosp")
param_trace_plot = params %>%
  pivot_longer(
    cols = -c(chain, iteration)
  ) %>%
  filter(
    name %in% params_of_inter
  ) %>%
  mutate(
    name = factor(
      name,
      levels = params_of_inter
    ),
    chain = factor(chain)
  ) %>%
  ggplot(
    aes(
      x = iteration,
      y = value,
      color = chain
    )
  ) +
  geom_line() +
  scale_color_paletteer_d("fishualize::Acanthurus_coeruleus") +
  facet_wrap(
    .~name,
    scales = "free"
  ) +
  labs(
    y = "Sample",
    x = "Iteration",
    color = "Chain"
  ) +
  custom_theme +
  theme(
    axis.text.y = element_text()
  )

param_effective_samp_size = params %>%
  pivot_longer(
    cols = -c(chain, iteration)
  ) %>%
  filter(
    name %in% params_of_inter
  ) %>%
  mutate(
    Parameter = factor(
      name,
      levels = params_of_inter
    ),
    Chain = factor(chain)
  ) %>%
  select(-chain, -name) %>%
  group_by(Chain, Parameter) %>%
  summarise(
    `Effective Sample Size` = ess(value)
  ) %>%
  formattable(
    list(
      `Effective Sample Size` = formatter(
        "span", 
        style = x ~ ifelse(
          x > 1000, #TODO: What should this number be?
          "color: green;", 
          "color: red;"
        )
      )
    )
  )


acf_plots = params %>%
  pivot_longer(
    cols = -c(chain, iteration)
  ) %>%
  filter(
    name %in% params_of_inter
  ) %>%
  mutate(
    name = factor(
      name,
      levels = params_of_inter
    ),
    chain = factor(chain)
  ) %>%
  group_by(chain, name) %>%
  nest() %>%
  mutate(
    acf_result = map(data, ~ acf(.x$value, plot = FALSE)),
    acf_data = map(acf_result, ~ data.frame(
      lag = .$lag, 
      acf = .$acf, 
      ci = qnorm(0.975) / sqrt(.$n.used)
    ))
  ) %>%
  unnest(acf_data) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = lag, 
      y = acf, 
      color = chain
    )
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0.1,
    linetype = "dashed",
    color = "red"
  ) +
  geom_hline(
    yintercept = -0.1,
    linetype = "dashed",
    color = "red"
  ) +
  geom_segment(
    aes(
      xend = lag, 
      yend = 0
    )
  ) +
  scale_color_paletteer_d("fishualize::Acanthurus_coeruleus") +
  facet_wrap(
    .~name, 
    scales = "free_y"
  ) +
  labs(
    x = "Lag", 
    y = "ACF"
  ) +
  custom_theme +
  theme(
    axis.text.y = element_text(),
  )


# gq = rbind(
#   raw_gq %>%
#     mutate(bayes = "Posterior"),
#   prior_gq %>%
#     mutate(bayes = "Prior")
# )
# params_of_inter = c("E_init", "I_init", "H_init", 
#                     "gamma", "nu", "epsilon", 
#                     "rho_gene", "sigma_ww", "sigma_hosp")
# hosp_ww_nt_plot = gq %>%
#   select(-chain, -iteration) %>%
#   pivot_longer(
#     cols = c(-bayes)
#   ) %>%
#   filter(
#     name %in% params_of_inter
#   ) %>%
#   mutate(
#     name = factor(
#       name,
#       levels = params_of_inter
#     )
#   ) %>%
#   ggplot(
#     aes(
#       x = value,
#       fill = bayes,
#       color = bayes
#     )
#   ) +
#   geom_density(
#     color = "black",
#     alpha = 0.6
#   ) +
#   geom_rug(
#     sides = "b"
#   ) +
#   scale_color_paletteer_d("fishualize::Acanthurus_coeruleus") +
#   scale_fill_paletteer_d("fishualize::Acanthurus_coeruleus") +
#   facet_wrap(
#     .~name,
#     scales = "free"
#   ) +
#   labs(
#     y = "Density",
#     x = "Parameter"
#   ) +
#   custom_theme
# hosp_ww_nt_plot



