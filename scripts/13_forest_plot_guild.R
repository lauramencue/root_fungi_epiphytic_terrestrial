

# Load data
load("analyses.Rdata")
load("alpha_diversity_results.Rdata")

library(broom.mixed)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

forest_plot_glmmTMB <- function(model,
                                xlab,
                                title,
                                exponentiate = TRUE) {
  
  td <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      estimate_exp = if (exponentiate) exp(estimate) else estimate,
      conf.low_exp = if (exponentiate) exp(conf.low) else conf.low,
      conf.high_exp = if (exponentiate) exp(conf.high) else conf.high
    )
  
  # Clean and harmonize labels
  td <- td %>%
    mutate(
      term = str_replace_all(term, "`", ""),
      term = str_replace_all(term, "host_family", ""),
      term = str_replace(term, "^habitatTerrestrial$", 
                         "Habitat: Terrestrial"),
      term = str_replace(term, "^epiphitic/terrestrialT$", 
                         "Growth form: Terrestrial"),
      term = str_replace(term, "^host_family", ""),
      term = str_replace_all(term, ":", " × "),
      term = str_replace(term, "habitatTerrestrial × ", 
                         "Habitat: Terrestrial × "),
      term = str_replace(term, "epiphitic/terrestrialT × ", 
                         "Growth form: Terrestrial × "),
      term = fct_rev(fct_inorder(term))
    )
  
  ggplot(td, aes(y = term, x = estimate_exp)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = conf.low_exp, xmax = conf.high_exp),
                   height = 0.25, linewidth = 0.7) +
    geom_point(size = 2.3) +
    scale_x_log10() +
    labs(x = xlab, y = NULL, title = title) +
    theme_classic(base_size = 12)
}

p_symb <- forest_plot_glmmTMB(
  m_symb,
  xlab  = "Odds Ratio (OR)",
  title = "Effects on proportion of symbiotic ASVs (binomial GLMM)"
)

p_symb


library(dplyr)

habitat_effects <- fixed_tbl_all %>%
  filter(term == "habitatTerrestrial") %>%
  mutate(
    lower = exp(estimate - 1.96 * std.error),
    upper = exp(estimate + 1.96 * std.error)
  ) %>%
  select(guild, rate_ratio, lower, upper, p.value)

habitat_effects <- habitat_effects %>%
  mutate(
    guild = factor(
      guild,
      levels = c(
        "arbuscular_mycorrhizal",
        "dark_septate_endophyte",
        "ectomycorrhizal",
        "root_endophyte_endomycorrhizal",
        "ericoid_mycorrhizal",
        "ericoid_orchid_mycorrhizal",
        "orchid_mycorrhizal"
      )
    )
  )

ggplot(habitat_effects,
       aes(x = rate_ratio, y = guild)) +
  
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    linewidth = 0.6
  ) +
  
  geom_point(size = 3) +
  
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4")
  ) +
  
  labs(
    x = "Rate ratio (Terrestrial / Epiphytic)",
    y = "Functional guild"
  ) +
  
  theme_classic(base_size = 12)





