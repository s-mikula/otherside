library(dtplyr)
library(tidyverse)
library(fixest)
library(Formula)
library(lfe)
library(modelsummary)

rm(list = ls())

# File with replication data / chnage path if needed.
load("DATA/data_replication.RData")

#### Baseline estimation sample ####
rdata_other <- rdata %>% 
  mutate(
    dff = PORCISLO - lastrow,
    otherside = PORCISLO > lastrow
  ) %>% 
  filter(
    between(dff,0,1)
  ) %>% 
  group_by(ballot) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n==2) %>% 
  mutate(
    lastrow_cat = cut(lastrow, breaks = c(0,23,30,40), labels = c("c1","c2","c3"))
  )

# Variables:
# year = "Election year/Election ID",
# ballot = "Ballot ID -- combination of KSTRANA (election-specific party ID), VOLKRAJ (constituency ID) and year",
# pref_hlasy = "Number of preferential votes (n)",
# otherside = "Indicator variable for candidates listed on reverse side",
# male = "Indicator variable for males",
# agecat = "Age category",
# VEK = "Age",
# ISCO = "ISCO 1 (other == not categorized in ISCO)",
# maxBallot = "Total number of candidates on the ballot",
# POC_HLASU = "Total number of votes casted for the party",
# dff = "Distance to break",
# ISCED* = "Indicator variables for ISCED categories"


#### Table 1: Balance and descriptives ####

# Part A: age, gender, education
rdata_other %>% 
  ungroup() %>% 
  distinct(ballot,dff, .keep_all = TRUE) %>% 
  select(dff,pref_hlasy,starts_with("ISCED"),male,VEK) %>% 
  mutate(
    VEK = as.integer(VEK)
  ) %>% 
  mutate(
    across(
      where(is.logical),
      as.integer
    )
  ) %>% 
  mutate(
    dff = ifelse(dff == 0, "front side", "reverse side")
  ) %>% 
  datasummary_balance(
    ~ dff,
    .,
    dinm_statistic = "p.value",
    fmt = 3
  )

# Part B: occupation
rdata_other %>% 
  ungroup() %>% 
  distinct(ballot,dff, .keep_all = TRUE) %>% 
  select(ballot,dff,ISCO) %>% 
  mutate(
    ISCO = str_c("ISCO_",ISCO),
    value = 1
  ) %>% 
  complete(ballot,dff,ISCO,fill = list(value = 0)) %>% 
  pivot_wider(names_from = "ISCO") %>% 
  select(-ballot) %>% 
  mutate(
    dff = ifelse(dff == 0, "front side", "reverse side")
  ) %>% 
  datasummary_balance(
    ~dff,
    .,
    dinm_statistic = "p.value",
    fmt = 3
  )

#### Table 2: Results ####
modell <- list(
  #pref_hlasy ~ otherside,
  pref_hlasy ~ otherside | ballot,
  pref_hlasy ~ otherside | male + ISCED5 + ISCED6 + ISCED7 + ISCED8 + ISCEDB + agecat + ballot + maxBallot,
  pref_hlasy ~ otherside | male + ISCED5 + ISCED6 + ISCED7 + ISCED8 + ISCEDB + agecat + ISCO  + ballot + maxBallot 
) %>% map(as.formula)

# Models (1)-(3)
modell %>% 
  map(
    feglm,
    data = rdata_other,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 8
  ) %>% 
  etable(digits = 3)

# Models (4)-(6)
rdata_other_p3 <- rdata %>% 
  filter(ballot %in% rdata_other$ballot) %>% 
  mutate(
    dff = PORCISLO - lastrow,
    otherside = PORCISLO > lastrow
  ) %>% 
  filter(
    dff %in% c(-1,1)
  ) %>% 
  group_by(ballot) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n == 2)

modell %>% 
  map(
    feglm,
    data = rdata_other_p3,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

#### Table 3: Placebo ####

# Placebo shift: Two positions upwards
placebo_m2 <- rdata %>% 
  filter(ballot %in% rdata_other$ballot) %>% 
  mutate(
    lastrow = lastrow - 2,
    dff = PORCISLO - lastrow,
    otherside = PORCISLO > lastrow
  ) %>% 
  filter(
    between(dff,0,1)
  ) %>% 
  group_by(ballot) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n == 2)

modell %>% 
  map(
    feglm,
    data = placebo_m2,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

# Placebo shift: One position upwards
placebo_m1 <- rdata %>% 
  filter(ballot %in% rdata_other$ballot) %>% 
  mutate(
    lastrow = lastrow - 1,
    dff = PORCISLO - lastrow,
    otherside = PORCISLO > lastrow
  ) %>% 
  filter(
    between(dff,0,1)
  ) %>% 
  group_by(ballot) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n == 2)

modell %>% 
  map(
    feglm,
    data = placebo_m1,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

# Placebo shift: One position downwards
placebo_p1 <- rdata %>% 
  filter(ballot %in% rdata_other$ballot) %>% 
  mutate(
    lastrow = lastrow + 1,
    dff = PORCISLO - lastrow,
    otherside = PORCISLO > lastrow
  ) %>% 
  filter(
    between(dff,0,1)
  ) %>% 
  group_by(ballot) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n == 2)

modell %>% 
  map(
    feglm,
    data = placebo_p1,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

# Placebo shift: Two positions downwards
placebo_p2 <- rdata %>% 
  filter(ballot %in% rdata_other$ballot) %>% 
  mutate(
    lastrow = lastrow + 2,
    dff = PORCISLO - lastrow,
    otherside = PORCISLO > lastrow
  ) %>% 
  filter(
    between(dff,0,1)
  ) %>% 
  group_by(ballot) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n == 2)

pp2 <- modell %>% 
  map(
    feglm,
    data = placebo_p2,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

#### Table 4 ####
keep <- rdata_other %>% 
  distinct(ballot,.keep_all = TRUE) %>% 
  select(year,KSTRANA,VOLKRAJ,lastrow) %>% 
  group_by(KSTRANA,year) %>% 
  summarise(
    lrs = length(unique(lastrow))
  ) %>% 
  ungroup() %>% 
  filter(lrs > 1) %>% 
  select(-lrs) %>% 
  mutate(
    year = as.character(year)
  )

rdata_other_rc2 <- rdata_other %>% 
  semi_join(.,keep)

modell %>% 
  map(
    feglm,
    data = rdata_other_rc2,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

#### Table 5 ####
keep <- rdata_other %>% 
  distinct(ballot,.keep_all = TRUE) %>% 
  select(year,KSTRANA,VOLKRAJ,lastrow) %>% 
  group_by(VOLKRAJ,year) %>% 
  summarise(
    lrs = length(unique(lastrow))
  ) %>% 
  ungroup() %>% 
  filter(lrs > 1) %>% 
  select(-lrs) %>% 
  mutate(
    year = as.character(year)
  )

rdata_other_rc2c <- rdata_other %>% 
  semi_join(.,keep)

modell %>% 
  map(
    feglm,
    data = rdata_other_rc2c,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

#### Table 6 ####
rdata_other_rc1 <- rdata_other %>% 
  filter(lastrow > 9) 

modell %>% 
  map(
    feglm,
    data = rdata_other_rc1,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)


#### Table 7 ####
rdata_other_rc1b <- rdata_other %>% 
  filter(lastrow > 23) 

modell %>% 
  map(
    feglm,
    data = rdata_other_rc1b,
    cluster = "ballot",
    offset = ~log(POC_HLASU),
    family = quasipoisson,
    nthreads = 12
  ) %>% 
  etable(digits = 3)

#### Figure 1 ####
ranking_df <- rdata %>% 
  filter(maxBallot > lastrow) %>%
  mutate(
    dff = PORCISLO - lastrow,
  ) %>% 
  group_by(PORCISLO) %>% 
  summarise(
    pref = sum(pref_hlasy),
    .groups = "drop"
  ) %>%
  mutate(
    pref = 100*pref/sum(pref)
  )

rdata %>% 
  distinct(ballot, .keep_all = TRUE) %>% 
  filter(!is.na(lastrow)) %>% 
  filter(maxBallot > lastrow) %>% 
  group_by(lastrow) %>% 
  summarise(
    lastr = n(),
    .groups = "drop"
  ) %>% 
  mutate(
    lastr = 100*lastr/sum(lastr)
  ) %>% 
  rename(
    PORCISLO = lastrow
  ) %>% 
  left_join(ranking_df,.) %>% 
  ggplot(
    aes(x = PORCISLO)
  ) +
  geom_col(
    aes(y = lastr, fill = "lastr")#,
    #,fill = "grey50"
  ) +
  geom_line(
    aes(y = pref, color = "pvotes")
  ) +
  geom_point(
    aes(y = pref, shape = "pvotes"),
    #shape = 21,
    fill = "white"
  ) +
  scale_fill_manual(
    values = c("lastr"="grey50"),
    labels = c("lastr"="Last row"),
  ) +
  scale_color_manual(
    "a",
    values = c("pvotes"="black"),
    labels = c("pvotes"="Preferential votes"),
  ) +
  scale_shape_manual(
    "a",
    values = c("pvotes"=21),
    labels = c("pvotes"="Preferential votes"),
  ) +
  scale_x_continuous(
    "Position on the ballot",
    breaks = seq(from =1, to = 36, by = 2)
  ) +
  scale_y_continuous(
    "Share (%)"
  ) +
  theme_classic(
    base_family = "Times"
  ) +
  theme(
    legend.title = element_blank()
  )

#### Figure 2 ####
rdata %>% 
  filter(maxBallot > lastrow) %>%
  mutate(
    dff = PORCISLO - lastrow,
  ) %>% 
  filter(
    between(dff,-2,3)
  ) %>% 
  group_by(dff) %>% 
  summarise(
    ll = 100*sum(pref_hlasy)/sum(POC_HLASU)
  ) %>% 
  unnest(ll) %>% 
  ggplot(
    aes(x = factor(dff), y = ll, fill = dff>0)
  ) +
  #geom_errorbar() +
  geom_col() +
  scale_y_continuous(
    "Preferential votes (% of all votes)"
  ) +
  scale_x_discrete(
    "Position on the ballot",
    labels = c(
      "b-3","b-2","b-1",
      "b+1","b+2","b+3"
    )
  ) +
  scale_fill_grey(
    "Side of the ballot",
    start = 0.35,
    end = 0.65,
    #palette = "Greys",
    labels = c(
      "Front","Reverse"
    )
  ) +
  theme_classic(
    base_family = "Times"
  )

#### Figure 5a ####
rdata %>% 
  filter(maxBallot > lastrow) %>%
  distinct(year,KSTRANA,VOLKRAJ,.keep_all = TRUE) %>% 
  select(year,KSTRANA,VOLKRAJ,lastrow) %>% 
  group_by(KSTRANA,year) %>% 
  summarise(
    lrs = length(unique(lastrow))
  ) %>% 
  ungroup() %>% 
  ggplot(
    aes(x = lrs)
  ) +
  geom_bar(
    fill = "grey50"
  ) +
  scale_y_continuous(
    "Count (n)"
  ) +
  scale_x_continuous(
    "Number of unique last rows on the front side\nof the ballot per party and election",
    breaks = 1:8
  ) +
  theme_classic(
    base_family = "Times"
  )

#### Figure 5b ####
rdata %>% 
  filter(maxBallot > lastrow) %>%
  distinct(year,KSTRANA,VOLKRAJ,.keep_all = TRUE) %>% 
  select(year,KSTRANA,VOLKRAJ,lastrow) %>% 
  group_by(VOLKRAJ,year) %>% 
  summarise(
    lrs = length(unique(lastrow))
  ) %>% 
  ungroup() %>% 
  ggplot(
    aes(x = lrs)
  ) +
  geom_bar(
    fill = "grey50"
  ) +
  scale_y_continuous(
    "Count (n)"
  ) +
  scale_x_continuous(
    "Number of unique last rows on the front side\nof the ballot per constituency and election",
    breaks = 1:8
  ) +
  theme_classic(
    base_family = "Times"
  )

#### Figure 6 ####

lastrows <- rdata_other %>% select(lastrow_pl = lastrow) %>% filter(lastrow_pl >= 10)

set.seed(4326)

placebo_data <- 1:1000 %>% 
  map(
    function(x){
      rdata %>% 
        distinct(ballot) %>% 
        filter(!(ballot %in% rdata_other$ballot)) %>% 
        slice_sample(n = nrow(lastrows), replace = TRUE) %>% 
        bind_cols(.,lastrows) %>% 
        inner_join(rdata,., , by = "ballot") %>% 
        mutate(
          dff = PORCISLO - lastrow_pl,
          otherside = PORCISLO > lastrow_pl
        ) %>% 
        filter(
          between(dff,0,1)
        ) %>% 
        group_by(ballot) %>% 
        add_tally() %>% 
        ungroup() %>% 
        filter(n==2)
    }
  )

placebo_est <- placebo_data %>% 
  map_dfr(
    function(d){
      feglm(
        pref_hlasy ~ otherside | male + ISCED5 + ISCED6 + ISCED7 + ISCED8 + ISCEDB + agecat + ISCO  + ballot + maxBallot,
        data = d,
        offset = ~log(POC_HLASU),
        family = quasipoisson
      ) %>% 
        tidy()
    }
  )

placebo_est %>% 
  ggplot(
    aes(x = estimate)
  ) +
  stat_ecdf(geom = "step") +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  scale_x_continuous(
    "Placebo estimates",
    breaks = c(-0.45,-0.35,-0.25,-0.15,0,0.15,0.25,0.35,0.45)
  ) +
  scale_y_continuous(
    "Cumulative distribution",
    breaks = seq(from = 0, to = 1, by = 0.1)
  ) +
  theme_classic(
    base_family = "Times"
  )