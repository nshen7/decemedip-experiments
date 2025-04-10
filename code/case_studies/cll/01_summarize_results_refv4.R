source('code/SETPATHS.R')
library(rstan)
library(bayesplot)
library(argparse)
library(Rsamtools)

read_dir <-  here('data', 'interim', 'case_studies', 'cll', 'from_ze')
write_dir <- here('data', 'interim', 'case_studies', 'cll', '01_summarize_results_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

md_samples <- fread(here(read_dir, 'AML_CLL_ctDNA.csv')) |>
  mutate(Group = gsub("[0-9]", "", Sample)) |>
  filter(Group == 'CLL') |>
  mutate(SampleFileName = paste0('mHS_', Sample, '_unique.sorted.dedup'))
fwrite(md_samples, here(write_dir, 'metadata_samples.csv'))

######################################
# ---- estimated proportions ----
######################################

all_samples_pi.df <- data.frame()
all_samples_pi_nonblood.df <- data.frame()

for (i in seq_len(nrow(md_samples))) {

  sample <- md_samples$Sample[i]
  group <- md_samples$Group[i]
  file <- md_samples$SampleFileName[i]

  sample_pi.df <- fread(here(read_dir, 'output', paste0('fitted_pi_', file,'.csv.gz'))) |>
    select(cell_type, mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group) |>
    relocate(sample, group)

  # sample_pi_nonblood.df <- sample_pi.df |>
  #   dplyr::slice(8:n()) |> # remove blood cells
  #   group_by(sample, group) |>
  #   summarise(total = sum(mean)) |>
  #   ungroup()

  all_samples_pi.df <- all_samples_pi.df |> rbind(sample_pi.df)
  # all_samples_pi_nonblood.df <- all_samples_pi_nonblood.df |> rbind(sample_pi_nonblood.df)

  cat(sample, ' ')
}
fwrite(all_samples_pi.df, here(write_dir, paste0('estimated_pi.csv.gz')))
# fwrite(all_samples_pi_nonblood.df, here(write_dir, paste0('estimated_pi_nonblood.csv.gz')))

######################################
# ---- estimated w_mu ----
######################################

### RUN ONCE
all_samples_w_mu.df <- data.frame(sample = NULL, group = NULL, mean = NULL,
                                  `2.5%` = NULL, `25%` = NULL, `50%` = NULL, `75%` = NULL, `97.5%` = NULL)

for (i in seq_len(nrow(md_samples))) {
  
  sample <- md_samples$Sample[i]
  group <- md_samples$Group[i]
  file <- md_samples$SampleFileName[i]

  sample_w_mu.df <- fread(here(read_dir, 'output', paste0('fitted_w_mu_', file,'.csv.gz'))) |>
    select(mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group, parameter = paste0('w_mu_', c(1:(n()-1), 0))) |>
    relocate(sample, group, parameter)

  all_samples_w_mu.df <- all_samples_w_mu.df |> rbind(sample_w_mu.df)

  cat(sample, ' ')
}
fwrite(all_samples_w_mu.df, here(write_dir, paste0('estimated_w_mu.csv.gz')))



######################################
# ---- estimated w_sigma ----
######################################

### RUN ONCE
all_samples_w_sigma.df <- data.frame(sample = NULL, group = NULL, mean = NULL,
                                  `2.5%` = NULL, `25%` = NULL, `50%` = NULL, `75%` = NULL, `97.5%` = NULL)

for (i in seq_len(nrow(md_samples))) {

  sample <- md_samples$Sample[i]
  group <- md_samples$Group[i]
  file <- md_samples$SampleFileName[i]
  
  sample_w_sigma.df <- fread(here(read_dir, 'output', paste0('fitted_w_sigma_', file,'.csv.gz'))) |>
    select(mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group, parameter = paste0('w_mu_', c(1:(n()-1), 0))) |>
    relocate(sample, group, parameter)

  all_samples_w_sigma.df <- all_samples_w_sigma.df |> rbind(sample_w_sigma.df)

  cat(sample, ' ')
}
fwrite(all_samples_w_sigma.df, here(write_dir, paste0('estimated_w_sigma.csv.gz')))


