source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions
library(argparse)
library(rstan)

write_dir <- here('data', 'interim', 'case_studies', 'cell_lines', '01_apply_deconv_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'case_studies', 'cell_lines', '01_apply_deconv_refv4')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

# ---- Load data ----
md_medip <- read_tsv(here('data', 'metadata', 'medip_vs_wgbs', 'metadata_medip.tsv')) |>
  dplyr::mutate(bam_file = here('data', 'raw', 'medip_vs_wgbs', 'medip', paste0(`File accession`, '.bam')))

# ---- Fit decemedip-model model ----

for (INDEX in 1:2) {
  sample_name <- md_medip$`Biosample term name`[INDEX]
  
  ## !! MCMC chains get stuck for some samples, in that case, change a seed
  output <- decemedip(sample_bam_file = md_medip$bam_file[INDEX],
                      paired_end = FALSE,
                      diagnostics = TRUE)
  fit <- output$posterior
  # plotDiagnostics(output, 'model_fit')
  
  smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
  fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_', sample_name,'.csv')))
  
  smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
  fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_', sample_name,'.csv')))
  
  smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
    as.data.frame() |>
    mutate(cell_type = factor(colnames(decemedip::hg19.ref.cts.se), levels = rev(colnames(decemedip::hg19.ref.cts.se)))) |>
    relocate(cell_type)
  fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_', sample_name,'.csv')))
  
  plotDiagnostics(decemedip_output = output, plot_type = 'model_fit')
  ggsave(here(plot_dir, paste0('bayesplot_diagnostics_', sample_name,'.png')), width = 12, height = 10)
  
  plotDiagnostics(decemedip_output = output, plot_type = 'y_fit')
  ggsave(here(plot_dir, paste0('point_pred_y_vs_x_', sample_name,'.png')), width = 12, height = 10)
}
