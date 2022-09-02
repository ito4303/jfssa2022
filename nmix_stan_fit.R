model <- cmdstan_model(file.path("models", "nmix.stan"))
fit6 <- model$sample(data = list(M = nmix_data$nsites,
                                 J = nmix_data$nvisits,
                                 Y = nmix_data$C,
                                 K = 100 + max(nmix_data$C),
                                 Xsite = nmix_data$site.cov[, 2],
                                 Xsurvey = nmix_data$survey.cov),
                     seed = 1, chains = 4, parallel_chains = 4,
                     iter_sampling = 1000, iter_warmup = 1000,
                     output_dir = "output_dir")
