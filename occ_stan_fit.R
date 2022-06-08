model <- cmdstan_model(file.path("models", "occ.stan"))
fit3 <- model$sample(data = list(M = occ_data$M,
                                 J = occ_data$J,
                                 Y = occ_data$y,
                                 Elev = occ_data$elev,
                                 Forest = occ_data$forest,
                                 Wind = occ_data$wind),
                     seed = 1,
                     chains = 4, parallel_chains = 4,
                     iter_sampling = 1000,
                     iter_warmup = 1000,
                     output_dir = output_dir)
