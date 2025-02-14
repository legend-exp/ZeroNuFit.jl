# API

## Functions and macros

```@index
Pages = ["internal_api.md"]
Order = [:macro, :function]
```

## Documentation

### analysis.jl
```@docs
ZeroNuFit.Analysis.retrieve_real_fit_results
ZeroNuFit.Analysis.run_analysis
ZeroNuFit.Analysis.save_outputs
```

### likelihood.jl
```@docs
ZeroNuFit.Likelihood.norm_uniform
ZeroNuFit.Likelihood.norm_linear
ZeroNuFit.Likelihood.norm_exponential
ZeroNuFit.Likelihood.exp_stable
ZeroNuFit.Likelihood.gaussian_plus_lowEtail
ZeroNuFit.Likelihood.get_mu_b
ZeroNuFit.Likelihood.get_mu_s
ZeroNuFit.Likelihood.get_mu_s_b
ZeroNuFit.Likelihood.build_likelihood_zero_obs_evts
ZeroNuFit.Likelihood.build_likelihood_per_partition
ZeroNuFit.Likelihood.build_likelihood_looping_partitions
ZeroNuFit.Likelihood.get_stat_blocks
ZeroNuFit.Likelihood.run_fit_over_partitions
ZeroNuFit.Likelihood.get_signal_prior_info
ZeroNuFit.Likelihood.get_signal_bkg_priors
ZeroNuFit.Likelihood.get_bkg_pdf
ZeroNuFit.Likelihood.get_signal_pdf
ZeroNuFit.Likelihood.build_prior
ZeroNuFit.Likelihood.generate_data
```

### plotting.jl
```@docs
ZeroNuFit.Plotting.plot_data
ZeroNuFit.Plotting.plot_fit_and_data
ZeroNuFit.Plotting.plot_correlation_matrix
ZeroNuFit.Plotting.plot_marginal_distr
ZeroNuFit.Plotting.plot_two_dim_posteriors
ZeroNuFit.Plotting.fit_model
```

### utils.jl
```@docs
ZeroNuFit.Utils.check_key
ZeroNuFit.Utils.get_settings
ZeroNuFit.Utils.get_events
ZeroNuFit.Utils.event_is_contained
ZeroNuFit.Utils.get_partitions
ZeroNuFit.Utils.get_partitions_new
ZeroNuFit.Utils.get_partitions_events
ZeroNuFit.Utils.get_partition_event_index
ZeroNuFit.Utils.get_corr_info
ZeroNuFit.Utils.get_energy_scale_pars
ZeroNuFit.Utils.get_efficiency
ZeroNuFit.Utils.get_deltaE
ZeroNuFit.Utils.get_range
ZeroNuFit.Utils.get_bkg_info
ZeroNuFit.Utils.get_global_mode
ZeroNuFit.Utils.get_marginalized_mode
ZeroNuFit.Utils.get_par_posterior
ZeroNuFit.Utils.inverse_uniform_cdf
ZeroNuFit.Utils.generate_disjoint_uniform_samples
ZeroNuFit.Utils.save_generated_samples
ZeroNuFit.Utils.save_results_into_json
```