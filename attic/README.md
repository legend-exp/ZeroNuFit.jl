Here we collected a series of example scripts that might be useful to retrieve information or to plot results:
- `mbb.jl`: script to evaluate the effective Majorana neutrino mass (mbb) and plot different distributions (signal, NME, mbb). The script contains also a series of tests performed to understand how the different probability distributions transform and impacts the final mbb limit
- `ovbb_plot.jl`: script to retrieve useful information (e.g. best fit, upper 90% CI signal limit, ...) from a performed fit, in order to be later plotted with `py_scripts/get_ovbb_results.py`
- `py_scripts/get_fit_results.py`: script to print on terminal the main results obtained via a fit (e.g. global mode, 68% CI interval, ...) for the free parameters
- `py_scripts/plot_marginalized_posteriors.py`: script to plot marginalized posterior pdfs for comparing different signal prior choices, different background shape modeling choices, etc.
- `py_scripts/plot_toys_results.py`: script to retrieve mean, median and 68% CI interval out of the distribution of 90% CI upper limits obtained through sensitivity studies over toy data generated under the no-signal hypothesis
- `py_scripts/read_bi.py`: script to retrieve the total average BI in a fit where multiple BIs were included. The average is taken as the exposure-weighted average of each single BIs.

> [!NOTE]
> These scripts are intended as a starting point for more advanced analyses. Some of these are pretty hardcoded and related to the type of fits you previously ran.