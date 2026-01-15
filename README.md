# Environmental variability drives carbon and nitrogen allocation in the coral-algal symbiosis

In this repository you will find data and `R` scripts used in the manuscript by Tilstra and Daraghmeh et al. (2026).

Data files and scripts can be found in [Environmental parameters](Environmental%20parameters) and [Biological parameters](Biological%20parameters).

[Environmental parameters](Environmental%20parameters) contains `.xlsx` files with input data and associated `R` scripts for the analysis and (supplementary) figures / tables as used in the manuscript for the following environmental parameters:
- water temperature
- photosynthetically active radiation
- salinity
- dissolved oxygen
- dissolved nutrients (nitrite, nitrate, ammonium, phosphate)
- dissolved carbon

[Biological parameters](Biological%20parameters) contains `.xlsx` files with input data and associated `R` scripts for the analysis and (supplementary) figures / tables as used in the manuscript for the following biological parameters:
- Density, mitotic index and chlorophyll *a* content of algal symbiont cells in coral host tissue
- Elemental and stable isotope composition of nitrogen and organic carbon of coral host tissue and algal symbionts
- relative *nifH* gene copy numbers (diazotrophy)

[Biological parameters](Biological%20parameters) also contains data and code used for combined biological-environmental correlation analysis.

For our analysis, we used the following `R` packages / collections:

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). 
Welcome to the tidyverse.
Journal of Open Source Software, 4(43), 1686. doi: 10.21105/joss.01686. 

Pedersen, T.L. (2025).
ggforce: Accelerating 'ggplot2'. R package version 0.5.0.9000. [https://github.com/thomasp85/ggforce/](https://ggforce.data-imaginist.com)

The stable isotope data used for the SIBER analysis were visualized using the siber_ggplot.R script with siber_data_plot.csv as input data. 

The SIBER analysis was run with the SIBER_script_total.R script. Input data are stored in the siber_data_total.csv and iso_niche_data.csv files. See the script for details. This includes the estimation of host tissue and Symbiodiniaceae standard ellipse areas corrected for sample size (SEAc) and the calculation of SEAc overlap. Further, the distance between SEAc centroids is calculated and the residual permutation test to assess significance of trophic niche separation is carried out. This script relies on a source code by Turner et al. (2010; see below), which is the Turner.et.al.ecology.source.r script. It is read into R in the SIBER_script_total.R script.

For the stable isotope visualization and analysis with SIBER, we used parts of the scripts developed by:

Conti-Jerpe, I.E., Thompson, P.D., Wong, C.W.M., Oliveira, N.L., Duprey, N.N., Moynihan, M.A., BAker, D.M. 2020. 
Trophic strategy and bleaching resistance in reef-building corals. 
Science Advances, 6: eaaz5443. 
DOI: https://doi.org/10.1126/sciadv.aaz5443

Their scripts in turn are based on scripts, an `R` package (SIBER) and theory by: 

Jackson, A.L., Parnell, A.C., Inger R., & Bearhop, S. 2011. 
Comparing isotopic niche widths among and within communities: SIBER - Stable Isotope Bayesian Ellipses in R. 
Journal of Animal Ecology, 80: 595-602. 
DOI: https://doi.org/10.1111/j.1365-2656.2011.01806.x 

Turner, T.F., Collyer, M.L. & Krabbenhoft, T.J. 2010. 
A general hypothesis-testing framework for stable isotope ratios in ecological studies. 
Ecology, 91: 2227-2233. 
DOI: https://doi.org/10.1890/09-1454.1
Journal of Open Source Software, 4(43), 1686. doi: 10.21105/joss.01686. 
