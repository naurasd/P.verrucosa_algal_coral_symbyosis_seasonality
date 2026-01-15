The stable isotope data used for the SIBER analysis were visualized using the siber_ggplot.R script with siber_data_plot.csv as input data. 

The SIBER analysis was run with the SIBER_script_total.R script. Input data are stored in the siber_data_total.csv and iso_niche_data.csv files. See the script for details. This includes the estimation of host tissue and Symbiodiniaceae standard ellipse areas corrected for sample size (SEAc) and the calculation of SEAc overlap. Further, the distance between SEAc centroids is calculated and the residual permutation test to assess significance of trophic niche separation is carried out. This script relies on a source code by Turner et al. (2010; see README.md of main repository page), which is the Turner.et.al.ecology.source.r script. It is read into R in the SIBER_script_total.R script.

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
