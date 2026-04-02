# Input data and code for statistical analysis and plotting of biological parameters and biological-environmental correlation analysis
The file [Biological_Parameters.xlsx](Biological_Parameters.xlsx) contains input data formatted for the statistical analysis of seasonal differentiation of biological parameters. It is used in the `R` script [Biological_parameters_ANOVA_KW_Welch.R](Biological_parameters_ANOVA_KW_Welch.R) which performs the following analyses for each parameter:
 - Shapiro-Wilk test to check for normality of data distribution (normally distributed parameters were tested for homogeneity of variances using Levene’s test)
 - Parameters that were normally distributed and homoscedastic, either in their raw form or after transformation (logit transformation for elemental percentages, and log(x+1) or log(-x+1) transformation for all other parameters), were analyzed using a one-way ANOVA followed by Tukey’s post-hoc test (for ANOVA p < 0.05)
 - When parameters were normally distributed but violated homogeneity of variances, Welch’s ANOVA was used instead, followed by Games–Howell post-hoc comparisons
 - Parameters that remained non-normal after transformation were analyzed using non-parametric Kruskal-Wallis rank-sum tests with Dunn’s post-hoc test and Bonferroni correction (for Kruskal-Wallis p < 0.05)

The file [Biological Data file Tilstra et al.xlsx](Biological%20Data%20file%20Tilstra%20et%20al.xlsx) is used as input data for the remaining two `R`  scripts:
- The sheet `Temperature` contains all water temperature measurements in °C recorded every 30 or 60 minutes from 1 January 2017 to 31 January 2018
  
- The sheet `all_monthly` contains all values of all parameters which were measured/calculated monthly on three consecutive days from January 2017 to January 2018:
  - Salinity in Practical Salinity Units (PSU)
  - Photosynthetically Active Radiation (PAR) in umol photons m<sup>-2</sup>s<sup>-1</sup>
  - Dissolved nitrate in uM NO<sub>3</sub><sup>-</sup>
  - Dissolved nitrite in uM NO<sub>2</sub><sup>-</sup>
  - Dissolved ammonium in uM NH<sub>4</sub><sup>+</sup>
  - Dissolved inorganic phosphorus (DIP) in uM PO<sub>4</sub><sup>3-</sup>
  - Dissolved inorganic nitrogen (DIN) in uM N calculated as the sum of dissolved nitrate, nitrite and ammonium
  - The ratio of DIN:DIP as the calculated values of DIN divided by the respective DIP value
    
- The sheet `all_bimonthly` contains all values of all parameters which were measured/calculated bimonthly on two consecutive days from January 2017 to January 2018:
  - Dissolved oxygen (DO) in mg L<sup>-1</sup> (eight measurements each day) 
  - Dissolved organic carbon (DOC) in uM C L<sup>-1</sup> (eight measurements each day)
  - Total dissolved nitrogen (TDN) in uM N L<sup>-1</sup> (eight measurements each day)
  - Dissolved organic nitrogen (DON) in uM N L<sup>-1</sup> calculated as TDN-DIN using the calculated DIN value of the same date
  - The ratio of DOC:DON as DOC divided by the respective calculated DON values
 
- The sheet `season_mean_sem_env` contains the seasonal mean and standard error of the mean (SEM) for all environmental parameters (except TDN, which was not further analysed):
  - For water temperature and monthly recorded parameters, all measurements of the months before and within the months of coral fragment sampling were taken into account (i.e., March and April 2017 for spring, July and August 2017 for summer, October and November 2017 for fall and December 2017 and January 2018 for winter)
  - For bimonthly recorded parameters, all measurements of the months before and after or within the month of coral fragment sampling were used (i.e., March and May 2017 for spring, July and September 2017 for summer, November 2017 for fall and January 2018 for winter)
 
- The sheet `month_mean_sem_env` contains the mean and standard error of the mean (SEM) for all environmental parameters (except TDN, which was not further analysed):
  - For water temperature and monthly recorded parameters, mean and SEM were calculated for all measurements of each month
  - For bimonthly recorded parameters, mean and SEM were calculated for all measurements of each of the months data was recorded in
 
The `R` file [Fig_1C-J_S1_S2.R](Fig_1C-J_S1_S2.R) is used to plot the environmental data as shown in Figure 1C-J and Supplementary Figures S1 and S2. It denotes statistical differences between seasons for each parameter as assessed in [Environmental_parameters_ANOVA_KW_Welch.R](Environmental_parameters_ANOVA_KW_Welch.R).
 
The `R` file [Fig_S3.R](Fig_S3.R) is used to analyze and plot the environmental data as shown in Supplementary Figure S3. It assesses collinearity among environmental paramaters that exhibited seasonal variability (see analysis above) through Pearson correlations and also plots the results. The script also performs hierarchical clustering and plots a corresponding dendrogram, further illustrating collinearity structure among environmental parameters.
 



The SIBER analysis was run with the SIBER_script_total.R script. Input data are stored in the siber_data_total.csv and iso_niche_data.csv files. See the script for details. This includes the estimation of host tissue and Symbiodiniaceae standard ellipse areas corrected for sample size (SEAc) and the calculation of SEAc overlap. Further, the distance between SEAc centroids is calculated and the residual permutation test to assess significance of trophic niche separation is carried out. This script relies on a source code by Turner et al. (2010; see README.md of main repository page), which is the Turner.et.al.ecology.source.r script. It is read into R in the SIBER_script_total.R script.
