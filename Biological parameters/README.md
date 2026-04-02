# Input data and code for statistical analysis and plotting of biological parameters and biological-environmental correlation analysis
The file [Biological_Parameters.xlsx](Biological_Parameters.xlsx) contains input data formatted for the statistical analysis of seasonal differentiation of biological parameters. It is used in the `R` script [Biological_parameters_ANOVA_KW_Welch.R](Biological_parameters_ANOVA_KW_Welch.R) which performs the following analyses for each parameter:
 - Shapiro-Wilk test to check for normality of data distribution (normally distributed parameters were tested for homogeneity of variances using Levene’s test)
 - Parameters that were normally distributed and homoscedastic, either in their raw form or after transformation (logit transformation for elemental percentages, and log(x+1) or log(-x+1) transformation for all other parameters), were analyzed using a one-way ANOVA followed by Tukey’s post-hoc test (for ANOVA p < 0.05)
 - When parameters were normally distributed but violated homogeneity of variances, Welch’s ANOVA was used instead, followed by Games–Howell post-hoc comparisons
 - Parameters that remained non-normal after transformation were analyzed using non-parametric Kruskal-Wallis rank-sum tests with Dunn’s post-hoc test and Bonferroni correction (for Kruskal-Wallis p < 0.05)

The file [Biological Data file Tilstra et al.xlsx](Biological%20Data%20file%20Tilstra%20et%20al.xlsx) is used as input data for the remaining `R` scripts:
- The sheet `symbionts_chla_MI_nifH` contains the following data:
  - Algal symbiont density in cells x10<sup>6</sup> cm<sup>-2</sup>
  - Mitotic index of algal symbiont cells in %
  - Chlorophyll _a_ content of algal symbionts in pg cell<sup>-1</sup>
  - Relative _nifH_ gene abundance as qPCR Ct values and as fold change relative to the sample with the lowest abundance
  
- The sheet `elemental_and_isotopic_data` contains the following data separately for coral host tissue and algal symbionts:
  - δ<sup>15</sup>N and δ<sup>13</sup>C in ‰
  - organic carbon and nitrogen content in % and the calculated carbon:nitrogen ratio
  - the calculated difference between coral host tissue and algal symbiont δ<sup>15</sup>N and δ<sup>13</sup>C in ‰ for each sample
    
- The sheet `bio + env data` contains all data of biological and environmental parameters formatted for Pearson correlation analysis
  
- The sheet `season_mean_sem_bio` contains the seasonal mean and standard error of the mean (SEM) for all biological parameters
 
The `R` file [Fig_2B-D.R](Fig_2B-D.R) is used to plot the biological data as shown in Figure 2B-D in the mansucript. It denotes statistical differences between seasons for each parameter as assessed in [Biological_parameters_ANOVA_KW_Welch.R](Biological_parameters_ANOVA_KW_Welch.R).
 
The `R` file [Fig_2E_S5.R](Fig_2E_S5.R) is used to analyze and plot the biolgical and environmental data as shown in Figure 2E and Supplementary Figure S5 in the manuscript. It assesses correlation between environmental and biological paramaters and also plots the results.

The `R` file [Fig_3A_D_E_S6.R](Fig_3A_D_E_S6.R) is used to plot the biological data as shown in Figure 3A, 3D and 3E and Supplementary Figure S6 in the mansucript. It denotes statistical differences between seasons for each parameter as assessed in [Biological_parameters_ANOVA_KW_Welch.R](Biological_parameters_ANOVA_KW_Welch.R).

The `R` file [Fig_3B.R](Fig_3B.R) is used to perform linear regression analysis and plot the data shown in Figure 3B in the manuscript.

The `R` file [Fig_3C.R](Fig_3C.R) is used to perform Canonical Analysis of Principal coordinates (CAP) and ANOSIM analysis and plot the data shown in Figure 3C in the manuscript.

The `R` file [Fig_4A_S8.R](Fig_4A_S8.R) is used to perform Stable Isotope Bayesian Ellipses (SIBER) Analysis for stable isotope data and to plot the data shown in Figure 4A and Supplemenmtary Figure S8 in the manuscript.

The `R` file [Fig_4B.R](Fig_4B.R) is used to perform linear regression analysis and plot the data shown in Figure 4B in the manuscript.

The `R` file [Fig_S4.R](Fig_S4.R) is used to plot the _nifH_ data as shown in Supplementary Figure S4 in the mansucript. It denotes statistical differences between seasons as assessed in [Biological_parameters_ANOVA_KW_Welch.R](Biological_parameters_ANOVA_KW_Welch.R).

The `R` file [Fig_S7.R](Fig_S7.R) is used to plot the isotopic data as shown in Supplementary Figure S7 in the mansucript. It denotes statistical differences between seasons for each parameter as assessed in [Biological_parameters_ANOVA_KW_Welch.R](Biological_parameters_ANOVA_KW_Welch.R).
