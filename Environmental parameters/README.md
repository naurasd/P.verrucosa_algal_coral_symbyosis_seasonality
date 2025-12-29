# Input data, analysis and figures of environmental parameters
## Input data
The file `Environmental Data file Tilstra et al.xlsx` is used for all of the here presented R scripts. It contains five sheets:
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
 
