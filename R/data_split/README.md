# RAATD
## data_split

This folder holds code, configuration files and dates for splitting tracks into breeding stages.

Stage dates are defined in configuration files ('config' folder). Individual-specific dates and times defining breeding stages are output into the 'stage_dates' folder by running 'trip_splitting_script.R' in the 'code' folder.

Based on inspection of the output by species experts, some changes are made using the file 'stageDatesManualCorrections.R' in the 'code' folder.

For GLS tracks of BBAL from Kerguelen:
Data around 21 days before and after the equinoxes are discarded, as these locations were deemed unreliable after visaul inspection by Karine Delord (who curates those data, and estimated the GLS locations orginally). This is implemented in 'trip_splitting_scriptBBAL.R'. Of particular concern were the apparent trips to Prydz bay just before returning to Kerguelen. In a further processing stage (implemented in 'stageDatesManualCorrections.R') , post-breeding trips which straddle the discarded period around the austral spring equinox were truncated at 1 September to avoid extreme extrapolations. These data around the equinoxes must again be filtered out prior to fitting stage-specific SSMs.

Stages are not defined for HUWH or CRAS. All WESE individuals are post-moult and there is no stage dates file. For each of these species there is a config file, but the 'stages' dataframe is NULL.
