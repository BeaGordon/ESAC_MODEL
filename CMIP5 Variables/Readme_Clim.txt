The data included in the CMIP5 processing folder include:
1. Data that are processed to run each of the scripts, raw data are available from sources provided in the main manuscript
 a. Data included are all SWE and Q data obtained from the codes in this repository labelled as such.
2. Lookup Tables required to run the codes included
3. Codes to process the reservoir inflow and unregulated streamflow for the main manuscript. Brief descriptions are provided
 a. Step 1: CMIP5_variable extraction can be used to access the raw data from the VIC CMIP5 hydrologic outputs including SWE and runoff
 b. Step 2: In the case that MizuRoute data are used, the MizuRoute_processing script can be used to extract data and prepare them for the model.
 c. Step 3: Observed climatology can be obtained for SWE and Q using the historical record using the codes 'Observed_climatology_X" 
 d. Step 5: Simulated historical climatology can be obtained for SWE and Q over the historical period using the codes "Simulated_climatology_X"
 e. Step 6: SWE and Q can be bias corrected following the methdology outlined in the main manuscript and supplemental information for Gordonetal using the code "Bias_correction"