# NatClimCh_EmergingNewArctic
scripts for paper figures
Figures and data for Nature Climate Change, “Extremes become routine in an emerging New Arctic”, Landrum & Holland, 2020.

Below are outlined the data and scripts used for each of the figures. All scripts are in ncl. To use ncl, note that it does have a license:
https://www.ncl.ucar.edu/Download/NCL_source_license.shtml

Each figure has a master script (with figure number in title). To run:
ncl script_name
You will need to change locations for the data files.

Also, each figure has multiple steps to create the necessary data file from the CMIP5-MMLE. These intermediary data files are not saved here, although the scripts to make them are saved here and the CMIP5-MMLE data are publicly available. 


Figure 1
Figure1_rev.NH_SIE.CMIP_MMLE.ncl

Data:
monthly_NH_SIE_min_max_GFDL-CM3_LE.20_sim.1920-2100.nc
monthly_NH_SIE_min_max_MPI-ESM_LE.99_sim.1920-2005.nc
monthly_NH_SIE_min_max_CanESM2_LE.50_sim.1950-2100.nc
monthly_NH_SIE_min_max_GFDL-ESM2M_LE.30_sim.1950-2100.nc
monthly_NH_SIE_min_max_MPI-ESM_LE.99_sim.1920-2099.nc
monthly_NH_SIE_min_max_CESM_LE.40_sim.1920-2100.nc

These are created using:
min_max_monthly_NH_SIE_LE.ncl
from intermediary files created using:
pre_proc_hem_MHolland_NHregions_CMIP5_LE.ncl
(there is one file per simulation from each model)

Obs
NH_SIE_min_max_nsidc.from_monthly_sea_ice_index.1979-2019.nc
 Created using: 
netCDF_NH_SIE_ssmi_monthly_index.ncl
along with downloaded NSIDC monthly sea ice index.

Figure 2
Figure2_rev.monthly_tas_emerg.ncl

Using TAS emergence files:
emergence.2_std.CanESM2.50sim.NH_10yr_stats.1950_base.OCT_monthly_tas.nc 
emergence.2_std.GFDL-CM3.20sim.NH_10yr_stats.1950_base.OCT_monthly_tas.nc
emergence.2_std.GFDL-ESM2M.30sim.NH_10yr_stats.1950_base.OCT_monthly_tas.nc
emergence.2_std.MPI-ESM.100sim.NH_10yr_stats.1950_base.OCT_monthly_tas.nc
emergence.2_std.LE.40sim.NH_10yr_stats.1950_base.OCT_monthly_TREFHT.nc

emergence.2_std.CanESM2.50sim.NH_10yr_stats.1950_base.OCT_monthly_tas.nc 
emergence.2_std.GFDL-CM3.20sim.NH_10yr_stats.1950_base.FEB_monthly_tas.nc
emergence.2_std.GFDL-ESM2M.30sim.NH_10yr_stats.1950_base.FEB_monthly_tas.nc
emergence.2_std.MPI-ESM.100sim.NH_10yr_stats.1950_base.FEB_monthly_tas.nc
emergence.2_std.LE.40sim.NH_10yr_stats.1950_base.FEB_monthly_TREFHT.nc

These are made using:
netCDF_monthly_stats_tas.CMIP5_MMLE.ncl

And Chuckchi, Fairbanks files:
CanESM2.Chukchi.10yr_stats.OCT_monthly_tas.1950-2091.nc
GFDL-CM3.Chukchi.10yr_stats.OCT_monthly_tas.1950-2091.nc
GFDL-ESM2M.Chukchi.10yr_stats.OCT_monthly_tas.1950-2091.nc
MPI-ESM.Chukchi.10yr_stats.OCT_monthly_tas.1920-2090.nc
LE.Chukchi.10yr_stats.OCT_monthly_TREFHT.1920-2091.nc

CanESM2.Chukchi.10yr_stats.FEB_monthly_tas.1950-2091.nc
GFDL-CM3.Chukchi.10yr_stats.FEB_monthly_tas.1950-2091.nc
GFDL-ESM2M.Chukchi.10yr_stats.FEB_monthly_tas.1950-2091.nc
MPI-ESM.Chukchi.10yr_stats.FEB_monthly_tas.1920-2090.nc
LE.Chukchi.10yr_stats.FEB_monthly_TREFHT.1920-2091.nc

Made with:
netCDF_monthly_point_var.CMIP5_MMLE.ncl


Figure 3
CanESM2_regional_areavg.nc
CESM1-CAM5_regional_areavg.nc
GFDL-CM3_regional_areavg.nc
GFDL-ESM2M_regional_areavg.nc
MPI-ESM_regional_areavg.nc

CanESM2.50sim.tas_wgtavg.195001-210012.nc
CESM1-CAM5.40sim.tas_wgtavg.192001-210012.nc
CSIRO-Mk3-6-0.29sim.tas_wgtavg.192001-210012.nc
GFDL-CM3.20sim.tas_wgtavg.192001-210012.nc
GFDL-ESM2M.30sim.tas_wgtavg.195001-210012.nc
MPI-ESM.100sim.tas_wgtavg.192001-209912.nc

Ice_vol files – one for each simulation in each ensemble. 


created with :
regional_areas_CMIP5_LE.ncl
global_nh_sh_avg_tas_LE.ncl
pre_proc_hem_MHolland_NHregions_CMIP5_LE.ncl (in Figure 1 directory).

Figure 4
CanESM2.50sim.NH.10yr_stats.freeze_thaw.min60percent.1950-2091.nc
LE.40sim.NH.10yr_stats.freeze_thaw.min60percent.1920-2091.nc

emergence.2_std.CanESM2.50sim.NH_10yr_stats.freeze_thaw.min60percent.rain_season_length.nc
emergence.-2_std.CanESM2.50sim.NH_10yr_stats.freeze_thaw.min60percent.first_rain.nc
emergence.2_std.CanESM2.50sim.NH_10yr_stats.freeze_thaw.min60percent.last_rain.nc

emergence.2_std.1950-1959_baseclim.LE.40sim.NH_10yr_stats.freeze_thaw.min60percent.rain_season_length.nc
emergence.-2_std.1950-1959_baseclim.LE.40sim.NH_10yr_stats.freeze_thaw.min60percent.first_rain.nc
emergence.2_std.1950-1959_baseclim.LE.40sim.NH_10yr_stats.freeze_thaw.min60percent.last_rain.nc

created with:
netCDF_freeze_thaw_stats.ncl
map_year_emergence_freeze_thaw_LE.ncl

using:
CanESM2.50sim.NH.freeze_thaw_day.min60percent.1950-2100.nc
LE.40sim.NH.freeze_thaw_day.min60percent.1920-2100

Which in turn are made with
netCDF_freeze_thaw_day.CMIP_MMLE.ncl

Supplementary Figures
Figure 1
SuppFigure1.CMIP5_MMLE.min_max_NH_SIE.histograms.ncl

Uses same files used in (main) Figure 1.

Figure 2
SuppFigure2.CMIP5_MMLE.NH_SIE.ndecades.ncl
monthly_clim.10_yrstats.CanESM2.50_sim.1950-2091.nc
monthly_clim.10_yrstats.GFDL-CM3.20_sim.1920-2091.nc
monthly_clim.10_yrstats.GFDL-ESM2M.30_sim.1950-2091.nc
monthly_clim.10_yrstats.MPI-ESM.99_sim.1920-2090.nc
monthly_clim.10_yrstats.LE.40_sim.1920-2091.nc

climatologies are made with:
netCDF_annual_NH_clim.ncl
this uses ice_vol files (see Figure 3).

Figure 3
SuppFigure3.min_max_allsims_SIE_MMLE.ncl

Uses same monthly data as Supplementary Figure 2.
Also uses Walsh data:
ice_ext_MHolland_regions_Walsh_extended_ArcticSIC.1850-2017.nc

made with:
pre_proc_Walsh_extendedSIC_MHolland_NHregions.ncl
along with downloaded Walsh extended data set (from NSIDC).

Figure 4.
SuppFigure4.CMIP5_MMLE.monthly_tas_emerg.ncl
uses tas emergence files used in main Figure 2.

Figures 5-6
SuppFigures5-6.map_ndiff_decadal_monthly_var_ice.ncl
Uses same files created for main Figure 2

Figure 7-9
SuppFigures7-9.map_CMIP5_ndecade_drain.ncl
Uses same files created for main Figure 4

