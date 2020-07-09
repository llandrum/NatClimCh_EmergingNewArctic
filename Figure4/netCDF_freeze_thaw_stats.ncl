; Calculates statistics (by decade) for large ensemble freeze-thaw
; file and saves to netCDF.
; LLandrum
; Apr 2020
;
;----------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
; ==============================================================
; User defined parameters that specify region of globe to consider
; ==============================================================
  mname = "CanESM2"
  dirmname = "canesm2"
  ncases = 50

  case_numbers = ispan(1,ncases,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  datadir = "/glade/p/cesm/pcwg/landrum/LE/precip/"
  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/precip/"
  outdir = "/glade/scratch/landrum/CMIP5/atm/precip/"+mname+"/"

  min_snow_percent = 55.                        ; percent of total precip that must be snow for "freeze up onset" definition

  yr1 = 1950
  yr2 = 2100

  nyr_stats = 10
  yr2_stats = yr2 - nyr_stats + 1
  nyr = yr2_stats - yr1 + 1
  year = ispan(yr1,yr2_stats,1)
  sstat_dates = sprinti("%0.4i",yr1) + "-"+sprinti("%0.4i",yr2_stats)

  ffreeze = datadir+mname+"."+ncases+"sim.NH.freeze_thaw_day.min"+min_snow_percent+"percent."+yr1+"-"+yr2+".nc"
  fstats = outdir+mname+"."+ncases+"sim.NH."+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+yr1+"-"+yr2_stats+".nc"

;------------------------------------------------------------------------------
; read in data
;------------------------------------------------------------------------------
  f      = addfile (ffreeze, "r")
  lat = f->lat
  lon = f->lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  first_rain = int2flt(f->first_rain)
  last_rain = int2flt(f->last_rain)
  first_snow = int2flt(f->first_snow)
  thaw_onset = int2flt(f->thaw_onset)
  freeze_onset = int2flt(f->freeze_onset)

  rain_season_length = last_rain
  rain_season_length = last_rain - first_rain

  first_rain_stats = new((/nyr,4,nlat,nlon/),float)
  first_rain_stats!0 = "year"
  first_rain_stats&year = year
  first_rain_stats!1 = "moment"
  first_rain_stats!2 = "lat"
  first_rain_stats&lat = lat
  first_rain_stats!3 = "lon"
  first_rain_stats&lon = lon
  first_rain_stats&moment = ispan(1,4,1)

  last_rain_stats = first_rain_stats
  first_snow_stats = first_rain_stats
  thaw_onset_stats = first_rain_stats
  freeze_onset_stats = first_rain_stats
  rain_season_length_stats = first_rain_stats

  first_rain_min_annual = new((/nyr,nlat,nlon/),float)
  first_rain_min_annual!0 = "year"
  first_rain_min_annual&year = year
  first_rain_min_annual!1 = "lat"
  first_rain_min_annual&lat = lat
  first_rain_min_annual!2 = "lon"
  first_rain_min_annual&lon = lon

  last_rain_min_annual = first_rain_min_annual
  rain_season_length_min_annual = first_rain_min_annual
  first_snow_min_annual = first_rain_min_annual
  thaw_onset_min_annual = first_rain_min_annual
  freeze_onset_min_annual = first_rain_min_annual

  first_rain_max_annual = first_rain_min_annual
  last_rain_max_annual = first_rain_min_annual
  rain_season_length_max_annual = first_rain_min_annual

  first_snow_max_annual = first_rain_min_annual
  thaw_onset_max_annual = first_rain_min_annual
  freeze_onset_max_annual = first_rain_min_annual

  first_rain_min_decadal = first_rain_min_annual
  last_rain_min_decadal = first_rain_min_annual
  rain_season_length_min_decadal = first_rain_min_annual

  first_snow_min_decadal = first_rain_min_annual
  thaw_onset_min_decadal = first_rain_min_decadal
  freeze_onset_min_decadal = first_rain_min_decadal

  first_rain_max_decadal = first_rain_min_decadal
  last_rain_max_decadal = first_rain_min_decadal
  rain_season_length_max_decadal = first_rain_min_annual

  first_snow_max_decadal = first_rain_min_decadal
  thaw_onset_max_decadal = first_rain_min_decadal
  freeze_onset_max_decadal = first_rain_min_decadal

  dimNames = (/"year","moment","lat","lon"/)
  dimSizes = (/nyr,4,nlat,nlon/)
  dimUnlim = (/True,False,False,False/)

  statLabel = (/"Mean","Var","Skew","Kurtosis"/)
  statLabel!0 = "moment"
  moment = ispan(1,4,1)
  statLabel&moment = moment

  do i=0,nyr-1
     i2 = i+nyr_stats-1
     first_rain3d = onedtond(ndtooned(first_rain(i:i2,:,:,:)),(/nyr_stats*ncases,nlat,nlon/))
     last_rain3d = onedtond(ndtooned(last_rain(i:i2,:,:,:)),(/nyr_stats*ncases,nlat,nlon/))
     first_snow3d = onedtond(ndtooned(first_snow(i:i2,:,:,:)),(/nyr_stats*ncases,nlat,nlon/))
     thaw_onset3d = onedtond(ndtooned(thaw_onset(i:i2,:,:,:)),(/nyr_stats*ncases,nlat,nlon/))
     freeze_onset3d = onedtond(ndtooned(freeze_onset(i:i2,:,:,:)),(/nyr_stats*ncases,nlat,nlon/))
     rain_season_length3d = onedtond(ndtooned(rain_season_length(i:i2,:,:,:)),(/nyr_stats*ncases,nlat,nlon/))

     first_rain_stats(i,:,:,:) = (/ dim_stat4_n(first_rain3d,0) /)
     first_rain_min_annual(i,:,:) = (/ dim_min_n(first_rain(i,:,:,:),0) /)
     first_rain_max_annual(i,:,:) = (/ dim_max_n(first_rain(i,:,:,:),0) /)
     first_rain_min_decadal(i,:,:) = (/ dim_min_n(first_rain3d,0) /)
     first_rain_max_decadal(i,:,:) = (/ dim_max_n(first_rain3d,0) /)

     last_rain_stats(i,:,:,:) = (/ dim_stat4_n(last_rain3d,0) /)
     last_rain_min_annual(i,:,:) = (/ dim_min_n(last_rain(i,:,:,:),0) /)
     last_rain_max_annual(i,:,:) = (/ dim_max_n(last_rain(i,:,:,:),0) /)
     last_rain_min_decadal(i,:,:) = (/ dim_min_n(last_rain3d,0) /)
     last_rain_max_decadal(i,:,:) = (/ dim_max_n(last_rain3d,0) /)

     rain_season_length_stats(i,:,:,:) = (/ dim_stat4_n(rain_season_length3d,0) /)
     rain_season_length_min_annual(i,:,:) = (/ dim_min_n(rain_season_length(i,:,:,:),0) /)
     rain_season_length_max_annual(i,:,:) = (/ dim_max_n(rain_season_length(i,:,:,:),0) /)
     rain_season_length_min_decadal(i,:,:) = (/ dim_min_n(rain_season_length3d,0) /)
     rain_season_length_max_decadal(i,:,:) = (/ dim_max_n(rain_season_length3d,0) /)

     first_snow_stats(i,:,:,:) = (/ dim_stat4_n(first_snow3d,0) /)
     first_snow_min_annual(i,:,:) = (/ dim_min_n(first_snow(i,:,:,:),0) /)
     first_snow_max_annual(i,:,:) = (/ dim_max_n(first_snow(i,:,:,:),0) /)
     first_snow_min_decadal(i,:,:) = (/ dim_min_n(first_snow3d,0) /)
     first_snow_max_decadal(i,:,:) = (/ dim_max_n(first_snow3d,0) /)

     thaw_onset_stats(i,:,:,:) = (/ dim_stat4_n(thaw_onset3d,0) /) 
     thaw_onset_min_annual(i,:,:) = (/ dim_min_n(thaw_onset(i,:,:,:),0) /) 
     thaw_onset_max_annual(i,:,:) = (/ dim_max_n(thaw_onset(i,:,:,:),0) /) 
     thaw_onset_min_decadal(i,:,:) = (/ dim_min_n(thaw_onset3d,0) /) 
     thaw_onset_max_decadal(i,:,:) = (/ dim_max_n(thaw_onset3d,0) /) 

     freeze_onset_stats(i,:,:,:) = (/ dim_stat4_n(freeze_onset3d,0) /)
     freeze_onset_min_annual(i,:,:) = (/ dim_min_n(freeze_onset(i,:,:,:),0) /)
     freeze_onset_max_annual(i,:,:) = (/ dim_max_n(freeze_onset(i,:,:,:),0) /)
     freeze_onset_min_decadal(i,:,:) = (/ dim_min_n(freeze_onset3d,0) /)
     freeze_onset_max_decadal(i,:,:) = (/ dim_max_n(freeze_onset3d,0) /)

  end do 						; end nstats
;-----------------------------------------------------------------------------------------
;    save netCDF
;-----------------------------------------------------------------------------------------
  system("/bin/rm -f "+fstats)
  setfileoption("nc","Format","NetCDF4")
  fout = addfile(fstats,"c")
  setfileoption(fout,"DefineMode",True)

  file_title = "LE ensemble: 4 moments of statistics and min, max for annual first measurable rain, snow, thaw onset and freeze onset"
  fileAtt = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
  fileAtt@title = file_title
  fileattdef(fout,fileAtt)

  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout,"lat",typeof(lat),getvardims(lat))
  filevardef(fout,"lon",typeof(lon),getvardims(lon))
  filevardef(fout,"first_rain_stats",typeof(first_rain_stats),getvardims(first_rain_stats))
  filevardef(fout,"first_rain_min_annual",typeof(first_rain_min_annual),getvardims(first_rain_min_annual))
  filevardef(fout,"first_rain_max_annual",typeof(first_rain_max_annual),getvardims(first_rain_max_annual))
  filevardef(fout,"first_rain_min_decadal",typeof(first_rain_min_decadal),getvardims(first_rain_min_decadal))
  filevardef(fout,"first_rain_max_decadal",typeof(first_rain_max_decadal),getvardims(first_rain_max_decadal))
  filevardef(fout,"last_rain_stats",typeof(last_rain_stats),getvardims(last_rain_stats))
  filevardef(fout,"last_rain_min_annual",typeof(last_rain_min_annual),getvardims(last_rain_min_annual))
  filevardef(fout,"last_rain_max_annual",typeof(last_rain_max_annual),getvardims(last_rain_max_annual))
  filevardef(fout,"last_rain_min_decadal",typeof(last_rain_min_decadal),getvardims(last_rain_min_decadal))
  filevardef(fout,"last_rain_max_decadal",typeof(last_rain_max_decadal),getvardims(last_rain_max_decadal))
  filevardef(fout,"rain_season_length_stats",typeof(rain_season_length_stats),getvardims(rain_season_length_stats))
  filevardef(fout,"rain_season_length_min_annual",typeof(rain_season_length_min_annual),getvardims(rain_season_length_min_annual))
  filevardef(fout,"rain_season_length_max_annual",typeof(rain_season_length_max_annual),getvardims(rain_season_length_max_annual))
  filevardef(fout,"rain_season_length_min_decadal",typeof(rain_season_length_min_decadal),getvardims(rain_season_length_min_decadal))
  filevardef(fout,"rain_season_length_max_decadal",typeof(rain_season_length_max_decadal),getvardims(rain_season_length_max_decadal))
  filevardef(fout,"first_snow_stats",typeof(first_snow_stats),getvardims(first_snow_stats))
  filevardef(fout,"first_snow_min_annual",typeof(first_snow_min_annual),getvardims(first_snow_min_annual))
  filevardef(fout,"first_snow_max_annual",typeof(first_snow_max_annual),getvardims(first_snow_max_annual))
  filevardef(fout,"first_snow_min_decadal",typeof(first_snow_min_decadal),getvardims(first_snow_min_decadal))
  filevardef(fout,"first_snow_max_decadal",typeof(first_snow_max_decadal),getvardims(first_snow_max_decadal))
  filevardef(fout,"thaw_onset_stats",typeof(thaw_onset_stats),getvardims(thaw_onset_stats))
  filevardef(fout,"thaw_onset_min_annual",typeof(thaw_onset_min_annual),getvardims(thaw_onset_min_annual))
  filevardef(fout,"thaw_onset_max_annual",typeof(thaw_onset_max_annual),getvardims(thaw_onset_max_annual))
  filevardef(fout,"thaw_onset_min_decadal",typeof(thaw_onset_min_decadal),getvardims(thaw_onset_min_decadal))
  filevardef(fout,"thaw_onset_max_decadal",typeof(thaw_onset_max_decadal),getvardims(thaw_onset_max_decadal))
  filevardef(fout,"freeze_onset_stats",typeof(freeze_onset_stats),getvardims(freeze_onset_stats))
  filevardef(fout,"freeze_onset_min_annual",typeof(freeze_onset_min_annual),getvardims(freeze_onset_min_annual))
  filevardef(fout,"freeze_onset_max_annual",typeof(freeze_onset_max_annual),getvardims(freeze_onset_max_annual))
  filevardef(fout,"freeze_onset_min_decadal",typeof(freeze_onset_min_decadal),getvardims(freeze_onset_min_decadal))
  filevardef(fout,"freeze_onset_max_decadal",typeof(freeze_onset_max_decadal),getvardims(freeze_onset_max_decadal))

  filevarattdef(fout,"lat",lat)
  filevarattdef(fout,"lon",lon)
  filevarattdef(fout,"first_rain_stats",first_rain_stats)
  filevarattdef(fout,"first_rain_min_annual",first_rain_min_annual)
  filevarattdef(fout,"first_rain_max_annual",first_rain_max_annual)
  filevarattdef(fout,"first_rain_min_decadal",first_rain_min_decadal)
  filevarattdef(fout,"first_rain_max_decadal",first_rain_max_decadal)

  filevarattdef(fout,"last_rain_stats",last_rain_stats)
  filevarattdef(fout,"last_rain_min_annual",last_rain_min_annual)
  filevarattdef(fout,"last_rain_max_annual",last_rain_max_annual)
  filevarattdef(fout,"last_rain_min_decadal",last_rain_min_decadal)
  filevarattdef(fout,"last_rain_max_decadal",last_rain_max_decadal)

  filevarattdef(fout,"rain_season_length_stats",rain_season_length_stats)
  filevarattdef(fout,"rain_season_length_min_annual",rain_season_length_min_annual)
  filevarattdef(fout,"rain_season_length_max_annual",rain_season_length_max_annual)
  filevarattdef(fout,"rain_season_length_min_decadal",rain_season_length_min_decadal)
  filevarattdef(fout,"rain_season_length_max_decadal",rain_season_length_max_decadal)

  filevarattdef(fout,"first_snow_stats",first_snow_stats)
  filevarattdef(fout,"first_snow_min_annual",first_snow_min_annual)
  filevarattdef(fout,"first_snow_max_annual",first_snow_max_annual)
  filevarattdef(fout,"first_snow_min_decadal",first_snow_min_decadal)
  filevarattdef(fout,"first_snow_max_decadal",first_snow_max_decadal)
  filevarattdef(fout,"thaw_onset_stats",thaw_onset_stats)
  filevarattdef(fout,"thaw_onset_min_annual",thaw_onset_min_annual)
  filevarattdef(fout,"thaw_onset_max_annual",thaw_onset_max_annual)
  filevarattdef(fout,"thaw_onset_min_decadal",thaw_onset_min_decadal)
  filevarattdef(fout,"thaw_onset_max_decadal",thaw_onset_max_decadal)
  filevarattdef(fout,"freeze_onset_stats",freeze_onset_stats)
  filevarattdef(fout,"freeze_onset_min_annual",freeze_onset_min_annual)
  filevarattdef(fout,"freeze_onset_max_annual",freeze_onset_max_annual)
  filevarattdef(fout,"freeze_onset_min_decadal",freeze_onset_min_decadal)
  filevarattdef(fout,"freeze_onset_max_decadal",freeze_onset_max_decadal)

  fout->lat = (/ lat /)
  fout->lon = (/ lon /)
  fout->first_rain_stats = (/ first_rain_stats /)
  fout->first_rain_min_annual = (/ first_rain_min_annual /)
  fout->first_rain_max_annual = (/ first_rain_max_annual /)
  fout->first_rain_min_decadal = (/ first_rain_min_decadal /)
  fout->first_rain_max_decadal = (/ first_rain_max_decadal /)
  fout->last_rain_stats = (/ last_rain_stats /)
  fout->last_rain_min_annual = (/ last_rain_min_annual /)
  fout->last_rain_max_annual = (/ last_rain_max_annual /)
  fout->last_rain_min_decadal = (/ last_rain_min_decadal /)
  fout->last_rain_max_decadal = (/ last_rain_max_decadal /)
  fout->rain_season_length_stats = (/ rain_season_length_stats /)
  fout->rain_season_length_min_annual = (/ rain_season_length_min_annual /)
  fout->rain_season_length_max_annual = (/ rain_season_length_max_annual /)
  fout->rain_season_length_min_decadal = (/ rain_season_length_min_decadal /)
  fout->rain_season_length_max_decadal = (/ rain_season_length_max_decadal /)
  fout->first_snow_stats = (/ first_snow_stats /)
  fout->first_snow_min_annual = (/ first_snow_min_annual /)
  fout->first_snow_max_annual = (/ first_snow_max_annual /)
  fout->first_snow_min_decadal = (/ first_snow_min_decadal /)
  fout->first_snow_max_decadal = (/ first_snow_max_decadal /)
  fout->thaw_onset_stats = (/ thaw_onset_stats /)
  fout->thaw_onset_min_annual = (/ thaw_onset_min_annual /)
  fout->thaw_onset_max_annual = (/ thaw_onset_max_annual /)
  fout->thaw_onset_min_decadal = (/ thaw_onset_min_decadal /)
  fout->thaw_onset_max_decadal = (/ thaw_onset_max_decadal /)
  fout->freeze_onset_stats = (/ freeze_onset_stats /)
  fout->freeze_onset_min_annual = (/ freeze_onset_min_annual /)
  fout->freeze_onset_max_annual = (/ freeze_onset_max_annual /)
  fout->freeze_onset_min_decadal = (/ freeze_onset_min_decadal /)
  fout->freeze_onset_max_decadal = (/ freeze_onset_max_decadal /)
  fout->statLabel = (/ statLabel /)
  delete(fout)
end
