; CMIP5 MMLE
; Creates netCDF file of running multi-year, regional (user-specified) timeseries
; of 4-moment statistics (mean, variance, skewness, kurtosis).
; This is for single point, monthly tas.
;
; LLandrum
; Feb, 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
; ================================================;
begin
  imonth_plot = 9                      ; 0-11
  var_name = "tas"

;  mname ="GFDL-CM3"
;  dirmname = "gfdl_cm3"
;  ncases = 20
;
  mname = "CSIRO-Mk3-6-0"
  dirmname = "csiro_mk36"
  ncases = 29
;
;  mname = "CanESM2"
;  dirmname = "canesm2"
;  ncases = 50
;
;  mname = "GFDL-ESM2M"
;  dirmname = "gfdl_esm2m"
;  ncases = 30
;
;  mname = "MPI-ESM"
;  dirmname = "mpi"
;  ncases = 100

  case_numbers = ispan(1,ncases,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  case_yr1 = 1920
  if (mname.eq."CSIRO-Mk3-6-0") then
     case_yr1 = 1850
  else 
     if ((mname.eq."CanESM2").or.(mname.eq."GFDL-ESM2M")) then
        case_yr1 = 1950
     end if
  end if
  case_yr2 = 2100

  nyr_stats = 10

;  point_name = "Chukchi Sea"
;  point_short = "Chukchi"
;  point_lat = (/75./)
;  point_lon = (/190./)                      ; deg E

  point_name = "Fairbanks, AK"
  point_short = "Fairbanks"
  point_lat = (/64.8378/)
  point_lon = (/147.7164/)                      ; deg W

  point_lat@units = "degN"			; Russian points (Lena, Yakutsk) are in degE, others degW
  if (point_short.eq."Chukchi")
     point_lon@units = "degE"
  else
     point_lon@units = "degW"
  end if

  if (point_lon@units.eq."degW") then
     point_lon = 360. - point_lon
     point_lon@units = "degE"
  end if
  if (point_lon.lt.0.) then
     point_lon = point_lon + 360.
  end if

  datadir = "/glade/collections/cdg/data/CLIVAR_LE/"+dirmname+"_lens/Amon/"+var_name+"/"
  outdir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"

  stats_yr1 = 1920				; beginning and end yrs for plotting LE histogram - so we can compare from different time periods
;  stats_yr1 = 1950
  stats_yr2 = 2091                              ; must be at least nyr_stats-1 less than case_yr2
;  stats_yr2 = 2090
  last_yr_data = stats_yr2 + nyr_stats - 1
  nyr_LE = (stats_yr2 - stats_yr1 + 1) + (nyr_stats - 1)
  nyr_LE_stat = stats_yr2 - stats_yr1 + 1

  fdata = datadir+var_name+"_Amon_"+mname+"_historical_rcp85_r"+case_numbers+"i1p1_"+case_yr1+"01-"+case_yr2+"12.nc"
  fgrid = fdata(0)
 
  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
  ndays_per_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  if (imonth_plot.gt.0) then
      iday1 = dim_sum(ndays_per_month(:imonth_plot-1))
  else
      iday1 = 0
  end if
  iday2 = iday1 + ndays_per_month(imonth_plot) - 1
  ndays = ndays_per_month(imonth_plot)
 
  fstats = outdir+mname+"."+point_short+"."+nyr_stats+"yr_stats."+month_names(imonth_plot)+"_monthly_"+var_name+"."+stats_yr1+"-"+stats_yr2+".nc"
;------------------------------------------------------------------------------
; get indices to limit amount of data loaded...
;------------------------------------------------------------------------------
  f      = addfile (fgrid, "r")              ; get grid info and create variable arrays
  tmplat = f->lat
  tmplon = f->lon
  ilat_pt = ind_nearest_coord(point_lat,tmplat,0)
  jlon_pt = ind_nearest_coord(point_lon,tmplon,0)
  lon = tmplon(jlon_pt)
  lat = tmplat(ilat_pt)
  delete(tmplat)
  delete(tmplon)

  lat_round = (round(100.*lat,0))/100.
  lon_round = (round(100.*lon,0))/100.

  point_loc_txt = lat_round+"~S~o~N~N, "+lon_round+"~S~o~N~E"

;________________START HERE TUES JAN 14_________________________________________

  var_LE = new((/nyr_LE,ncases/),float)
  stat_LE = new((/nyr_LE_stat,4/),float)
  min_LE_annual = new(nyr_LE_stat,float)
  max_LE_annual = new(nyr_LE_stat,float)
  min_LE_decadal = new(nyr_LE_stat,float)
  max_LE_decadal = new(nyr_LE_stat,float)
;---------------------------------------------------------------------------------
; SKIPPING FOR NOW....
; if ICEFRAC then create a mask to only count grid cells that contain ice at 
; some time in comparison record (LE1)
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
; read in LE data
;---------------------------------------------------------------------------------
  iy1 = stats_yr1 - case_yr1
  j1 = iy1*12 + imonth_plot
  do i = 0,ncases-1
     f = addfile(fdata(i),"r")
     var_LE(:,i) = (/ f->$var_name$(j1::12,ilat_pt,jlon_pt) /)
  end do

  print("Done reading in data, now calculating statistics")

  do i=0,nyr_LE_stat - 1
     i2 = i + nyr_stats - 1 
     LE_var_1d = ndtooned(var_LE(i:i2,:))
     stat_LE(i,:) = dim_stat4(LE_var_1d)
     min_LE_decadal(i) = min(var_LE(i:i2,:))				; save min and max by year....
     max_LE_decadal(i) = max(var_LE(i:i2,:))
     min_LE_annual(i) = min(var_LE(i,:))                          ; save min and max by year....
     max_LE_annual(i) = max(var_LE(i,:))
  end do

  statLabel = (/"Mean","Var","Skew","Kurtosis"/)
  statLabel!0 = "moment"
  moment = ispan(1,4,1)
  statLabel&moment = moment

  f      = addfile(fdata(0), "r")
  tmp = f->$var_name$(0,0,0)
  stat_LE@units = tmp@units
  min_LE_decadal@units = tmp@units
  max_LE_decadal@units = tmp@units
  min_LE_annual@units = tmp@units
  max_LE_annual@units = tmp@units
  delete(tmp)

  stat_LE!0 = "year"
  year = ispan(stats_yr1,stats_yr2,1)
  stat_LE&year = year
  stat_LE@long_name = "4 moments of statistics for 10 years of monthly "+month_names(imonth_plot)+" "+var_name
  stat_LE!1 = "moment"
  stat_LE&moment = moment

  min_LE_decadal!0 = "year"
  min_LE_decadal&year = year
  max_LE_decadal!0 = "year"
  max_LE_decadal&year = year
  min_LE_decadal@long_name = "decadal min monthly "+month_names(imonth_plot)+" "+var_name
  max_LE_decadal@long_name = "decadal max monthly "+month_names(imonth_plot)+" "+var_name

  min_LE_annual!0 = "year"
  min_LE_annual&year = year
  max_LE_annual!0 = "year"
  max_LE_annual&year = year
  min_LE_annual@long_name = "annual min monthly "+month_names(imonth_plot)+" "+var_name
  max_LE_annual@long_name = "annual max monthly "+month_names(imonth_plot)+" "+var_name

   system("/bin/rm -f "+fstats)
   setfileoption("nc","Format","NetCDF4")
   fout = addfile(fstats,"c")
   setfileoption(fout,"DefineMode",True)

   file_title = "LE ensemble: 4 moments of statistics and yearly min, max for "+point_name+" "+nyr_stats+" monthly "+month_names(imonth_plot)+" "+var_name
   fileAtt = True
   fileAtt@creation_date = systemfunc("date")
   fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
   fileAtt@title = file_title
   fileattdef(fout,fileAtt)

   dimNames = (/"year","moment"/)
   dimSizes = (/nyr_LE_stat,4/)
   dimUnlim = (/True,False/)
   filedimdef(fout,dimNames,dimSizes,dimUnlim)

   filevardef(fout,"stat_LE",typeof(stat_LE),getvardims(stat_LE))
   filevardef(fout,"min_LE_annual",typeof(min_LE_annual),getvardims(min_LE_annual))
   filevardef(fout,"max_LE_annual",typeof(max_LE_annual),getvardims(max_LE_annual))
   filevardef(fout,"min_LE_decadal",typeof(min_LE_decadal),getvardims(min_LE_decadal))
   filevardef(fout,"max_LE_decadal",typeof(max_LE_decadal),getvardims(max_LE_decadal))
   filevardef(fout,"statLabel",typeof(statLabel),getvardims(statLabel))

   filevarattdef(fout,"stat_LE",stat_LE)
   filevarattdef(fout,"min_LE_annual",min_LE_annual)
   filevarattdef(fout,"max_LE_annual",max_LE_annual)
   filevarattdef(fout,"min_LE_decadal",min_LE_decadal)
   filevarattdef(fout,"max_LE_decadal",max_LE_decadal)

   fout->stat_LE = (/ stat_LE /)
   fout->min_LE_annual = (/ min_LE_annual /)
   fout->max_LE_annual = (/ max_LE_annual /)
   fout->min_LE_decadal = (/ min_LE_decadal /)
   fout->max_LE_decadal = (/ max_LE_decadal /)
   fout->statLabel = (/ statLabel /)

   fout->point_name = point_name
   fout->point_short = point_short
   fout->point_lat = point_lat
   fout->point_lon = point_lon
   fout->lat = lat
   fout->lon = lon
end

