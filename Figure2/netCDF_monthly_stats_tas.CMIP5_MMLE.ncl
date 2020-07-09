; March, 2018
; CMIP5 MMLE monthly data
; Calculates and saves to netCDF file (for plotting) 4-moment statistics
; (mean, variance, skewness and kurtosis) for atm variables.
;
; LLandrum
; February 2020
;----------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
; ==============================================================
; User defined parameters that specify region of globe to consider
; ==============================================================
  imonth_stats = 11                      ; 0-11
  var_name = "tas"

  mname = "MPI-ESM"
  model_names = (/"CanESM2","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2M","MPI-ESM"/)
  dir_names = (/"canesm2","csiro_mk36","gfdl_cm3","gfdl_esm2m","mpi"/)
  model_case_yr1 = (/1950,1850,1920,1950,1850/)
  model_case_yr2 = (/2100,2100,2100,2100,2099/)
  model_nsim = (/50,29,20,30,100/)

  imodel = ind(model_names.eq.mname)
  dirmname = dir_names(imodel)
  ncases = model_nsim(imodel)
  case_yr1 = model_case_yr1(imodel)
  case_yr2 = model_case_yr2(imodel)

  case_numbers = ispan(1,ncases,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  LEdir = "/glade/collections/cdg/data/CLIVAR_LE/"+dirmname+"_lens/Amon/"+var_name+"/"
;  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"
  datadir = "/glade/scratch/landrum/CMIP5/atm/"+dirmname+"/"

  nyr_stats = 10
  last_year_stats = case_yr2 - nyr_stats + 1
  if (case_yr1.le.1920) then
     yr1_stats = ispan(1920,last_year_stats,1)
  else
     yr1_stats = ispan(1950,last_year_stats,1)
  end if
  yr2_stats = yr1_stats+nyr_stats-1
  nstats = dimsizes(yr1_stats)
  sstat_dates = sprinti("%0.4i",yr1_stats) + "-"+sprinti("%0.4i",yr2_stats)

  lat_min = 50.				; looking at Arctic so retrieve only data north of this latitude

  fdata = LEdir+var_name+"_Amon_"+mname+"_historical_rcp85_r"+case_numbers+"i1p1_"+case_yr1+"01-"+case_yr2+"12.nc"
  fgrid = fdata(0)                                                                                     ; to get grid info and limit to nh alone

  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
  fstats = datadir+mname+"."+ncases+"sim.NH."+nyr_stats+"yr_stats."+month_names(imonth_stats)+"_monthly_"+var_name+"."+sstat_dates+".nc"

  if (.not.(isfilepresent(datadir))) then
     system("mkdir "+datadir)
  end if
 
  ntot = nyr_stats*ncases

;------------------------------------------------------------------------------
; read in data
;------------------------------------------------------------------------------
  f      = addfile (fgrid, "r")		; get grid info and create variable arrays
  tmplat = f->lat
  iNH = ind(tmplat.ge.lat_min)
  inh_lat = iNH(0)				; lat range -90 to +90 - retrieve only data north of lat_min
  delete(tmplat)
  delete(iNH)

  lon = f->lon
  lat = f->lat(inh_lat:)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  var = new((/ntot,nlat,nlon/),float)
  var!0 = "year"
  var!1 = "lat"
  var!2 = "lon"

  var_stat = new((/1,4,nlat,nlon/),float)
  var_stat!0 = "year"
  var_stat!1 = "moment"
  var_stat!2 = "lat"
  var_stat&lat = lat
  var_stat!3 = "lon"
  var_stat&lon = lon
  moment = ispan(1,4,1)
  var_stat&moment = moment

  var_min = new((/1,nlat,nlon/),float)
  var_min!0 = "year"
  var_min!1 = "lat"
  var_min&lat = lat
  var_min!2 = "lon"
  var_min&lon = lon

  var_max = var_min

  statLabel = (/"Mean","Var","Skew","Kurtosis"/)
  statLabel!0 = "moment"
  statLabel&moment = moment

  dimNames = (/"year","moment","lat","lon"/)
  dimSizes = (/1,4,nlat,nlon/)
  dimUnlim = (/True,False,False,False/)

  do is=0,nstats-1
     var_stat&year = yr1_stats(is)
     var_min&year = yr1_stats(is)
     var_max&year = yr1_stats(is)
     print("Calculating statistics for year: "+yr1_stats(is))
     do i = 0,ncases-1
        i1 = nyr_stats*i
        i2 = i1 + nyr_stats - 1
        f      = addfile (fdata(i), "r")
        iLE1 = (yr1_stats(is) - case_yr1)*12 + imonth_stats
        iLE2 = (yr2_stats(is) - case_yr1)*12 + imonth_stats
        var(i1:i2,:,:) = (/ f->$var_name$(iLE1:iLE2:12,inh_lat:,:) /)
     end do								; end ncases

     if (is.eq.0) then
        tmp = f->$var_name$(0,0,0)
        var_stat@units = tmp@units
        var_min@units = tmp@units
        var_max@units = tmp@units
        delete(tmp)
     end if

     var_stat(0,:,:,:) = (/ dim_stat4_n(var,0) /)
     var_min(0,:,:) = (/ dim_min_n(var,0) /)
     var_max(0,:,:) = (/ dim_max_n(var,0) /)
;-----------------------------------------------------------------------------------------
;    save monthly netCDF - in case takes too long to do all years
;-----------------------------------------------------------------------------------------
      system("/bin/rm -f "+fstats(is))
      setfileoption("nc","Format","NetCDF4")
      fout = addfile(fstats(is),"c")
      setfileoption(fout,"DefineMode",True)

      file_title = "LE ensemble: 4 moments of statistics and min, max for "+nyr_stats+" of daily "+month_names(imonth_stats)+" "+var_name
      fileAtt = True
      fileAtt@creation_date = systemfunc("date")
      fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
      fileAtt@title = file_title
      fileattdef(fout,fileAtt)

      filedimdef(fout,dimNames,dimSizes,dimUnlim)

      filevardef(fout,"lat",typeof(lat),getvardims(lat))
      filevardef(fout,"lon",typeof(lon),getvardims(lon))
      filevardef(fout,"var_stat",typeof(var_stat),getvardims(var_stat))
      filevardef(fout,"var_min",typeof(var_min),getvardims(var_min))
      filevardef(fout,"var_max",typeof(var_max),getvardims(var_max))
      filevardef(fout,"statLabel",typeof(statLabel),getvardims(statLabel))

      filevarattdef(fout,"lat",lat)
      filevarattdef(fout,"lon",lon)
      filevarattdef(fout,"var_stat",var_stat)
      filevarattdef(fout,"var_min",var_min)
      filevarattdef(fout,"var_max",var_max)

      fout->lat = (/ lat /)
      fout->lon = (/ lon /)
      fout->var_stat = (/ var_stat /)
      fout->var_min = (/ var_min /)
      fout->var_max = (/ var_max /)
      fout->statLabel = (/ statLabel /)
      delete(fout)
   end do 					; end nstats
end
