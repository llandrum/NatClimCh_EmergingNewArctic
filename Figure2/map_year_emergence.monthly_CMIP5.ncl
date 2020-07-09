; Calculates and saves year that the mean variable exceeds threshold for given month of daily statistics.
; User defines variable, month and threshold. Currently thresholds include above # of STDS, 
; or exceeds maxima, etc. Thresholds are by comparison with user-defined time period.
;
; This script also saves mean and variance at the time of emergence for given variable.
; User can plot changes in mean and variance as well.
;
; assumes netCDF files of gridded statistics for specified variable exist - script does
; save a netCDF of "emergence" calculation (or uses already saved one)..
;
; LLandrum
; March 2020
;----------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
; ==============================================================
; User defined parameters that specify region of globe to consider
; ==============================================================
  imonth = 11
  var_name = "tas"

  mname = "MPI-ESM"
  model_names = (/"CanESM2","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2M","MPI-ESM"/)
  dir_names = (/"canesm2","csiro_mk36","gfdl_cm3","gfdl_esm2m","mpi"/)
  model_stats_yr1 = (/1950,1920,1920,1950,1920/)
  model_stats_yr2 = (/2091,2091,2091,2091,2090/)
  model_nsim = (/50,29,20,30,100/)

  imodel = ind(model_names.eq.mname)
  dirmname = dir_names(imodel)
  nsim = model_nsim(imodel)
  stats_yr1 = model_stats_yr1(imodel)
  stats_yr2 = model_stats_yr2(imodel)

  nyr_stats = 10 

  base_clim_yr1 = 1950
  base_clim_yr2 = base_clim_yr1 + nyr_stats - 1
  plot_change = True					; True to plot changes in mean, variance (or min, max if thresh_type.eq. min or max)

  netCDF = True						; True to save netCDF, False if it exists

  thresh_type = "std"					; if "max" then maximum value seen in base climate
  thresh = 2.						; x thresh_type

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"

  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

  fstats = datadir+mname+"."+nsim+"sim.NH."+nyr_stats+"yr_stats."+month_names(imonth)+"_monthly_"+var_name+"."+stats_yr1+"-"+stats_yr2+".nc"
  femerg = datadir+"emergence."+thresh+"_"+thresh_type+"."+mname+"."+nsim+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(imonth)+"_monthly_"+var_name+".nc"
 
;------------------------------------------------------------------------------
; read in data - both comparison decade and future...
;------------------------------------------------------------------------------
  if (netCDF) then
     f = addfile(fstats,"r")
     lon = f->lon
     lat = f->lat
     nlat = dimsizes(lat)
     nlon = dimsizes(lon)

     i1 = base_clim_yr1 - stats_yr1
 
     base_mean = f->var_stat(i1,0,:,:)
     base_var = f->var_stat(i1,1,:,:)
     base_min = f->var_min(i1,:,:)
     base_max = f->var_max(i1,:,:)
     if (thresh_type.eq."std") then
        std = sqrt(base_var)
        base_thresh = base_mean + thresh*std
     else
        if (thresh_type.eq."min") then
           base_thresh = thresh*base_min
        else						; thresh_type "max"
           base_thresh = thresh*base_max
        end if
     end if

     var_mean = f->var_stat(i1:,0,:,:)		; mean	
     var_var = f->var_stat(i1:,1,:,:)            ; variance
     var_min = f->var_min(i1:,:,:)	               ; min
     var_max = f->var_max(i1:,:,:)	               ; max
     dims = dimsizes(var_mean)
     nyr = dims(0)
     year = ispan(base_clim_yr1,base_clim_yr1+nyr-1,1)            
;------------------------------------------------------------------------------
; calculate year of emergence
;------------------------------------------------------------------------------
     yr_emergence = new((/nlat,nlon/),integer)
     yr_emergence!0 = "lat"
     yr_emergence&lat = lat
     yr_emergence!1 = "lon"
     yr_emergence&lon = lon
     yr_emergence@long_name = "Year at which mean "+month_names(imonth)+" monthly "+var_name+" exceeds "+thresh+"*("+base_clim_yr1+"-"+base_clim_yr2+") "+thresh_type

     delta_mean= new((/nlat,nlon/),float)
     delta_mean!0 = "lat"
     delta_mean&lat = lat
     delta_mean!1 = "lon"
     delta_mean&lon = lon
     delta_mean@long_name = "Year at which mean "+month_names(imonth)+" monthly "+var_name+" exceeds "+thresh+"*("+base_clim_yr1+"-"+base_clim_yr2+") "+thresh_type

     delta_var = delta_mean
     delta_max = delta_mean
     delta_min = delta_mean

     delta_mean@long_name = "Change in mean at year of emergence"
     delta_min@long_name = "Change in minimum at year of emergence" 
     delta_max@long_name = "Change in maximum at year of emergence" 

     do i=0,nlat-1
        do j=0,nlon-1
           if ((thresh_type.eq."max").or.(thresh_type.eq."std")) then
              ind_emerg = ind(var_mean(:,i,j).gt.base_thresh(i,j))
           else							; for min
              ind_emerg = ind(var_mean(:,i,j).lt.base_thresh(i,j))
           end if
           if (.not.(ismissing(ind_emerg(0)))) then
              yr_emergence(i,j) = (/ year(ind_emerg(0)) /)
              delta_mean(i,j) = (/ var_mean(ind_emerg(0),i,j) - base_mean(i,j) /)
              delta_var(i,j) = (/ var_var(ind_emerg(0),i,j) - base_var(i,j) /)
              delta_max(i,j) = (/ var_max(ind_emerg(0),i,j) - base_max(i,j) /)
              delta_min(i,j) = (/ var_min(ind_emerg(0),i,j) - base_min(i,j) /)
           end if
           delete(ind_emerg)
        end do
     end do
;------------------------------------------------------------------------------
; save netCDF file
;------------------------------------------------------------------------------
     dimNames = (/"lat","lon"/)
     dimSizes = (/nlat,nlon/)
     dimUnlim = (/False,False/)

     system("/bin/rm -f "+femerg)
     setfileoption("nc","Format","NetCDF4")
     fout = addfile(femerg,"c")
     setfileoption(fout,"DefineMode",True)

     file_title = "Computed using "+fstats
     fileAtt = True
     fileAtt@creation_date = systemfunc("date")
     fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
     fileAtt@title = file_title
     fileattdef(fout,fileAtt)

     filedimdef(fout,dimNames,dimSizes,dimUnlim)

     filevardef(fout,"lat",typeof(lat),getvardims(lat))
     filevardef(fout,"lon",typeof(lon),getvardims(lon))
     filevardef(fout,"yr_emergence",typeof(yr_emergence),getvardims(yr_emergence))
     filevardef(fout,"delta_mean",typeof(delta_mean),getvardims(delta_mean))
     filevardef(fout,"delta_var",typeof(delta_var),getvardims(delta_var))
     filevardef(fout,"delta_min",typeof(delta_min),getvardims(delta_min))
     filevardef(fout,"delta_max",typeof(delta_max),getvardims(delta_max))

     filevarattdef(fout,"lat",lat)
     filevarattdef(fout,"lon",lon)
     filevarattdef(fout,"yr_emergence",yr_emergence)
     filevarattdef(fout,"delta_mean",delta_mean)
     filevarattdef(fout,"delta_var",delta_var)
     filevarattdef(fout,"delta_min",delta_min)
     filevarattdef(fout,"delta_max",delta_max)

     fout->lat = (/ lat /)
     fout->lon = (/ lon /)
     fout->yr_emergence = (/ yr_emergence /)
     fout->delta_mean = (/ delta_mean /)
     fout->delta_var = (/ delta_var /)
     fout->delta_min = (/ delta_min /)
     fout->delta_max = (/ delta_max /)

   else
      f = addfile(femerg,"r")
      yr_emergence = f->yr_emergence
      lat = f->lat
      lon = f->lon
      delta_mean = f->delta_mean
      delta_var = f->delta_var
   end if
end
