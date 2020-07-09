; Calculates and saves mean and variance along with time of emergence
; for user specified variable and emergence criteria.
; User defines variable, month and threshold. Currently thresholds include above # of STDS, 
; or exceeds maxima, etc. Thresholds are by comparison with user-defined time period.
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
  mname = "CanESM2"
  dirmname = "canesm2"
  nsim = 50

;  var_name = "first_rain"			; choice include first_rain, last_rain, rain_season_length
  var_name = "last_rain"
  min_snow_percent = 65.     ; percent of total precip that must be snow for "freeze up onset" definition

  base_clim_yr1 = 1950
  base_clim_yr2 = 1959
  plot_change = True					; True to plot changes in mean, variance (or min, max if thresh_type.eq. min or max)

  netCDF = True						; True to save netCDF, False if it exists

  thresh_type = "std"					; if "max" then maximum value seen in base climate
  if ((var_name.eq."first_rain").or.(var_name.eq."thaw_onset")) then
     thresh = -2.
  else
     thresh = 2.						; x thresh_type
  end if

  nyr_stats = 10						; # yrs used in statistics
  stats_yr1 = 1950
  stats_yr2 = 2091					; end of statistics year - although no stats for yrs stats_yr2-nyrs through stats_yr2

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/precip/"

  fstats = datadir+mname+"."+nsim+"sim.NH."+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+stats_yr1+"-"+stats_yr2+".nc"
  femerg = datadir+"emergence."+thresh+"_"+thresh_type+"."+mname+"."+nsim+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_name+".nc"
;------------------------------------------------------------------------------
; read in data - both comparison decade and future...
;------------------------------------------------------------------------------
  if (netCDF) then
     var_stat_name = var_name+"_stats"
     var_min_name = var_name+"_min_annual"
     var_max_name = var_name+"_max_annual"

     f = addfile(fstats,"r")
     lon = f->lon
     lat = f->lat
     nlat = dimsizes(lat)
     nlon = dimsizes(lon)
 
     f = addfile(fstats,"r")
     var_mean = f->$var_stat_name$(:,0,:,:)		; mean	
     var_var = f->$var_stat_name$(:,1,:,:)            ; variance
     var_min = f->$var_min_name$	               ; min
     var_max = f->$var_max_name$	               ; max

     if (var_name.eq."rain_season_length") then				; need to set to FillValue points where first and last rain days are the same - i.e. no rain season
        var_mean = where(var_mean.lt.1.,var_mean@_FillValue,var_mean)
     end if
     ib = base_clim_yr1 - stats_yr1
     base_mean = var_mean(ib,:,:)
     base_var = var_var(ib,:,:)
     base_min = var_min(ib,:,:)
     base_max = var_max(ib,:,:)

     if (thresh_type.eq."std") then
        std = sqrt(base_var)
        base_thresh = base_mean + thresh*std
     else
        if (thresh_type.eq."min") then
           base_thresh = thresh*base_min
        else                                            ; thresh_type "max"
           base_thresh = thresh*base_max
        end if
     end if

     dims = dimsizes(var_mean)
     nyr = dims(0)
     year = ispan(stats_yr1,stats_yr1+nyr-1,1)            
;------------------------------------------------------------------------------
; calculate year of emergence
;------------------------------------------------------------------------------
     yr_emergence = new((/nlat,nlon/),integer)
     yr_emergence!0 = "lat"
     yr_emergence&lat = lat
     yr_emergence!1 = "lon"
     yr_emergence&lon = lon

     delta_mean= new((/nlat,nlon/),float)
     delta_mean!0 = "lat"
     delta_mean&lat = lat
     delta_mean!1 = "lon"
     delta_mean&lon = lon

     if (thresh.gt.0) then
        delta_mean@long_name = "Year at which mean "+var_name+" exceeds "+thresh+"*("+base_clim_yr1+"-"+base_clim_yr2+") "+thresh_type
        yr_emergence@long_name = "Year at which mean "+var_name+" exceeds "+thresh+"*("+base_clim_yr1+"-"+base_clim_yr2+") "+thresh_type
     else
        delta_mean@long_name = "Year at which mean "+var_name+" falls below "+thresh+"*("+base_clim_yr1+"-"+base_clim_yr2+") "+thresh_type
        yr_emergence@long_name = "Year at which mean "+var_name+" falls below "+thresh+"*("+base_clim_yr1+"-"+base_clim_yr2+") "+thresh_type
     end if

     delta_std = delta_mean
     delta_max = delta_mean
     delta_min = delta_mean

     delta_mean@long_name = "Change in mean at year of emergence"
     delta_min@long_name = "Change in minimum at year of emergence" 
     delta_max@long_name = "Change in maximum at year of emergence" 

     do i=0,nlat-1
        do j=0,nlon-1
           if ((thresh_type.eq."max").or.((thresh_type.eq."std").and.(thresh.gt.0.))) then
              ind_emerg = ind(var_mean(:,i,j).gt.base_thresh(i,j))
           else							; for min
              ind_emerg = ind(var_mean(:,i,j).lt.base_thresh(i,j))
           end if
           if (.not.(ismissing(ind_emerg(0)))) then
              yr_emergence(i,j) = (/ year(ind_emerg(0)) /)
              delta_mean(i,j) = (/ var_mean(ind_emerg(0),i,j) - base_mean(i,j) /)
              delta_std(i,j) = (/ sqrt(var_var(ind_emerg(0),i,j)) - sqrt(base_var(i,j)) /)
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
     filevardef(fout,"delta_std",typeof(delta_std),getvardims(delta_std))
     filevardef(fout,"delta_min",typeof(delta_min),getvardims(delta_min))
     filevardef(fout,"delta_max",typeof(delta_max),getvardims(delta_max))

     filevarattdef(fout,"lat",lat)
     filevarattdef(fout,"lon",lon)
     filevarattdef(fout,"yr_emergence",yr_emergence)
     filevarattdef(fout,"delta_mean",delta_mean)
     filevarattdef(fout,"delta_std",delta_std)
     filevarattdef(fout,"delta_min",delta_min)
     filevarattdef(fout,"delta_max",delta_max)

     fout->lat = (/ lat /)
     fout->lon = (/ lon /)
     fout->yr_emergence = (/ yr_emergence /)
     fout->delta_mean = (/ delta_mean /)
     fout->delta_std = (/ delta_std /)
     fout->delta_min = (/ delta_min /)
     fout->delta_max = (/ delta_max /)

   else
      f = addfile(femerg,"r")
      yr_emergence = f->yr_emergence
      lat = f->lat
      lon = f->lon
      delta_mean = f->delta_mean
      delta_std = f->delta_std
   end if
end
