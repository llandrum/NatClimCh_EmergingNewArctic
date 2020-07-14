; Maps changes in user-picked rain variable for CMIP5 MMLE 
; (two ensemles) for user-specified decades. One figure per model
; (ensemble mean).
; LLandrum
; Apr 2020
;----------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
begin
; ==============================================================
; User defined parameters that specify region of globe to consider
; ==============================================================
  base_yr1 = 1950
  base_yr2 = 1959

;  var_name = "rain_season_length"                       ; choice are first_rain, thaw_onset and freeze_onset
;  var_label = "Rain Season Duration"
;  var_name = "last_rain"
;  var_label = "Last Rain Day"
  var_name = "first_rain"
  var_label = "First Rain Day"

  min_snow_percent = 60.     ; percent of total precip that must be snow for "freeze up onset" definition

  nyr_stats = 10                                                ; # yrs used in statistics
  plot_yr1 = (/2015,2050,2085/)
  plot_yr2 = plot_yr1 + nyr_stats - 1

  mname = (/"CanESM2","CESM"/)
  dirmname = (/"canesm2","LE"/)
  ncases = (/50,40/)
  stats_yr1 = (/1950,1920/)                           ; some data for gridded data different than point data
  stats_yr2 = (/2091,2091/)
  nens = dimsizes(mname)

  ndecades = dimsizes(plot_yr1)
  nplot = dimsizes(ndecades*nens)

  sstat_base_dates = base_yr1+"-"+base_yr2
  sstat_dates = stats_yr1+"-"+stats_yr2

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/precip/"
  LEdatadir = "/glade/p/cesm/pcwg/landrum/LE/precip/"
  iCESM = ind(mname.eq."CESM")
  mtitle = mname
  mtitle(iCESM) = "CESM1-CAM5"
  datadir(iCESM) = LEdatadir

  max_drain = 120.				; scaling for variable changes (2015 12; 2050 24; 2085 28)
  min_drain = -1.*max_drain

  sstat2_dates = plot_yr1+"-"+plot_yr2

  fstats = datadir+mname+"."+ncases+"sim.NH."+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+sstat_dates+".nc"
  fstats(iCESM) = datadir(iCESM)+"LE."+ncases(iCESM)+"sim.NH."+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+sstat_dates(iCESM)+".nc"

  psout = "SupplementaryFigure7.d"+var_name+"."+ndecades+"_decades.CMIP5_"+nens+"MMLE"

  pi=4.*atan(1.0)
  rad=(pi/180.)
;------------------------------------------------------------------------------
;  need to read in and map at same time - different grids
;-----------------------------------------------------------------------------------------
  wks = gsn_open_wks("eps",psout)             ; send graphics to PNG file
;  wks = gsn_open_wks("pdf",psout)             ; send graphics to PNG file

  gsn_define_colormap(wks,"NCV_blu_red")

  nplot = ndecades*nens
  plot = new(nplot,graphic)

  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@gsnFrame = False
  res@gsnDraw = False
  res@trGridType = "TriangularMesh"
  res@gsnPolar  = "NH"
  res@mpMinLatF =   50.
  res@gsnAddCyclic = True
  res@gsnPolarLabelSpacing = 30.

  res@cnLinesOn       = False                   ; turn off contour lines
  res@cnLineLabelsOn       = False                   ; turn off contour lines
 
  res@lbLabelBarOn        = False
  res@gsnSpreadColors      = True                 ; use full colormap

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnLeftStringFontHeightF   = 0.038
  res@gsnRightStringFontHeightF   = 0.03
  res@gsnCenterStringFontHeightF   = 0.03

  res@gsnRightString = " "
  res@gsnLeftString = ""
  res@gsnCenterString = " "

  ncntrs = 28
  minXY   = min_drain
  maxXY   = max_drain
  mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)
;  print("Contouring from "+mnmxint(0)+" to "+ mnmxint(1)+" by "+mnmxint(2))

  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = mnmxint(0)
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)

  print("Mean "+var_name+" change from ("+base_yr1+"-"+base_yr2+") by model,  area averaged 70-90 N")
  avg_drain = new(ndecades,float)

  stat_var_name = var_name+"_stats"
  do im=0,nens-1
     i1 = base_yr1 - stats_yr1(im)
     f = addfile(fstats(im),"r")
     lon = f->lon
     lat = f->lat
     nlat = dimsizes(lat)
     nlon = dimsizes(lon)

     coswgt=cos(rad*lat)
     coswgt!0 = "lat"
     coswgt&lat= lat
     
     var_base = f->$stat_var_name$(i1,0,:,:)

     var_diff = new((/nlat,nlon/),float)
     var_diff!0 = "lat"
     var_diff&lat = lat
     var_diff!1 = "lon"
     var_diff&lon = lon

     do i=0,ndecades-1
        i1 = plot_yr1(i) - stats_yr1(im)
        var_mean =  f->$stat_var_name$(i1,0,:,:) 
        var_diff = (/ var_mean - var_base /)				; no need to check and convert from deg Kelvin to deg C as taking a difference
        if (i.eq.0) then
           res@gsnLeftString = mtitle(im)
        else
           res@gsnLeftString = " "
        end if
        if (im.eq.0) then
           res@gsnRightString = plot_yr1(i)+"-"+plot_yr2(i)
        else
           res@gsnRightString = " "
        end if
        plot(im*ndecades+i) = gsn_csm_contour_map_polar(wks,var_diff,res)     ; create plot

        avg_drain(i) = wgt_areaave(var_diff({lat|70.:90.},lon|:),coswgt({70.:90.}),1.0,0)
     end do
     print(mname(im)+":        "+RoundPastDecimal(avg_drain(0),3)+"             "+RoundPastDecimal(avg_drain(1),3)+"             "+RoundPastDecimal(avg_drain(2),3))
     delete(var_base)
     delete(var_mean)
     delete(var_diff)
     delete(lat)
     delete(lon)
     delete(coswgt)
  end do

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar  = False
  pres@lbTitleString    = "days"
  pres@lbTitlePosition  = "Right"
  pres@lbTitleOffsetF  = 0.001
  pres@lbTitleFontHeightF= .015
  pres@lbTitleDirection = "Across"

  pres@gsnPanelTop = 0.9
     
  pres@gsnPanelLabelBar  = True
  pres@gsnMaximize     = True                 ; fill up the page
  pres@gsnPanelMainString = var_label+" changes ~C~"
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.02
;
  gsn_panel(wks,plot,(/nens,ndecades/),pres)
  frame(wks)

end
