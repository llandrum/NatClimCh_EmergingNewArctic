; Maps changes in monthly tas from CMIP5 MMLE for
; user specified month and decades. One figure per model
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
  imonth = 1
  imonth_ice = 1

  base_yr1 = 1950
  base_yr2 = 1959

  var_name = "tas"
  var_label = "Surface Air Temperature"

  nyr_stats = 10                                                ; # yrs used in statistics
  plot_yr1 = (/2015,2050,2085/)
  plot_yr2 = plot_yr1 + nyr_stats - 1
  sstat2_dates = plot_yr1+"-"+plot_yr2

  mname = (/"MPI-ESM","CanESM2","CESM","GFDL-CM3","GFDL-ESM2M"/)
  dirmname = (/"mpi","canesm2","LE","gfdl_cm3","gfdl_esm2m"/)
  ncases = (/100,50,40,20,30/)
  stats_yr1 = (/1920,1950,1950,1920,1950/)                           ; some data for gridded data different than point data
  stats_yr2 = (/2090,2091,2091,2091,2091/)
  nens = dimsizes(mname)

  ndecades = dimsizes(plot_yr1)
  nplot = dimsizes(ndecades*nens)

  overlay_15per = True
  overlay_85per = True
  sstat_base_dates = base_yr1+"-"+base_yr2
  sstat_dates = stats_yr1+"-"+stats_yr2

  atm_datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"
  ice_datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"
  LEdatadir = "/glade/p/cesm/pcwg/landrum/LE/"
  iCESM = ind(mname.eq."CESM")
  mtitle = mname
  mtitle(iCESM) = "CESM1-CAM5"
  atm_datadir(iCESM) = LEdatadir+"TREFHT/monthly/"

  max_dT = 24.				; scaling for variable changes (2015 12; 2050 24; 2085 28)
  min_dT = -1.*max_dT

  sstat2_dates = plot_yr1+"-"+plot_yr2

  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

  fbase = atm_datadir+mname+"."+ncases+"sim.NH."+nyr_stats+"yr_stats."+month_names(imonth)+"_monthly_"+var_name+"."+sstat_base_dates+".nc"
  fstats = atm_datadir+mname+"."+ncases+"sim.NH."+nyr_stats+"yr_stats."+month_names(imonth)+"_monthly_"+var_name+"."+sstat_dates+".nc"
  fbase(iCESM) = atm_datadir(iCESM)+"/LE."+ncases(iCESM)+"sim.NH."+nyr_stats+"yr_stats."+month_names(imonth)+"_monthly_TREFHT."+sstat_base_dates+".nc"
  fstats(iCESM) = atm_datadir(iCESM)+"/LE."+ncases(iCESM)+"sim.NH."+nyr_stats+"yr_stats."+month_names(imonth)+"_monthly_TREFHT."+sstat_dates(iCESM)+".nc"

  fbase2 = ice_datadir+mname+"."+ncases+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_sic."+sstat_base_dates+".nc"
  fbase2(iCESM) = LEdatadir+"aice/LE."+ncases(iCESM)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_aice_nh."+sstat_base_dates+".nc"

  psout = "SupplementaryFigure6.dTAS."+ndecades+"_decades.CMIP5_"+nens+"MMLE."+month_names(imonth)

  pi=4.*atan(1.0)
  rad=(pi/180.)
;------------------------------------------------------------------------------
;  need to read in and map at same time - different grids
;-----------------------------------------------------------------------------------------
  wks = gsn_open_wks("eps",psout)             ; send graphics to PNG file
;  wks = gsn_open_wks("pdf",psout)             ; send graphics to PNG file

  gsn_define_colormap(wks,"temp_diff_18lev")

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

  res@lbLabelBarOn        = False
  res@gsnSpreadColors      = True                 ; use full colormap

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnLeftStringFontHeightF   = 0.045
  res@gsnRightStringFontHeightF   = 0.04
  res@gsnCenterStringFontHeightF   = 0.03

  res@gsnRightString = " "
  res@gsnLeftString = ""
  res@gsnCenterString = " "

;  color_start              = n_colors_one + 1
;  color_end                = n_colors_one+n_colors_two-1
;  color_start              = n_colors_one + 2
;  color_end                = n_colors_one+n_colors_two
;  res@gsnSpreadColorStart  = color_start
;  res@gsnSpreadColorEnd    = color_end-1

  ncntrs = 18
  minXY   = min_dT
  maxXY   = max_dT
  mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)
;  print("Contouring from "+mnmxint(0)+" to "+ mnmxint(1)+" by "+mnmxint(2))

  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = mnmxint(0)
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)

  print("Mean "+month_names(imonth)+" temperature change from ("+base_yr1+"-"+base_yr2+") by model,  area averaged 70-90 N")
  avg_dtas = new(ndecades,float)

  if (overlay_15per) then
     res@cnLinesOn       = False                   ; turn off contour lines

     res15per = True
     res15per@gsnDraw          = False     ; don't draw
     res15per@gsnFrame         = False     ; don't advance frame
     res15per@cnLevelSelectionMode =  "ExplicitLevels"
     res15per@cnInfoLabelOn    =  False
     res15per@gsnLeftString    =  ""
     res15per@cnLevels         =  (/15./)
     res15per@tiXAxisString    = ""
     res15per@tiYAxisString    = ""
     res15per@cnLineThicknessF = 3.0       ; thicker contours
     res15per@cnLineLabelsOn   = False     ; no line labels
     res15per@cnMonoLineDashPattern = True
     res15per@trGridType = "TriangularMesh"

     res15per@gsnRightString = " "

     res15per@cnMonoLineColor = True
     res15per@cnLineColor = "black"

     plot_15per = new(2*nplot,graphic)
  end if
  if (overlay_85per) then
     res@cnLinesOn       = False                   ; turn off contour lines

     res85per = True
     res85per@gsnDraw          = False     ; don't draw
     res85per@gsnFrame         = False     ; don't advance frame
     res85per@cnLevelSelectionMode =  "ExplicitLevels"
     res85per@cnInfoLabelOn    =  False
     res85per@gsnLeftString    =  ""
     res85per@cnLevels         =  (/85./)
     res85per@tiXAxisString    = ""
     res85per@tiYAxisString    = ""
     res85per@cnLineThicknessF = 3.0       ; thicker contours
     res85per@cnLineLabelsOn   = False     ; no line labels
     res85per@cnMonoLineDashPattern = True
     res85per@trGridType = "TriangularMesh"

     res85per@gsnRightString = " "

     res85per@cnMonoLineColor = True
     res85per@cnLineColor = "black"

     plot_85per = new(2*nplot,graphic)
  end if

  do im=0,nens-1
     f = addfile(fbase(im),"r")
     lon = f->lon
     lat = f->lat
     nlat = dimsizes(lat)
     nlon = dimsizes(lon)

     coswgt=cos(rad*lat)
     coswgt!0 = "lat"
     coswgt&lat= lat

     var_base = f->var_stat(0,0,:,:)

     var_diff = new((/nlat,nlon/),float)
     var_diff!0 = "lat"
     var_diff&lat = lat
     var_diff!1 = "lon"
     var_diff&lon = lon

     f1 = addfile(fstats(im),"r")

     if ((overlay_15per).or.(overlay_85per)) then
        f2 = addfile(fbase2(im),"r")
        if (im.ne.iCESM) then
           fstats2 = ice_datadir(im)+mname(im)+"."+ncases(im)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_sic."+sstat2_dates+".nc"
           ocn_lon = f2->rlon
           ocn_lat = f2->rlat
        else
           fstats2 = LEdatadir+"aice/LE."+ncases(im)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_aice_nh."+sstat2_dates+".nc"
           ocn_lon = f2->tlon
           ocn_lat = f2->tlat
        end if
        dim_lat_ocn = dimsizes(ocn_lat)
        dim_lon_ocn = dimsizes(ocn_lon)

        dimslat = dimsizes(dim_lat_ocn)
        aice_base = f2->var_stat(0,0,:,:)
        if (dimslat.eq.2) then
           ocn_grid_type = "2D"
           nj = dim_lat_ocn(0)
           ni = dim_lat_ocn(1)
        else
           ocn_grid_type = "1D"
           nj = dim_lat_ocn
           ni = dim_lon_ocn
        end if
        delete(dim_lat_ocn)
        delete(dim_lon_ocn)

        aice_mean = new((/nj,ni/),float)
        if (ocn_grid_type.eq."2D") then
           aice_mean@lat2d = ocn_lat
           aice_mean@lon2d = ocn_lon

           aice_base@lat2d = ocn_lat
           aice_base@lon2d = ocn_lon
        else
           aice_base!0 = "lat"
           aice_base!1 = "lon"
           aice_base&lat = ocn_lat
           aice_base&lon = ocn_lon

           aice_mean!0 = "lat"
           aice_mean!1 = "lon"
           aice_mean&lat = ocn_lat
           aice_mean&lon = ocn_lon
        end if

        if (max(aice_base).lt.1.1) then                            ; convert from fraction to percent
           aice_base = 100.*aice_base
        end if

        if (im.ne.iCESM) then
           fstats2 = ice_datadir(im)+mname(im)+"."+ncases(im)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_sic."+sstat2_dates+".nc"
        else
           fstats2 = LEdatadir+"aice/LE."+ncases(im)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_aice_nh."+sstat2_dates+".nc"
        end if
     end if

     do i=0,ndecades-1
        i1 = plot_yr1(i) - stats_yr1(im)
        var_mean =  f1->var_stat(i1,0,:,:) 
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

        avg_dtas(i) = wgt_areaave(var_diff({lat|70.:90.},lon|:),coswgt({70.:90.}),1.0,0)
        if ((overlay_15per).or.(overlay_85per)) then
           f2 = addfile(fstats2(i),"r")
           aice_mean = (/ f2->var_stat(0,0,:,:) /)
           if (max(aice_mean).lt.1.1) then                            ; convert from fraction to percent
              aice_mean = 100.*aice_mean
           end if
           if (overlay_15per) then
              res15per@cnLineColor = "black"
              plot_15per(im*2*ndecades+i) = gsn_csm_contour(wks,aice_base,res15per)
              overlay(plot(im*ndecades+i), plot_15per(im*2*ndecades+i))
              res15per@cnLineColor = "white"
              plot_15per(im*2*ndecades+i+ndecades) = gsn_csm_contour(wks,aice_mean,res15per)
              overlay(plot(im*ndecades+i), plot_15per(im*2*ndecades+i+ndecades))
           end if
           if (overlay_85per) then
              res85per@cnLineColor = "grey"
              plot_85per(im*2*ndecades+i+ndecades) = gsn_csm_contour(wks,aice_mean,res85per)
              overlay(plot(im*ndecades+i), plot_85per(im*2*ndecades+i+ndecades))
           end if
        end if
     end do
     print(mname(im)+":        "+RoundPastDecimal(avg_dtas(0),3)+"             "+RoundPastDecimal(avg_dtas(1),3)+"             "+RoundPastDecimal(avg_dtas(2),3))
     delete(var_base)
     delete(var_mean)
     delete(var_diff)
     delete(lat)
     delete(lon)
     delete(coswgt)

     if (overlay_15per) then
        delete(fstats2)
        delete(aice_base)
        delete(aice_mean)
        delete(ocn_lat)
        delete(ocn_lon)
     end if

  end do

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar  = False
  pres@lbTitleString    = "~S~o~N~C"
  pres@lbTitlePosition  = "Right"
  pres@lbTitleOffsetF  = 0.001
  pres@lbTitleFontHeightF= .015
  pres@lbTitleDirection = "Across"
;  pres@lbBoxEndCapStyle = "TriangleBothEnds"

  pres@gsnPanelTop = 0.9
     
  pres@gsnPanelLabelBar  = True
  pres@gsnMaximize     = True                 ; fill up the page
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.016
  pres@gsnPanelMainString = month_names(imonth)+" "+var_label+" changes"
;
  gsn_panel(wks,plot,(/nens,ndecades/),pres)
  frame(wks)

end
