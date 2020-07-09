; CMIP5 MMLE (only CanESM2 and CESM1-CAM5 have daily precip including snow)
; Feb 2020
;
; Maps year that the mean freeze-thaw variables exceed threshold
; along with timeseries of mean, std one land gridpoint (Fairbanks, AK) 
; and one ocean (East Siberian Shelf) gridpoint.
; User defines threshold. Currently thresholds include above # of STDS, 
; or exceeds maxima, etc. Thresholds are by comparison with user-defined time period.
;
; assumes netCDF files of gridded statistics for specified variable exist - script does
; save a netCDF of "emergence" calculation (or uses already saved one)..
;
; LLandrum
; Feb 2020
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

  mname = (/"CanESM2","CESM1-CAM5"/)
  dirmname = (/"canesm2","cesm"/)
  stats_yr1 = (/1950,1920/)
;  stats_yr1 = (/1950,1950/)
  stats_yr2 = (/2091,2091/)
  nsim = (/50,40/)

  nmodel = dimsizes(mname)
  iCESM = ind(mname.eq."CESM1-CAM5")


  var_names = (/"first_rain","last_rain","rain_season_length"/)
  long_var_names = (/"First Rain Day","Last Rain Day","Rain Season Length"/)
  thresh_type = "std"                                   ; if "max" then maximum value seen in base climate
  thresh = (/-2.,2.,2./)                                 ; x thresh_type
  xy2 = True                                   ; True to plot mn and std on same plot with different y axis, false on 2 different plots

  mean_min = (/20.,240.,80./)
  mean_max = (/160.,360.,360./)
  std_min = (/10.,16.,20./)
  std_max = (/60.,30.,80./)

  nvar = dimsizes(var_names)
  min_snow_percent = 60.     ; percent of total precip that must be snow for "freeze up onset" definition

  nyr_stats = 10                                                ; # yrs used in statistics

  base_clim_yr1 = 1950
  base_clim_yr2 = base_clim_yr1 + nyr_stats - 1

  plot_points = True                                    ; True to plot polymarkers for timeseries points

  point_short = (/"Chukchi","Fairbanks"/)
  point_long = (/"Chukchi Sea","Fairbanks, AK"/)
  point_lat = (/75.,64.8378/)
  point_lon = (/190.,147.7164/)
  point_lon(1) = 360. - point_lon(1)                    ; now Fairbanks in degE

  npnt = dimsizes(point_lat)
  nyr_stats = 10						; # yrs used in statistics

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/precip/"
  datadir(iCESM) = "/glade/p/cesm/pcwg/landrum/LE/precip/"

  fstats = datadir+mname+"."+nsim+"sim.NH."+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+stats_yr1+"-"+stats_yr2+".nc"
  fstats(iCESM) = datadir(iCESM)+"LE."+nsim(iCESM)+"sim.NH."+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+stats_yr1(iCESM)+"-"+stats_yr2(iCESM)+".nc"

  psout = "Figure4_rev.freeze_thaw."+thresh_type+"_emergence_from_"+base_clim_yr1+"-"+base_clim_yr2+".min"+min_snow_percent+"percent"

  pi=4.*atan(1.0)
  rad=(pi/180.)
;------------------------------------------------------------------------------
; read in point data for both models
;------------------------------------------------------------------------------
  syr1 = min(stats_yr1)
  syr2 = max(stats_yr2)
  year = ispan(syr1,syr2,1)
  nyr = syr2 - syr1 + 1

  stats = new((/nmodel,npnt,nvar,nyr,4/),float)
  stats!0 = "model"
  stats!1 = "point"
  stats!2 = "variable"
  stats!3 = "year"
  stats&year = ispan(syr1,syr2,1)
  stats!4 = "moment"

  pnt_yr_emerg = new((/nmodel,npnt,nvar/),integer)
  pnt_yr_emerg!0 = "model"
  pnt_yr_emerg!1 = "point"
  pnt_yr_emerg!2 = "variable"

  mean_emerge = new((/nmodel,npnt,nvar,nyr/),float)
  mean_emerge!0 = "model"
  mean_emerge!1 = "point"
  mean_emerge2 = "variable"
  mean_emerge!3 = "year"
  mean_emerge&year = year

  do im=0,nmodel-1
     f = addfile(fstats(im),"r")
     lat = f->lat
     lon = f->lon
     i1 = stats_yr1(im) - syr1
     i2 = stats_yr2(im) - syr1

     if (im.eq.iCESM) then
       femerg = datadir(iCESM)+"emergence."+thresh+"_"+thresh_type+"."+base_clim_yr1+"-"+base_clim_yr2+"_baseclim.LE."+nsim(iCESM)+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_names+".nc"
     else
       femerg = datadir(im)+"emergence."+thresh+"_"+thresh_type+"."+mname(im)+"."+nsim(im)+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_names+".nc"
     end if

     do i=0,npnt-1
        ilat_pt = ind_nearest_coord(point_lat(i),lat,0)
        jlon_pt = ind_nearest_coord(point_lon(i),lon,0)
        do iv=0,nvar-1
           tmpvar_name = var_names(iv)+"_stats"     
           stats(im,i,iv,i1:i2,:) = (/ f->$tmpvar_name$(:,:,ilat_pt,jlon_pt) /)

           f1 = addfile(femerg(iv),"r")

           ToE = f1->yr_emergence(ilat_pt,jlon_pt)
           pnt_yr_emerg(im,i,iv) = ToE
           if (.not.ismissing(ToE)) then
              iyr = ind(year.eq.ToE)
              mean_emerge(im,i,iv,iyr:i2) = stats(im,i,iv,iyr:i2,0)
           end if
        end do
     end do 						; end i
     delete(lat)
     delete(lon)
  end do 						; end im

;------------------------------------------------------------------------------
; read in emergence data for both models
;------------------------------------------------------------------------------

  print("Mean year of emergence by model for "+var_names(0)+", "+var_names(1)+" and "+var_names(2)+" area averaged 70-90 N")
  avg_ToE = new(nvar,float)

  if (iCESM.eq.0) then
     femerg = datadir(iCESM)+"emergence."+thresh+"_"+thresh_type+"."+base_clim_yr1+"-"+base_clim_yr2+"_baseclim.LE."+nsim(iCESM)+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_names+".nc"
  else
     femerg = datadir(0)+"emergence."+thresh+"_"+thresh_type+"."+mname(0)+"."+nsim(0)+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_names+".nc"
  end if
  do i=0,nvar-1
     f = addfile(femerg(i),"r")
     if (i.eq.0) then
        lat=f->lat
        lon=f->lon
        nlat = dimsizes(lat)
        nlon = dimsizes(lon)
        yr_emergence1 = new((/nvar,nlat,nlon/),integer)
        yr_emergence1!0 = "variable"
        yr_emergence1!1 = "lat"
        yr_emergence1&lat = lat
        yr_emergence1!2 = "lon"
        yr_emergence1&lon = lon

        coswgt=cos(rad*lat)
        coswgt!0 = "lat"
        coswgt&lat= lat
     end if
     yr_emergence1(i,:,:) = f->yr_emergence

     avg_ToE(i) = wgt_areaave(yr_emergence1(variable|i,{lat|70.:90.},lon|:),coswgt({70.:90.}),1.0,0)
  end do

  print(mname(0)+":        "+RoundPastDecimal(avg_ToE(0),3)+"             "+RoundPastDecimal(avg_ToE(1),3)+"             "+RoundPastDecimal(avg_ToE(2),3))

  delete(lat)
  delete(lon)
  delete(coswgt)

  if (iCESM.eq.1) then
     femerg = datadir(iCESM)+"emergence."+thresh+"_"+thresh_type+"."+base_clim_yr1+"-"+base_clim_yr2+"_baseclim.LE."+nsim(iCESM)+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_names+".nc"
  else
     femerg = datadir(1)+"emergence."+thresh+"_"+thresh_type+"."+mname(1)+"."+nsim(1)+"sim.NH_"+nyr_stats+"yr_stats.freeze_thaw.min"+min_snow_percent+"percent."+var_names+".nc"
  end if

  do i=0,nvar-1
     f = addfile(femerg(i),"r")
     if (i.eq.0) then     
        lat=f->lat     
        lon=f->lon     
        nlat = dimsizes(lat)
        nlon = dimsizes(lon)
        yr_emergence2 = new((/nvar,nlat,nlon/),integer)
        yr_emergence2!0 = "variable"
        yr_emergence2!1 = "lat"
        yr_emergence2&lat = lat
        yr_emergence2!2 = "lon"
        yr_emergence2&lon = lon

        coswgt=cos(rad*lat)
        coswgt!0 = "lat"
        coswgt&lat= lat
     end if
     yr_emergence2(i,:,:) = f->yr_emergence
     avg_ToE(i) = wgt_areaave(yr_emergence2(variable|i,{lat|70.:90.},lon|:),coswgt({70.:90.}),1.0,0)
  end do
  print(mname(1)+":        "+RoundPastDecimal(avg_ToE(0),3)+"             "+RoundPastDecimal(avg_ToE(1),3)+"             "+RoundPastDecimal(avg_ToE(2),3))
  delete(lat)
  delete(lon)
  delete(coswgt)
;-----------------------------------------------------------------------------------------
; plot - first timeseries
;-----------------------------------------------------------------------------------------
;  wks = gsn_open_wks("pdf",psout)             ; send graphics to PNG file
  wks = gsn_open_wks("eps",psout)             ; send graphics to PNG file
  gsn_define_colormap(wks,"posneg_1")
  colors = (/8,17/)				; colors for CESM1 and CanESM in other plots

  plot_ts = new((/npnt*nvar/),graphic)
  plot_ts_emerge = new((/npnt*nvar/),graphic)

  res1 = True                           ; res for stats plot
  res1@gsnDraw = False
  res1@gsnFrame = False

  res1@vpWidthF    = 0.6      ; Make plots wider than
  res1@vpHeightF   = 0.2      ; they are high.

  res1@tiXAxisString = "year"
  res1@trXMinF = base_clim_yr1
  res1@trXMaxF = syr2

  res1@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font

  res2 = True

  res1@xyLineColors = colors   ; mean
  res2@xyLineColors = colors     ; variance or std

  res1@tiYAxisString    = "Mean (day)"
;  res2@tiYAxisString    = "Variance (day~S~2~N~)"
  res2@tiYAxisString    = "STD (day)"

  res1@xyDashPattern = 0
  res1@xyLineThicknessF = 20

  res2@xyDashPattern = 1
  res1@xyLineThicknessF = 4

  res_emerge = True
  res_emerge = res1
  res1@xyLineOpacityF = 0.3
  res_emerge@xyLineOpacityF = 1.

  do ip=0,npnt-1
     res1@gsnLeftString    =  point_short(ip)
     do ivar=0,nvar-1
        res1@trYMinF = mean_min(ivar)
        res1@trYMaxF = mean_max(ivar)
        res2@trYMinF = std_min(ivar)
        res2@trYMaxF = std_max(ivar)

        res1@gsnRightString    =  long_var_names(ivar)
        if (xy2) then
           plot_ts(ip*nvar+ivar)  = gsn_csm_xy2(wks,stats&year,stats(:,ip,ivar,:,0),sqrt(stats(:,ip,ivar,:,1)),res1,res2)
        else
           plot_ts(ip*nvar+ivar)  = gsn_csm_xy(wks,stats&year,stats(:,ip,ivar,:,0),res1)
        end if
        plot_ts_emerge(ip*nvar+ivar)  = gsn_csm_xy(wks,mean_emerge&year,mean_emerge(:,ip,ivar,:),res_emerge)
        overlay(plot_ts(ip*nvar+ivar),plot_ts_emerge(ip*nvar+ivar))
      end do                            ; end ivar
  end do                                ; end ip

  pres = True
  pres@gsnFrame = False
  pres@txFontHeightF = 0.012
  pres@gsnMaximize     = True                 ; fill up the page
  pres@gsnPanelMainString = "CMIP5 MMLE First and Last rain days"
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.02

  pres@gsnPanelTop = 0.95
  pres@gsnPanelBottom = 0.6

  lgres                    = True
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgDashIndex        = 0
  lgres@lgLineColors       = colors
  labels = mname
  lgres@lgLineThicknessF = 3.
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.05           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter

  nitems = dimsizes(labels)                 ; number of legend items

  gsn_panel(wks,(/plot_ts(0),plot_ts(1),plot_ts(3),plot_ts(4)/),(/npnt,2/),pres)
  gsn_legend_ndc(wks,nitems,labels,0.75,0.6,lgres)

;-----------------------------------now plot contour maps-------------------------------------

  gsn_define_colormap(wks,"GMT_panoply")
  gsn_reverse_colormap(wks)

  nplot = nmodel*nvar
  plot = new(nplot,graphic)
  dum = new(npnt*nplot,graphic)

  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@gsnFrame = False
  res@gsnDraw = False
  res@trGridType = "TriangularMesh"
  res@gsnPolar  = "NH"
  res@mpMinLatF =   45.
  res@gsnAddCyclic = True
  res@gsnPolarLabelSpacing = 30.

  res@lbLabelBarOn        = False       ; Turn off labelbar

  res@gsnSpreadColors      = True                 ; use full colormap

  ncntrs = 16
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  minXY   = 2006
  maxXY   = 2080
  mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)
  print("Contouring from "+mnmxint(0)+" to "+ mnmxint(1)+" by "+mnmxint(2))

  res@cnLinesOn       = False                   ; turn off contour lines
  res@cnLineLabelsOn       = False                   ; turn off contour lines

  res@cnMinLevelValF       = mnmxint(0)
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnCenterString = ""
  res@gsnRightStringFontHeightF   = 0.02
  res@gsnLeftStringFontHeightF   = 0.03
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@gsnCenterStringFontHeightF   = 0.02

  do i=0,nvar-1
     left_strings = (/mname(0),mname(1)/)
     res@gsnCenterString = left_strings(0)
     plot(i) = gsn_csm_contour_map_polar(wks,yr_emergence1(i,:,:),res)     ; create plot
     res@gsnCenterString = left_strings(1)
     plot(nvar+i) = gsn_csm_contour_map_polar(wks,yr_emergence2(i,:,:),res)     ; create plot
  end do

  point_colors = (/"black","red"/)

  if (plot_points) then
     polyres                   = True
     polyres@gsMarkerIndex     = 16          ; polymarker style
     polyres@gsMarkerSizeF     = 10.          ; polymarker size
     polyres@gsMarkerColor = point_colors(0)
     do i=0,npnt-1
        polyres@gsMarkerColor = point_colors(i)
        do j=0,nplot-1
           dum(i*nplot+j) = gsn_add_polymarker(wks,plot(j),point_lon(i),point_lat(i),polyres)
        end do
     end do
  end if

  pres@gsnPanelLabelBar  = True
  pres@gsnPanelMainString = "ToE: Rain Season Length"
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.02

  pres@lbTitleOn        =  True               ; turn on title
  pres@lbTitleString    = "yr"                ; title string
  pres@lbTitleFontHeightF= .013                ; make title smaller
  pres@lbTitlePosition  = "Right"              ; title position
  pres@lbTitleDirection = "Across"             ; title direction
  pres@lbLabelFontHeightF = 0.01

  pres@gsnPanelTop = pres@gsnPanelBottom - 0.1
  pres@gsnPanelBottom = 0.08
  y1txt = 0.27
  y2txt = 0.253
  gsn_panel(wks,(/plot(2),plot(5)/),(/1,2/),pres)

  tres               = True                    ; text res
  tres@txFontHeightF = 0.010
  tres@txFontColor = point_colors(0)
  gsn_text_ndc(wks,point_long(0),0.1,y1txt, tres)
  tres@txFontColor = point_colors(1)
  gsn_text_ndc(wks,point_long(1),0.1,y2txt, tres)

  frame(wks)
end
