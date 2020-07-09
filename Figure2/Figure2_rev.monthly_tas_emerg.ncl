; Modified to plot example timeseries for tas from all
; CMIP5 MMLE and from one example for maps of ToE for Oct, 
; Allows user to plot change in tas and up to two sea ice 
; variables from the example as well.
; LLandrum
; Mar 2020
;
; Maps year that the mean variable exceeds threshold for OCT, FEB 
; and plots timeseries of mean, std for TREFHT for OCT, FEB at one land
; gridpoint (Fairbanks, AK) and one ocean (East Siberian Shelf) gridpoint.
;----------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
; ==============================================================
; User defined parameters that specify region of globe to consider
; ==============================================================
  mname = (/"MPI-ESM","CanESM2","CESM","GFDL-CM3","GFDL-ESM2M"/)
  dirmname = (/"mpi","canesm2","LE","gfdl_cm3","gfdl_esm2m"/)
  ncases = (/100,50,40,20,30/)
  stats_yr1 = (/1920,1950,1950,1950,1950/)				; some data for gridded data different than point data
  stats_yr2 = (/2090,2091,2091,2091,2091/)
  pnt_stats_yr1 = (/1920,1950,1920,1950,1950/)
  pnt_stats_yr2 = (/2090,2091,2091,2091,2091/)
  nens = dimsizes(mname)

  ex_mname = "CESM"
  iexample = ind(mname.eq.ex_mname)
  ex_dirmname = dirmname(iexample)

  var_name = "tas"
  base_clim_yr1 = 1950
  base_clim_yr2 = 1959
  base_clim_txt = base_clim_yr1+"-"+base_clim_yr2

  thresh_type = "std"					; if "max" then maximum value seen in base climate
  thresh = 2.						; x thresh_type
  xy2 = True                                   ; True to plot mn and std on same plot with different y axis, false on 2 different plots

  plot_points = True					; True to plot polymarkers for timeseries points
  overlay_obs = True					; True to overlay a point on timeseries plots for 1951-1980 mean obs value from Berkeley Earth Surface Temperature data

  point_short = (/"Chukchi","Fairbanks"/)
  point_long = (/"Chukchi Sea","Fairbanks, AK"/)
  point_lat = (/75.,64.8378/)
  point_lon = (/190.,147.7164/)
  point_lon(1) = 360. - point_lon(1)			; now Fairbanks in degE
  npnt = dimsizes(point_lat)
 
  nyr_stats = 10						; # yrs used in statistics

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"
  LEdatadir = "/glade/p/cesm/pcwg/landrum/LE/"
  iCESM = ind(mname.eq."CESM")
  mtitle = mname
  mtitle(iCESM) = "CESM1-CAM5"
  datadir(iCESM) = LEdatadir+"TREFHT/monthly/"
  ex_dirmname = dirmname(iexample)
  ex_mtitle = mtitle(iexample)

  obs_dir = "/glade/p/cesm/pcwg/landrum/obs/"
  obs_yr1 = 1979
  obs_yr2 = 1985
  nyr_stats_obs = 30
  obs_mn_yr = 1980					; which year to take nyr_stats_obs-yr-mean for plot

  plot_months = (/9,1/)
  plot_delta_tas = False					; True to plot change in tas, False otherwise
  plot_ice_var = False					; True to plot ice variables, false otherwise
  imonth_ice = (/1,1/)                                  ; must have same dimensions as plot_months
  plot_ice_yr = (/2015,2050/)
  var2_name = (/"sic","sit"/)                          ; choices are sit and sic
  var2_label = (/"SIC","SIT"/)
  var2_CESM = where(var2_name.eq."sic","aice","hi")
  var2_plot_type = (/"diff","mean"/)               ; mean or diff
  overlay_15per = (/True,False/)                  ; True to overlay 15% contours

  max_dT = 12.                          ; scaling for variable changes (2015 12; 2050 24; 2085 28)
  min_dT = -1.*max_dT

  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
  nmnth = dimsizes(plot_months)

  psout = "Figure2_rev.CMIP5_"+nens+"MMLE.monthly_tas_emergence"
;------------------------------------------------------------------------------
; read in grid point statistics timeseries
;------------------------------------------------------------------------------
  syr1 = min(pnt_stats_yr1)
  syr2 = max(pnt_stats_yr2)
  year = ispan(syr1,syr2,1)
  nyr = syr2 - syr1 + 1

  stats = new((/nens,npnt,nmnth,nyr,4/),float)
  stats!0 = "model"
  stats!1 = "point"
  stats!2 = "month"
  stats!3 = "year"
  stats&year = year
  stats!4 = "moment"

  pnt_yr_emerg = new((/nens,npnt,nmnth/),integer)
  pnt_yr_emerg!0 = "model"
  pnt_yr_emerg!1 = "point"
  pnt_yr_emerg!2 = "month"

  mean_emerge = new((/nens,npnt,nmnth,nyr/),float)
  mean_emerge!0 = "model"
  mean_emerge!1 = "point"
  mean_emerge2 = "month"
  mean_emerge!3 = "year"
  mean_emerge&year = year

  do im=0,nens-1
     if (im.eq.iCESM) then
       fstats1 = datadir(im)+"LE."+point_short+"."+nyr_stats+"yr_stats."+month_names(plot_months(0))+"_monthly_TREFHT."+pnt_stats_yr1(im)+"-"+pnt_stats_yr2(im)+".nc"
       fstats2 = datadir(im)+"LE."+point_short+"."+nyr_stats+"yr_stats."+month_names(plot_months(1))+"_monthly_TREFHT."+pnt_stats_yr1(im)+"-"+pnt_stats_yr2(im)+".nc"
       femerg = datadir(im)+"emergence."+thresh+"_"+thresh_type+".LE."+ncases(im)+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(plot_months)+"_monthly_TREFHT.nc"
     else
       femerg = datadir(im)+"emergence."+thresh+"_"+thresh_type+"."+mname(im)+"."+ncases(im)+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(plot_months)+"_monthly_"+var_name+".nc"
       fstats1 = datadir(im)+mname(im)+"."+point_short+"."+nyr_stats+"yr_stats."+month_names(plot_months(0))+"_monthly_"+var_name+"."+pnt_stats_yr1(im)+"-"+pnt_stats_yr2(im)+".nc"
       fstats2 = datadir(im)+mname(im)+"."+point_short+"."+nyr_stats+"yr_stats."+month_names(plot_months(1))+"_monthly_"+var_name+"."+pnt_stats_yr1(im)+"-"+pnt_stats_yr2(im)+".nc"
     end if

     f1 = addfile(femerg(0),"r")
     tmplat = f1->lat
     tmplon = f1->lon

     i1 = pnt_stats_yr1(im) - syr1
     i2 = pnt_stats_yr2(im) - syr1
     do ip=0,npnt-1
        f = addfile(fstats1(ip),"r")
        stats(im,ip,0,i1:i2,:) = (/ f->stat_LE /)
        f = addfile(fstats2(ip),"r")
        stats(im,ip,1,i1:i2,:) = (/ f->stat_LE /)
        do i=0,nmnth-1
           f1 = addfile(femerg(i),"r")
           ilat_pt = ind_nearest_coord(point_lat(ip),tmplat,0)
           jlon_pt = ind_nearest_coord(point_lon(ip),tmplon,0)

           ToE = f1->yr_emergence(ilat_pt,jlon_pt)
           pnt_yr_emerg(im,ip,i) = ToE
           if (.not.ismissing(ToE)) then
              iyr = ind(year.eq.ToE)
              mean_emerge(im,ip,i,iyr:i2) = stats(im,ip,i,iyr:i2,0)
           end if
        end do                                                  ; end i
     end do                                                     ; end ip
     delete(tmplat)
     delete(tmplon)
  end do                                                        ; end im

  stats(:,:,:,:,0) = stats(:,:,:,:,0) - 273.15                  ; convert to deg C
  mean_emerge = mean_emerge - 273.15
;------------------------------------------------------------------------------
; read in observational estimates
;------------------------------------------------------------------------------
  if (overlay_obs) then
     obs_mean = new((/npnt,nmnth/),float)
     obs_std = new((/npnt,nmnth/),float)
     i1 = obs_mn_yr - obs_yr1
     do i=0,npnt-1
;        fdata = obs_dir+"best.tas."+point_short(i)+"."+nyr_stats_obs+"yr_stats."+month_names(plot_months)+"."+obs_yr1+"-"+obs_yr2+".nc"
        fdata = obs_dir+"erai.t2m."+point_short(i)+"."+nyr_stats_obs+"yr_stats."+month_names(plot_months)+"."+obs_yr1+"-"+obs_yr2+".nc"
        do j=0,nmnth-1
           f = addfile(fdata(j),"r")
           obs_mean(i,j) = (/ f->stat_tas(i1,0) /)
           var = (/ f->stat_tas(i1,1) /)
           obs_std(i,j) = sqrt(var)
        end do
     end do
  end if
;------------------------------------------------------------------------------
; read in example data
;------------------------------------------------------------------------------
  if (ex_mname.eq."CESM") then
     femerg = datadir(iexample)+"emergence."+thresh+"_"+thresh_type+".LE."+ncases(iexample)+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(plot_months)+"_monthly_TREFHT.nc"
  else
     femerg = datadir(iexample)+"emergence."+thresh+"_"+thresh_type+"."+mname(iexample)+"."+ncases(iexample)+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(plot_months)+"_monthly_"+var_name+".nc"
  end if
  f = addfile(femerg(0),"r")
  lat = f->lat
  lon = f->lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  yr_emergence = new((/nmnth,nlat,nlon/),integer)
  yr_emergence!0 = "month"
  yr_emergence!1 = "lat"
  yr_emergence&lat = lat
  yr_emergence!2 = "lon"
  yr_emergence&lon = lon

  delta_mean = new((/nmnth,nlat,nlon/),float)
  delta_mean!0 = "month"
  delta_mean!1 = "lat"
  delta_mean&lat = lat
  delta_mean!2 = "lon"
  delta_mean&lon = lon
  delta_var = delta_mean

  do i=0,nmnth-1
     f = addfile(femerg(i),"r")
     yr_emergence(i,:,:) = (/ f->yr_emergence /)
     delta_mean(i,:,:) = (/ f->delta_mean /)
     delta_var(i,:,:) = (/ f->delta_var /)
  end do

  if (plot_delta_tas) then
     if (iexample.eq.iCESM) then
        fbase1 = datadir(iexample)+"LE."+ncases(iexample)+"sim.NH."+nyr_stats+"yr_stats."+month_names(plot_months)+"_monthly_TREFHT."+base_clim_yr1+"-"+base_clim_yr2+".nc"
        fstats1 = datadir(iexample)+"LE."+ncases(iexample)+"sim.NH."+nyr_stats+"yr_stats."+month_names(plot_months)+"_monthly_TREFHT."+stats_yr1(iexample)+"-"+stats_yr2(iexample)+".nc"
     else
        fbase1 = datadir(iexample)+mname(iexample)+"."+ncases(iexample)+"sim.NH."+nyr_stats+"yr_stats."+month_names(plot_months)+"_monthly_tas."+base_clim_yr1+"-"+base_clim_yr2+".nc"
        fstats1 = datadir(iexample)+mname(iexample)+"."+ncases(iexample)+"sim.NH."+nyr_stats+"yr_stats."+month_names(plot_months)+"_monthly_tas."+stats_yr1(iexample)+"-"+stats_yr2(iexample)+".nc"
     end if

     diff_tas = new((/nmnth,nlat,nlon/),float)
     diff_tas!0 = "month"
     diff_tas!1 = "lat"
     diff_tas&lat = lat
     diff_tas!2 = "lon"
     diff_tas&lon = lon

     do ii = 0,nmnth-1
        f = addfile(fbase1(ii),"r")
        base_tas = f->var_stat(0,0,:,:)

        f1 = addfile(fstats1(ii),"r")
        i1 = plot_ice_yr(ii) - stats_yr1(iexample)
        tas_mean = (/ f1->var_stat(i1,0,:,:) /)       ; mean atm variable

        diff_tas(ii,:,:) = (/ tas_mean - base_tas /)
     end do 					; end ii
     diff_tas@units = "deg C"
  end if
  if (plot_ice_var) then
     sstat_base_dates = base_clim_yr1+"-"+base_clim_yr2
     plot_ice_yr2 = plot_ice_yr + nyr_stats - 1
     sstat2_dates = plot_ice_yr+"-"+plot_ice_yr2
     LEdatadir = "/glade/p/cesm/pcwg/landrum/LE/"

     if (iexample.eq.iCESM) then
        fbase2 = LEdatadir+var2_CESM+"/LE."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_"+var2_CESM+"_nh."+sstat_base_dates+".nc"
        fstats2 = LEdatadir+var2_CESM+"/LE."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_"+var2_CESM+"_nh."+sstat2_dates+".nc"

        faice_base = LEdatadir+"aice/LE."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_aice_nh."+sstat_base_dates+".nc"
        fstats_aice = LEdatadir+"aice/LE."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_aice_nh."+sstat2_dates+".nc"
     else
        fbase2 = datadir+"ice/"+dirmname(iexample)+"/"+mname(iexample)+"."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_"+var2_name+"."+sstat_base_dates+".nc"
        fstats2 = datadir+"ice/"+dirmname(iexample)+"/"+mname(iexample)+"."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_"+var2_name+"."+sstat2_dates+".nc"

        faice_base = datadir+"ice/"+dirmname(iexample)+"/"+mname(iexample)+"."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_sic."+sstat_base_dates+".nc"
        fstats_aice = datadir+"ice/"+dirmname(iexample)+"/"+mname(iexample)+"."+ncases(iexample)+"sim."+nyr_stats+"yr_stats."+month_names(imonth_ice)+"_sic."+sstat2_dates+".nc"
     end if

     f = addfile(fbase2(0),"r")
     if (iexample.eq.iCESM) then
        tlat = f->tlat
        tlon = f->tlon
     else
        tlat = f->rlat
        tlon = f->rlon
     end if
     dim_ocn_lat = dimsizes(tlat)
     dimslat = dimsizes(dim_ocn_lat)
     dim_ocn_lon = dimsizes(tlon)
     if (dimslat.eq.2) then
        ocn_grid_type = "2D"
        nj = dim_ocn_lat(0)
        ni = dim_ocn_lat(1)
     else
        ocn_grid_type = "1D"
        nj = dim_ocn_lat
        ni = dim_ocn_lon
     end if

     base_mn_var2 = new((/nmnth,nj,ni/),float)
     var2_mean = new((/nmnth,nj,ni/),float)
     base_aice = base_mn_var2
     aice_mn = base_mn_var2
 
     do ii=0,nmnth-1
        f = addfile(fbase2(ii),"r")
        base_mn_var2(ii,:,:) = (/ f->var_stat(0,0,:,:) /)
        f2 = addfile(fstats2(ii),"r")
        var2_mean(ii,:,:) = (/ f2->var_stat(0,0,:,:) /)     ; mean ice variable

        if ((var2_name(ii).eq."sic").and.(max(base_mn_var2).lt.1.1)) then                          ; convert from fraction to percent
           base_mn_var2(ii,:,:) = 100.*base_mn_var2(ii,:,:)
           var2_mean(ii,:,:) = 100.*var2_mean(ii,:,:)
        end if

        if (var2_name(ii).eq."sic") then             ;mask out 0 areas for plotting purposes
           var2_mean(ii,:,:) = where(var2_mean(ii,:,:).lt.1.,var2_mean@_FillValue,var2_mean(ii,:,:))
           base_mn_var2(ii,:,:) = where(base_mn_var2(ii,:,:).lt.1.,base_mn_var2@_FillValue,base_mn_var2(ii,:,:))

           base_aice(ii,:,:) = base_mn_var2(ii,:,:)
           aice_mn(ii,:,:) = var2_mean(ii,:,:)
        else
           f = addfile(faice_base(ii),"r")
           base_aice(ii,:,:) = (/ f->var_stat(0,0,:,:) /)

           f = addfile(fstats_aice(ii),"r")
           aice_mn(ii,:,:) = (/ f->var_stat(0,0,:,:) /)

           if (max(base_aice(ii,:,:)).lt.1.1) then                            ; convert from fraction to percent
              base_aice(ii,:,:) = 100.*base_aice(ii,:,:)
           end if
           var2_mean(ii,:,:) = where(base_aice(ii,:,:).lt.1.,var2_mean@_FillValue,var2_mean(ii,:,:))           ; in case plotting mean
        end if
     end do
     diff_var2 = var2_mean
     diff_var2 = var2_mean - base_mn_var2
  end if
;-----------------------------------------------------------------------------------------
; plot
;-----------------------------------------------------------------------------------------
;  wks = gsn_open_wks("pdf",psout)             ; send graphics to PNG file
  wks = gsn_open_wks("eps",psout)             ; send graphics to PNG file
  gsn_define_colormap(wks,"posneg_1")
  colors = (/5,8,17,15,10/)

  plot_ts = new((/npnt*nmnth/),graphic)
  plot_ts2 = new((/npnt*nmnth/),graphic)			; for xy2 attribute
  plot_ts_emerge = new((/npnt*nmnth/),graphic)

  res1 = True                           ; res for stats plot
  res1@gsnDraw = False
  res1@gsnFrame = False

  res1@vpWidthF    = 0.6      ; Make plots wider than
  res1@vpHeightF   = 0.2      ; they are high.

  res1@tiXAxisString = "year"
  res1@trXMinF = base_clim_yr1
  res1@trXMaxF = syr2

  res1@trYMinF = -36.
  res1@trYMaxF = 8.

  res1@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font

  res2 = True
  res2@trYMinF = 0.
  res2@trYMaxF = 6.

  res1@xyLineColors = colors   ; mean
  res2@xyLineColors = colors     ; variance or std

  res1@tiYAxisString    = "Mean (~S~o~N~C)"
;  res2@tiYAxisString    = "Variance (~S~o~N~C~S~2~N~)"
  res2@tiYAxisString    = "STD (~S~o~N~C)"

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
     do im=0,nmnth-1
        res1@gsnRightString    =  month_names(plot_months(im))
        if (xy2) then
;           plot_ts(ip*nmnth+im)  = gsn_csm_xy2(wks,stats&year,stats(:,ip,im,:,0),sqrt(stats(:,ip,im,:,1)),res1,res2)
           tmpplot = gsn_csm_xy2(wks,stats&year,stats(:,ip,im,:,0),sqrt(stats(:,ip,im,:,1)),res1,res2)
           plot_ts(ip*nmnth+im) = tmpplot
           plot_ts2(ip*nmnth+im)  = tmpplot@xy2
        else
           plot_ts(ip*nmnth+im)  = gsn_csm_xy(wks,stats&year,stats(:,ip,im,:,0),res1)
        end if
        plot_ts_emerge(ip*nmnth+im)  = gsn_csm_xy(wks,mean_emerge&year,mean_emerge(:,ip,im,:),res_emerge)
        overlay(plot_ts(ip*nmnth+im),plot_ts_emerge(ip*nmnth+im))
      end do                            ; end im
  end do                                ; end ip
  if (overlay_obs) then
     mkplot = new(npnt*nmnth,graphic)
     mkres = True
     mkres@gsMarkerIndex = 16
     mkres@gsMarkerSizeF = 0.008

     lnplot = new(npnt*nmnth,graphic)
     lnstdplot = new(npnt*nmnth,graphic)

     lnres = True
     lnres@gsLineThicknessF = 3
     lnres@gsLineColor = "black"
     lnres@tfPolyDrawOrder = "PostDraw"
 
     do i=0,npnt-1
        do j=0,nmnth-1
           lnplot(i*nmnth+j) = gsn_add_polyline(wks,plot_ts(i*nmnth+j),(/obs_mn_yr,obs_mn_yr+nyr_stats_obs/),(/obs_mean(i,j),obs_mean(i,j)/),lnres)
           if (xy2) then
              lnres@gsLineDashPattern = 1
              lnstdplot(i*nmnth+j) = gsn_add_polyline(wks,plot_ts2(i*nmnth+j),(/obs_mn_yr,obs_mn_yr+nyr_stats_obs/),(/obs_std(i,j),obs_std(i,j)/),lnres)
              lnres@gsLineDashPattern = 0
           end if
        end do
     end do
  end if

  pres = True
  pres@gsnFrame = False
  pres@gsnMaximize     = True                 ; fill up the page
  pres@gsnPanelMainString = "CMIP5 MMLE tas "
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.02
;  pres@gsnPanelRowSpec  = True
  pres@gsnPanelTop = 0.95
  if ((plot_delta_tas).or.(plot_ice_var)) then
     pres@gsnPanelBottom = 0.7
  else
     pres@gsnPanelBottom = 0.6
  end if

  lgres                    = True
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgDashIndex        = 0
  lgres@lgLineColors       = colors
  labels = mtitle
  lgres@lgLineThicknessF = 3.
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.08           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter

  nitems = dimsizes(labels)                 ; number of legend items

  gsn_panel(wks,plot_ts,(/npnt,nmnth/),pres)
  gsn_legend_ndc(wks,nitems,labels,0.75,0.6,lgres)

  obs_label = "ERAI 1980-2009"
  delete(lgres@lgLineColors)
  lgres@vpHeightF          = 0.04           ; height of legend (NDC)
  lgres@lgLineColor = lnres@gsLineColor
  gsn_legend_ndc(wks,1,obs_label,0.59,0.6,lgres)

;-----------------------------------now plot example contour maps-------------------------------------

  gsn_define_colormap(wks,"GMT_panoply")
  gsn_reverse_colormap(wks)

  plot = new(nmnth,graphic)
  plot_pnt2 = new(2,graphic)

  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@gsnFrame = False
  res@gsnDraw = False
  res@trGridType = "TriangularMesh"
  res@gsnPolar  = "NH"
  res@mpMinLatF =   45.
  res@gsnAddCyclic = True
  res@gsnPolarLabelSpacing = 30.
  res@cnLineLabelsOn       = False              ; turn the line labels off
  res@cnLinesOn       = False              ; turn the line labels off

  res@lbLabelBarOn        = False       ; Turn off labelbar

  res@gsnSpreadColors      = True                 ; use full colormap

  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  minXY   = 2006
  maxXY   = 2080
  ncntrs = 16
  mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)
  print("Contouring from "+mnmxint(0)+" to "+ mnmxint(1)+" by "+mnmxint(2))

  res@cnMinLevelValF       = mnmxint(0)
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnRightStringFontHeightF   = 0.03
  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnCenterStringFontHeightF   = 0.035
  res@gsnCenterString = ""
  res@gsnRightString = ""
  res@gsnLeftString = ""
  do i=0,nmnth-1
     res@gsnCenterString = month_names(plot_months(i))
     plot(i) = gsn_csm_contour_map_polar(wks,yr_emergence(i,:,:),res)     ; create plot
     res@gsnLeftString = ""
  end do

  pres@gsnPanelLabelBar  = True
  pres@gsnPanelMainString = "ToE: "+ex_mtitle+" tas"
  pres@lbTitleOn        =  True               ; turn on title
  pres@lbTitleString    = "yr"                ; title string
  pres@lbTitleFontHeightF= .013                ; make title smaller
  pres@lbTitlePosition  = "Right"              ; title position
  pres@lbTitleDirection = "Across"             ; title direction
  pres@lbLabelFontHeightF = 0.01
;
  point_colors = (/"black","red"/)
  if (plot_points) then
     polyres                   = True
     polyres@gsMarkerIndex     = 16          ; polymarker style
     polyres@gsMarkerSizeF     = 10.          ; polymarker size
     polyres@gsMarkerColor = point_colors(0)
     dum1 = gsn_add_polymarker(wks,plot(0),point_lon(0),point_lat(0),polyres)
     dum2 = gsn_add_polymarker(wks,plot(1),point_lon(0),point_lat(0),polyres)
     polyres@gsMarkerColor = point_colors(1)
     dum3 = gsn_add_polymarker(wks,plot(0),point_lon(1),point_lat(1),polyres)
     dum4 = gsn_add_polymarker(wks,plot(1),point_lon(1),point_lat(1),polyres)
  end if

  pres@gsnPanelTop = pres@gsnPanelBottom - 0.1
  if (plot_delta_tas) then
     pres@gsnPanelBottom = 0.5
     y1txt = 0.57
     y2txt = 0.553
  else
     pres@gsnPanelBottom = 0.08
     y1txt = 0.37 
     y2txt = 0.353 
  end if
  gsn_panel(wks,plot,(/1,2/),pres)

  tres               = True                    ; text res
  tres@txFontHeightF = 0.010
  tres@txFontColor = point_colors(0)
  gsn_text_ndc(wks,point_long(0),0.1,y1txt, tres)
  tres@txFontColor = point_colors(1)
  gsn_text_ndc(wks,point_long(1),0.1,y2txt, tres)

  if (plot_delta_tas) then
     delete(minXY)
     delete(maxXY)
     delete(mnmxint)

     n_colors_two = 19
     n_colors_one = 16
 
     gsn_define_colormap(wks,"GMT_panoply")
     gsn_merge_colormaps(wks,"GMT_panoply", "temp_diff_18lev") ; merge 2 colormaps

     plot_tas = new((nmnth)*nens,graphic)

     res@cnLinesOn       = True              ; turn the line labels off

     res@cnLevels = ispan(-24,20,4)
     ncntrs = 18

     color_start              = n_colors_one + 2
     color_end                = n_colors_one+n_colors_two
     res@gsnSpreadColors      = True                 ; use full colormap
     res@gsnSpreadColorStart  = color_start
     res@gsnSpreadColorEnd    = color_end-1

     minXY   = min_dT
     maxXY   = max_dT
     mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       = mnmxint(0)
     res@cnMaxLevelValF       = mnmxint(1)
     res@cnLevelSpacingF      = mnmxint(2)

     res@gsnRightString = "tas"

     do i=0,nmnth-1
        plot_tas(i) = gsn_csm_contour_map_polar(wks,diff_tas(i,:,:),res)     ; create plot
     end do

     pres@lbTitleString    = "~S~o~N~C"                ; title string
     pres@gsnPanelTop = pres@gsnPanelBottom - 0.1
     if (plot_ice_var) then
        pres@gsnPanelBottom = 0.3
     else
        pres@gsnPanelBottom = 0.08
     end if
     gsn_panel(wks,plot_tas,(/1,2/),pres)
  end if

  if (plot_ice_var) then
     plot_ice = new((nmnth)*nens,graphic)
     plot_15per = new((nmnth)*nens,graphic)

     res@lbLabelBarOn        = True       ; Turn off labelbar

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
     res15per@trGridType = "TriangularMesh"

     res15per@cnMonoLineColor = True
     res15per@cnLineColor = "white"

     do i=0,nmnth-1
        if (var2_name(i).eq."sic") then
           max_dSIC = 100.                    ; scaling for sea ice concentration changes
        else
           max_dSIC = 2.                       ; scaling for sea ice thickness changes
        end if
        min_dSIC = -1.*max_dSIC

        if (var2_plot_type(i).eq."mean") then
           minXY   = 0.
           maxXY   = max_dSIC
           res@gsnSpreadColorStart  = 2
           res@gsnSpreadColorEnd    = n_colors_one
           data = (/ var2_mean(i,:,:) /)
        else
           res@gsnSpreadColorStart  = color_start
           res@gsnSpreadColorEnd    = color_end-1

           minXY   = min_dSIC
           maxXY   = max_dSIC

           data = (/ diff_var2(i,:,:) /)
        end if

        mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)
        res@cnMinLevelValF       = mnmxint(0)
        res@cnMaxLevelValF       = mnmxint(1)
        res@cnLevelSpacingF      = mnmxint(2)

        if (ocn_grid_type.eq."2D") then
           data@lat2d = tlat
           data@lon2d = tlon
           aice_mn@lat2d = tlat
           aice_mn@lon2d = tlon
        else
           data!0 = "lat"
           data!1 = "lon"
           data&lat = tlat
           data&lon = tlon

           aice_mn!0 = "lat"
           aice_mn!1 = "lon"
           aice_mn&lat = tlat
           aice_mn&lon = tlon
        end if

        plot_ice(i) = gsn_csm_contour_map_polar(wks,data,res)     ; create plot

        if (overlay_15per(i)) then
           plot_15per(i) = gsn_csm_contour(wks,aice_mn(i,:,:),res15per)
           overlay(plot_tas(i),plot_15per(i))
           plot_15per(i+1) = gsn_csm_contour(wks,aice_mn(i,:,:),res15per)
           overlay(plot_ice(i), plot_15per(i+1))
        end if
     end do
     pres@gsnPanelLabelBar  = False
     pres@gsnPanelTop = pres@gsnPanelBottom - 0.08
     pres@gsnPanelBottom = 0.08
     gsn_panel(wks,plot_ice,(/1,2/),pres)
  end if

  frame(wks)

end
