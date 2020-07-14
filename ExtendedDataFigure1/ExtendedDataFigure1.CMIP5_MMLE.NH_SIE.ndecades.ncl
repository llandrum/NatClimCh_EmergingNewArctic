; Plots CMIP5 MMLE annual monthly climatology
; for SIE (or vhi or vai - user specified) for LE for user-specified
; time periods. Assumes 10-yr running mean statistics have already
; been created. Use netCDF_annual_NH_clim.ncl if not.
; 
; LLandrum
; Apr 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin
;
  include_CSIRO = False
  if (include_CSIRO) then
     mname = (/"MPI-ESM","CanESM2","CESM","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2M"/)
     dirmname = (/"mpi","canesm2","LE","csiro_mk36","gfdl_cm3","gfdl_esm2m"/)
     ncases = (/99,50,40,29,20,30/)
     stats_yr1 = (/1920,1950,1920,1920,1920,1950/)                              ; some dates for gridded data different than point data
     stats_yr2 = (/2090,2091,2091,2091,2091,2091/)
  else
     mname = (/"MPI-ESM","CanESM2","CESM","GFDL-CM3","GFDL-ESM2M"/)
     dirmname = (/"mpi","canesm2","LE","gfdl_cm3","gfdl_esm2m"/)
     ncases = (/99,50,40,20,30/)
     stats_yr1 = (/1920,1950,1920,1920,1950/)                           ; some data for gridded data different than point data
     stats_yr2 = (/2090,2091,2091,2091,2091/)
  end if
  nens = dimsizes(mname)

  iCESM = ind(mname.eq."CESM")
  mtitle = mname
  mtitle(iCESM) = "CESM1-CAM5"

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"
  datadir(iCESM) = "/glade/p/cesm/pcwg/landrum/LE/ice/pre_process/monthly/"

  nyr_stats = 10

  grid_month = True					; True to overlay grid of months 
  plot_yr1 = (/1980,2015,2050,2085/)		; each yr is beginning of nyr_stats running statistics
  ndec = dimsizes(plot_yr1)
  plot_yr2 = plot_yr1 + nyr_stats - 1

  plot_var = "ext"		; choices are ext, vhi, vai
  plot_obs = True		; True to overlay NSIDC sea ice index obs. only for extent!

  plot_std = False		; True to plot shaded std

  fclim = datadir+"monthly_clim."+nyr_stats+"_yrstats."+mname+"."+ncases+"_sim."+stats_yr1+"-"+stats_yr2+".nc"
  fclim(iCESM) = datadir(iCESM)+"monthly_clim."+nyr_stats+"_yrstats.LE."+ncases(iCESM)+"_sim."+stats_yr1(iCESM)+"-"+stats_yr2(iCESM)+".nc"

  yr1_ssmi = 1978					; both of these start/end years are incomplete so will use 1979-2018
  yr2_ssmi = 2020
  fssmi_monthly = "/glade/u/home/landrum/SSMI/nsidc_monthly_ice_index."+yr1_ssmi+"-"+yr2_ssmi+".nc"

  psout = "ExtendedDataFigure1.CMIP5_"+nens+"MMLE.NH_SIE."+ndec+"dec."+plot_var+"_clim"

;---------------------------------------------------------
; for making nice x axis labels
;---------------------------------------------------------
  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
  month_init = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  ndays_per_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  nday_month_start = ndays_per_month
  nday_month_start(0) = 1
  do i=1,11
     nday_month_start(i) = nday_month_start(i-1)+ndays_per_month(i-1)
  end do
;---------------------------------------------------------
; get data...
;---------------------------------------------------------
  cmip_ice_data = new((/nens,ndec,12/),float)

  var_name = plot_var+"_stats"

  do im=0,nens-1
     f = addfile(fclim(im),"r")
     var_name = plot_var+"_stats"
     mn_data = f->$var_name$(:,:,0)

     do i=0,ndec-1
        i1 = plot_yr1(i) - stats_yr1(im)
        cmip_ice_data(im,i,:) = mn_data(i1,:)
     end do
     delete(mn_data)
  end do
  
  cmip_ice_data!2 = "month"
  month = ispan(1,12,1)
  cmip_ice_data&month = month

  if (plot_obs) then
     obs_yr1 = yr1_ssmi + 1                      ; first and last years are not full years of data - so do not use
     obs_yr2 = yr2_ssmi - 1                      ; first yr = 1979, last year = 2017

     f = addfile(fssmi_monthly,"r")
     nh_ext = f->nh_ext

     i1 = (obs_yr1 - yr1_ssmi)*12
     i2 = (obs_yr1 + nyr_stats - yr1_ssmi)*12 - 1
     nh_ext1 = nh_ext(i1:i2)

     i1 = (obs_yr2 - nyr_stats - yr1_ssmi)*12
     i2 = (obs_yr2 - yr1_ssmi)*12 - 1
     nh_ext2 = nh_ext(i1:i2)

     ssmi_clim = new((/2,12/),float)
     ssmi_clim!0 = "decade"
     ssmi_clim!1 = "month"
     ssmi_clim&month = month
     do i=0,11
        ssmi_clim(0,i) = (/ dim_avg(nh_ext1(i::12)) /)
        ssmi_clim(1,i) = (/ dim_avg(nh_ext2(i::12)) /)
     end do
  end if
;===================================================================================================;
; plot
;===================================================================================================;
;  wks  = gsn_open_wks("pdf",psout)             ; open a ps file
  wks  = gsn_open_wks("eps",psout)             ; open a ps file
  gsn_define_colormap(wks,"posneg_1")
  if (include_CSIRO) then
     colors = (/5,8,10,13,15,17/)
  else
     colors = (/5,8,17,15,10/)
  end if

  plot = new(ndec,graphic)
  dum_icefree = new(ndec,graphic)
  ssmi_plot = new(ndec,graphic)
  dum = new(12*ndec,graphic)

  res = True                            ; res for stats plot
  res@gsnDraw = False
  res@gsnFrame = False

  res@trXMinF = 1
  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues    = month
  res@tmXBLabels    = month_init

  res@trYMinF = 0
  if (plot_var.eq."ext") then
     res@trYMaxF = 20.
     res@tiYAxisString= "Sea Ice Extent (10~S~12~N~ m~S~2~N~)"
  else
     res@trYMaxF = 4.5
     res@tiYAxisString= "Sea Ice Volume (10~S~13~N~ m~S~3~N~)"
  end if

  res@xyLineColors = colors
  res@xyMonoDashPattern = True                         ; all lines solid
  res@xyDashPattern = 0.
  res@xyLineThicknessF = 3. 

  res@tiMainString = ""

  res@tiXAxisString = "Month"
  res@trXMaxF = 12

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnCenterStringFontHeightF   = 0.035
  res@gsnCenterString = ""

  lres = True
  lres@gsLineDashPattern = 0
  lres@gsLineColor = "grey"
  lres@gsLineThicknessF = 1.

  if (plot_obs)
     res_obs = True
     res_obs@gsnDraw = False
     res_obs@gsnFrame = False
     res_obs@xyLineThicknessF = 3.

     res_obs@xyDashPattern = 1

;     res_obs@pmLegendDisplayMode    = "Always"   ; Display a legend.
;     res_obs@lgPerimOn              = False      ; No legend perimeter.
;     res_obs@lgJustification        = "TopRight"
;     res_obs@pmLegendOrthogonalPosF = -0.45     ; Move into plot
;     res_obs@pmLegendWidthF         = 0.4       ; Decrease width
;     res_obs@pmLegendHeightF        = 0.15       ; Decrease height
  end if

  do i=0,ndec-1
     res@gsnCenterString = plot_yr1(i)+"-"+plot_yr2(i)
     plot(i) = gsn_csm_xy(wks,cmip_ice_data&month,cmip_ice_data(:,i,:),res)

     if (plot_var.eq."ext") then
        dum_icefree(i) = gsn_add_polyline(wks,plot(i),(/1,12/),(/1.,1./),lres)
     end if

     if (grid_month) then
        do j=0,10
           dum(11*i+j) = gsn_add_polyline(wks,plot(i),(/j+1,j+1/),(/0,20/),lres)
        end do
     end if

     if (plot_obs) then
        if (plot_yr1(i).lt.1990) then
           res_obs@xyDashPattern = 1
           res_obs@xyExplicitLegendLabels = "Sea Ice Index"+obs_yr1+"-"+(obs_yr1 + 9)
           ssmi_plot(i) = gsn_csm_xy(wks,ssmi_clim&month,ssmi_clim(0,:),res_obs)
        else
           res_obs@xyDashPattern = 2
           res_obs@xyExplicitLegendLabels = "Sea Ice Index"+(obs_yr2 - 9)+"-"+obs_yr2
           ssmi_plot(i) = gsn_csm_xy(wks,ssmi_clim&month,ssmi_clim(1,:),res_obs)
        end if
        overlay(plot(i),ssmi_plot(i))
     end if
  end do

  pres = True
  pres@gsnFrame = False
  pres@gsnMaximize     = True                 ; fill up the page
  pres@gsnPanelMainString = "CMIP5 MMLE NH SIE"
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.02
;  pres@gsnPanelDebug = True

  pres@gsnPanelTop = 0.9
  pres@gsnPanelBottom = 0.53
  gsn_panel(wks,plot(:1),(/1,2/),pres)
; 
  pres@gsnPanelMainString = ""
;
  pres@gsnPanelTop = 0.47
  pres@gsnPanelBottom = 0.1
  gsn_panel(wks,plot(2:),(/1,2/),pres)

;  gsn_panel(wks,plot,(/2,2/),pres)
 
  lgres                    = True
  lgres@lgMonoDashIndex    = True
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

  gsn_legend_ndc(wks,nitems,labels,0.45,0.53,lgres)

  lgres2                    = True
  lgres2@lgMonoDashIndex    = True
  lgres2@lgDashIndex        = 1
  lgres2@lgLineThicknessF = 2.
;  lgres2@lgItemType         = "Lines"        ; show lines only (default)
  lgres2@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres2@vpWidthF           = 0.15           ; width of legend (NDC)
  lgres2@vpHeightF          = 0.02           ; height of legend (NDC)
  lgres2@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  ssmi_label = "Sea Ice Index "+obs_yr1+"-"+(obs_yr1+nyr_stats-1)
  gsn_legend_ndc(wks,1,ssmi_label,0.2,0.5,lgres2)

  lgres2@lgDashIndex        = 2
  ssmi_label = "Sea Ice Index "+(obs_yr2 - nyr_stats + 1)+"-"+obs_yr2
  gsn_legend_ndc(wks,1,ssmi_label,0.7,0.5,lgres2)

  frame(wks)

end

