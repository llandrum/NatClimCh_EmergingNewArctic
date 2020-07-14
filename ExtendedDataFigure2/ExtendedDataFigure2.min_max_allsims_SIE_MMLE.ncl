; Plots ensemble mean min, max SIE, and also
; polygons showing variability for each CMIP5 ensemble
; and overlays the Walsh et al extended sea ice min, max
;
; LLandrum
; April 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
; ================================================;
begin
  mname = (/"MPI-ESM","CanESM2","CESM","GFDL-CM3","GFDL-ESM2M"/)
  dirmname = (/"mpi","canesm2","LE","gfdl_cm3","gfdl_esm2m"/)
  ncases = (/99,50,40,20,30/)
  case_yr1 = (/1920,1950,1920,1920,1950/)
  case_yr2 = (/2099,2100,2100,2100,2100/)
  nens = dimsizes(mname)

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"
  iCESM = ind(mname.eq."CESM")
  datadir(iCESM) = "/glade/p/cesm/pcwg/landrum/LE/ice/pre_process/monthly/"

  remove_mean = False				; True to plot after removing mean, False to plot absolute
  mn_yr1 = 1950
  mn_yr2 = 1990
  plot_obs = True                              ; True to overlay mn obs for first and last decade of obs on min/max plots
  add_PDF = True				; True to add plots of PDFs of trends
  nyr_trnds = 20
  trnds_yr2 = case_yr2 - nyr_trnds + 1

  pdf_yr1 = 1920                                ; yrs over which to calculate pdfs
  pdf_yr2 = 1980                                ; note this is last yr of nyr_trnds trends...
  nyr = pdf_yr2 - pdf_yr1 + 1                   ; some cases begin later (e.g. 1950) than others (e.g. 1920)
  ftrnds = datadir+"monthly_NH_SIE_min_max_"+mname+"_LE."+ncases+"_sim."+nyr_trnds+"yr_trnds."+case_yr1+"-"+trnds_yr2+".nc"

  obs_yr1 = 1979
  obs_yr2 = 2019
  nyr_obs = obs_yr2 - obs_yr1 + 1
  fobs = "/glade/p/cesm/pcwg/landrum/SSMI/NH/NH_SIE_min_max_nsidc.from_monthly_sea_ice_index."+obs_yr1+"-"+obs_yr2+".nc"

  Walshdir = "/glade/p/cesm/pcwg/landrum/NSIDC/Arctic_extended_SeaIce_Walsh/"
  Walsh_yr1 = 1850                                                    ; start and end yrs for time series
  Walsh_yr2 = 2017
  nyrs_Walsh = Walsh_yr2 - Walsh_yr1 + 1
  styrs_proc = sprinti("%0.4i",Walsh_yr1) + "-"+sprinti("%0.4i",Walsh_yr2)
  fice_ext = Walshdir+"/ice_ext_MHolland_regions_Walsh_extended_ArcticSIC."+styrs_proc+".nc"
 

  overlay_emerg = False				; True to overlay line at emergence year, false otherwise

  plot_yr1 = 1850
  plot_yr2 = 2100 
  nyr = plot_yr2 - plot_yr1 + 1
  year = ispan(plot_yr1,plot_yr2,1)

  fLE_min_max = datadir+"monthly_NH_SIE_min_max_"+mname+"_LE."+ncases+"_sim."+case_yr1+"-"+case_yr2+".nc"
  min_var = "min_ext"
  max_var = "max_ext"

  ice_free_min = 1.					; definition of ice-free minimum (in units of 10^6 km^2)
  psout = "SupplementaryFigure3.min_max_SIE_allsims_ts.CMIP5_MMLE"
  if (plot_obs) then
     psout = psout+".NSIDC"
  end if
  if (remove_mean) then
     psout = psout + "."+mn_yr1+"-"+mn_yr2+".mn_removed"
  end if

;---------------------------------------------------------------------------------
; read in LE data
;---------------------------------------------------------------------------------
  CMIP5_avg_min = new((/nens+2,nyr/),float)						; each ensemble member, Walsh extended and NSIDC ice index
  CMIP5_avg_min!0 = "model"
  CMIP5_avg_min!1 = "year"
  CMIP5_avg_min&year = year

  CMIP5_avg_max = CMIP5_avg_min

  min_MIN = new((/nens,nyr/),float)				; to show ensemble spread in 10-yr running mean MIN
  min_MIN!0 = "model"
  min_MIN!1 = "year"
  min_MIN&year = year
  max_MIN = min_MIN

  min_MAX = min_MIN						 ; to show ensemble spread in 10-yr running mean MAX
  max_MAX = min_MIN

  do ie = 0,nens-1
     i1 = case_yr1(ie) - plot_yr1
     i2 = case_yr2(ie) - plot_yr1

     nyr_data = case_yr2(ie) - case_yr1(ie) + 1

     max_ext_LE = new((/nyr_data,ncases(ie)/),float)
     min_ext_LE = new((/nyr_data,ncases(ie)/),float)

     f1 = addfile(fLE_min_max(ie),"r")
     max_ext_LE = (/ f1->$max_var$ /) 
     min_ext_LE = (/ f1->$min_var$ /) 

     CMIP5_avg_min(ie,i1:i2) = (/ dim_avg_n(min_ext_LE,1) /)
     CMIP5_avg_max(ie,i1:i2) = (/ dim_avg_n(max_ext_LE,1) /)

     min_MIN(ie,i1:i2) = (/ dim_min_n(min_ext_LE,1) /)						; min 10-yr mean from ensemble
     max_MIN(ie,i1:i2) = (/ dim_max_n(min_ext_LE,1) /)

     min_MAX(ie,i1:i2) = (/ dim_min_n(max_ext_LE,1) /)
     max_MAX(ie,i1:i2) = (/ dim_max_n(max_ext_LE,1) /)

     tmp = f1->$max_var$
     max_ext_LE@units = tmp@units
     min_ext_LE@units = tmp@units
     delete(tmp)

     delete(max_ext_LE)
     delete(min_ext_LE)
  end do					; end ie
  tmp = f1->$max_var$
  CMIP5_avg_max@units = tmp@units
  CMIP5_avg_min@units = tmp@units
  min_MIN@units = tmp@units
  max_MIN@units = tmp@units
  min_MAX@units = tmp@units 
  max_MAX@units = tmp@units 
  delete(tmp)

;---------------------------------------------------------------------------------
; Walsh data
;---------------------------------------------------------------------------------
  f = addfile(fice_ext,"r")
  ext = f->ext_mo_nh

  i1 = Walsh_yr1 - plot_yr1

  do i=0,nyrs_Walsh-1
     ii1 = 12*i
     ii2 = ii1+11
     CMIP5_avg_min(nens,i1+i) = min(ext(ii1:ii2))
     CMIP5_avg_max(nens,i1+i) = max(ext(ii1:ii2))
  end do
  delete(ext)

  if (remove_mean) then
     i1 = mn_yr1 - plot_yr1
     i2 = mn_yr2 - plot_yr1
     min_mean = dim_avg_n(CMIP5_avg_min(:,i1:i2),1)
     max_mean = dim_avg_n(CMIP5_avg_max(:,i1:i2),1)

     CMIP5_avg_min = CMIP5_avg_min - conform(CMIP5_avg_min,min_mean,0)
     CMIP5_avg_max = CMIP5_avg_max - conform(CMIP5_avg_max,max_mean,0)

     min_MIN = min_MIN - conform(min_MIN,min_mean(:nens-1),0)
     max_MIN = max_MIN - conform(max_MIN,min_mean(:nens-1),0)
     min_MAX = min_MAX - conform(min_MAX,max_mean(:nens-1),0)
     max_MAX = max_MAX - conform(max_MAX,max_mean(:nens-1),0)

  end if
;---------------------------------------------------------------------------------
; trend data if plotting pdfs
;---------------------------------------------------------------------------------
  if (add_PDF) then
     max_ncase = max(ncases)
     LE_yr1 = min(case_yr1)
     LE_yr2 = max(case_yr2)
     nyr_LE = LE_yr2 - LE_yr1 + 1

     min_trnds = new((/nens,nyr,max_ncase/),float)
     max_trnds = new((/nens,nyr,max_ncase/),float)

     do im=0,nens-1
        ic2 = ncases(im) - 1
        if (case_yr1(im).le.pdf_yr1) then
           ii1 = pdf_yr1 - case_yr1(im)                    ; for reading case data
           i1 = 0                                          ; for saving to data
        else
           i1 = case_yr1(im) - pdf_yr1
           ii1 = 0
        end if
        i2 = pdf_yr2 - pdf_yr1
        ii2 = pdf_yr2 - case_yr1(im)

        f1 = addfile(ftrnds(im),"r")
        max_trnds(im,i1:i2,:ic2) = f1->max_ext_trnds(ii1:ii2,:)
        min_trnds(im,i1:i2,:ic2) = f1->min_ext_trnds(ii1:ii2,:)
     end do
  end if
;---------------------------------------------------------------------------------
; plot
;---------------------------------------------------------------------------------
;  wks  = gsn_open_wks("pdf",psout)
  wks  = gsn_open_wks("eps",psout)

  gsn_define_colormap(wks,"posneg_1")
  colors = (/5,8,17,15,10,3/)

  plot = new(2,graphic)
  poly_min = new(nens,graphic)
  poly_max = new(nens,graphic)

  res = True				; res for stats plot
  res@gsnDraw = False
  res@gsnFrame = False

;  res@vpWidthF    = 0.5      ; Make plots wider than
;  res@vpHeightF   = 0.2      ; they are high.

  res@vpWidthF    = 0.4      ; Make plots wider than
  res@vpHeightF   = 0.19      ; they are high.

  res@tiXAxisString = "year"
  res@trXMinF = plot_yr1
  res@trXMaxF = plot_yr2

  res@xyMonoDashPattern = True
  res@xyLineColors = colors

  res@xyLineThicknessF = 5.
  res@tiYAxisString    = "10~S~6~N~ km~S~2~N~"
  res@tfPolyDrawOrder    = "Predraw"               ; put line on top 

  res@gsnLeftString    =  "MIN SIE"
  if (remove_mean) then
     res@gsnRightString    =  "("+mn_yr1+"-"+mn_yr2+") mean removed"
  end if
  plot(0)  = gsn_csm_xy(wks,CMIP5_avg_min&year,CMIP5_avg_min,res)  ; Create filled XY plot.

  res@gsnLeftString    =  "MAX SIE"
  plot(1)  = gsn_csm_xy(wks,CMIP5_avg_max&year,CMIP5_avg_max,res)  ; Create filled XY plot.

  res@gsnLeftString = ""
  res@gsnRightString = ""

  pres                   = True                        ; poly res
  pres@gsFillOpacityF = 0.25

  do i=0,nens-1
     pres@gsFillColor       = colors(i)
     nyr_case = case_yr2(i) - case_yr1(i) + 1
     xp = new(2*nyr_case,integer)
     yp = new(2*nyr_case,float)
     i1 = case_yr1(i) - plot_yr1
     i2 = case_yr2(i) - plot_yr1

     xp(:nyr_case-1) = year(i1:i2)
     xp(nyr_case:) = year(i1:i2:-1)
     yp(:nyr_case-1) = max_MIN(i,i1:i2)
     yp(nyr_case:) = min_MIN(i,i1:i2:-1)
     poly_min(i) = gsn_add_polygon (wks,plot(0),xp,yp,pres)

     yp(:nyr_case-1) = max_MAX(i,i1:i2:)
     yp(nyr_case:) = min_MAX(i,i1:i2:-1)
     poly_max(i) = gsn_add_polygon (wks,plot(1),xp,yp,pres)

     delete(xp)
     delete(yp)
  end do

  if (plot_obs) then
     f1 = addfile(fobs,"r")
     min_SSMI = f1->min_ext
     max_SSMI = f1->max_ext

     print("NSIDC 1979 mean       min     max   2009 mean     max     diff_min     diff_max")
     y2 = nyr_obs - 10
     obs1 = dim_avg(min_SSMI(:0))
     obs2 = dim_avg(min_SSMI(y2:))
     obs_diff = obs2 - obs1
     obs3 = dim_avg(max_SSMI(:0))
     obs4 = dim_avg(max_SSMI(y2:))
     obs_diff2 = obs4 - obs3

     print(obs1+"       "+obs3+"     "+obs2+"      "+obs4+"        "+obs_diff+"        "+obs_diff2)

     if (remove_mean) then
        i1 = mn_yr2 - obs_yr1
        min_SSMI = min_SSMI - dim_avg(min_SSMI(:i1))
        max_SSMI = max_SSMI - dim_avg(max_SSMI(:i1))
     end if

     obs_data = new(nyr,float)
     obs_data!0 = "year"
     obs_data&year = year
     i1 = obs_yr1 - plot_yr1
     i2 = obs_yr2 - plot_yr1

     delete(res@xyLineColors)
;     res@xyLineColors = "grey"
     res@xyLineColors = 12

     res@tiYAxisString    = " "
     res@gsnLeftString    =  " "
;     res@xyDashPattern = 1

     obs_data(i1:i2) = (/ min_SSMI /)
     plot_obs_min = gsn_csm_xy(wks,obs_data&year,obs_data,res)
     overlay(plot(0),plot_obs_min)
 
     obs_data(i1:i2) = (/ max_SSMI /)
     plot_obs_max = gsn_csm_xy(wks,obs_data&year,obs_data,res)
     overlay(plot(1),plot_obs_max)
  end if

  if (add_PDF) then
     plot_PDF = new(2,graphic)
     dum0 = new(2,graphic)
     dum_minus = new(2,graphic)
     dum_plus = new(2,graphic)
     nbin = 13

     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@tiYAxisString = "PDF (%)"
     xyres@tiXAxisString = "10~S~6~N~ km~S~2~N~/yr"

     xyres@tiYAxisFontHeightF = 0.005
     xyres@tiXAxisFontHeightF = 0.005
     xyres@tiYAxisFont = "helvetica-bold"
     xyres@tiXAxisFont = "helvetica-bold"                  

     xyres@xyMonoDashPattern = True
     xyres@xyLineColors = colors
     xyres@xyLineThicknessF = 4.

     xyres@tiMainString = " "
     xyres@vpWidthF           = 0.09
     xyres@vpHeightF          = 0.09
;     xyres@vpWidthF           = 0.15
;     xyres@vpHeightF          = 0.15

     xyres@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
     xyres@gsnStringFontHeightF   = 0.005

     do im=0,nens-1
        ic2 = ncases(im) - 1
        if (case_yr1(im).le.pdf_yr1) then
           i1 = 0                                          ; for saving to data
        else
           i1 = case_yr1(im) - pdf_yr1
        end if
        tmp_min = min_trnds(im,i1:,:ic2)
        pdf_min = pdfx(tmp_min,nbin,False)

        tmp_max = max_trnds(im,i1:,:ic2)
        pdf_max = pdfx(tmp_max,nbin,False)

        if (im.eq.0) then
           xx_min = new((/nens,nbin/),typeof(pdf_min))
           yy_min = new((/nens,nbin/),typeof(pdf_min))

           xx_max = new((/nens,nbin/),typeof(pdf_min))
           yy_max = new((/nens,nbin/),typeof(pdf_min))
        end if
        xx_min(im,:) = pdf_min@bin_center
        yy_min(im,:) = (/ pdf_min /)
        xx_max(im,:) = pdf_max@bin_center
        yy_max(im,:) = (/ pdf_max /)
        delete(tmp_min)
        delete(tmp_max)
     end do                                ; end im
     xyres@trXMaxF = 0.12
     xyres@trXMinF = -0.12
     xyres@vpXF = 0.2
     xyres@vpYF = 0.7
     xyres@gsnCenterString = "20~S~th~N~C "+nyr_trnds+"-yr trends"
     plot_PDF(0) = gsn_csm_xy (wks, xx_min, yy_min, xyres)

     xyres@vpXF = 0.2
     xyres@vpYF = 0.21
     plot_PDF(1) = gsn_csm_xy (wks, xx_max, yy_max, xyres)

     polyres = True
     polyres@gsLineColor = "grey"
     polyres@gsLineThicknessF = 2.
;     polyres@tfPolyDrawOrder = "PostDraw"

     xminus = -1.*1./20.				;-1 million km^2 in 20 yrs
     xplus = 1.*1./20.				; +1 million km^2 in 20 yrs

;     draw(plot_PDF(0))
;     draw(plot_PDF(1))

     do i=0,1
        polyres@gsLineDashPattern = 0
        dum0(i) = gsn_add_polyline(wks,plot(i),(/0,0/),(/0,30/),polyres)
        polyres@gsLineDashPattern = 1
;        gsn_polyline(wks,plot(i),(/xminus,xminus/),(/0,30/),polyres)
;        gsn_polyline(wks,plot(i),(/xplus,xplus/),(/0,30/),polyres)
        dum_minus(i) = gsn_add_polyline(wks,plot(i),(/xminus,xminus/),(/0,30/),polyres)
        dum_plus(i) = gsn_add_polyline(wks,plot(i),(/xplus,xplus/),(/0,30/),polyres)
     end do

     amres                  = True
     amres@amParallelPosF   = -0.218    ; -0.5 is the left edge of the plot.
     amres@amOrthogonalPosF = -0.08   ; -0.5 is the top edge of the plot.
     amres@amJust           = "TopRight"

     min_anno = gsn_add_annotation(plot(0), plot_PDF(0), amres)  ; Attach PDF plot to time series plot
     max_anno = gsn_add_annotation(plot(1), plot_PDF(1), amres)  ; Attach PDF plot to time series plot

  end if

  lgres                    = True
  if (plot_obs) then
     ncolors = dimsizes(colors)
     colors2 = new(ncolors+1,typeof(colors))
     colors2(2:) = colors(:ncolors-2)
     colors2(0) = 12
     colors2(1) = colors(ncolors-1)				; Walsh is stored after CMIP5 MMLE as last in array
     lgres@lgLineColors       = colors2
     lgres@lgMonoDashIndex    = True
     lgres@lgDashIndex = 0
     labels = new((nens+2),typeof(mname))
     labels(0) = "NSIDC Sea Ice Index"
     labels(1) = "Walsh et al"
     labels(2:) = mname
  else
     lgres@lgMonoDashIndex    = True
     lgres@lgDashIndex        = 0
     lgres@lgLineColors       = colors
     labels = new((nens+1),typeof(mname))
     labels(0) = "Walsh et al"
     labels(1:) = mname
  end if
  lgres@lgLineThicknessF = 3.
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.18           ; width of legend (NDC)
  lgres@vpHeightF          = 0.12           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter

  nitems = dimsizes(labels)                 ; number of legend items

  pres = True
  pres@gsnFrame = False
  gsn_panel(wks,plot,(/2,1/),pres)
  if (add_PDF) then
     gsn_legend_ndc(wks,nitems,labels,0.72,0.94,lgres)
  else
     gsn_legend_ndc(wks,nitems,labels,0.25,0.25,lgres)
  end if
  frame(wks)
end

