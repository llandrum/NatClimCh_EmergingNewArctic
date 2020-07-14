; Plots histograms for 4 of the CMIP5 MMLE
; from monthly data
; LLandrum
; Apr 2020
;
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
; ================================================;
begin
  mname = (/"MPI-ESM","CanESM2","GFDL-CM3","GFDL-ESM2M"/)
  dirmname = (/"mpi","canesm2","gfdl_cm3","gfdl_esm2m"/)
  ncases = (/99,50,20,30/)
  case_yr1 = (/1920,1950,1920,1950/)                           ; some data for gridded data different than point data
  case_yr2 = (/2099,2100,2100,2100/)

  nens = dimsizes(mname)

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"

  plot_obs = True				; True to overlay mn obs for first and last decade of obs on min/max plots
  base_yr1 = 1950                       ; base yrs of LE for comparison
  base_yr2 = 1959

  obs_yr1 = 1979
  obs_yr2 = 2019
  fobs = "/glade/p/cesm/pcwg/landrum/SSMI/NH/NH_SIE_min_max_nsidc.from_monthly_sea_ice_index."+obs_yr1+"-"+obs_yr2+".nc"

  LE_plot_yr1 = (/1980,2010,2030,2050/)
  LE_plot_yr2 = LE_plot_yr1 + 9
  nplot = dimsizes(LE_plot_yr1)

  fLE_min_max = datadir+"monthly_NH_SIE_min_max_"+mname+"_LE."+ncases+"_sim."+case_yr1+"-"+case_yr2+".nc"

  psout = "SupplementaryFigure1.CMIP5_"+nens+"MMLE.min_max_NH_SIE."+nplot+"histograms"
;---------------------------------------------------------------------------------
; read in min/max data
;---------------------------------------------------------------------------------
  max_ncase = max(ncases)
  LE_yr1 = min(case_yr1)
  LE_yr2 = max(case_yr2)
  nyr_LE = LE_yr2 - LE_yr1 + 1

  min_ext_base = new((/nens,10*max_ncase/),float)
  max_ext_base = new((/nens,10*max_ncase/),float)

  max_ext_MMLE = new((/nens,max_ncase,nyr_LE/),float)
  min_ext_MMLE = new((/nens,max_ncase,nyr_LE/),float)

  ndata_max = 10*max_ncase
  min_data = new((/nens,nplot+1,ndata_max/),typeof(min_ext_base))
  max_data = new((/nens,nplot+1,ndata_max/),typeof(max_ext_base))

  do im=0,nens-1
     ic2 = 10*ncases(im) - 1
     i1 = base_yr1 - case_yr1(im)
     i2 = base_yr2 - case_yr1(im)
     f1 = addfile(fLE_min_max(im),"r")
     max_data(im,0,:ic2) = ndtooned(f1->max_ext(i1:i2,:))			; base data for each model
     min_data(im,0,:ic2) = ndtooned(f1->min_ext(i1:i2,:))

     do i=0,nplot-1
        i1 = LE_plot_yr1(i) - case_yr1(im)
        i2 = LE_plot_yr2(i) - case_yr1(im)
        LE_1d = ndtooned(f1->min_ext(i1:i2,:))
        min_data(im,i+1,:ic2) = (/ LE_1d /)
        LE_1d = ndtooned(f1->max_ext(i1:i2,:))
        max_data(im,i+1,:ic2) = (/ LE_1d /)
        delete(LE_1d)
     end do
  end do
;---------------------------------------------------------------------------------
; get obs if needed 
;---------------------------------------------------------------------------------
  if (plot_obs) then
     f = addfile(fobs,"r")
     obs_min_ext = f->min_ext
     obs_max_ext = f->max_ext
     year_obs = f->year
     nyr_obs = dimsizes(year_obs)

     obs_min_stats = new((/2,4/),float)
     obs_min_stats(0,:) = dim_stat4(obs_min_ext(:9))		; first decade
     obs_min_stats(1,:) = dim_stat4(obs_min_ext(nyr_obs-10:))	; last decade

     obs_max_stats = new((/2,4/),float)
     obs_max_stats(0,:) = dim_stat4(obs_max_ext(:9))               ; first decade
     obs_max_stats(1,:) = dim_stat4(obs_max_ext(nyr_obs-10:))      ; last decade

  end if
;---------------------------------------------------------------------------------
; plot
;---------------------------------------------------------------------------------
;  wks  = gsn_open_wks("pdf",psout)
  wks  = gsn_open_wks("eps",psout)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  ltgry    = NhlNewColor(wks,0.7,0.7,0.7)      ; add gray to map
  dkgry    = NhlNewColor(wks,0.4,0.4,0.4)
  blk    = NhlNewColor(wks,0.,0.,0.)      ; add gray to map

  nplots_per_hist = 2
  plot = new((/nens,2*nplot,nplots_per_hist/),graphic)

  ext_min_levels = fspan(0.,9.5,39)
  dext_min = ext_min_levels(1) - ext_min_levels(0)                  ; diff in levels same for min, max
  plot_colors = (/"black","blue","grey"/)
  hres = True
  hres@gsnDraw            = False
  hres@gsnFrame           = False
;  hres@gsnHistogramNumberOfBins = 20
;  hres@gsnHistogramSelectNiceIntervals = True

  xyres = True
  xyres@gsnDraw = False
  xyres@gsnFrame = False
  xyres@tiYAxisString = "Frequency (%)"
;  xyres@tiXAxisString = max_ext_LE@units
  xyres@tiXAxisString = "million sq km"

  xyres@gsnXYBarChart         = True                 ; turn on bar chart
;  xyres@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  xyres@gsnXYBarChartColors = "blue"
  xyres@gsnXYBarChartBarWidth = 0.4*(ext_min_levels(1)-ext_min_levels(0))  ; 40% of the width

  xyres@tmXBMode = "Explicit"
;  xyres@trYMaxF = 0.5
  xyres@trYMaxF = 0.35
  xyres@trYMinF = 0.
  xyres@tmYLMode = "Explicit"
  per_levels = fspan(0.,1.,21)
;  per_levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1./)
  xyres@tmYLValues = per_levels
  xyres@tmYLLabels = ""+(100.*per_levels)

;  res@gsEdgeColor = "black"
;  res@gsFillColor = (/"transparent"/)

;   res@tiYAxisString = "Frequency of occurence"

  nlevels_min = dimsizes(ext_min_levels)

;  xyres@tmXBLabelDirection = "Down"

  xyres@tfPolyDrawOrder    = "Predraw"               ; put line on top

  xyres@tiMainString = " "

  if (plot_obs) then
     gsres                   = True                        ; poly res
     gsres@gsFillColor       = "Grey"                 ; color chosen
     gsres@gsFillOpacityF = 0.25
     dum = new((/nens,8*nplot/),graphic)
     dumtxt = new((/nens,4*nplot/),graphic)

     lres = True
     lres@gsLineDashPattern = 1
;     lres@gsLineColor = "grey"
     lres@gsLineColor = "black"
     lres@gsLineThicknessF = 1.

     txres               = True                     ; text mods desired
     txres@txFontHeightF = 0.015                     ; font smaller. default big

     yp= (/0.,0.,.35,.35/)
  end if

  do im=0,nens-1
     x = (ext_min_levels(1:)+ext_min_levels(0:nlevels_min-2))*0.5         ; halfway to mid point of each range
     x2 = x - dext_min/8.
     x3 = x + dext_min/8.
     hres@gsnHistogramBinIntervals = ext_min_levels

     xyres@trXMaxF = ext_min_levels(nlevels_min-2)
     xyres@trXMinF = ext_min_levels(0)
     xyres@tmXBValues = ext_min_levels(0::2)
     xyres@tmXBLabels = ""+ext_min_levels(0::2)

     if (mname(im).eq."GFDL-ESM2M") then
        ext_max_levels = fspan(13.0,21.0,33)                                     ; was 9.5,17,5,33 for LENS
     else
        ext_max_levels = fspan(9.0,17.0,33)                                        ; was 9.5,17,5,33 for LENS
     end if
     dext_max = ext_max_levels(1) - ext_max_levels(0)                  ; diff in levels same for min, max

     ic2 = 10*ncases(im) - 1
     data = min_data(im,:1,:ic2)
     do i=0,nplot-1					; min SIE
        data(1,:) = min_data(im,i+1,:ic2)
        hist_plot = gsn_histogram(wks,data,hres)
        tmpy = hist_plot@NumInBins
        ybinned = int2flt(tmpy)
        ybinned(0,:) = ybinned(0,:)/dim_sum(ybinned(0,:))
        ybinned(1,:) = ybinned(1,:)/dim_sum(ybinned(1,:))
;        ybinned = ybinned/dim_sum(ybinned)

        xyres@tiMainString = LE_plot_yr1(i)+"-"+LE_plot_yr2(i)
        xyres@gsnXYBarChartColors = "blue"
        xyres@xyLineColors = "blue"
        plot(im,i,0) = gsn_csm_xy(wks,x2,ybinned(0,:),xyres)
        xyres@gsnXYBarChartColors = "red"
        xyres@xyLineColors = "red"
        plot(im,i,1) = gsn_csm_xy(wks,x3,ybinned(1,:),xyres)
        overlay(plot(im,i,0),plot(im,i,1))

        if (plot_obs) then
           lres@gsLineDashPattern = 1
           lres@gsLineColor = "black"
           dum(im,4*i) = gsn_add_polyline(wks,plot(im,i,0),(/obs_min_stats(0,0),obs_min_stats(0,0)/),(/0,.35/),lres)
           lres@gsLineDashPattern = 0
           dum(im,4*i+1) = gsn_add_polyline(wks,plot(im,i,0),(/obs_min_stats(1,0),obs_min_stats(1,0)/),(/0,.35/),lres)

           dumtxt(im,2*i) = gsn_add_text(wks,plot(im,i,0),"Sea Ice Index ("+year_obs(0)+"-"+year_obs(9)+")",obs_min_stats(0,0),0.33,txres)
           dumtxt(im,2*i+1) = gsn_add_text(wks,plot(im,i,0),"Sea Ice Index ("+year_obs(nyr_obs-10)+"-"+year_obs(nyr_obs-1)+")",obs_min_stats(1,0),0.30,txres)

           obs_min = min(obs_min_ext(:9))
           obs_max = max(obs_min_ext(:9))
           xp_min = (/obs_min,obs_max,obs_max,obs_min/)
           dum(im,4*i+2) = gsn_add_polygon (wks,plot(im,i,0),xp_min,yp,gsres)

           obs_min = min(obs_min_ext(nyr_obs-10:))
           obs_max = max(obs_min_ext(nyr_obs-10:))
           xp_min = (/obs_min,obs_max,obs_max,obs_min/)
           dum(im,4*i+3) = gsn_add_polygon (wks,plot(im,i,0),xp_min,yp,gsres)
        end if
     end do
     delete(tmpy)
     delete(ybinned)
     delete(hres@gsnHistogramBinIntervals)
     delete(x)
     delete(x2)
     delete(x3)
     delete(xyres@tmXBValues)
     delete(xyres@tmXBLabels)

     nlevels = dimsizes(ext_max_levels)
     hres@gsnHistogramBinIntervals = ext_max_levels

     x = (ext_max_levels(1:)+ext_max_levels(0:nlevels-2))*0.5         ; halfway to mid point of each range
     x2 = x - dext_max/8.
     x3 = x + dext_max/8.

     xyres@trXMaxF = ext_max_levels(nlevels-2)
     xyres@trXMinF = ext_max_levels(0)
     xyres@tmXBValues = ext_max_levels(0::2)
     xyres@tmXBLabels = ""+ext_max_levels(0::2)
     xyres@gsnXYBarChartBarWidth = 0.4*(ext_max_levels(1)-ext_max_levels(0))  ; 40% of the width

     data = max_data(im,:1,:ic2)
     do i=0,nplot-1                                        ; max SIE
        data(1,:) = max_data(im,i+1,:ic2)
        hist_plot = gsn_histogram(wks,data,hres)
        tmpy = hist_plot@NumInBins
        ybinned = int2flt(tmpy)
        ysum0 = dim_sum(ybinned(0,:))
        ysum1 = dim_sum(ybinned(1,:))

        if ((ysum0.gt.0.).and.(ysum1.gt.0.)) then
           ybinned(0,:) = ybinned(0,:)/ysum0
           ybinned(1,:) = ybinned(1,:)/ysum1
;           ybinned = ybinned/dim_sum(ybinned)

           xyres@tiMainString = LE_plot_yr1(i)+"-"+LE_plot_yr2(i)
           xyres@gsnXYBarChartColors = "blue"
           xyres@xyLineColors = "blue"
           plot(im,nplot+i,0) = gsn_csm_xy(wks,x2,ybinned(0,:),xyres)
           xyres@gsnXYBarChartColors = "red"
           xyres@xyLineColors = "red"
           plot(im,nplot+i,1) = gsn_csm_xy(wks,x3,ybinned(1,:),xyres)
           overlay(plot(im,nplot+i,0),plot(im,nplot+i,1))
        end if
        if (plot_obs) then
           lres@gsLineDashPattern = 1
           lres@gsLineColor = "black"
           dum(im,4*(i+nplot)) = gsn_add_polyline(wks,plot(im,nplot+i,0),(/obs_max_stats(0,0),obs_max_stats(0,0)/),(/0,.35/),lres)
           lres@gsLineDashPattern = 0
           dum(im,4*(i+nplot)+1) = gsn_add_polyline(wks,plot(im,nplot+i,0),(/obs_max_stats(1,0),obs_max_stats(1,0)/),(/0,.35/),lres)

           dumtxt(im,2*(i+nplot)) = gsn_add_text(wks,plot(im,nplot+i,0),"Sea Ice Index ("+year_obs(0)+"-"+year_obs(9)+")",obs_max_stats(0,0),0.33,txres)
           dumtxt(im,2*(i+nplot)+1) = gsn_add_text(wks,plot(im,nplot+i,0),"Sea Ice Index ("+year_obs(nyr_obs-10)+"-"+year_obs(nyr_obs-1)+")",obs_max_stats(1,0),0.30,txres)
   
           obs_min = min(obs_max_ext(:9))
           obs_max = max(obs_max_ext(:9))
           xp_min = (/obs_min,obs_max,obs_max,obs_min/)
           dum(im,4*(i+nplot)+2) = gsn_add_polygon (wks,plot(im,nplot+i,0),xp_min,yp,gsres)

           obs_min = min(obs_max_ext(nyr_obs-10:))
           obs_max = max(obs_max_ext(nyr_obs-10:))
           xp_min = (/obs_min,obs_max,obs_max,obs_min/)
           dum(im,4*(i+nplot)+3) = gsn_add_polygon (wks,plot(im,nplot+i,0),xp_min,yp,gsres)
        end if
     end do
     delete(tmpy)
     delete(ybinned)
     delete(hres@gsnHistogramBinIntervals)
     delete(x)
     delete(x2)
     delete(x3)
     delete(xyres@tmXBValues)
     delete(xyres@tmXBLabels)

     delete(data)
  end do				; end im

  pres = True
  pres@gsnFrame = False
  panel_diff = ((0.95 - 0.05)/2) - 0.16
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.015

  pres@gsnPanelTop = 0.95
  pres@gsnPanelBottom = (pres@gsnPanelTop - 0.5*panel_diff)
  pres@gsnPanelRight = 0.5
  pres@gsnPanelMainPosXF = 0.12
  do im=0,1
     pres@gsnPanelMainString = mname(im)+"~C~ ~Z50~MIN NH SIE"
     gsn_panel(wks,plot(im,:nplot-1,0),(/1,nplot/),pres)
     pres@gsnPanelTop = pres@gsnPanelBottom - 0.04
     pres@gsnPanelBottom = pres@gsnPanelTop - 0.5*panel_diff
     pres@gsnPanelMainString = "~Z50~MAX NH SIE"
     gsn_panel(wks,plot(im,nplot:,0),(/1,nplot/),pres)
     pres@gsnPanelTop = pres@gsnPanelBottom - 0.08
     pres@gsnPanelBottom = pres@gsnPanelTop - 0.5*panel_diff
  end do
  pres@gsnPanelMainPosXF = 0.61
  pres@gsnPanelTop = 0.95
  pres@gsnPanelBottom = (pres@gsnPanelTop - 0.5*panel_diff)
  pres@gsnPanelRight = 1.
  pres@gsnPanelLeft = 0.5
  do im=2,3
     pres@gsnPanelMainString = mname(im)+"~C~ ~Z50~MIN NH SIE"
     gsn_panel(wks,plot(im,:nplot-1,0),(/1,nplot/),pres)
     pres@gsnPanelTop = pres@gsnPanelBottom - 0.04
     pres@gsnPanelBottom = pres@gsnPanelTop - 0.5*panel_diff
     pres@gsnPanelMainString = "~Z50~MAX NH SIE"
     gsn_panel(wks,plot(im,nplot:,0),(/1,nplot/),pres)
     pres@gsnPanelTop = pres@gsnPanelBottom - 0.08
     pres@gsnPanelBottom = pres@gsnPanelTop - 0.5*panel_diff
  end do

  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.18           ; labelbar width
  lbres@vpHeightF          = 0.08           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbPerimOn          = False

  colors = (/"blue","red"/)
  labels = (/"("+base_yr1+"-"+base_yr2+")", "MMLE future yrs"/)
;  xpos = (/0.1, 0.35/)
  xpos = (/0.3, 0.55/)
  do i=0,1
     lbres@lbFillColors = colors(i)
;     gsn_labelbar_ndc(wks,1,labels(i),xpos(i),0.9,lbres)
     gsn_labelbar_ndc(wks,1,labels(i),xpos(i),0.18,lbres)
  end do	
  frame(wks)
end

