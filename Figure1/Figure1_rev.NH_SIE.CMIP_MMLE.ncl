; Nature Climate Change
; Emerging New Arctic
; Figure1 for New Arctic paper.
; LLandrum
; July 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
; ================================================;
begin
  mname = (/"MPI-ESM","CanESM2","CESM","GFDL-CM3","GFDL-ESM2M"/)
  dirmname = (/"mpi","canesm2","LE","gfdl_cm3","gfdl_esm2m"/)
  ncases = (/99,50,40,20,30/)
  case2ID = (/"034","44","103","3","11"/)
  case2_plot_yr = (/2012,2012,2017,2017,2018/)
  case3ID = (/"055","45","28","18","30"/)
  case3_plot_yr = (/2018,2011,2013,2015,2014/)
  case2_plot_month = 9
  case3_plot_month = 9
  nens = dimsizes(mname)

  ex_mname = "CESM"

  sea_ice_var = "ext"
  remove_mean = True                            ; True to plot after removing mean, False to plot absolute
  min_var = "min_"+sea_ice_var
  max_var = "max_"+sea_ice_var

  iexample = ind(mname.eq.ex_mname)
  ex_dirmname = dirmname(iexample)
  ex_ncases = ncases(iexample)

  impi = ind(dirmname.eq."mpi")
  igfdl_cm3 = ind(dirmname.eq."gfdl_cm3")
  igfdl_esm2m = ind(dirmname.eq."gfdl_esm2m")
  iCESM = ind(mname.eq."CESM")
  mtitle = mname
  mtitle(iCESM) = "CESM1-CAM5"
  ex_title = mtitle(iexample)

  case_yr1 = new(nens,integer)
  case_yr1 = 1950
  case_yr1(iCESM) = 1920
  case_yr1(igfdl_cm3) = 1920
  case_yr1(impi) = 1920
  case_yr2 = case_yr1
  case_yr2 = 2100
  case_yr2(impi) = 2099
  ex_case_yr1 = case_yr1(iexample)
  ex_case_yr2 = case_yr2(iexample)
  ex_case2_plot_yr = case2_plot_yr(iexample)
  ex_case3_plot_yr = case3_plot_yr(iexample)
  ex_case2ID = case2ID(iexample)
  ex_case3ID = case3ID(iexample)

  if (ex_mname.eq."MPI-ESM") then
     ex_sic_case_yr1 = 2006				; MPI historical ends 2005 then continues in different rcp scenarios
     case_post = "_rcp85_r"
     case2 = ex_mname+case_post+ex_case2ID+"i2005p3"
     case3 = ex_mname+case_post+ex_case3ID+"i2005p3"
  else
     ex_sic_case_yr1 = ex_case_yr1
     case_post = "_historical_rcp85_r"
     case2 = ex_mname+case_post+ex_case2ID+"i1p1"
     case3 = ex_mname+case_post+ex_case3ID+"i1p1"
  end if
  if (ex_mname.eq."CESM") then
     sic_name = "aice"
     ex_datadir = "/glade/p/cesm/pcwg/landrum/LE/ice/pre_process/monthly/"

     case2 = "b.e11.BRCP85C5CNBDRD.f09_g16.103"
     case2_label = "LE103"
     ex_sic_case_yr1 = 2006
     ex_sic_case2_yr2 = 2100
     ex_case2_plot_yr = 2017

     case3 = "b.e11.BRCP85C5CNBDRD.f09_g16.028"
     case3_label = "LE028"
     ex_sic_case3_yr2 = 2080

     LE_datadir = "/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/aice/"
     fsic2 = LE_datadir+case2+".cice.h.aice_nh."+ex_sic_case_yr1+"01-"+ex_sic_case2_yr2+"12.nc"
     fsic3 = LE_datadir+case3+".cice.h.aice_nh."+ex_sic_case_yr1+"01-"+ex_sic_case3_yr2+"12.nc"
  else
     sic_name = "sic"
     ex_datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+ex_dirmname+"/"
     case2_label = ex_mname+"_"+ex_case2ID
     case3_label = ex_mname+"_"+ex_case3ID

     LE_datadir = "/glade/collections/cdg/data/CLIVAR_LE/"+ex_dirmname+"_lens/OImon/sic"
     fsic2 = LE_datadir+"/sic_OImon_"+case2+"_"+ex_sic_case_yr1+"01-"+ex_case_yr2+"12.nc"
     fsic3 = LE_datadir+"/sic_OImon_"+case3+"_"+ex_sic_case_yr1+"01-"+ex_case_yr2+"12.nc"
  end if

  plot_obs = True				; True to overlay mn obs for first and last decade of obs on min/max plots
  base_yr1 = 1950                       ; base yrs of LE for comparison
  base_yr2 = 1959

  obs_yr1 = 1979
  obs_yr2 = 2019
  fobs = "/glade/p/cesm/pcwg/landrum/SSMI/NH/NH_SIE_min_max_nsidc.from_monthly_sea_ice_index."+obs_yr1+"-"+obs_yr2+".nc"

  nyr_LE_data1 = ex_case_yr2 - ex_case_yr1 + 1
  nyr_LE_data_all = ex_case_yr2 - ex_case_yr1 + 1

  LE_plot_yr1 = (/2010,2030,2050/)
  LE_plot_yr2 = LE_plot_yr1 + 9
  nplot = dimsizes(LE_plot_yr1)

  fex_min_max = ex_datadir+"monthly_NH_SIE_min_max_"+ex_mname+"_LE."+ex_ncases+"_sim."+ex_case_yr1+"-"+ex_case_yr2+".nc"

  LE_min_yr1 = 2010				; years from which to pick min and max SIE
  LE_min_yr2 = 2019

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"
  datadir(iCESM) = "/glade/p/cesm/pcwg/landrum/LE/ice/pre_process/monthly/"
  fLE_min_max = datadir+"monthly_NH_SIE_min_max_"+mname+"_LE."+ncases+"_sim."+case_yr1+"-"+case_yr2+".nc"

  stats_yr1 = 1950
  stats_yr2 = 2091                              ; must be at least nyr_stats-1 less than case2_yr2
;  stats_yr2 = 2090
  nyr_stats = 10
  nyr = stats_yr2 - stats_yr1 + 1

  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

  ice_free_min = 1.                                     ; definition of ice-free minimum (in units of 10^6 km^2)

  psout = "Figure1_rev.monthly_NH_SIE_min_max.CMIP5_"+nens+"MMLE"
;---------------------------------------------------------------------------------
; read in MMLE data for ts plots
;---------------------------------------------------------------------------------
  stats_min = new((/nens,nyr,4/),float)
  stats_min!0 = "model"
  stats_min!1 = "year"
  stats_min&year = ispan(stats_yr1,stats_yr2,1)
  stats_min!2 = "moment"

  stats_max = stats_min
  stats_min_month = stats_min
  stats_max_month = stats_max
  stats_min_nmonths = stats_max
  stats_icefree_month = stats_min
  stats_icefree_nmonths = stats_max
  do ie = 0,nens-1
     i1 = stats_yr1 - case_yr1(ie)
     if (i1.ge.0) then                                          ; possible to use stat_yr1 before an ensemble begins. Added this to check if 1920 vs 1950
                                                                ; makes a significant difference in ToE for those ensembles that start 1920
        nyr_data = case_yr2(ie) - stats_yr1 + 1

        max_ext_LE = new((/nyr_data,ncases(ie)/),float)
        min_ext_LE = new((/nyr_data,ncases(ie)/),float)
        max_month_LE = new((/nyr_data,ncases(ie)/),float)
        min_month_LE = new((/nyr_data,ncases(ie)/),float)
        min_nmonths_LE = new((/nyr_data,ncases(ie)/),float)
        icefree_month_LE = new((/nyr_data,ncases(ie)/),float)
        icefree_nmonths_LE = new((/nyr_data,ncases(ie)/),float)

        f1 = addfile(fLE_min_max(ie),"r")
        max_ext_LE = (/ f1->$max_var$(i1:,:) /)
        min_ext_LE = (/ f1->$min_var$(i1:,:) /)
        max_month_LE = (/ int2flt(f1->max_month(i1:,:)) /)
        min_month_LE = (/ int2flt(f1->min_month(i1:,:)) /)
        min_nmonths_LE = (/ int2flt(f1->min_nmonths(i1:,:)) /)
        icefree_month_LE = (/ int2flt(f1->icefree_month(i1:,:)) /)
        icefree_nmonths_LE = (/ int2flt(f1->icefree_nmonths(i1:,:)) /)

        tmp = f1->$max_var$
        max_ext_LE@units = tmp@units
        min_ext_LE@units = tmp@units
        delete(tmp)
;------------------------------------------------------------------------------
; calculate statistics
;------------------------------------------------------------------------------
        nyr_last = nyr_data - nyr_stats + 1
        do i=0,nyr_last-1
           i2 = i+nyr_stats-1
           stats_min(ie,i,:) = (/ dim_stat4_n(ndtooned(min_ext_LE(i:i2,:)),0) /)
           stats_min_month(ie,i,:) = (/ dim_stat4_n(ndtooned(min_month_LE(i:i2,:)),0) /)
           stats_max(ie,i,:) = (/ dim_stat4_n(ndtooned(max_ext_LE(i:i2,:)),0) /)
           stats_max_month(ie,i,:) = (/ dim_stat4_n(ndtooned(max_month_LE(i:i2,:)),0) /)
           stats_min_nmonths(ie,i,:) = (/ dim_stat4_n(ndtooned(min_nmonths_LE(i:i2,:)),0) /)
           if (stats_min(ie,i,0).lt.ice_free_min) then
              stats_icefree_month(ie,i,:) = (/ dim_stat4_n(ndtooned(icefree_month_LE(i:i2,:)),0) /)
              stats_icefree_nmonths(ie,i,:) = (/ dim_stat4_n(ndtooned(icefree_nmonths_LE(i:i2,:)),0) /)
           end if
        end do
        delete(max_ext_LE)
        delete(min_ext_LE)
        delete(max_month_LE)
        delete(min_month_LE)
        delete(min_nmonths_LE)
        delete(icefree_month_LE)
        delete(icefree_nmonths_LE)
      end if
  end do                                        ; end ie

  stats = new((/nens,5,nyr,4/),float)
  stats(:,0,:,:) = stats_min
  stats(:,1,:,:) = stats_max
  stats(:,2,:,:) = stats_min_month
  stats(:,3,:,:) = stats_max_month
  stats(:,4,:,:) = stats_icefree_month

  if (sea_ice_var.eq."ext") then
     stats_title = (/"Minimum Sea Ice Extent","Maximum Sea Ice Extent","month of min Sea Ice Extent","month of max Sea Ice Extent","first Ice-Free month (< 1 mill sq km)"/)
  else
    stats_title = (/"Minimum Ice Volume","Maximum Ice Volume","month of min vhi","month of max vhi","IGNORE"/)
  end if

  min_emerg_crit = stats(:,0,0,0)-2.*sqrt(stats(:,0,0,1))
  max_emerg_crit = stats(:,1,0,0)-2.*sqrt(stats(:,1,0,1))
  imin = new(nens,integer)
  imax = new(nens,integer)
  i_ice_free = new(nens,integer)
  print("model           mn MIN  (std)  ToE min SIE              mn MAX  (std)  ToE max SIE        first mean ice free min")
  do i=0,nens-1
     tmpmin = min(ind(stats(i,0,:,0).lt.min_emerg_crit(i)))
     tmpmax = min(ind(stats(i,1,:,0).lt.max_emerg_crit(i)))
     tmpfree = min(ind(stats(i,0,:,0).lt.1.0))
     if (.not.ismissing(tmpmin)) then
        imin(i) = tmpmin(0)
        yrmin = stats&year(imin(i))
     else
        yrmin = -999
     end if
     if (.not.ismissing(tmpmax)) then
        imax(i) = tmpmax(0)
        yrmax = stats&year(imax(i))
     else
        yrmax = -999
     end if
     if (.not.ismissing(tmpfree)) then
        i_ice_free(i) = tmpfree(0)
        yrfree = stats&year(i_ice_free(i))
     else
        yrfree = -999
     end if
     mn1 = stats(i,0,0,0)
     mn2 = stats(i,1,0,0)
     std1 = sqrt(stats(i,0,0,1))
     std2 = sqrt(stats(i,1,0,1))
     print(mname(i)+"           "+mn1+" ("+std1+")     "+yrmin+"          "+mn2+" ("+std2+")     "+yrmax+"             "+yrfree)
     delete(tmpmin)
     delete(tmpmax)
     delete(tmpfree)
  end do

  if (remove_mean) then
     mnmn2d = conform(stats(:,0,:,0),stats(:,0,0,0),0)
     mnmx2d = conform(stats(:,1,:,0),stats(:,1,0,0),0)

     stats(:,0,:,0) = stats(:,0,:,0) - mnmn2d
     stats(:,1,:,0) = stats(:,1,:,0) - mnmx2d
  end if

;---------------------------------------------------------------------------------
; read in min/max data from example for histogram
;---------------------------------------------------------------------------------
  i1 = base_yr1 - ex_case_yr1
  i2 = base_yr2 - ex_case_yr1
  f1 = addfile(fex_min_max,"r")
  max_ext_base = f1->max_ext(i1:i2,:)
  min_ext_base = f1->min_ext(i1:i2,:)

  max_ext_LE_ex = f1->max_ext
  min_ext_LE_ex = f1->min_ext
;---------------------------------------------------------------------------------
; calculate year and month for sic contour plots
;---------------------------------------------------------------------------------
  y1 = LE_min_yr1 - ex_case_yr1
  y2 = LE_min_yr2 - ex_case_yr1
  min_min = min(min_ext_LE_ex(y1:y2,:))
  max_min = max(min_ext_LE_ex(y1:y2,:))
  do i=0,ex_ncases-1
     case_number = i+1
     case_min_min = min(min_ext_LE_ex(y1:y2,i))
     case_max_min = max(min_ext_LE_ex(y1:y2,i))
     if (case_min_min.eq.min_min) then
        year_min_min = LE_min_yr1+ind(min_ext_LE_ex(y1:y2,i).eq.case_min_min)
        print("Case: "+case_number+" year: "+year_min_min+" min: "+min_min)
     end if
     if (case_max_min.eq.max_min) then
        year_max_min = LE_min_yr1+ind(min_ext_LE_ex(y1:y2,i).eq.case_max_min)
        print("Case: "+case_number+" year: "+year_max_min+" max: "+max_min)
     end if
  end do
;---------------------------------------------------------------------------------
; put data from both PI and LE into single array
; for plotting purposes
;---------------------------------------------------------------------------------
  nyrs_LE = LE_plot_yr2(0) - LE_plot_yr1(0) + 1
  nLE = nyrs_LE*ex_ncases

  ndata = nLE
  min_data = new((/nplot+1,ndata/),typeof(min_ext_base))
  min_data(0,:) = ndtooned(min_ext_base)

  max_data = new((/nplot+1,ndata/),typeof(max_ext_base))
  max_data(0,:) = ndtooned(max_ext_base)
  do i=0,nplot-1
     y1 = LE_plot_yr1(i) - ex_case_yr1
     y2 = LE_plot_yr2(i) - ex_case_yr1
     LE_1d = ndtooned(min_ext_LE_ex(y1:y2,:))
     min_data(i+1,:nLE-1) = (/ LE_1d /)
     LE_1d = ndtooned(max_ext_LE_ex(y1:y2,:))
     max_data(i+1,:nLE-1) = (/ LE_1d /)
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
; get sic for contouring 
;---------------------------------------------------------------------------------
  i1 = (ex_case2_plot_yr - ex_sic_case_yr1)*12 + case2_plot_month - 1
  f = addfile(fsic2,"r")
  tmp = f->$sic_name$
  sic2 = tmp(i1,:,:)
  delete(tmp)
;  sic2 = f->sic(i1,:,:)			; for some reason getting a "subscript out of range" error with this
  if (ex_mname.ne."CESM") then
     tlat = f->lat
     tlon = f->lon
  else
     tlat = f->TLAT
     tlon = f->TLON
  end if

  i1 = (ex_case3_plot_yr - ex_sic_case_yr1)*12 + case3_plot_month - 1
  f = addfile(fsic3,"r")
  tmp = f->$sic_name$
  sic3 = tmp(i1,:,:)
  delete(tmp)
;  sic3 = f->sic(i1,:,:)
  if ((ex_mname.ne."CanESM2").and.(ex_mname.ne."GFDL-ESM2M")) then
     sic2@lat2d = tlat
     sic2@lon2d = tlon

     sic3@lat2d = tlat
     sic3@lon2d = tlon
  end if
  if (max(sic2).lt.1.1) then                    ; sic is fraction rather than %
     sic2 = 100.*sic2
     sic3 = 100.*sic3
  end if
;---------------------------------------------------------------------------------
; plot
;---------------------------------------------------------------------------------
  wks  = gsn_open_wks("pdf",psout)
;  wks  = gsn_open_wks("eps",psout)
;---------------------------------------------------------------------------------
; timeseries for all CMIP5 MMLE
;---------------------------------------------------------------------------------
  gsn_define_colormap(wks,"posneg_1")						;NOTE: figure out colors
  colors = (/5,8,17,15,10/)

  plot_ts = new(2,graphic)

  res1 = True                           ; res for stats plot
  res1@gsnDraw = False
  res1@gsnFrame = False

  res1@vpWidthF    = 0.4      ; Make plots wider than
  res1@vpHeightF   = 0.19      ; they are high.

  res1@tiXAxisString = "year"
  res1@trXMinF = stats_yr1
  res1@trXMaxF = stats_yr2

  res1@xyMonoDashPattern = True
  res1@xyLineColors = colors

  res2 = True
  res2 = res1
  res2@xyDashPattern = 1

  res1@xyLineThicknessF = 5.

  res1@tiYAxisString    = "Mean (month)"
  res2@tiYAxisString    = "STD (month)"

  res1@tiYAxisString    = "Mean (10~S~6~N~ km~S~2~N~)"
  res2@tiYAxisString    = "STD (10~S~6~N~ km~S~2~N~)"

  res1@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res1@gsnCenterStringFontHeightF   = 0.02

  res1@xyLineOpacityF = 0.3
  do i=0,1
     res1@gsnCenterString    =  stats_title(i)
     plot_ts(i)  = gsn_csm_xy(wks,stats&year,stats(:,i,:,0),res1)  ; Create filled XY plot.
  end do
  res1@gsnCenterString    = ""
  res1@xyLineOpacityF= 1.

  sie_emerge = new((/nens,2,nyr/),float)
  sie_emerge!2 = "year"
  sie_emerge&year = ispan(stats_yr1,stats_yr2,1)
  do i=0,nens-1
     if (.not.ismissing(imin(i))) then
        sie_emerge(i,0,imin(i):) = stats(i,0,imin(i):,0)
     end if
     if (.not.ismissing(imax(i))) then
        sie_emerge(i,1,imax(i):) = stats(i,1,imax(i):,0)
     end if
  end do
  plot_min_emerge  = gsn_csm_xy(wks,sie_emerge&year,sie_emerge(:,0,:),res1)  ; Create filled XY plot.
  overlay(plot_ts(0),plot_min_emerge)
  plot_max_emerge  = gsn_csm_xy(wks,sie_emerge&year,sie_emerge(:,1,:),res1)  ; Create filled XY plot.
  overlay(plot_ts(1),plot_max_emerge)

  if (plot_obs) then
     min_SSMI = obs_min_ext
     max_SSMI = obs_max_ext

     if (remove_mean) then
        min_SSMI = min_SSMI - dim_avg(min_SSMI(:9))
        max_SSMI = max_SSMI - dim_avg(max_SSMI(:9))
     end if

     obs_data = new(nyr,float)
     obs_data!0 = "year"
     obs_data&year = ispan(stats_yr1,stats_yr2,1)
     i1 = obs_yr1 - stats_yr1
     i2 = obs_yr2 - stats_yr1

     delete(res1@xyLineColors)
     res1@xyLineColors = "black"

     res1@tiYAxisString    = " "
     res1@gsnLeftString    =  " "
     res1@xyDashPattern = 1

     obs_data(i1:i2) = (/ min_SSMI /)
     plot_obs_min = gsn_csm_xy(wks,obs_data&year,obs_data,res1)
     overlay(plot_ts(0),plot_obs_min)

     obs_data(i1:i2) = (/ max_SSMI /)
     plot_obs_max = gsn_csm_xy(wks,obs_data&year,obs_data,res1)
     overlay(plot_ts(1),plot_obs_max)

  end if
  lgres                    = True
  if (plot_obs) then
     ncolors = dimsizes(colors)
     colors2 = new(ncolors+1,typeof(colors))
     colors2(1:) = colors
     colors2(0) = 1
     lgres@lgLineColors       = colors2
     lgres@lgMonoDashIndex    = True
     lgres@lgDashIndex = 0
     lgres@lgDashIndexes = (/1,0,0/)
     labels = new((nens+1),typeof(mtitle))
     labels(0) = "Sea Ice Index"
     labels(1:) = mtitle
  else
     lgres@lgMonoDashIndex    = True
     lgres@lgDashIndex        = 0
     lgres@lgLineColors       = colors
     labels = mtitle
  end if
  lgres@lgLineThicknessF = 3.
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .07            ; legend label font thickness
  lgres@vpWidthF           = 0.17           ; width of legend (NDC)
  lgres@vpHeightF          = 0.11           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter

  nitems = dimsizes(labels)                 ; number of legend items

;---------------------------------------------------------------------------------
; panel these figures first to preserve colormap
;---------------------------------------------------------------------------------

  pres = True
  pres@gsnFrame = False
  pdiff = 0.18
  pres@gsnPanelTop = 0.98
;  pres@gsnPanelBottom   = 0.8
  pres@gsnPanelBottom   = pres@gsnPanelTop - pdiff

  gsn_panel(wks,plot_ts,(/1,2/),pres)                                           ; MMLE timeseries
  gsn_legend_ndc(wks,nitems,labels,0.8,0.82,lgres)

;---------------------------------------------------------------------------------
; histograms and example sic maps 
;---------------------------------------------------------------------------------

  nplots_per_hist = 2
  plot = new((/2*nplot,nplots_per_hist/),graphic)

  ext_min_levels = fspan(0.,9.5,39)
  if (ex_mname.eq."GFDL-ESM2M") then
     ext_max_levels = fspan(13.0,21.0,33)                                     ; was 9.5,17,5,33 for LENS
  else
     ext_max_levels = fspan(9.0,17.0,33)					; was 9.5,17,5,33 for LENS
  end if
  dext_min = ext_min_levels(1) - ext_min_levels(0)                  ; diff in levels same for min, max
  dext_max = ext_max_levels(1) - ext_max_levels(0)                  ; diff in levels same for min, max
  plot_colors = (/"black","blue","grey"/)
  hres = True
  hres@gsnDraw            = False
  hres@gsnFrame           = False

  xyres = True
  xyres@gsnDraw = False
  xyres@gsnFrame = False
  xyres@tiYAxisString = "Frequency (%)"
  xyres@tiXAxisString = "million sq km"

  xyres@tmXTOn = False

  xyres@gsnLeftStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  xyres@gsnStringFontHeightF   = 0.035

  xyres@gsnXYBarChart         = True                 ; turn on bar chart
  xyres@gsnXYBarChartColors = "blue"
  xyres@gsnXYBarChartBarWidth = 0.4*(ext_min_levels(1)-ext_min_levels(0))  ; 40% of the width

  xyres@tmXBMode = "Explicit"
;  xyres@trYMaxF = 0.5
  xyres@trYMaxF = 0.35
  xyres@trYMinF = 0.
  xyres@tmYLMode = "Explicit"
  per_levels = fspan(0.,1.,21)
  xyres@tmYLValues = per_levels
  xyres@tmYLLabels = ""+(100.*per_levels)

  nlevels = dimsizes(ext_min_levels)
  hres@gsnHistogramBinIntervals = ext_min_levels

  x = (ext_min_levels(1:)+ext_min_levels(0:nlevels-2))*0.5         ; halfway to mid point of each range
  x2 = x - dext_min/8.
  x3 = x + dext_min/8.

  xyres@trXMaxF = ext_min_levels(nlevels-2)
  xyres@trXMinF = ext_min_levels(0)
  xyres@tmXBValues = ext_min_levels(0::2)
  xyres@tmXBLabels = ""+ext_min_levels(0::2)

  xyres@tfPolyDrawOrder    = "Predraw"               ; put line on top

  xyres@tiMainString = " "

  if (plot_obs) then
     gsres                   = True                        ; poly res
     gsres@gsFillColor       = "Grey"                 ; color chosen
     gsres@gsFillOpacityF = 0.25
     dum = new(8*nplot,graphic)
     dumtxt = new(4*nplot,graphic)

     lres = True
     lres@gsLineDashPattern = 1
     lres@gsLineColor = "black"
     lres@gsLineThicknessF = 1.

     txres               = True                     ; text mods desired
     txres@txFontHeightF = 0.015                     ; font smaller. default big

     yp= (/0.,0.,.35,.35/)
  end if

  xyres@gsnLeftString = "MIN NH SIE"
  data = min_data(:1,:)
  do i=0,nplot-1					; min SIE
     data(1,:) = min_data(i+1,:)
     hist_plot = gsn_histogram(wks,data,hres)
     tmpy = hist_plot@NumInBins
     ybinned = int2flt(tmpy)
     ybinned(0,:) = ybinned(0,:)/dim_sum(ybinned(0,:))
     ybinned(1,:) = ybinned(1,:)/dim_sum(ybinned(1,:))

     xyres@gsnRightString = LE_plot_yr1(i)+"-"+LE_plot_yr2(i)
     xyres@gsnXYBarChartColors = "blue"
     xyres@xyLineColors = "blue"
     plot(i,0) = gsn_csm_xy(wks,x2,ybinned(0,:),xyres)
     xyres@gsnXYBarChartColors = "red"
     xyres@xyLineColors = "red"
     plot(i,1) = gsn_csm_xy(wks,x3,ybinned(1,:),xyres)
     overlay(plot(i,0),plot(i,1))

     if (plot_obs) then
        lres@gsLineDashPattern = 1
        lres@gsLineColor = "black"
        dum(4*i) = gsn_add_polyline(wks,plot(i,0),(/obs_min_stats(0,0),obs_min_stats(0,0)/),(/0,.35/),lres)
        lres@gsLineDashPattern = 0
        dum(4*i+1) = gsn_add_polyline(wks,plot(i,0),(/obs_min_stats(1,0),obs_min_stats(1,0)/),(/0,.35/),lres)

        dumtxt(2*i) = gsn_add_text(wks,plot(i,0),"Sea Ice Index ("+year_obs(0)+"-"+year_obs(9)+")",obs_min_stats(0,0),0.33,txres)
        dumtxt(2*i+1) = gsn_add_text(wks,plot(i,0),"Sea Ice Indes ("+year_obs(nyr_obs-10)+"-"+year_obs(nyr_obs-1)+")",obs_min_stats(1,0),0.30,txres)

        obs_min = min(obs_min_ext(:9))
        obs_max = max(obs_min_ext(:9))
        xp_min = (/obs_min,obs_max,obs_max,obs_min/)
        dum(4*i+2) = gsn_add_polygon (wks,plot(i,0),xp_min,yp,gsres)

        obs_min = min(obs_min_ext(nyr_obs-10:))
        obs_max = max(obs_min_ext(nyr_obs-10:))
        xp_min = (/obs_min,obs_max,obs_max,obs_min/)
        dum(4*i+3) = gsn_add_polygon (wks,plot(i,0),xp_min,yp,gsres)
     end if
     xyres@gsnLeftString = ""
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

  xyres@gsnLeftString = "MAX NH SIE"
  data = max_data(:1,:)
  do i=0,nplot-1                                        ; max SIE
     data(1,:) = max_data(i+1,:)
     hist_plot = gsn_histogram(wks,data,hres)
     tmpy = hist_plot@NumInBins
     ybinned = int2flt(tmpy)
     ysum0 = dim_sum(ybinned(0,:))
     ysum1 = dim_sum(ybinned(1,:))
     xyres@gsnRightString = LE_plot_yr1(i)+"-"+LE_plot_yr2(i)

     if ((ysum0.gt.0.).and.(ysum1.gt.0.)) then
        ybinned(0,:) = ybinned(0,:)/ysum0
        ybinned(1,:) = ybinned(1,:)/ysum1

        xyres@gsnXYBarChartColors = "blue"
        xyres@xyLineColors = "blue"
        plot(nplot+i,0) = gsn_csm_xy(wks,x2,ybinned(0,:),xyres)
        xyres@gsnXYBarChartColors = "red"
        xyres@xyLineColors = "red"
        plot(nplot+i,1) = gsn_csm_xy(wks,x3,ybinned(1,:),xyres)
        overlay(plot(nplot+i,0),plot(nplot+i,1))
     end if
     if (plot_obs) then
        lres@gsLineDashPattern = 1
        lres@gsLineColor = "black"
        dum(4*(i+nplot)) = gsn_add_polyline(wks,plot(nplot+i,0),(/obs_max_stats(0,0),obs_max_stats(0,0)/),(/0,.35/),lres)
        lres@gsLineDashPattern = 0
        dum(4*(i+nplot)+1) = gsn_add_polyline(wks,plot(nplot+i,0),(/obs_max_stats(1,0),obs_max_stats(1,0)/),(/0,.35/),lres)

        dumtxt(2*(i+nplot)) = gsn_add_text(wks,plot(nplot+i,0),"Sea Ice Index ("+year_obs(0)+"-"+year_obs(9)+")",obs_max_stats(0,0),0.33,txres)
        dumtxt(2*(i+nplot)+1) = gsn_add_text(wks,plot(nplot+i,0),"Sea Ice Index ("+year_obs(nyr_obs-10)+"-"+year_obs(nyr_obs-1)+")",obs_max_stats(1,0),0.30,txres)

        obs_min = min(obs_max_ext(:9))
        obs_max = max(obs_max_ext(:9))
        xp_min = (/obs_min,obs_max,obs_max,obs_min/)
        dum(4*(i+nplot)+2) = gsn_add_polygon (wks,plot(nplot+i,0),xp_min,yp,gsres)

        obs_min = min(obs_max_ext(nyr_obs-10:))
        obs_max = max(obs_max_ext(nyr_obs-10:))
        xp_min = (/obs_min,obs_max,obs_max,obs_min/)
        dum(4*(i+nplot)+3) = gsn_add_polygon (wks,plot(nplot+i,0),xp_min,yp,gsres)
     end if
     xyres@gsnLeftString = ""
  end do
;-------------------------
; resources for contouring
;------------------------
  sicplot = new(2,graphic)
  sic15per = new(2,graphic)
  cmap = RGBtoCmap("/glade/u/home/landrum/mydiags/ice_diag_swift/rgb/julies.rgb")
;  gsn_define_colormap(wks,cmap)
  gsn_define_colormap(wks,"cmocean_ice")
  res                      = True                 ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnFillOn             = True                 ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels" ; set manual contour levels
;  res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)				; cmap
  res@cnFillColors = (/27,43,59,75,91,107,123,139,155,171,187,203,219,235,251/)				; cmocean_ice

  res@lbLabelAutoStride    = True                 ; control labelbar labels
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn the line labels off

  res@lbLabelBarOn        = False       ; Turn off labelbar
  res@lbTitleOn        = True                  ; turn on title
  res@lbTitleString    = "m"
  res@lbTitleFontHeightF= .015                 ; make title smaller

  res@cnInfoLabelOn       = False       ; Turn off informational label
  res@trGridType = "TriangularMesh"
  res@gsnStringFontHeightF = 0.015
  res@cnMonoLineThickness = True
  res@cnMonoLineDashPattern = True

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnLeftStringFontHeightF   = 0.035
  res@gsnLeftStringOrthogonalPosF = 0.03

  res@gsnAddCyclic = True

  res@gsnCenterString = ""
  res@gsnRightString = ""

  res@gsnPolar  = "NH"
  res@mpLandFillColor      = "grey"           ; choose color of continents.
  res@mpFillDrawOrder      = "PostDraw"       ; Fill continents later to cover up lines
  res@mpPerimDrawOrder     = "PostDraw"       ; Perimeter gets covered up by land
  res@mpMinLatF = 55.
  res@gsnPolarLabelSpacing = 30.
  ncntrs = 14

  sic_cntrs=(/1,5,10,15,20,30,40,50, 60, 70, 80,90,95, 99/) ; sic
  res@cnLevels = sic_cntrs
  res@gsnLeftString    =  ex_case2_plot_yr
  sicplot(0) = gsn_csm_contour_map_polar(wks,sic2,res)

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
  res15per@cnLineColor = "white"
  sic15per(0) = gsn_csm_contour(wks,sic2,res15per)
  overlay(sicplot(0),sic15per(0))

  res@gsnLeftString    =  ex_case3_plot_yr

  sicplot(1) = gsn_csm_contour_map_polar(wks,sic3,res)
  sic15per(1) = gsn_csm_contour(wks,sic3,res15per)
  overlay(sicplot(1),sic15per(1))

;------------------------------------------------------------------------------------------------------------
; now panel remaining plots
;------------------------------------------------------------------------------------------------------------

  pres@gsnPanelTop = pres@gsnPanelBottom - 0.05
  pres@gsnPanelBottom   = pres@gsnPanelTop - 2*pdiff

  pres@gsnPanelMainString = ex_title
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.015

  gsn_panel(wks,plot(:,0),(/2,nplot/),pres)                              ; first set of histograms (min SIE)

  pres@gsnPanelTop = pres@gsnPanelBottom - 0.07
  pres@gsnPanelBottom   = 0.05
  pres@gsnPanelLabelBar  = True
  pres@lbTitleString    = "%"
  pres@lbTitlePosition  = "Right"
  pres@lbTitleFontHeightF= .015
  pres@gsnPanelMainString = " Sea Ice Concentrations during~C~ ~Z85~ (2010-2019) min SIE extremes"
  gsn_panel(wks,sicplot,(/1,2/),pres)						; example maps
 
  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.13           ; labelbar width
  lbres@vpHeightF          = 0.07           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.014         ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbPerimOn          = False

  label_colors = (/"blue","red"/)
  ex_labels = (/"("+base_yr1+"-"+base_yr2+")", "future decades"/)
  xpos = (/0.78, 0.78/)
  ypos = (/0.49,0.47/)
  do i=0,1
     lbres@lbFillColors = label_colors(i)
     gsn_labelbar_ndc(wks,1,ex_labels(i),xpos(i),ypos(i),lbres)
  end do	
  frame(wks)
end

