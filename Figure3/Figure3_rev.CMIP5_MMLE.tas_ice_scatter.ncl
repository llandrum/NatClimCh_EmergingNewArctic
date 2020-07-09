; Figure 3: CMIP5-MMLE scatter plots
; OCT dtas: Sept SIE
; FEB dtas: FEB VHI
; SEPT SIE: FEB VHI
;
; Modified to plot all CMIP5 MMLE together....
;
; Plots scatter plot of area weighted average tas vs
; user-defined area averaged sea ice variable (e.g. sie, vai, vhi, etc). User
; also defines if scatter is vs change in ice variable
; or mean ice variable (e.g. mean hi).
;
; LLandrum
; March, 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
; ================================================;
begin
  imonth = (/9,1/)
  var1_name = "tas"
  var1_plot_type = "diff"		; for this figure use diff

  imonth_ice = (/8,1/)                        ; needs to be same dimension as ice_var_name
  nice_var = dimsizes(imonth_ice)
  ice_varname = (/"ext","vhi"/)           ; ext, vai, vhi and same dims as imonth_ice
  ice_varname_abbrev = (/"SIE","SIT"/)
  plot_fitline = (/True,True,True/)                           ; True to plot regression line, False otherwise
  fit_type = "poly"                             ; poly (polynomial fit) or regline (regression line)
  npoly_coef = (/3,3,3/)                                ; number of coefficients for polynomial fit

  ice_region = "ArcOc"                             ; eg sh, nh, ArcOc, CArc
  ice_region_longname = "Arctic Ocean"
  var2_name = ice_varname+"_mo_"+ice_region
  var2_plot_type = (/"mean","mean"/)
  avg_vhi = True				; True to plot as area average ice thickness, false to plot as area volume
  mask_icefree = True                          ; True to mask out points in sea ice variable when NH SIE is "ice free" (<1x10^12 m2)
                                               ; for calculating polynomial fit (otherwise polynomial fit can be overly weighted at low ice values and 
					       ; give inaccurate and misleading curves for high ice values)

  model_names = (/"MPI-ESM","CanESM2","CESM1-CAM5","GFDL-CM3","GFDL-ESM2M"/)
  dir_names = (/"mpi","canesm2","cesm","gfdl_cm3","gfdl_esm2m"/)
  model_case_yr1 = (/1920,1950,1920,1920,1950/)
  model_case_yr2 = (/2099,2100,2100,2100,2100/)
  model_nsim = (/100,50,40,20,30/)
  nmodel = dimsizes(model_names)
  nsim_max = dim_max(model_nsim)

  base_yr1 = 1950
  base_yr2 = 1959
  nyr_base = base_yr2 - base_yr1 + 1

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/"
  nyr_stats = 10                                                ; # yrs used in statistics

  nh_reg = 70
  
  plot_yr1 = 1980
  plot_yr2 = dim_min(model_case_yr2)
  nyr_scatter = plot_yr2 - plot_yr1 + 1

  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
  max_ext = (/19.,19.,16.,14.,14.,12.,12.,10.,10.,12.,14.,16./)
  if (ice_region.eq."ArcOc") then
     max_ext = 8.
  end if
  min_ext = (/2.,2.,1.5,1.5,1.,1.,1.,0.,0.,0.,1.,1./)
  max_vol = (/4.,4.,4.,3.,3.,3.,2.,2.,1.,1.,2.,3./)
  max_vol = 3.5

  ice_free_min = 1.					; definition of ice-free minimum (in units of 10^6 km^2)

  psout = "Figure3_rev.CMIP5_"+nmodel+"MMLE.tas_ice_scatter."+ice_region+"."+plot_yr1+"-"+plot_yr2
;---------------------------------------------------------------------------------
  var1 = new((/nmodel,2,nyr_scatter,nsim_max/),float)
  var1!0 = "model"
  var1!1 = "region"
  var1!2 = "year"
  var1!3 = "case"
  var1_diff = var1
  var1_base_clim = new((/nmodel,2/),float)

  var2 = new((/nmodel,nice_var,nyr_scatter,nsim_max/),float)
  var2!1 = "ice_var"
  var2!2 = "year"
  var2!3 = "case"
  nh_ext = var2                                ; for use masking out ice free points

  var2_base = new((/nmodel,nice_var,nyr_base,nsim_max/),float)
  var2_base!0 = "ice_var"
  var2_base!1 = "year"
  var2_base!2 = "case"
;---------------------------------------------------------------------------------
; read in LE data
;---------------------------------------------------------------------------------
  do im=0,nmodel-1
     mname = model_names(im)
     dirmname = dir_names(im)
     nsim = model_nsim(im)
     case_yr1 = model_case_yr1(im)
     case_yr2 = model_case_yr2(im)

     case_numbers = ispan(1,nsim,1)
     styrs_proc = sprinti("%0.4i",case_yr1) + "01-"+sprinti("%0.4i",case_yr2)+"12"

     farea = datadir+"ice/"+dirmname+"/"+mname+"_regional_areavg.nc"
     ftas = datadir+"atm/"+dirmname+"/"+mname+"."+nsim+"sim.tas_wgtavg."+styrs_proc+".nc"

     f = addfile(ftas,"r")
 
;------------------------------------------------------------------------------
; now ice data 
;------------------------------------------------------------------------------
     if ((mname.eq."GFDL-CM3").or.(mname.eq."GFDL-ESM2M")) then
        fice = datadir+"ice/"+dirmname+"/ice_vol_MHolland_regions_"+mname+"_historical_rcp85_r"+case_numbers+"i1p1_sitXsic_"+case_yr1+"-"+case_yr2+".nc"
     else
        if (mname.eq."MPI-ESM") then
           scnumbers = sprinti("%0.3i",case_numbers)
           fice1 = datadir+"ice/"+dirmname+"/ice_vol_MHolland_regions_"+mname+"_historical_r"+scnumbers+"i1850p3_"+case_yr1+"-2005.nc"
           fice2 = datadir+"ice/"+dirmname+"/ice_vol_MHolland_regions_"+mname+"_rcp85_r"+scnumbers+"i2005p3_2006-"+case_yr2+".nc"
           delete(scnumbers)
        else
           if (mname.eq."CESM1-CAM5") then
              case_numbers(35:) = ispan(101,105,1)
              scnumbers = sprinti("%0.3i",case_numbers)

              datadirLE = "/glade/p/cesm/pcwg/landrum/LE/ice/pre_process/monthly/"
              fice1 = datadirLE+"ice_vol_MHolland_regions_b.e11.B20TRC5CNBDRD.f09_g16."+scnumbers+"_"+case_yr1+"-2005.nc"
              fice2 = datadirLE+"ice_vol_MHolland_regions_b.e11.BRCP85C5CNBDRD.f09_g16."+scnumbers+"_2006-"+case_yr2+".nc"
              delete(scnumbers)
           else
              fice = datadir+"ice/"+dirmname+"/ice_vol_MHolland_regions_"+mname+"_historical_rcp85_r"+case_numbers+"i1p1_"+case_yr1+"-"+case_yr2+".nc"
           end if
        end if
     end if

     do ii=0,nice_var-1
        i1 = (base_yr1 - case_yr1)*12 + imonth(ii)
        i2 = i1 + nyr_base*12 - 1
        ii1 = (plot_yr1 - case_yr1)*12 + imonth(ii)
        ii2 = (plot_yr2 - case_yr1)*12 + imonth(ii)

        var_name = "tas_mo_nh"+nh_reg
        var1_base_clim(im,ii) = dim_avg(dim_avg(f->$var_name$(i1:i2:12,:)))
        var1(im,ii,:,:nsim-1) = (/ f->$var_name$(ii1:ii2:12,:) /)    ; base climatological ensemble global  mean
        var1_diff(im,ii,:,:nsim-1) = var1(im,ii,:,:nsim-1) - var1_base_clim(im,ii)

        i1 = (base_yr1 - case_yr1)*12 + imonth_ice(ii)
        i2 = i1 + nyr_base*12 - 1
        ii1 = (plot_yr1 - case_yr1)*12 + imonth_ice(ii)
        ii2 = (plot_yr2 - case_yr1)*12 + imonth_ice(ii)

        vname = var2_name(ii)
        do is=0,nsim-1
           if ((mname.eq."MPI-ESM").or.(mname.eq."CESM1-CAM5")) then
              iy2 = 2005 - plot_yr1
              f1 = addfile(fice1(is),"r")
              var2(im,ii,:iy2,is) = (/ f1->$vname$(ii1::12) /)
              nh_ext(im,ii,:iy2,is) = (/ f1->ext_mo_nh(ii1::12) /)

              f2 = addfile(fice2(is),"r")
              ii2 = (plot_yr2 - 2006)*12 + imonth_ice(ii)
              var2(im,ii,iy2+1:,is) = (/ f2->$vname$(imonth_ice(ii):ii2:12) /)
              var2_base(im,ii,:,is) = (/ f1->$vname$(i1:i2:12) /)
              nh_ext(im,ii,iy2+1:,is) = (/ f2->ext_mo_nh(imonth_ice(ii):ii2:12) /)
           else
              f1 = addfile(fice(is),"r")
              var2(im,ii,:,is) = (/ f1->$vname$(ii1:ii2:12) /)
              var2_base(im,ii,:,is) = (/ f1->$vname$(i1:i2:12) /)
              nh_ext(im,ii,:,is) = (/ f1->ext_mo_nh(ii1:ii2:12) /)
           end if
        end do                                          ; end nsim

        if ((ice_varname(ii).eq."vhi").and.(avg_vhi)) then
           f3 = addfile(farea,"r")
           avg_name = "area_"+ice_region
           reg_areaavg = f3->$avg_name$
           var2_base(im,ii,:,:) = var2_base(im,ii,:,:)/reg_areaavg
           var2(im,ii,:,:) = var2(im,ii,:,:)/reg_areaavg
        end if
        if (ice_varname(ii).eq."ext") then			; mask out low areas
           var2_base(im,ii,:,:) = where(var2_base(im,ii,:,:).gt.0.05,var2_base(im,ii,:,:),var2_base@_FillValue)
           var2(im,ii,:,:) = where(var2(im,ii,:,:).gt.0.05,var2(im,ii,:,:),var2_base@_FillValue)
        end if
     end do						; end ii

     delete(case_numbers)
     delete(styrs_proc)
     if ((mname.eq."MPI-ESM").or.(mname.eq."CESM1-CAM5")) then
        delete(fice1)
        delete(fice2)
     else
        delete(fice)
     end if
  end do						; end im

  var1 = var1 - 273.15
  var1@units = "deg C" 
  var1_base_clim = var1_base_clim - 273.15
  var1_base_clim@units = "deg C"
  var1_diff@units = "deg C"

  var2_base_clim = dim_avg(dim_avg(var2_base))		; averaged over ensemble then years....var2_base_clim has 2 base values

  var2_base_clim4d = conform(var2,var2_base_clim,(/0,1/))
  var2_diff = var2
  var2_diff = var2 - var2_base_clim4d

  if (mask_icefree) then
     var2_icefree = where(nh_ext.le.1.,var2@_FillValue,var2)
     var2_diff_icefree = var2_icefree
     var2_diff_icefree = var2_icefree - var2_base_clim4d
  end if
;---------------------------------------------------------------------------------
; plot
;---------------------------------------------------------------------------------
;  wks  = gsn_open_wks("pdf",psout)
  wks  = gsn_open_wks("eps",psout)

  plot = new((/3,nmodel/),graphic)
  plot_fit = new((/3,nmodel/),graphic)

  gsn_define_colormap(wks,"posneg_1")
  colors = (/5,8,17,15,10/)
  opacities = (/0.07,0.1,0.05,0.1,0.1/)                                      ; red needs lower opacity, and also MPI w 100 simulations!

  hres = True
  hres@gsnDraw            = False
  hres@gsnFrame           = False

  res = True				; res for stats plot
  res@gsnDraw = False
  res@gsnFrame = False

;  res@gsnMaximize         = True                   ; maximize plot in frame
  res@xyMarkers           = 16                     ; choose type of marker 
  res@xyMarkerSizeF       = 0.005                  ; Marker size (default 0.01)

  res@xyMarkerOpacityF = 0.1
  res@xyLineOpacityF = 1.
  res@xyLineThicknessF = 4.

  res@xyMarkLineModes     = "Markers"

  res@tiYAxisString = "~S~o~N~C"

  res@trYMinF = -5.
  res@trYMaxF = 25.

  res@gsnStringFont = "helvetica-bold"                  ; set gsnLeft/Center/RightString font
  res@gsnStringFontHeightF   = 0.027

;---------------------------------------------------------------------------------
; scatter plots for tas vs ice
;---------------------------------------------------------------------------------

  do ii=0,nice_var-1

     if ((ice_varname(ii).eq."ext").or.(ice_varname(ii).eq."vai")) then                                   ; area
        res@tiXAxisString    = "10~S~6~N~ km~S~2~N~"
        if (var2_plot_type(ii).eq."mean") then
           res@trXMinF = 0.
           res@trXMaxF = max_ext(imonth_ice(ii))
        end if
     else    
        if (avg_vhi) then
           res@tiXAxisString    = "m"
        else
           res@tiXAxisString    = "10~S~13~N~ km~S~3~N~"
        end if
        if (var2_plot_type(ii).eq."diff") then
           res@trXMinF = -10.
           res@trXMaxF = 0.
        else
           res@trXMinF = 0.
           res@trXMaxF = max_vol(imonth_ice(ii))
        end if
     end if

     if (var2_plot_type(ii).eq."diff") then
        text2 = "~F8~D~F21~"+ice_varname_abbrev(ii)
     else
        text2 = ice_varname_abbrev(ii)
     end if

     if (imonth(ii).eq.imonth_ice(ii)) then
        res@gsnCenterString = month_names(imonth(ii))+" ~F8~D~F22~TAS:"+text2
     else
        res@gsnCenterString = month_names(imonth(ii))+" ~F8~D~F22~TAS:"+month_names(imonth_ice(ii))+" "+text2
     end if

     do i=0,nmodel-1
        res@xyMarkerOpacityF = opacities(i)
        nsim = model_nsim(i)
        if (var2_plot_type(ii).eq."diff") then
           var2_plot = var2_diff(i,ii,:,:nsim-1)
        else
           var2_plot = var2(i,ii,:,:nsim-1)
        end if

        x = ndtooned(var2_plot)

        y = ndtooned(var1_diff(i,ii,:,:nsim-1))
        res@xyMarkerColor       = colors(i)

        plot(ii,i)  = gsn_csm_xy(wks,x,y,res)

        if (plot_fitline(ii)) then
           if (mask_icefree) then
              if (var2_plot_type(ii).eq."diff") then
                 var2_icefree_plot = var2_icefree_diff(i,ii,:,:nsim-1)
              else
                 var2_icefree_plot = var2_icefree(i,ii,:,:nsim-1)
              end if
              x = ndtooned(var2_icefree_plot)
              delete(var2_icefree_plot)
           end if
           res@xyMarkLineModes     = "Lines"

           res@xyLineColor       = colors(i)
           x1 = x
           jj = dim_pqsort_n(x1,1,0)
           xx = x1(jj)
           yy = y(jj)

           if (num(.not.ismissing(x1)).gt.100) then			; requiring a minimum of 100 points 
              if (fit_type.eq."regline") then
                 rc = regline_stats(xx,yy)
                 yfit = rc*(xx-rc@xave) + rc@yave
; y = mx+b 
; m is the slope:       rc      returned from regline
; b is the y intercept: rc@yave attribute of rc returned from regline
;                 diff_b = rc@yintercept - rc@y95(0)
;                 diff_m = rc - rc@b95(0)

;        print("+/- "+diff_b+" (intercept) and +/- "+diff_m+" (slope)")

;                 print("("+plot_yr1+"-"+plot_yr2+") (npnts: "+rc@N+") slope: "+RoundPastDecimal(rc,3)+" ("+RoundPastDecimal(diff_m,3)+") intercept: "+RoundPastDecimal(rc@yintercept,3)+" ("+RoundPastDecimal(diff_b,3)+")")
              else							; polynomial fit
                 c = lspoly(x1,y,1,npoly_coef(ii))
 
                 yfit = c(0) + c(1)*xx
                 if (npoly_coef(ii).gt.2) then
                    do icoef=2,npoly_coef(ii)-1
                       yfit = yfit + c(icoef)*(xx^icoef)
                    end do
                 end if 
                 delete(c)
              end if
              plot_fit(ii,i)  = gsn_csm_xy(wks,xx,yfit,res)

;              overlay(plot(ii,i),plot_fit(ii,i))

              res@xyMarkLineModes     = "Markers"

              delete(yfit)
           end if
           delete(x1)
           delete(jj)
           delete(xx)
           delete(yy)
        end if
        if (i.gt.0) then
           overlay(plot(ii,0),plot(ii,i))
        end if
        delete(x)
        delete(var2_plot)
        delete(y)
     end do							; end nmodel
  end do 							; end nice_var
  
;---------------------------------------------------------------------------------
; scatter plots for ice vs ice: ice_varname(0) on x, ice_varname(1) on y
;---------------------------------------------------------------------------------
  if ((ice_varname(0).eq."ext").or.(ice_varname(0).eq."vai")) then                                   ; area
     res@tiXAxisString    = "10~S~6~N~ km~S~2~N~"
     res@trXMinF = 0.
     res@trXMaxF = max_ext(imonth_ice(0))
  else
     if (avg_vhi) then
        res@tiXAxisString    = "m"
     else
        res@tiXAxisString    = "10~S~13~N~ km~S~3~N~"
     end if
     if (var2_plot_type(0).eq."diff") then
        res@trXMinF = -10.
        res@trXMaxF = 0.
     else
        res@trXMinF = 0.
        res@trXMaxF = max_vol(imonth_ice(0))
     end if
  end if

  if ((ice_varname(1).eq."ext").or.(ice_varname(1).eq."vai")) then                                   ; area
     res@tiYAxisString    = "10~S~6~N~ km~S~2~N~"
     res@trYMinF = 0.
     res@trYMaxF = max_ext(imonth_ice(0))
  else
     if (avg_vhi) then
        res@tiYAxisString    = "m"
     else
        res@tiYAxisString    = "10~S~13~N~ km~S~3~N~"
     end if
     if (var2_plot_type(0).eq."diff") then
        res@trYMinF = -10.
        res@trYMaxF = 0.
     else
        res@trYMinF = 0.
        res@trYMaxF = max_vol(imonth_ice(0))
     end if
  end if

  if (var2_plot_type(0).eq."diff") then
     text1 = "~F8~D~F21~"+ice_varname_abbrev(0)
  else
     text1 = ice_varname_abbrev(0)
  end if

  if (imonth_ice(0).eq.imonth_ice(1)) then
     res@gsnCenterString = month_names(imonth_ice(0))+" "+text1+":"+text2
  else
     res@gsnCenterString = month_names(imonth_ice(0))+" "+text1+":"+month_names(imonth_ice(1))+" "+text2
  end if

  do i=0,nmodel-1
     res@xyMarkerOpacityF = opacities(i)
     nsim = model_nsim(i)
     if (var2_plot_type(1).eq."diff") then
        y = ndtooned(var2_diff(i,1,:,:nsim-1))
     else
        y = ndtooned(var2(i,1,:,:nsim-1))
     end if

     if (var2_plot_type(1).eq."diff") then
        x = ndtooned(var2_diff(i,0,:,:nsim-1))
     else
        x = ndtooned(var2(i,0,:,:nsim-1))
     end if
     res@xyMarkerColor       = colors(i)

     plot(2,i)  = gsn_csm_xy(wks,x,y,res)

     if (plot_fitline(2)) then
        if (mask_icefree) then
           if (var2_plot_type(1).eq."diff") then
              y = ndtooned(var2_icefree_diff(i,1,:,:nsim-1))
           else
              y = ndtooned(var2_icefree(i,1,:,:nsim-1))
           end if
           if (var2_plot_type(0).eq."diff") then
              x = ndtooned(var2_icefree_diff(i,0,:,:nsim-1))
           else
              x = ndtooned(var2_icefree(i,0,:,:nsim-1))
           end if 
        end if
        res@xyMarkLineModes     = "Lines"

        res@xyLineColor       = colors(i)
        x1 = x
        jj = dim_pqsort_n(x1,1,0)
        xx = x1(jj)
        yy = y(jj)

        if (num(.not.ismissing(x1)).gt.100) then                     ; requiring a minimum of 100 points
           if (fit_type.eq."regline") then
              rc = regline_stats(xx,yy)
              yfit = rc*(xx-rc@xave) + rc@yave
           else                                                      ; polynomial fit
              c = lspoly(x1,y,1,npoly_coef(2))

              yfit = c(0) + c(1)*xx
              if (npoly_coef(2).gt.2) then
                 do icoef=2,npoly_coef(2)-1
                    yfit = yfit + c(icoef)*(xx^icoef)
                 end do
              end if
              delete(c)
           end if
           plot_fit(2,i)  = gsn_csm_xy(wks,xx,yfit,res)

;           overlay(plot(ii,i),plot_fit(ii,i))

           res@xyMarkLineModes     = "Markers"

           delete(yfit)
        end if
        delete(x1)
        delete(jj)
        delete(xx)
        delete(yy)
     end if
     if (i.gt.0) then
        overlay(plot(2,0),plot(2,i))
     end if
     delete(x)
     delete(y)
 end do                                                     ; end nmodel

   do ii=0,2							; overlay last for better figure
      if (plot_fitline(ii)) then
         do i=0,nmodel-1
            overlay(plot(ii,i),plot_fit(ii,i))
         end do
      end if
   end do

   resM                      = True            ; resources for map
   resM@mpOutlineOn = True
   resM@mpFillOn = True
   resM@gsnDraw              = False        ; don't draw yet
   resM@gsnFrame             = False        ; don't advance frame yet
   resM@mpLandFillColor      = "grey"
   resM@mpFillDrawOrder = "PreDraw"
   resM@mpGridAndLimbDrawOrder = "PreDraw"
;   resM@gsnSpreadColors     = True
   resM@mpMinLatF            = 60
   resM@gsnPolar             = "NH"               ; choose hemisphere
   resM@vpWidthF           = 0.1
   resM@vpHeightF          = 0.1
   resM@vpXF = 0.55
   resM@vpYF = 0.61

   resP = True
   resP@gsFillOpacityF = 0.65
   resP@gsFillColor = "navy"
   resP@tfPolyDrawOrder = "PostDraw"

   map = gsn_csm_map_polar(wks,resM)
   ArcOc_lat = (/90., 80., 80., 68., 68., 90. /)
   ArcOc_long = (/-117., -117., 100., 100., 243., 243./)
   gsn_polygon(wks,map,ArcOc_long,ArcOc_lat,resP)

   res_lat = True
   res_lat@gsLineThicknessF = 3.
   draw(map)
   gsn_polyline(wks,map,(/0.,180.,360./),(/70.,70.,70./),res_lat)

;----------------------------------------------------------------------
; panel plots
;----------------------------------------------------------------------
 
   pres = True
   pres@gsnFrame = False
   pres@txString    =  "CMIP5 MMLE ("+plot_yr1+"-"+plot_yr2+")"

   gsn_panel(wks,(/plot(:,0)/),(/1,3/),pres)
;----------------------------------------------------------------------
; add a legend
;----------------------------------------------------------------------
  lgres                    = True
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgLineColors       = colors
  labels = model_names
  lgres@lgLineThicknessF = 3.
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.15           ; width of legend (NDC)
  lgres@vpHeightF          = 0.1           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter

  nitems = dimsizes(labels)                 ; number of legend items

  gsn_legend_ndc(wks,nitems,labels,0.14,0.6,lgres)
  frame(wks)

 end

