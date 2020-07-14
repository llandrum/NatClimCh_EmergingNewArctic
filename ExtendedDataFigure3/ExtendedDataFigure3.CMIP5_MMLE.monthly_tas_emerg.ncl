; Modified to print out mean ToE for 70N-90N in tas
; LLandrum
; April 2020
;
; Maps tas ToE for Oct, Feb for CMIP5 MMLE
; LLandrum
; March 2020
;
; Maps year that the mean variable exceeds threshold for OCT, FEB 
; and plots timeseries of mean, std for TREFHT for OCT, FEB at one land
; gridpoint (Fairbanks, AK) and one ocean (East Siberian Shelf) gridpoint.
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
  mname = (/"MPI-ESM","CanESM2","CESM","GFDL-CM3","GFDL-ESM2M"/)
  dirmname = (/"mpi","canesm2","LE","gfdl_cm3","gfdl_esm2m"/)
  ncases = (/100,50,40,20,30/)
  stats_yr1 = (/1920,1950,1920,1950,1920/)
  stats_yr2 = (/2090,2091,2091,2091,2091/)
  nens = dimsizes(mname)

  var_name = "tas"
  base_clim_yr1 = 1950
  base_clim_yr2 = 1959
  base_clim_txt = base_clim_yr1+"-"+base_clim_yr2

  thresh_type = "std"					; if "max" then maximum value seen in base climate
  thresh = 2.						; x thresh_type

  nyr_stats = 10						; # yrs used in statistics
  nh_reg = 70.						        ; minimum Latitude for calculatin area averaged ToE

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"
  iCESM = ind(mname.eq."CESM")
  mtitle = mname
  mtitle(iCESM) = "CESM1-CAM5"
  datadir(iCESM) = "/glade/p/cesm/pcwg/landrum/LE/TREFHT/monthly/"

  plot_months = (/9,1/)
  nmonths = dimsizes(plot_months)
  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

  psout = "SupplementaryFigure4.CMIP5_"+nens+"MMLE.monthly_tas_emergence"

  pi=4.*atan(1.0)
  rad=(pi/180.)

;-----------------------------------------------------------------------------------------
; plot
;-----------------------------------------------------------------------------------------
;  wks = gsn_open_wks("pdf",psout)             ; send graphics to PNG file
  wks = gsn_open_wks("eps",psout)             ; send graphics to PNG file
  cmap = RGBtoCmap("/glade/u/home/landrum/mydiags/ice_diag_swift/rgb/julies.rgb")
  gsn_define_colormap(wks,"GMT_panoply")
  gsn_reverse_colormap(wks)

  plot = new(nens*nmonths,graphic)

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
;  res@cnMinLevelValF       = 2006
;  res@cnMaxLevelValF       = 2080

  res@gsnStringFont = "helvetica-bold"			; set gsnLeft/Center/RightString font
  res@gsnRightStringFontHeightF   = 0.04
  res@gsnLeftStringFontHeightF   = 0.045
  res@gsnCenterString = " "
  res@gsnRightString = ""
;------------------------------------------------------------------------------
; read in data - within plotting as grids differ for each model
;------------------------------------------------------------------------------
  print("Mean year of emergence by model for "+month_names(plot_months(0))+" and "+month_names(plot_months(1))+" tas, area averaged "+nh_reg+"-90 N")
  avg_ToE = new(nmonths,float)

  do im=0,nens-1
     if (im.eq.iCESM) then
        femerg = datadir(im)+"emergence."+thresh+"_"+thresh_type+".LE."+ncases(im)+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(plot_months)+"_monthly_TREFHT.nc"
     else
        femerg = datadir(im)+"emergence."+thresh+"_"+thresh_type+"."+mname(im)+"."+ncases(im)+"sim.NH_"+nyr_stats+"yr_stats."+base_clim_yr1+"_base."+month_names(plot_months)+"_monthly_"+var_name+".nc"
     end if
     do i=0,nmonths-1
        f = addfile(femerg(i),"r")
        if (i.eq.0) then
           lat = f->lat
           lon = f->lon
           nlat = dimsizes(lat)
           nlon = dimsizes(lon)

           yr_emergence = new((/nmonths,nlat,nlon/),integer)
           yr_emergence!0 = "month"
           yr_emergence!1 = "lat"
           yr_emergence&lat = lat
           yr_emergence!2 = "lon"
           yr_emergence&lon = lon

           delta_mean = new((/nmonths,nlat,nlon/),float)
           delta_mean!0 = "month"
           delta_mean!1 = "lat"
           delta_mean&lat = lat
           delta_mean!2 = "lon"
           delta_mean&lon = lon
           delta_var = delta_mean
           res@gsnLeftString = mtitle(im)

           coswgt=cos(rad*lat)
           coswgt!0 = "lat"
           coswgt&lat= lat
        else
           res@gsnLeftString = ""
        end if
        if (im.eq.0) then
           res@gsnRightString = month_names(plot_months(i))
        else
           res@gsnRightString = " "
        end if

        yr_emergence(i,:,:) = (/ f->yr_emergence /)
        delta_mean(i,:,:) = (/ f->delta_mean /)
        delta_var(i,:,:) = (/ f->delta_var /)
;        res@gsnLeftString = month_names(plot_months(i))
        plot(i+nmonths*im) = gsn_csm_contour_map_polar(wks,yr_emergence(i,:,:),res)     ; create plot

        avg_ToE(i) = wgt_areaave(yr_emergence(month|i,{lat|70.:90.},lon|:),coswgt({70.:90.}),1.0,0)
     end do
     delete(femerg)
     delete(lat)
     delete(lon)
     delete(yr_emergence)
     delete(delta_mean)
     delete(delta_var)
     delete(coswgt)					; each model different grid

     print(mname(im)+":        "+RoundPastDecimal(avg_ToE(0),3)+"             "+RoundPastDecimal(avg_ToE(1),3))                         
   end do

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar  = True
  pres@gsnMaximize     = True                 ; fill up the page
  pres@gsnPanelMainString = "Time of Emergence"
  pres@gsnPanelMainFont = "helvetica-bold"
  pres@gsnPanelMainFontHeightF = 0.016
  pres@lbTitleOn        =  True               ; turn on title
  pres@lbTitleString    = "yr"                ; title string
  pres@lbTitleFontHeightF= .013                ; make title smaller
  pres@lbTitlePosition  = "Right"              ; title position
  pres@lbTitleDirection = "Across"             ; title direction
  pres@lbLabelFontHeightF = 0.01
;  pres@gsnPanelRowSpec  = True
;

  gsn_panel(wks,plot,(/nens,nmonths/),pres)
  frame(wks)

end
