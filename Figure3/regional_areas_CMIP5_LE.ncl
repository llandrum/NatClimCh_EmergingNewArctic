; Calculates total area by region for different CMIP5 models. Will be used to plot
; regional sea ice thickness as area average thickeness (vhi/tarea)
; Area grid files (areacello, areacella, 1x1) do not consistently indicate land vs 
; ocean grid cells. For consistency then I define area grid cells to be used in area
; averaging as the grid cells that contain greater than 0.00 sea ice concentration during 
; January(NH)/August(SH) of the first year of the sic file. For some regions (Arctic Basin,
; for example), this should be fine. For others this may result in different area calculations
; for different models based on initial sic (e.g. Lab sea). Not sure how else to do it...
;
; LLandrum
; Mar 2020
;***************************************************************
; This script is a copy of pre_process.ncl from the ice diagnostics
; package modified to use already processed monthly aice, hi, hs 
; files.
; Laura Landrum
; December, 2010
;
; modified to create pre processed ice file from single-variable
; files that are saved by hemisphere (e.g. Large Ensemble).
; September, 2014
;***************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  mname = "MPI-ESM"
  model_names = (/"CanESM2","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2M","MPI-ESM","CESM1-CAM5"/)
  dir_names = (/"canesm2","csiro_mk36","gfdl_cm3","gfdl_esm2m","mpi","cesm"/)
  model_case_yr1 = (/1950,1850,1920,1950,1850,1920/)
  model_case_yr2 = (/2100,2100,2100,2100,2099,2100/)
  model_nsim = (/50,29,20,30,100,40/)
  grids = (/"atm","atm","ocn","1x1","ocn","ocn"/)
  grid_dirs = "/glade/p/cesm/pcwg/"+(/"landrum/CMIP5","mholland/cmip5","mholland/cmip5","landrum/CMIP5","mholland/cmip5","mholland/cmip5"/)+"/fx/"

  imodel = ind(model_names.eq.mname)
  dirmname = dir_names(imodel)
  ncase = model_nsim(imodel)
  grid_type = grids(imodel)
  grid_dir = grid_dirs(imodel)
  yr1 = model_case_yr1(imodel)
  yr2 = model_case_yr2(imodel)

  syr0_mod = yr1+"01"
  syr1_mod = yr2+"12"

;  case_numbers = ispan(2,ncase,1)
  case_numbers = 2
  scase_numbers = sprinti("%0.3i",case_numbers)
  ncase = dimsizes(case_numbers)

  if (mname.eq."MPI-ESM") then
     case_names = mname+"_historical_r"+scase_numbers+"i1850p3"
     syr1_mod = "200512"
  else
     case_names = mname+"_historical_rcp85_r"+case_numbers+"i1p1"
  end if
;  case_names = mname+"_rcp85_r"+scase_numbers+"i2005p3"

  data_dir = "/glade/collections/cdg/data/CLIVAR_LE/"+dirmname+"_lens/OImon/"
  out_dir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"

  yr1_area = 1950                                                  ; start and end yrs for pre-processed time series
  yr1_proc = 2006
;  yr2_proc = 1925
;  yr2_proc = 2005
  yr2_proc = 2100
  yr2_proc = 2099
  nyears = yr2_proc - yr1_proc + 1

  offset = (yr1_proc - yr1)*12

  styrs_proc = sprinti("%0.4i",yr1_proc) + "-"+sprinti("%0.4i",yr2_proc)

  fsic = data_dir+"sic/sic_OImon_"+case_names+"_"+syr0_mod+"-"+syr1_mod+".nc"
  farea = out_dir+"/"+mname+"_regional_areavg.nc"

  if (grid_type.eq."ocn") then
     fgrid = grid_dir+"areacello_fx_"+mname+"_historical_r0i0p0.nc"
  else
     if (grid_type.eq."atm") then
        if (mname.eq."CSIRO-Mk3-6-0") then
           fgrid = grid_dir+"areacella_fx_"+mname+"_historical_r0i0p0.nc"
        else
           fgrid = grid_dir+"areacella_fx_"+mname+"_rcp85_r0i0p0.nc"
        end if
     else
        fgrid = grid_dir+"area_1x1.nc"
     end if
  end if

  if (.not.(isfilepresent(out_dir))) then
     system("mkdir "+out_dir)
  end if

;****************************************************************
; Boundary areas for regional masks:
; Northern Hemisphere:
;       Lab            Lab Sea
;       GIN            GIN Sea
;       Bar            Barents/Kara Sea
;       ArcOc          Arctic Ocean
;       Sib            Siberian Sea
;       Beau           Beaufort Sea
;       CArc           Central Arctic
;       Bering         Bering Sea
;       Okhotsk        Sea of Okhotsk
;       Hudson         Hudson Bay
;       CAArch         Canadian Archipelago
;       Chukchi        Chukchi Sea
;       Kara           Kara Sea
;       Laptev         Laptev Sea
;
;Southern Hemisphere:
;       Wed            Weddell Sea
;       Ross           Ross Sea
;       Ind            Indian sector
;       Pac            Pacific sector
;       BAm            Bellingshausen-Amundsen Sea
;
;******************************************************************

  Lab_lat = (/50.,78./)
  Lab_long = (/284.,315./)
  GIN_lat = (/65.,80./)
  GIN_long1 = (/315.,360./)
  GIN_long2 = (/0.,15./)
  Bar_lat = (/66.,80./)
  Bar_long = (/15.,55./)
  ArcOc_lat1 = (/80.,90./)
  ArcOc_lat2 = ArcOc_lat1
  ArcOc_lat3 = (/68.,90./)
  ArcOc_long1 = (/243.,360./)
  ArcOc_long2 = (/0.,100./)
  ArcOc_long3 = (/100.,243./)
  Kara_lat = (/68.,80./)
  Kara_long = (/55.,100./)
  Laptev_lat = (/68.,80./)
  Laptev_long = (/100.,145./)

  Sib_lat = (/68.,80./)
  Sib_long = (/145.,185./)
  Beau_lat = (/68.,80./)
  Beau_long = (/205.,240./)
  Chukchi_lat = (/68.,80./)
  Chukchi_long = (/185.,205./)
;  Okhotsk_lat = (/45.,63./)			; have below - more complicated - but roughly 44-64N, 133-168E
;  Okhotsk_long = (/135.,160./)

  CArc_lat = (/80.,90./)			; Central Arctic
  CArc_long = (/0.,360./)

;  Bering_lat = (/50.,65./)			; have below - more complicated - but roughly 52.5-68N, 162-204E
;  Bering_long = (/163.,205./)

;---------------------------------------------------------------
; These Arctic regions are more complicated geographically
; and will be defined by polygons rather than min,max lat,lon
;--------------------------------------------------------------
  Lab_poly_lat = (/50., 64., 74., 78., 78., 50./)
  Lab_poly_long = (/295., 292., 280., 280., 315., 315./)
  Bering_poly_lat = (/56., 60., 68., 68., 58., 52.5, 52.5/)
  Bering_poly_long = (/162., 162., 178., 200., 204., 190., 180./)
  Okhotsk_poly_lat = (/44., 44., 46., 51., 60., 64., 60., 55., 44./)
  Okhotsk_poly_long = (/142., 145., 150., 157., 162., 168., 140., 133., 142./)
  Hudson_poly_lat = (/50., 66., 66., 70., 70., 68., 62., 60., 50./)
  Hudson_poly_long = (/265., 265., 274., 274., 290., 290., 292., 295., 295./)
  CAArch_poly_lat = (/67., 80., 80., 70., 67./)
  CAArch_poly_long = (/243., 243., 281., 281., 275./)
; ------SH regions------------------
  Wed_lat = (/-90.,-50./)
  Wed_long1 = (/300.,360./)
  Wed_long2 = (/0.,20./)
  Ross_lat = (/-90.,-50./)
  Ross_long = (/160.,220./)
  Ind_lat = (/-90.,-50./)
  Ind_long = (/20.,90./)
  Pac_lat = (/-90.,-50./)
  Pac_long = (/90.,160./)
  BAm_lat = (/-90.,-50./)
  BAm_long = (/220.,300./)
;*****************************
;
;  time = new((/nyears*12/),float)

f = addfile(fgrid,"r")
gdata = addfile(fgrid,"r")
if (grid_type.eq."ocn") then
   tarea = gdata->areacello
tlat = gdata->lat
tlon = gdata->lon
else
   if (grid_type.eq."atm") then
      tarea = gdata->areacella
      tlat = gdata->lat
      tlon = gdata->lon
   else
      tarea = dble2flt(gdata->AREA)
      f = addfile(fsic(0),"r")			; get lat, lon from data file - long's are slightly different than in area file
						; and calculating regional hs, hi, etc so need correct longitudes. 
						; areas will still be the same (regardless of longitude...lats are same in both)
      tlat = f->lat
      tlon = f->lon
   end if
end if
dimlat = dimsizes(tlat)
ndimlat = dimsizes(dimlat)
if ((tlon@units.ne."degrees_east").and.(tlon@units.ne."degrees_E")) then
   print("check longitude units")
   return
end if
if (min(tlon).lt.0.) then				; need all positive
   tlon = where(tlon.lt.0.,tlon+360.,tlon)
end if

print("calculating area totals for "+case_names(0))
f1 = addfile(fsic,"r")

ai_tmp = f1->sic(:11,:,:)				; necessary to prevent wierd dimension error for ai_aug
ai_jan = ai_tmp(0,:,:)
ai_aug = ai_tmp(7,:,:)
delete(ai_tmp)

if (max(ai_jan).gt.90.) then
   ai_jan = ai_jan/100.
   ai_aug = ai_aug/100.                                        ; convert from % to 0-1.
end if
if (ndimlat.eq.2) then
   tlat2d = tlat
   tlon2d = tlon
else
   tlat2d = conform(tarea,tlat,(/0/))
   tlon2d = conform(tarea,tlon,(/1/))
end if

tarea_jan = mask(tarea,(ai_jan.ge.0.00),True)
tarea_aug = mask(tarea,(ai_aug.ge.0.00),True)

tarea_nh = mask(tarea_jan,(tlat2d.gt.0.),True)
tarea_sh = mask(tarea_aug,(tlat2d.lt.0.),True)
;-------------------------------------------------------------
; create masks for Labrador Sea, Bering Sea, Sea of Okhotsk,
; Hudson Bay and Canadian Archipelago here. They are not simply
; rectangles...........
;------------------------------------------------------------
Lab_mask = new(dimsizes(tarea),typeof(tarea),tarea@_FillValue)
Bering_mask = new(dimsizes(tarea),typeof(tarea),tarea@_FillValue)
Okhotsk_mask = new(dimsizes(tarea),typeof(tarea),tarea@_FillValue)
CAArch_mask = new(dimsizes(tarea),typeof(tarea),tarea@_FillValue)
Hudson_mask = new(dimsizes(tarea),typeof(tarea),tarea@_FillValue)
if (ndimlat.eq.2) then
   nlat = dimlat(0)
   nlon = dimlat(1)
else
   nlat = dimsizes(tlat)
   nlon = dimsizes(tlon)
end if
do i=0,nlat-1
   do j=0,nlon-1
      if (ndimlat.eq.2) then
         plat = tlat(i,j)
         plon = tlon(i,j)
      else
         plat = tlat(i)
         plon = tlon(j)
      end if
      if ((.not.(ismissing(plat))).and.(plat.ge.min(Okhotsk_poly_lat)).and.(plon.ge.min(Okhotsk_poly_long))) then
         Lab_inout = gc_inout(plat,plon,Lab_poly_lat,Lab_poly_long)
         Bering_inout = gc_inout(plat,plon,Bering_poly_lat,Bering_poly_long)
         Okhotsk_inout = gc_inout(plat,plon,Okhotsk_poly_lat,Okhotsk_poly_long)
         Hudson_inout = gc_inout(plat,plon,Hudson_poly_lat,Hudson_poly_long)
         CAArch_inout = gc_inout(plat,plon,CAArch_poly_lat,CAArch_poly_long)
         if (Lab_inout .eq. True) then
	    Lab_mask(i,j) = 1.
	 end if
         if (Bering_inout .eq. True) then
            Bering_mask(i,j) = 1.
         end if
         if (Okhotsk_inout .eq. True) then
            Okhotsk_mask(i,j) = 1.
         end if
         if (Hudson_inout .eq. True) then
            Hudson_mask(i,j) = 1.
         end if
         if (CAArch_inout .eq. True) then
            CAArch_mask(i,j) = 1.
         end if
      end if
   end do
end do

tarea_Lab = mask(tarea_jan,(Lab_mask .eq. 1.), True)
 
tarea_GIN = mask(tarea_jan,(tlat2d.ge.GIN_lat(0) .and. tlat2d.le.GIN_lat(1) .and.\
               ((tlon2d.ge.GIN_long1(0) .and. tlon2d.le.GIN_long1(1)) .or.\
	       (tlon2d.ge.GIN_long2(0) .and. tlon2d.le.GIN_long2(1)))), True)
tarea_Bar = mask(tarea_jan,(tlat2d.ge.Bar_lat(0) .and. tlat2d.le.Bar_lat(1) .and.\
               tlon2d.ge.Bar_long(0) .and. tlon2d.le.Bar_long(1)), True)
tarea_ArcOc = mask(tarea_jan,((tlat2d.ge.ArcOc_lat1(0) .and. (tlon2d.ge.ArcOc_long1(0) .or.\
                  tlon2d.le.ArcOc_long2(1))) .or. (tlat2d.ge.ArcOc_lat3(0) .and.\
                  (tlon2d.ge.ArcOc_long3(0) .and. tlon2d.le.ArcOc_long3(1)))), True)
tarea_Kara = mask(tarea_jan,(tlat2d.ge.Kara_lat(0) .and. tlat2d.le.Kara_lat(1) .and.\
                 tlon2d.ge.Kara_long(0) .and. tlon2d.le.Kara_long(1)), True)
tarea_Laptev = mask(tarea_jan,(tlat2d.ge.Laptev_lat(0) .and. tlat2d.le.Laptev_lat(1) .and.\
                  tlon2d.ge.Laptev_long(0) .and. tlon2d.le.Laptev_long(1)), True)
tarea_Chukchi = mask(tarea_jan,(tlat2d.ge.Chukchi_lat(0) .and. tlat2d.le.Chukchi_lat(1) .and.\
                  tlon2d.ge.Chukchi_long(0) .and. tlon2d.le.Chukchi_long(1)), True)

tarea_Sib = mask(tarea_jan,(tlat2d.ge.Sib_lat(0) .and. tlat2d.le.Sib_lat(1) .and.\
                  tlon2d.ge.Sib_long(0) .and. tlon2d.le.Sib_long(1)), True)

tarea_Beau = mask(tarea_jan,(tlat2d.ge.Beau_lat(0) .and. tlat2d.le.Beau_lat(1) .and.\
                  tlon2d.ge.Beau_long(0) .and. tlon2d.le.Beau_long(1)), True)

tarea_CArc = mask(tarea_jan,(tlat2d.ge.CArc_lat(0) .and. tlat2d.le.CArc_lat(1) .and.\
                  tlon2d.ge.CArc_long(0) .and. tlon2d.le.CArc_long(1)), True)

tarea_Bering = mask(tarea_jan,(Bering_mask .eq. 1.), True)
tarea_Okhotsk = mask(tarea_jan,(Okhotsk_mask .eq. 1.), True)
tarea_Hudson = mask(tarea_jan,(Hudson_mask .eq. 1.), True)
tarea_CAArch = mask(tarea_jan,(CAArch_mask .eq. 1.), True)

tarea_Wed = mask(tarea_aug,(tlat2d.le.Wed_lat(1) .and. ((tlon2d.ge.Wed_long1(0) \
               .and. tlon2d.le.Wed_long1(1)) .or. (tlon2d.ge.Wed_long2(0) \
               .and. tlon2d.le.Wed_long2(1)))), True)

tarea_Ross = mask(tarea_aug,(tlat2d.le.Ross_lat(1) .and. tlon2d.ge.Ross_long(0) \
		.and. tlon2d.le.Ross_long(1)), True)

tarea_Ind = mask(tarea_aug,(tlat2d.le.Ind_lat(1) .and. tlon2d.ge.Ind_long(0) \
               .and. tlon2d.le.Ind_long(1)), True)

tarea_Pac = mask(tarea_aug,(tlat2d.le.Pac_lat(1) .and. tlon2d.ge.Pac_long(0) \
               .and. tlon2d.le.Pac_long(1)), True)

tarea_BAm = mask(tarea_aug,(tlat2d.le.BAm_lat(1) .and. tlon2d.ge.BAm_long(0) \
               .and. tlon2d.le.BAm_long(1)), True)
;********************************
; northern hemisphere areas ...
;********************************

area_nh = dim_sum(dim_sum(tarea_nh))*1.0e-13   
area_Lab = dim_sum(dim_sum(tarea_Lab))*1.0e-13 
area_GIN = dim_sum(dim_sum(tarea_GIN))*1.0e-13 
 
area_Bar = dim_sum(dim_sum(tarea_Bar))*1.0e-13 
area_ArcOc = dim_sum(dim_sum(tarea_ArcOc))*1.0e-13 

area_Kara = dim_sum(dim_sum(tarea_Kara))*1.0e-13

area_Laptev = dim_sum(dim_sum(tarea_Laptev))*1.0e-13
area_Chukchi = dim_sum(dim_sum(tarea_Chukchi))*1.0e-13
area_Sib = dim_sum(dim_sum(tarea_Sib))*1.0e-13 
area_Beau = dim_sum(dim_sum(tarea_Beau))*1.0e-13 
area_CArc = dim_sum(dim_sum(tarea_CArc))*1.0e-13 
area_Bering = dim_sum(dim_sum(tarea_Bering))*1.0e-13 

area_Okhotsk = dim_sum(dim_sum(tarea_Okhotsk))*1.0e-13 
area_Hudson = dim_sum(dim_sum(tarea_Hudson))*1.0e-13 
area_CAArch = dim_sum(dim_sum(tarea_CAArch))*1.0e-13 

;********************************
; southern hemisphere ...
;********************************

area_sh = dim_sum(dim_sum(tarea_sh))*1.0e-13 
area_Wed = dim_sum(dim_sum(tarea_Wed))*1.0e-13 
area_Ross = dim_sum(dim_sum(tarea_Ross))*1.0e-13 
area_Ind = dim_sum(dim_sum(tarea_Ind))*1.0e-13 
area_Pac = dim_sum(dim_sum(tarea_Pac))*1.0e-13 
area_BAm = dim_sum(dim_sum(tarea_BAm))*1.0e-13 

;********************************
; units and save...
;********************************

area_nh@units = "1.E+13 m2"
area_Lab@units = "1.E+13 m2"
area_GIN@units = "1.E+13 m2"
area_Bar@units = "1.E+13 m2"
area_ArcOc@units = "1.E+13 m2"
area_Kara@units = "1.E+13 m2"
area_Laptev@units = "1.E+13 m2"
area_Chukchi@units = "1.E+13 m2"
area_Sib@units = "1.E+13 m2"
area_Beau@units = "1.E+13 m2"
area_CArc@units = "1.E+13 m2"
area_Bering@units = "1.E+13 m2"
area_Okhotsk@units = "1.E+13 m2"
area_Hudson@units = "1.E+13 m2"
area_CAArch@units = "1.E+13 m2"

area_sh@units = "1.E+13 m2"
area_Wed@units = "1.E+13 m2"
area_Ross@units = "1.E+13 m2"
area_Ind@units = "1.E+13 m2"
area_Pac@units = "1.E+13 m2"
area_BAm@units = "1.E+13 m2"

system("/bin/rm -f "+farea)
fout = addfile(farea,"c")
setfileoption(fout,"DefineMode",True)

fileAtt = True
fileAtt@creation_date = systemfunc("date")
fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
fileattdef(fout,fileAtt)

fout->area_nh = area_nh
fout->area_Lab = area_Lab
fout->area_GIN = area_GIN
fout->area_Bar = area_Bar
fout->area_ArcOc = area_ArcOc
fout->area_Kara = area_Kara
fout->area_Laptev = area_Laptev
fout->area_Chukchi = area_Chukchi
fout->area_Sib = area_Sib
fout->area_Beau = area_Beau
fout->area_CArc = area_CArc
fout->area_Bering = area_Bering
fout->area_Okhotsk = area_Okhotsk
fout->area_Hudson = area_Hudson
fout->area_CAArch = area_CAArch

fout->area_sh = area_sh
fout->area_Wed = area_Wed
fout->area_Ross = area_Ross
fout->area_Ind = area_Ind
fout->area_Pac = area_Pac
fout->area_BAm = area_BAm

end
