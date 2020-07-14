; Modified to calculate time series from Walsh et al
; extended Arctic sea ice concentrations database.
; LLandrum
; April 2020
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

  case = "Walsh_extended_ArcticSIC"
  datadir = "/glade/p/cesm/pcwg/landrum/NSIDC/Arctic_extended_SeaIce_Walsh/"
  out_dir = datadir

  yr1 = 1850							; start and end yrs for time series
  yr2 = 2017
  nyears = yr2 - yr1 + 1

  styrs_proc = sprinti("%0.4i",yr1) + "-"+sprinti("%0.4i",yr2)

  data_filename = datadir+"G10010_sibt1850_v2.0.nc"
  data = addfile(data_filename,"r")

  ice_vol_filename = out_dir+"/ice_ext_MHolland_regions_"+case+"."+styrs_proc+".nc"

  R = 6371221.3d						; earth's radius in m. Need double precision to keep all
  pi = 4.*atan(1.0)
  rad=(pi/180.)

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
;--------------------------------------------------------------
 
  vai_mo_nh = new((/nyears*12/),float)
  ext_mo_nh = new((/nyears*12/),float)

;*****************************
; regional arrays: Arctic
;*****************************
  vai_mo_Lab = new((/nyears*12/),float)
  ext_mo_Lab = new((/nyears*12/),float)
  vai_mo_GIN = new((/nyears*12/),float)
  ext_mo_GIN = new((/nyears*12/),float)
  vai_mo_Bar = new((/nyears*12/),float)
  ext_mo_Bar = new((/nyears*12/),float)
  vai_mo_ArcOc = new((/nyears*12/),float)
  ext_mo_ArcOc = new((/nyears*12/),float)

  vai_mo_Kara = new((/nyears*12/),float)
  ext_mo_Kara = new((/nyears*12/),float)
  vai_mo_Laptev = new((/nyears*12/),float)
  ext_mo_Laptev = new((/nyears*12/),float)
  vai_mo_Chukchi = new((/nyears*12/),float)
  ext_mo_Chukchi = new((/nyears*12/),float)

  vai_mo_Sib = new((/nyears*12/),float)
  ext_mo_Sib = new((/nyears*12/),float)
  vai_mo_Beau = new((/nyears*12/),float)
  ext_mo_Beau = new((/nyears*12/),float)
  vai_mo_CArc = new((/nyears*12/),float)
  ext_mo_CArc = new((/nyears*12/),float)
  vai_mo_Bering = new((/nyears*12/),float)
  ext_mo_Bering = new((/nyears*12/),float)
  vai_mo_Okhotsk = new((/nyears*12/),float)
  ext_mo_Okhotsk = new((/nyears*12/),float)
  vai_mo_Hudson = new((/nyears*12/),float)
  ext_mo_Hudson = new((/nyears*12/),float)
  vai_mo_CAArch = new((/nyears*12/),float)
  ext_mo_CAArch = new((/nyears*12/),float)
;*****************************
   lat = data->latitude
   lon = data->longitude
   it2 = (nyears*12) - 1
   time = data->time(:it2)
   nlat = dimsizes(lat)
   nlon = dimsizes(lon)

;   coswgt=cos(rad*lat)
;   delta_lat = 0.25
;   delta_lon = 0.25
;   delta_lat_rad = delta_lat*rad
;   delta_lon_rad = delta_lon*rad

;   dx = R*coswgt*(delta_lon_rad)				; calculating grid cell area as dx*dy where dx, dy are
;   dy = R*(delta_lat_rad)				        ; length of segment of a circle. R = earth's radius (constant
;   area1d = dx*dy						; for dlat calculations), and R cos(lat) is radius for dlon

   area1d = data->Gridcell_Area
   area1d = area1d*1.e+6					; Gridcell_Area is in km2. Convert to m2

   ai = new((/12,nlat,nlon/),float)
   lat3d = conform(ai,lat,(/1/))
   lon3d = conform(ai,lon,(/2/))
   area3d = conform(ai,area1d,(/1/))
;   f_area3d = doubletofloat(area3d) 
   f_area3d = area3d
;-------------------------------------------------------------
; create masks for Labrador Sea, Bering Sea, Sea of Okhotsk,
; Hudson Bay and Canadian Archipelago here. They are not simply
; rectangles...........
; Only need to do this once (first year of nyears)
;------------------------------------------------------------
   Lab_mask = ai
   Bering_mask = ai
   Okhotsk_mask = ai
   CAArch_mask = ai
   Hudson_mask = ai
   do i=0,nlat-1
      do j=0,nlon-1
         plat = lat(i)
         plon = lon(j)
         if ((.not.(ismissing(plat))).and.(plat.ge.min(Okhotsk_poly_lat)).and.(plon.ge.min(Okhotsk_poly_long))) then
            Lab_inout = gc_inout(plat,plon,Lab_poly_lat,Lab_poly_long)
            Bering_inout = gc_inout(plat,plon,Bering_poly_lat,Bering_poly_long)
            Okhotsk_inout = gc_inout(plat,plon,Okhotsk_poly_lat,Okhotsk_poly_long)
            Hudson_inout = gc_inout(plat,plon,Hudson_poly_lat,Hudson_poly_long)
            CAArch_inout = gc_inout(plat,plon,CAArch_poly_lat,CAArch_poly_long)
            if (Lab_inout .eq. True) then
               Lab_mask(:,i,j) = 1.
            end if
            if (Bering_inout .eq. True) then
               Bering_mask(:,i,j) = 1.
            end if
            if (Okhotsk_inout .eq. True) then
               Okhotsk_mask(:,i,j) = 1.
            end if
            if (Hudson_inout .eq. True) then
               Hudson_mask(:,i,j) = 1.
            end if
            if (CAArch_inout .eq. True) then
               CAArch_mask(:,i,j) = 1.
            end if
         end if
      end do
   end do

   do nn = 0, nyears-1
      year = yr1 + nn

      print("processing year = "+ year)

      imo1 = 12*nn
      imo2 = 12*nn+11

      imo1_datain = 12*nn
      imo2_datain = 12*nn+11

      ai = (/ tofloat(data->seaice_conc(imo1_datain:imo2_datain,:,:)) /)
      ai = where(ai.gt.100.,ai@_FillValue,ai)				; this data set uses 120 for land grid

      vai = ai*f_area3d							; units are 10^-2 m^2
      ext = where(ai.ge.15.,1.,0.)*f_area3d

      vai_nh = vai
      ext_nh = ext

      vai_Lab = mask(vai,(Lab_mask .eq. 1.), True)
      ext_Lab = mask(ext,(Lab_mask .eq. 1.), True)
 
      vai_GIN = mask(vai,(lat3d.ge.GIN_lat(0) .and. lat3d.le.GIN_lat(1) .and.\
                 ((lon3d.ge.GIN_long1(0) .and. lon3d.le.GIN_long1(1)) .or.\
                 (lon3d.ge.GIN_long2(0) .and. lon3d.le.GIN_long2(1)))), True)
      ext_GIN = mask(ext,(lat3d.ge.GIN_lat(0) .and. lat3d.le.GIN_lat(1) .and.\
                 ((lon3d.ge.GIN_long1(0) .and. lon3d.le.GIN_long1(1)) .or.\
                 (lon3d.ge.GIN_long2(0) .and. lon3d.le.GIN_long2(1)))), True)
      vai_Bar = mask(vai,(lat3d.ge.Bar_lat(0) .and. lat3d.le.Bar_lat(1) .and.\
                  lon3d.ge.Bar_long(0) .and. lon3d.le.Bar_long(1)), True)
      ext_Bar = mask(ext,(lat3d.ge.Bar_lat(0) .and. lat3d.le.Bar_lat(1) .and.\
                 lon3d.ge.Bar_long(0) .and. lon3d.le.Bar_long(1)), True)

      vai_ArcOc = mask(vai,((lat3d.ge.ArcOc_lat1(0) .and. (lon3d.ge.ArcOc_long1(0) .or.\
                  lon3d.le.ArcOc_long2(1))) .or. (lat3d.ge.ArcOc_lat3(0) .and.\
                  (lon3d.ge.ArcOc_long3(0) .and. lon3d.le.ArcOc_long3(1)))), True)
      ext_ArcOc = mask(ext,((lat3d.ge.ArcOc_lat1(0) .and. (lon3d.ge.ArcOc_long1(0) .or.\
                  lon3d.le.ArcOc_long2(1))) .or. (lat3d.ge.ArcOc_lat3(0) .and.\
                  (lon3d.ge.ArcOc_long3(0) .and. lon3d.le.ArcOc_long3(1)))), True)

      vai_Kara = mask(vai,(lat3d.ge.Kara_lat(0) .and. lat3d.le.Kara_lat(1) .and.\
                  lon3d.ge.Kara_long(0) .and. lon3d.le.Kara_long(1)), True)
      ext_Kara = mask(ext,(lat3d.ge.Kara_lat(0) .and. lat3d.le.Kara_lat(1) .and.\
                  lon3d.ge.Kara_long(0) .and. lon3d.le.Kara_long(1)), True)

      vai_Laptev = mask(vai,(lat3d.ge.Laptev_lat(0) .and. lat3d.le.Laptev_lat(1) .and.\
                  lon3d.ge.Laptev_long(0) .and. lon3d.le.Laptev_long(1)), True)
      ext_Laptev = mask(ext,(lat3d.ge.Laptev_lat(0) .and. lat3d.le.Laptev_lat(1) .and.\
                  lon3d.ge.Laptev_long(0) .and. lon3d.le.Laptev_long(1)), True)

      vai_Chukchi = mask(vai,(lat3d.ge.Chukchi_lat(0) .and. lat3d.le.Chukchi_lat(1) .and.\
                  lon3d.ge.Chukchi_long(0) .and. lon3d.le.Chukchi_long(1)), True)
      ext_Chukchi = mask(ext,(lat3d.ge.Chukchi_lat(0) .and. lat3d.le.Chukchi_lat(1) .and.\
                  lon3d.ge.Chukchi_long(0) .and. lon3d.le.Chukchi_long(1)), True)

      vai_Sib = mask(vai,(lat3d.ge.Sib_lat(0) .and. lat3d.le.Sib_lat(1) .and.\
                  lon3d.ge.Sib_long(0) .and. lon3d.le.Sib_long(1)), True)
      ext_Sib = mask(ext,(lat3d.ge.Sib_lat(0) .and. lat3d.le.Sib_lat(1) .and.\
                  lon3d.ge.Sib_long(0) .and. lon3d.le.Sib_long(1)), True)

      vai_Beau = mask(vai,(lat3d.ge.Beau_lat(0) .and. lat3d.le.Beau_lat(1) .and.\
                  lon3d.ge.Beau_long(0) .and. lon3d.le.Beau_long(1)), True)
      ext_Beau = mask(ext,(lat3d.ge.Beau_lat(0) .and. lat3d.le.Beau_lat(1) .and.\
                  lon3d.ge.Beau_long(0) .and. lon3d.le.Beau_long(1)), True)

      vai_CArc = mask(vai,(lat3d.ge.CArc_lat(0) .and. lat3d.le.CArc_lat(1) .and.\
                  lon3d.ge.CArc_long(0) .and. lon3d.le.CArc_long(1)), True)
      ext_CArc = mask(ext,(lat3d.ge.CArc_lat(0) .and. lat3d.le.CArc_lat(1) .and.\
                  lon3d.ge.CArc_long(0) .and. lon3d.le.CArc_long(1)), True)

      vai_Bering = mask(vai,(Bering_mask .eq. 1.), True)
      ext_Bering = mask(ext,(Bering_mask .eq. 1.), True)

      vai_Okhotsk = mask(vai,(Okhotsk_mask .eq. 1.), True)
      ext_Okhotsk = mask(ext,(Okhotsk_mask .eq. 1.), True)

      vai_Hudson = mask(vai,(Hudson_mask .eq. 1.), True)
      ext_Hudson = mask(ext,(Hudson_mask .eq. 1.), True)

      vai_CAArch = mask(vai,(CAArch_mask .eq. 1.), True)
      ext_CAArch = mask(ext,(CAArch_mask .eq. 1.), True)

;********************************
; monthly data ...
;********************************

      vai_mo_nh(imo1:imo2) = dim_sum(dim_sum(vai_nh))*1.0e-14 
      ext_mo_nh(imo1:imo2) = dim_sum(dim_sum(ext_nh))*1.0e-12

      vai_mo_Lab(imo1:imo2) = dim_sum(dim_sum(vai_Lab))*1.0e-14 
      ext_mo_Lab(imo1:imo2) = dim_sum(dim_sum(ext_Lab))*1.0e-12
 
      vai_mo_GIN(imo1:imo2) = dim_sum(dim_sum(vai_GIN))*1.0e-14 
      ext_mo_GIN(imo1:imo2) = dim_sum(dim_sum(ext_GIN))*1.0e-12
 
      vai_mo_Bar(imo1:imo2) = dim_sum(dim_sum(vai_Bar))*1.0e-14 
      ext_mo_Bar(imo1:imo2) = dim_sum(dim_sum(ext_Bar))*1.0e-12
 
      vai_mo_ArcOc(imo1:imo2) = dim_sum(dim_sum(vai_ArcOc))*1.0e-14 
      ext_mo_ArcOc(imo1:imo2) = dim_sum(dim_sum(ext_ArcOc))*1.0e-12

      vai_mo_Kara(imo1:imo2) = dim_sum(dim_sum(vai_Kara))*1.0e-14
      ext_mo_Kara(imo1:imo2) = dim_sum(dim_sum(ext_Kara))*1.0e-12

      vai_mo_Laptev(imo1:imo2) = dim_sum(dim_sum(vai_Laptev))*1.0e-14
      ext_mo_Laptev(imo1:imo2) = dim_sum(dim_sum(ext_Laptev))*1.0e-12

      vai_mo_Chukchi(imo1:imo2) = dim_sum(dim_sum(vai_Chukchi))*1.0e-14
      ext_mo_Chukchi(imo1:imo2) = dim_sum(dim_sum(ext_Chukchi))*1.0e-12

      vai_mo_Sib(imo1:imo2) = dim_sum(dim_sum(vai_Sib))*1.0e-14 
      ext_mo_Sib(imo1:imo2) = dim_sum(dim_sum(ext_Sib))*1.0e-12

      vai_mo_Beau(imo1:imo2) = dim_sum(dim_sum(vai_Beau))*1.0e-14 
      ext_mo_Beau(imo1:imo2) = dim_sum(dim_sum(ext_Beau))*1.0e-12
 
      vai_mo_CArc(imo1:imo2) = dim_sum(dim_sum(vai_CArc))*1.0e-14 
      ext_mo_CArc(imo1:imo2) = dim_sum(dim_sum(ext_CArc))*1.0e-12

      vai_mo_Bering(imo1:imo2) = dim_sum(dim_sum(vai_Bering))*1.0e-14 
      ext_mo_Bering(imo1:imo2) = dim_sum(dim_sum(ext_Bering))*1.0e-12

      vai_mo_Okhotsk(imo1:imo2) = dim_sum(dim_sum(vai_Okhotsk))*1.0e-14 
      ext_mo_Okhotsk(imo1:imo2) = dim_sum(dim_sum(ext_Okhotsk))*1.0e-12

      vai_mo_Hudson(imo1:imo2) = dim_sum(dim_sum(vai_Hudson))*1.0e-14 
      ext_mo_Hudson(imo1:imo2) = dim_sum(dim_sum(ext_Hudson))*1.0e-12

      vai_mo_CAArch(imo1:imo2) = dim_sum(dim_sum(vai_CAArch))*1.0e-14 
      ext_mo_CAArch(imo1:imo2) = dim_sum(dim_sum(ext_CAArch))*1.0e-12

   end do ; nn


system("/bin/rm -f "+ice_vol_filename)
fout = addfile(ice_vol_filename,"c")
setfileoption(fout,"DefineMode",True)

fileAtt = True
fileAtt@creation_date = systemfunc("date")
fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
fileattdef(fout,fileAtt)

dimNames = (/"time"/)
dimSizes = (/-1/)
dimUnlim = (/True/)
filedimdef(fout,dimNames,dimSizes,dimUnlim)

filevardef(fout,"time",typeof(time),getvardims(time))
filevarattdef(fout,"time",time)

vai_mo_nh!0 = "time"
vai_mo_nh&time = time
vai_mo_nh@units = "1.E+12 m2"
vai_mo_nh@coordinates = "time"
vai_mo_nh@missing_value = 1.e+30
vai_mo_nh@_FillValue = 1.e+30
vai_mo_nh@cell_methods = "time:mean"

ext_mo_nh!0 = "time"
ext_mo_nh&time = time
ext_mo_nh@units = "1.E+12 m2"
ext_mo_nh@coordinates = "time"
ext_mo_nh@missing_value = 1.e+30
ext_mo_nh@_FillValue = 1.e+30
ext_mo_nh@cell_methods = "time:mean"

copy_VarMeta(vai_mo_nh,vai_mo_Lab)
copy_VarMeta(ext_mo_nh,ext_mo_Lab)

copy_VarMeta(vai_mo_nh,vai_mo_GIN)
copy_VarMeta(ext_mo_nh,ext_mo_GIN)

copy_VarMeta(vai_mo_nh,vai_mo_Bar)
copy_VarMeta(ext_mo_nh,ext_mo_Bar)

copy_VarMeta(vai_mo_nh,vai_mo_ArcOc)
copy_VarMeta(ext_mo_nh,ext_mo_ArcOc)

copy_VarMeta(vai_mo_nh,vai_mo_Kara)
copy_VarMeta(ext_mo_nh,ext_mo_Kara)

copy_VarMeta(vai_mo_nh,vai_mo_Laptev)
copy_VarMeta(ext_mo_nh,ext_mo_Laptev)

copy_VarMeta(vai_mo_nh,vai_mo_Chukchi)
copy_VarMeta(ext_mo_nh,ext_mo_Chukchi)

copy_VarMeta(vai_mo_nh,vai_mo_Sib)
copy_VarMeta(ext_mo_nh,ext_mo_Sib)

copy_VarMeta(vai_mo_nh,vai_mo_Beau)
copy_VarMeta(ext_mo_nh,ext_mo_Beau)

copy_VarMeta(vai_mo_nh,vai_mo_CArc)
copy_VarMeta(ext_mo_nh,ext_mo_CArc)

copy_VarMeta(vai_mo_nh,vai_mo_Bering)
copy_VarMeta(ext_mo_nh,ext_mo_Bering)

copy_VarMeta(vai_mo_nh,vai_mo_Okhotsk)
copy_VarMeta(ext_mo_nh,ext_mo_Okhotsk)

copy_VarMeta(vai_mo_nh,vai_mo_Hudson)
copy_VarMeta(ext_mo_nh,ext_mo_Hudson)

copy_VarMeta(vai_mo_nh,vai_mo_CAArch)
copy_VarMeta(ext_mo_nh,ext_mo_CAArch)

filevardef(fout,"vai_mo_nh",typeof(vai_mo_nh),getvardims(vai_mo_nh))
filevardef(fout,"ext_mo_nh",typeof(ext_mo_nh),getvardims(ext_mo_nh))

filevardef(fout,"vai_mo_Lab",typeof(vai_mo_Lab),getvardims(vai_mo_Lab))
filevardef(fout,"ext_mo_Lab",typeof(ext_mo_Lab),getvardims(ext_mo_Lab))

filevardef(fout,"vai_mo_GIN",typeof(vai_mo_GIN),getvardims(vai_mo_GIN))
filevardef(fout,"ext_mo_GIN",typeof(ext_mo_GIN),getvardims(ext_mo_GIN))

filevardef(fout,"vai_mo_Bar",typeof(vai_mo_Bar),getvardims(vai_mo_Bar))
filevardef(fout,"ext_mo_Bar",typeof(ext_mo_Bar),getvardims(ext_mo_Bar))

filevardef(fout,"vai_mo_ArcOc",typeof(vai_mo_ArcOc),getvardims(vai_mo_ArcOc))
filevardef(fout,"ext_mo_ArcOc",typeof(ext_mo_ArcOc),getvardims(ext_mo_ArcOc))

filevardef(fout,"vai_mo_Kara",typeof(vai_mo_Kara),getvardims(vai_mo_Kara))
filevardef(fout,"ext_mo_Kara",typeof(ext_mo_Kara),getvardims(ext_mo_Kara))

filevardef(fout,"vai_mo_Laptev",typeof(vai_mo_Laptev),getvardims(vai_mo_Laptev))
filevardef(fout,"ext_mo_Laptev",typeof(ext_mo_Laptev),getvardims(ext_mo_Laptev))

filevardef(fout,"vai_mo_Chukchi",typeof(vai_mo_Chukchi),getvardims(vai_mo_Chukchi))
filevardef(fout,"ext_mo_Chukchi",typeof(ext_mo_Chukchi),getvardims(ext_mo_Chukchi))

filevardef(fout,"vai_mo_Sib",typeof(vai_mo_Sib),getvardims(vai_mo_Sib))
filevardef(fout,"ext_mo_Sib",typeof(ext_mo_Sib),getvardims(ext_mo_Sib))

filevardef(fout,"vai_mo_Beau",typeof(vai_mo_Beau),getvardims(vai_mo_Beau))
filevardef(fout,"ext_mo_Beau",typeof(ext_mo_Beau),getvardims(ext_mo_Beau))

filevardef(fout,"vai_mo_CArc",typeof(vai_mo_CArc),getvardims(vai_mo_CArc))
filevardef(fout,"ext_mo_CArc",typeof(ext_mo_CArc),getvardims(ext_mo_CArc))

filevardef(fout,"vai_mo_Bering",typeof(vai_mo_Bering),getvardims(vai_mo_Bering))
filevardef(fout,"ext_mo_Bering",typeof(ext_mo_Bering),getvardims(ext_mo_Bering))

filevardef(fout,"vai_mo_Okhotsk",typeof(vai_mo_Okhotsk),getvardims(vai_mo_Okhotsk))
filevardef(fout,"ext_mo_Okhotsk",typeof(ext_mo_Okhotsk),getvardims(ext_mo_Okhotsk))

filevardef(fout,"vai_mo_Hudson",typeof(vai_mo_Hudson),getvardims(vai_mo_Hudson))
filevardef(fout,"ext_mo_Hudson",typeof(ext_mo_Hudson),getvardims(ext_mo_Hudson))

filevardef(fout,"vai_mo_CAArch",typeof(vai_mo_CAArch),getvardims(vai_mo_CAArch))
filevardef(fout,"ext_mo_CAArch",typeof(ext_mo_CAArch),getvardims(ext_mo_CAArch))

filevarattdef(fout,"vai_mo_nh",vai_mo_nh)
filevarattdef(fout,"ext_mo_nh",ext_mo_nh)

filevarattdef(fout,"vai_mo_Lab",vai_mo_Lab)
filevarattdef(fout,"ext_mo_Lab",vai_mo_Lab)

filevarattdef(fout,"vai_mo_GIN",vai_mo_GIN)
filevarattdef(fout,"ext_mo_GIN",ext_mo_GIN)

filevarattdef(fout,"vai_mo_Bar",vai_mo_Bar)
filevarattdef(fout,"ext_mo_Bar",ext_mo_Bar)

filevarattdef(fout,"vai_mo_ArcOc",vai_mo_ArcOc)
filevarattdef(fout,"ext_mo_ArcOc",ext_mo_ArcOc)

filevarattdef(fout,"vai_mo_Kara",vai_mo_Kara)
filevarattdef(fout,"ext_mo_Kara",ext_mo_Kara)

filevarattdef(fout,"vai_mo_Laptev",vai_mo_Laptev)
filevarattdef(fout,"ext_mo_Laptev",ext_mo_Laptev)

filevarattdef(fout,"vai_mo_Chukchi",vai_mo_Chukchi)
filevarattdef(fout,"ext_mo_Chukchi",ext_mo_Chukchi)

filevarattdef(fout,"vai_mo_Sib",vai_mo_Sib)
filevarattdef(fout,"ext_mo_Sib",ext_mo_Sib)

filevarattdef(fout,"vai_mo_Beau",vai_mo_Beau)
filevarattdef(fout,"ext_mo_Beau",ext_mo_Beau)

filevarattdef(fout,"vai_mo_CArc",vai_mo_CArc)
filevarattdef(fout,"ext_mo_CArc",ext_mo_CArc)

filevarattdef(fout,"vai_mo_Bering",vai_mo_Bering)
filevarattdef(fout,"ext_mo_Bering",ext_mo_Bering)
filevarattdef(fout,"vai_mo_Okhotsk",vai_mo_Okhotsk)
filevarattdef(fout,"ext_mo_Okhotsk",ext_mo_Okhotsk)
filevarattdef(fout,"vai_mo_Hudson",vai_mo_Hudson)
filevarattdef(fout,"ext_mo_Hudson",vai_mo_Hudson)
filevarattdef(fout,"vai_mo_CAArch",vai_mo_CAArch)
filevarattdef(fout,"ext_mo_CAArch",ext_mo_CAArch)

fout->time = time
;fout->time = (/ time /)
fout->vai_mo_nh = vai_mo_nh
fout->ext_mo_nh = ext_mo_nh

fout->vai_mo_Lab = vai_mo_Lab
fout->ext_mo_Lab = ext_mo_Lab
fout->vai_mo_GIN = vai_mo_GIN
fout->ext_mo_GIN = ext_mo_GIN
fout->vai_mo_Bar = vai_mo_Bar
fout->ext_mo_Bar = ext_mo_Bar

fout->vai_mo_ArcOc = vai_mo_ArcOc
fout->ext_mo_ArcOc = ext_mo_ArcOc

fout->vai_mo_Kara = vai_mo_Kara
fout->ext_mo_Kara = ext_mo_Kara

fout->vai_mo_Laptev = vai_mo_Laptev
fout->ext_mo_Laptev = ext_mo_Laptev

fout->vai_mo_Chukchi = vai_mo_Chukchi
fout->ext_mo_Chukchi = ext_mo_Chukchi

fout->vai_mo_Sib = vai_mo_Sib
fout->ext_mo_Sib = ext_mo_Sib
fout->vai_mo_Beau = vai_mo_Beau
fout->ext_mo_Beau = ext_mo_Beau
fout->vai_mo_CArc = vai_mo_CArc
fout->ext_mo_CArc = ext_mo_CArc

fout->vai_mo_Bering = vai_mo_Bering
fout->ext_mo_Bering = ext_mo_Bering
fout->vai_mo_Okhotsk = vai_mo_Okhotsk
fout->ext_mo_Okhotsk = ext_mo_Okhotsk
fout->vai_mo_Hudson = vai_mo_Hudson
fout->ext_mo_Hudson = ext_mo_Hudson
fout->vai_mo_CAArch = vai_mo_CAArch
fout->ext_mo_CAArch = ext_mo_CAArch


end
