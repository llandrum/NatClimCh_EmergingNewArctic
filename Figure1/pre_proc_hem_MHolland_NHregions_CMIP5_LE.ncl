; Modified for multi-model large ensemble
; LLandrum
; Jan 2020
;
; CMIP5 MMLE notes
;
; MPI has historical (1850-1920) and rcp8.5 stored separately
;
; mname=(/"CanESM2", "CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2M","MPI-ESM","MPI-ESM"/)
; dirname=(/"canesm2", "csiro_mk36","gfdl_cm3","gfdl_esm2m","mpi","mpi"/)
; n_ens=(/50,30,20,30,100,100/)
; syr0_mod=(/195001,185001,192001,195001,185001,200601/)
; syr1_mod=(/210012,210012,210012,210012,200512,209912/)

;***************************************************************
; This script is a copy of pre_process.ncl from the ice diagnostics
; package modified to use already processed monthly aice, hi, hs 
; files.
; Laura Landrum
; December, 2010
;***************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  mname = "CESM1-CAM5"
  model_names = (/"CanESM2","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2M","MPI-ESM","CESM1-CAM5"/)
  dir_names = (/"canesm2","csiro_mk36","gfdl_cm3","gfdl_esm2m","mpi","cesm"/)
  model_case_yr1 = (/1950,1850,1920,1950,1850,1920/)
  model_case_yr2 = (/2100,2100,2100,2100,2099,2100/)
  model_nsim = (/50,29,20,30,100,40/)
  model_snow = (/True,False,False,False,False,False/)

  imodel = ind(model_names.eq.mname)
  dirmname = dir_names(imodel)
  ncase = model_nsim(imodel)
  snow_exists = model_snow(imodel)
  yr1 = model_case_yr1(imodel)
  yr2 = model_case_yr2(imodel)

  syr0_mod = yr1+"01"
  syr1_mod = yr2+"12"

;  syr0_mod = 192001
;  syr0_mod = 185001
;  syr0_mod = 200601
;  syr1_mod = 200512
;  syr1_mod = 210012
;  syr1_mod = 209912
  case_numbers = ispan(2,ncase,1)
;  case_numbers = 1
;  case_numbers = ispan(3,ncase,1)
  scase_numbers = sprinti("%0.3i",case_numbers)
  ncase = dimsizes(case_numbers)

  case_names = mname+"_historical_rcp85_r"+case_numbers+"i1p1"
;  case_names = mname+"_historical_r"+scase_numbers+"i1850p3"
;  case_names = mname+"_rcp85_r"+scase_numbers+"i2005p3"

  data_dir = "/glade/collections/cdg/data/CLIVAR_LE/"+dirmname+"_lens/OImon/"
  out_dir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"
  grid_dir = "/glade/p/cesm/pcwg/mholland/cmip5/fx/"
;  grid_dir = "/glade/p/cesm/pcwg/landrum/CMIP5/fx/"
;  grid_type = "atm"					; atm for areacella ocn for areacello
;  grid_type = "1x1"					; or 1x1 for calculated from PI*r2*cos for regular 1 deg grid
  grid_type = "ocn"

;  yr1 = floattoint(floor(syr0_mod/100))
;  yr2 = floattoint(floor(syr1_mod/100))

  yr1_proc = 1920                                                  ; start and end yrs for pre-processed time series
;  yr1_proc = 2006
;  yr2_proc = 1925
;  yr2_proc = 2005
  yr2_proc = 2100
  yr2_proc = 2099
  nyears = yr2_proc - yr1_proc + 1

  offset = (yr1_proc - yr1)*12

  styrs_proc = sprinti("%0.4i",yr1_proc) + "-"+sprinti("%0.4i",yr2_proc)

  fsic = data_dir+"sic/sic_OImon_"+case_names+"_"+syr0_mod+"-"+syr1_mod+".nc"
  fsit = data_dir+"sit/sit_OImon_"+case_names+"_"+syr0_mod+"-"+syr1_mod+".nc"
  fsnow = data_dir+"snd/snd_OImon_"+case_names+"_"+syr0_mod+"-"+syr1_mod+".nc"

  if (grid_type.eq."ocn") then
     fgrid = grid_dir+"areacello_fx_"+mname+"_historical_r0i0p0.nc"
  else
     if (grid_type.eq."atm") then
;        fgrid = grid_dir+"areacella_fx_"+mname+"_historical_r0i0p0.nc"
        fgrid = grid_dir+"areacella_fx_"+mname+"_rcp85_r0i0p0.nc"
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

  vai_mo_nh = new((/nyears*12/),float)
  ext_mo_nh = new((/nyears*12/),float)
  vhi_mo_nh = new((/nyears*12/),float)
  vhs_mo_nh = new((/nyears*12/),float)
  vai_mo_sh = new((/nyears*12/),float)
  ext_mo_sh = new((/nyears*12/),float)
  vhi_mo_sh = new((/nyears*12/),float)
  vhs_mo_sh = new((/nyears*12/),float)

;*****************************
; regional arrays: Arctic
;*****************************
  vai_mo_Lab = new((/nyears*12/),float)
  ext_mo_Lab = new((/nyears*12/),float)
  vhi_mo_Lab = new((/nyears*12/),float)
  vhs_mo_Lab = new((/nyears*12/),float)
  vai_mo_GIN = new((/nyears*12/),float)
  ext_mo_GIN = new((/nyears*12/),float)
  vhi_mo_GIN = new((/nyears*12/),float)
  vhs_mo_GIN = new((/nyears*12/),float)
  vai_mo_Bar = new((/nyears*12/),float)
  ext_mo_Bar = new((/nyears*12/),float)
  vhi_mo_Bar = new((/nyears*12/),float)
  vhs_mo_Bar = new((/nyears*12/),float)
  vai_mo_ArcOc = new((/nyears*12/),float)
  ext_mo_ArcOc = new((/nyears*12/),float)
  vhi_mo_ArcOc = new((/nyears*12/),float)
  vhs_mo_ArcOc = new((/nyears*12/),float)

  vai_mo_Kara = new((/nyears*12/),float)
  ext_mo_Kara = new((/nyears*12/),float)
  vhi_mo_Kara = new((/nyears*12/),float)
  vhs_mo_Kara = new((/nyears*12/),float)
  vai_mo_Laptev = new((/nyears*12/),float)
  ext_mo_Laptev = new((/nyears*12/),float)
  vhi_mo_Laptev = new((/nyears*12/),float)
  vhs_mo_Laptev = new((/nyears*12/),float)
  vai_mo_Chukchi = new((/nyears*12/),float)
  ext_mo_Chukchi = new((/nyears*12/),float)
  vhi_mo_Chukchi = new((/nyears*12/),float)
  vhs_mo_Chukchi = new((/nyears*12/),float)

  vai_mo_Sib = new((/nyears*12/),float)
  ext_mo_Sib = new((/nyears*12/),float)
  vhi_mo_Sib = new((/nyears*12/),float)
  vhs_mo_Sib = new((/nyears*12/),float)
  vai_mo_Beau = new((/nyears*12/),float)
  ext_mo_Beau = new((/nyears*12/),float)
  vhi_mo_Beau = new((/nyears*12/),float)
  vhs_mo_Beau = new((/nyears*12/),float)
  vai_mo_CArc = new((/nyears*12/),float)
  ext_mo_CArc = new((/nyears*12/),float)
  vhi_mo_CArc = new((/nyears*12/),float)
  vhs_mo_CArc = new((/nyears*12/),float)
  vai_mo_Bering = new((/nyears*12/),float)
  ext_mo_Bering = new((/nyears*12/),float)
  vhi_mo_Bering = new((/nyears*12/),float)
  vhs_mo_Bering = new((/nyears*12/),float)
  vai_mo_Okhotsk = new((/nyears*12/),float)
  ext_mo_Okhotsk = new((/nyears*12/),float)
  vhi_mo_Okhotsk = new((/nyears*12/),float)
  vhs_mo_Okhotsk = new((/nyears*12/),float)
  vai_mo_Hudson = new((/nyears*12/),float)
  ext_mo_Hudson = new((/nyears*12/),float)
  vhi_mo_Hudson = new((/nyears*12/),float)
  vhs_mo_Hudson = new((/nyears*12/),float)
  vai_mo_CAArch = new((/nyears*12/),float)
  ext_mo_CAArch = new((/nyears*12/),float)
  vhi_mo_CAArch = new((/nyears*12/),float)
  vhs_mo_CAArch = new((/nyears*12/),float)
;*****************************
; regional arrays: Antarctic
;*****************************
  vai_mo_Wed = new((/nyears*12/),float)
  ext_mo_Wed = new((/nyears*12/),float)
  vhi_mo_Wed = new((/nyears*12/),float)
  vhs_mo_Wed = new((/nyears*12/),float)
  vai_mo_Ross = new((/nyears*12/),float)
  ext_mo_Ross = new((/nyears*12/),float)
  vhi_mo_Ross = new((/nyears*12/),float)
  vhs_mo_Ross = new((/nyears*12/),float)
  vai_mo_Ind = new((/nyears*12/),float)
  ext_mo_Ind = new((/nyears*12/),float)
  vhi_mo_Ind = new((/nyears*12/),float)
  vhs_mo_Ind = new((/nyears*12/),float)
  vai_mo_Pac = new((/nyears*12/),float)
  ext_mo_Pac = new((/nyears*12/),float)
  vhi_mo_Pac = new((/nyears*12/),float)
  vhs_mo_Pac = new((/nyears*12/),float)
  vai_mo_BAm = new((/nyears*12/),float)
  ext_mo_BAm = new((/nyears*12/),float)
  vhi_mo_BAm = new((/nyears*12/),float)
  vhs_mo_BAm = new((/nyears*12/),float)

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
do ii = 0,ncase-1
   print("calculating area totals for "+case_names(ii))
   f1 = addfile(fsic(ii),"r")
   f2 = addfile(fsit(ii),"r")
   if (snow_exists) then
      f3 = addfile(fsnow(ii),"r")
   end if

   do nn = 0, nyears-1

      imo1 = 12*nn
      imo2 = 12*nn+11

      imo1_datain = offset + 12*nn
      imo2_datain = offset + 12*nn+11

      hi = f2->sit(imo1_datain:imo2_datain,:,:)
      ai = f1->sic(imo1_datain:imo2_datain,:,:)
      if (snow_exists) then
         hs = f3->snd(imo1_datain:imo2_datain,:,:)
      end if

      if (max(ai(0,:,:)).gt.90.) then
         ai = ai/100.                                        ; convert from % to 0-1.
      end if
      if ((ii.eq.0).and.(nn.eq.0)) then
         it2 = offset + (nyears*12) - 1
         time = f1->time(imo1_datain:it2)
         tarea3d = conform(ai,tarea,(/1,2/))
         if (ndimlat.eq.2) then
            tlat3d = conform(ai,tlat,(/1,2/))
            tlon3d = conform(ai,tlon,(/1,2/))
         else
            tlat3d = conform(ai,tlat,(/1/))
            tlon3d = conform(ai,tlon,(/2/))
         end if
      end if

      vai = ai*tarea3d
      vhi = hi*tarea3d
      ext = where(ai.ge.0.15,1.,0.)*tarea3d

      vai_nh = mask(vai,(tlat3d.ge.0.),True)
      vhi_nh = mask(vhi,(tlat3d.ge.0.),True)
      ext_nh = mask(ext,(tlat3d.ge.0.),True)

      vai_sh = mask(vai,(tlat3d.ge.0.),False)
      vhi_sh = mask(vhi,(tlat3d.ge.0.),False)
      ext_sh = mask(ext,(tlat3d.ge.0.),False)

      if (snow_exists) then
         hs = f3->snd(imo1_datain:imo2_datain,:,:)
         vhs = hs*tarea3d
         vhs_nh = mask(vhs,(tlat3d.ge.0.),True)
         vhs_sh = mask(vhs,(tlat3d.ge.0.),False)
      end if

;-------------------------------------------------------------
; create masks for Labrador Sea, Bering Sea, Sea of Okhotsk,
; Hudson Bay and Canadian Archipelago here. They are not simply
; rectangles...........
; Only need to do this once (first year of nyears)
;------------------------------------------------------------
      if ((ii.eq.0).and.(nn .eq. 0)) then
         Lab_mask = new(dimsizes(vhi_nh),typeof(vhi_nh),vhi_nh@_FillValue)
         Bering_mask = new(dimsizes(vhi_nh),typeof(vhi_nh),vhi_nh@_FillValue)
         Okhotsk_mask = new(dimsizes(vhi_nh),typeof(vhi_nh),vhi_nh@_FillValue)
         CAArch_mask = new(dimsizes(vhi_nh),typeof(vhi_nh),vhi_nh@_FillValue)
         Hudson_mask = new(dimsizes(vhi_nh),typeof(vhi_nh),vhi_nh@_FillValue)
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
      end if

      vhi_Lab = mask(vhi_nh,(Lab_mask .eq. 1.), True)
      vai_Lab = mask(vai_nh,(Lab_mask .eq. 1.), True)
      ext_Lab = mask(ext_nh,(Lab_mask .eq. 1.), True)
 
      vhi_GIN = mask(vhi_nh,(tlat3d.ge.GIN_lat(0) .and. tlat3d.le.GIN_lat(1) .and.\
                  ((tlon3d.ge.GIN_long1(0) .and. tlon3d.le.GIN_long1(1)) .or.\
		   (tlon3d.ge.GIN_long2(0) .and. tlon3d.le.GIN_long2(1)))), True)
      vai_GIN = mask(vai_nh,(tlat3d.ge.GIN_lat(0) .and. tlat3d.le.GIN_lat(1) .and.\
                  ((tlon3d.ge.GIN_long1(0) .and. tlon3d.le.GIN_long1(1)) .or.\
                   (tlon3d.ge.GIN_long2(0) .and. tlon3d.le.GIN_long2(1)))), True)
      ext_GIN = mask(ext_nh,(tlat3d.ge.GIN_lat(0) .and. tlat3d.le.GIN_lat(1) .and.\
                  ((tlon3d.ge.GIN_long1(0) .and. tlon3d.le.GIN_long1(1)) .or.\
                   (tlon3d.ge.GIN_long2(0) .and. tlon3d.le.GIN_long2(1)))), True)
      vhi_Bar = mask(vhi_nh,(tlat3d.ge.Bar_lat(0) .and. tlat3d.le.Bar_lat(1) .and.\
                  tlon3d.ge.Bar_long(0) .and. tlon3d.le.Bar_long(1)), True)
      vai_Bar = mask(vai_nh,(tlat3d.ge.Bar_lat(0) .and. tlat3d.le.Bar_lat(1) .and.\
                  tlon3d.ge.Bar_long(0) .and. tlon3d.le.Bar_long(1)), True)
      ext_Bar = mask(ext_nh,(tlat3d.ge.Bar_lat(0) .and. tlat3d.le.Bar_lat(1) .and.\
                  tlon3d.ge.Bar_long(0) .and. tlon3d.le.Bar_long(1)), True)

      vhi_ArcOc = mask(vhi_nh,((tlat3d.ge.ArcOc_lat1(0) .and. (tlon3d.ge.ArcOc_long1(0) .or.\
                  tlon3d.le.ArcOc_long2(1))) .or. (tlat3d.ge.ArcOc_lat3(0) .and.\
                   (tlon3d.ge.ArcOc_long3(0) .and. tlon3d.le.ArcOc_long3(1)))), True)
      vai_ArcOc = mask(vai_nh,((tlat3d.ge.ArcOc_lat1(0) .and. (tlon3d.ge.ArcOc_long1(0) .or.\
                  tlon3d.le.ArcOc_long2(1))) .or. (tlat3d.ge.ArcOc_lat3(0) .and.\
                   (tlon3d.ge.ArcOc_long3(0) .and. tlon3d.le.ArcOc_long3(1)))), True)
      ext_ArcOc = mask(ext_nh,((tlat3d.ge.ArcOc_lat1(0) .and. (tlon3d.ge.ArcOc_long1(0) .or.\
                  tlon3d.le.ArcOc_long2(1))) .or. (tlat3d.ge.ArcOc_lat3(0) .and.\
                   (tlon3d.ge.ArcOc_long3(0) .and. tlon3d.le.ArcOc_long3(1)))), True)

      vhi_Kara = mask(vhi_nh,(tlat3d.ge.Kara_lat(0) .and. tlat3d.le.Kara_lat(1) .and.\
                  tlon3d.ge.Kara_long(0) .and. tlon3d.le.Kara_long(1)), True)
      vai_Kara = mask(vai_nh,(tlat3d.ge.Kara_lat(0) .and. tlat3d.le.Kara_lat(1) .and.\
                  tlon3d.ge.Kara_long(0) .and. tlon3d.le.Kara_long(1)), True)
      ext_Kara = mask(ext_nh,(tlat3d.ge.Kara_lat(0) .and. tlat3d.le.Kara_lat(1) .and.\
                  tlon3d.ge.Kara_long(0) .and. tlon3d.le.Kara_long(1)), True)

      vhi_Laptev = mask(vhi_nh,(tlat3d.ge.Laptev_lat(0) .and. tlat3d.le.Laptev_lat(1) .and.\
                  tlon3d.ge.Laptev_long(0) .and. tlon3d.le.Laptev_long(1)), True)
      vai_Laptev = mask(vai_nh,(tlat3d.ge.Laptev_lat(0) .and. tlat3d.le.Laptev_lat(1) .and.\
                  tlon3d.ge.Laptev_long(0) .and. tlon3d.le.Laptev_long(1)), True)
      ext_Laptev = mask(ext_nh,(tlat3d.ge.Laptev_lat(0) .and. tlat3d.le.Laptev_lat(1) .and.\
                  tlon3d.ge.Laptev_long(0) .and. tlon3d.le.Laptev_long(1)), True)

      vhi_Chukchi = mask(vhi_nh,(tlat3d.ge.Chukchi_lat(0) .and. tlat3d.le.Chukchi_lat(1) .and.\
                  tlon3d.ge.Chukchi_long(0) .and. tlon3d.le.Chukchi_long(1)), True)
      vai_Chukchi = mask(vai_nh,(tlat3d.ge.Chukchi_lat(0) .and. tlat3d.le.Chukchi_lat(1) .and.\
                  tlon3d.ge.Chukchi_long(0) .and. tlon3d.le.Chukchi_long(1)), True)
      ext_Chukchi = mask(ext_nh,(tlat3d.ge.Chukchi_lat(0) .and. tlat3d.le.Chukchi_lat(1) .and.\
                  tlon3d.ge.Chukchi_long(0) .and. tlon3d.le.Chukchi_long(1)), True)

      vhi_Sib = mask(vhi_nh,(tlat3d.ge.Sib_lat(0) .and. tlat3d.le.Sib_lat(1) .and.\
                  tlon3d.ge.Sib_long(0) .and. tlon3d.le.Sib_long(1)), True)
      vai_Sib = mask(vai_nh,(tlat3d.ge.Sib_lat(0) .and. tlat3d.le.Sib_lat(1) .and.\
                  tlon3d.ge.Sib_long(0) .and. tlon3d.le.Sib_long(1)), True)
      ext_Sib = mask(ext_nh,(tlat3d.ge.Sib_lat(0) .and. tlat3d.le.Sib_lat(1) .and.\
                  tlon3d.ge.Sib_long(0) .and. tlon3d.le.Sib_long(1)), True)

      vhi_Beau = mask(vhi_nh,(tlat3d.ge.Beau_lat(0) .and. tlat3d.le.Beau_lat(1) .and.\
                  tlon3d.ge.Beau_long(0) .and. tlon3d.le.Beau_long(1)), True)
      vai_Beau = mask(vai_nh,(tlat3d.ge.Beau_lat(0) .and. tlat3d.le.Beau_lat(1) .and.\
                  tlon3d.ge.Beau_long(0) .and. tlon3d.le.Beau_long(1)), True)
      ext_Beau = mask(ext_nh,(tlat3d.ge.Beau_lat(0) .and. tlat3d.le.Beau_lat(1) .and.\
                  tlon3d.ge.Beau_long(0) .and. tlon3d.le.Beau_long(1)), True)

      vhi_CArc = mask(vhi_nh,(tlat3d.ge.CArc_lat(0) .and. tlat3d.le.CArc_lat(1) .and.\
                  tlon3d.ge.CArc_long(0) .and. tlon3d.le.CArc_long(1)), True)
      vai_CArc = mask(vai_nh,(tlat3d.ge.CArc_lat(0) .and. tlat3d.le.CArc_lat(1) .and.\
                  tlon3d.ge.CArc_long(0) .and. tlon3d.le.CArc_long(1)), True)
      ext_CArc = mask(ext_nh,(tlat3d.ge.CArc_lat(0) .and. tlat3d.le.CArc_lat(1) .and.\
                  tlon3d.ge.CArc_long(0) .and. tlon3d.le.CArc_long(1)), True)

      vhi_Bering = mask(vhi_nh,(Bering_mask .eq. 1.), True)
      vai_Bering = mask(vai_nh,(Bering_mask .eq. 1.), True)
      ext_Bering = mask(ext_nh,(Bering_mask .eq. 1.), True)

      vhi_Okhotsk = mask(vhi_nh,(Okhotsk_mask .eq. 1.), True)
      vai_Okhotsk = mask(vai_nh,(Okhotsk_mask .eq. 1.), True)
      ext_Okhotsk = mask(ext_nh,(Okhotsk_mask .eq. 1.), True)

      vhi_Hudson = mask(vhi_nh,(Hudson_mask .eq. 1.), True)
      vai_Hudson = mask(vai_nh,(Hudson_mask .eq. 1.), True)
      ext_Hudson = mask(ext_nh,(Hudson_mask .eq. 1.), True)

      vhi_CAArch = mask(vhi_nh,(CAArch_mask .eq. 1.), True)
      vai_CAArch = mask(vai_nh,(CAArch_mask .eq. 1.), True)
      ext_CAArch = mask(ext_nh,(CAArch_mask .eq. 1.), True)

      if (snow_exists) then
         vhs_Lab = mask(vhs_nh,(Lab_mask .eq. 1.), True)
         vhs_GIN = mask(vhs_nh,(tlat3d.ge.GIN_lat(0) .and. tlat3d.le.GIN_lat(1) .and.\
                  ((tlon3d.ge.GIN_long1(0) .and. tlon3d.le.GIN_long1(1)) .or.\
                   (tlon3d.ge.GIN_long2(0) .and. tlon3d.le.GIN_long2(1)))), True)
         vhs_Bar = mask(vhs_nh,(tlat3d.ge.Bar_lat(0) .and. tlat3d.le.Bar_lat(1) .and.\
                  tlon3d.ge.Bar_long(0) .and. tlon3d.le.Bar_long(1)), True)
         vhs_ArcOc = mask(vhs_nh,((tlat3d.ge.ArcOc_lat1(0) .and. (tlon3d.ge.ArcOc_long1(0) .or.\
                  tlon3d.le.ArcOc_long2(1))) .or. (tlat3d.ge.ArcOc_lat3(0) .and.\
                   (tlon3d.ge.ArcOc_long3(0) .and. tlon3d.le.ArcOc_long3(1)))), True)
         vhs_Kara = mask(vhs_nh,(tlat3d.ge.Kara_lat(0) .and. tlat3d.le.Kara_lat(1) .and.\
                  tlon3d.ge.Kara_long(0) .and. tlon3d.le.Kara_long(1)), True)
         vhs_Laptev = mask(vhs_nh,(tlat3d.ge.Laptev_lat(0) .and. tlat3d.le.Laptev_lat(1) .and.\
                  tlon3d.ge.Laptev_long(0) .and. tlon3d.le.Laptev_long(1)), True)
         vhs_Chukchi = mask(vhs_nh,(tlat3d.ge.Chukchi_lat(0) .and. tlat3d.le.Chukchi_lat(1) .and.\
                  tlon3d.ge.Chukchi_long(0) .and. tlon3d.le.Chukchi_long(1)), True)

         vhs_Sib = mask(vhs_nh,(tlat3d.ge.Sib_lat(0) .and. tlat3d.le.Sib_lat(1) .and.\
                  tlon3d.ge.Sib_long(0) .and. tlon3d.le.Sib_long(1)), True)
         vhs_Beau = mask(vhs_nh,(tlat3d.ge.Beau_lat(0) .and. tlat3d.le.Beau_lat(1) .and.\
                  tlon3d.ge.Beau_long(0) .and. tlon3d.le.Beau_long(1)), True)
         vhs_CArc = mask(vhs_nh,(tlat3d.ge.CArc_lat(0) .and. tlat3d.le.CArc_lat(1) .and.\
                  tlon3d.ge.CArc_long(0) .and. tlon3d.le.CArc_long(1)), True)
         vhs_Bering = mask(vhs_nh,(Bering_mask .eq. 1.), True)
         vhs_Okhotsk = mask(vhs_nh,(Okhotsk_mask .eq. 1.), True)
         vhs_Hudson = mask(vhs_nh,(Hudson_mask .eq. 1.), True)
         vhs_CAArch = mask(vhs_nh,(CAArch_mask .eq. 1.), True)
      end if

      vhi_Wed = mask(vhi_sh,(tlat3d.le.Wed_lat(1) .and. ((tlon3d.ge.Wed_long1(0) \
                  .and. tlon3d.le.Wed_long1(1)) .or. (tlon3d.ge.Wed_long2(0) \
                   .and. tlon3d.le.Wed_long2(1)))), True)
      vai_Wed = mask(vai_sh,(tlat3d.le.Wed_lat(1) .and.((tlon3d.ge.Wed_long1(0) \
                  .and. tlon3d.le.Wed_long1(1)) .or. (tlon3d.ge.Wed_long2(0) \
                  .and. tlon3d.le.Wed_long2(1)))), True)
      ext_Wed = mask(ext_sh,(tlat3d.le.Wed_lat(1) .and.((tlon3d.ge.Wed_long1(0) \
                  .and. tlon3d.le.Wed_long1(1)) .or. (tlon3d.ge.Wed_long2(0) \
                  .and. tlon3d.le.Wed_long2(1)))), True)

      vhi_Ross = mask(vhi_sh,(tlat3d.le.Ross_lat(1) .and. tlon3d.ge.Ross_long(0) \
		   .and. tlon3d.le.Ross_long(1)), True)
      vai_Ross = mask(vai_sh,(tlat3d.le.Ross_lat(1) .and. tlon3d.ge.Ross_long(0) \
		   .and. tlon3d.le.Ross_long(1)), True)
      ext_Ross = mask(ext_sh,(tlat3d.le.Ross_lat(1) .and. tlon3d.ge.Ross_long(0) \
                   .and. tlon3d.le.Ross_long(1)), True)

      vhi_Ind = mask(vhi_sh,(tlat3d.le.Ind_lat(1) .and. tlon3d.ge.Ind_long(0) \
                   .and. tlon3d.le.Ind_long(1)), True)
      vai_Ind = mask(vai_sh,(tlat3d.le.Ind_lat(1) .and. tlon3d.ge.Ind_long(0) \
                   .and. tlon3d.le.Ind_long(1)), True)
      ext_Ind = mask(ext_sh,(tlat3d.le.Ind_lat(1) .and. tlon3d.ge.Ind_long(0) \
                   .and. tlon3d.le.Ind_long(1)), True)

      vhi_Pac = mask(vhi_sh,(tlat3d.le.Pac_lat(1) .and. tlon3d.ge.Pac_long(0) \
                   .and. tlon3d.le.Pac_long(1)), True)
      vai_Pac = mask(vai_sh,(tlat3d.le.Pac_lat(1) .and. tlon3d.ge.Pac_long(0) \
                   .and. tlon3d.le.Pac_long(1)), True)
      ext_Pac = mask(ext_sh,(tlat3d.le.Pac_lat(1) .and. tlon3d.ge.Pac_long(0) \
                   .and. tlon3d.le.Pac_long(1)), True)

      vhi_BAm = mask(vhi_sh,(tlat3d.le.BAm_lat(1) .and. tlon3d.ge.BAm_long(0) \
                   .and. tlon3d.le.BAm_long(1)), True)
      vai_BAm = mask(vai_sh,(tlat3d.le.BAm_lat(1) .and. tlon3d.ge.BAm_long(0) \
                   .and. tlon3d.le.BAm_long(1)), True)
      ext_BAm = mask(ext_sh,(tlat3d.le.BAm_lat(1) .and. tlon3d.ge.BAm_long(0) \
                   .and. tlon3d.le.BAm_long(1)), True)
      if (snow_exists) then
         vhs_Wed = mask(vhs_sh,(tlat3d.le.Wed_lat(1) .and. ((tlon3d.ge.Wed_long1(0) \
                  .and. tlon3d.le.Wed_long1(1)) .or. (tlon3d.ge.Wed_long2(0) \
                  .and. tlon3d.le.Wed_long2(1)))), True)
         vhs_Ross = mask(vhs_sh,(tlat3d.le.Ross_lat(1) .and. tlon3d.ge.Ross_long(0) \
                   .and. tlon3d.le.Ross_long(1)), True)
         vhs_Ind = mask(vhs_sh,(tlat3d.le.Ind_lat(1) .and. tlon3d.ge.Ind_long(0) \
                   .and. tlon3d.le.Ind_long(1)), True)
         vhs_Pac = mask(vhs_sh,(tlat3d.le.Pac_lat(1) .and. tlon3d.ge.Pac_long(0) \
                   .and. tlon3d.le.Pac_long(1)), True)
         vhs_BAm = mask(vhs_sh,(tlat3d.le.BAm_lat(1) .and. tlon3d.ge.BAm_long(0) \
                   .and. tlon3d.le.BAm_long(1)), True)
      end if
;********************************
; northern hemisphere monthly ...
;********************************

      vhi_mo_nh(imo1:imo2) = dim_sum(dim_sum(vhi_nh))*1.0e-13 
      vai_mo_nh(imo1:imo2) = dim_sum(dim_sum(vai_nh))*1.0e-14 
      ext_mo_nh(imo1:imo2) = dim_sum(dim_sum(ext_nh))*1.0e-12

      vhi_mo_Lab(imo1:imo2) = dim_sum(dim_sum(vhi_Lab))*1.0e-13 
      vai_mo_Lab(imo1:imo2) = dim_sum(dim_sum(vai_Lab))*1.0e-14 
      ext_mo_Lab(imo1:imo2) = dim_sum(dim_sum(ext_Lab))*1.0e-12
 
      vhi_mo_GIN(imo1:imo2) = dim_sum(dim_sum(vhi_GIN))*1.0e-13 
      vai_mo_GIN(imo1:imo2) = dim_sum(dim_sum(vai_GIN))*1.0e-14 
      ext_mo_GIN(imo1:imo2) = dim_sum(dim_sum(ext_GIN))*1.0e-12
 
      vhi_mo_Bar(imo1:imo2) = dim_sum(dim_sum(vhi_Bar))*1.0e-13 
      vai_mo_Bar(imo1:imo2) = dim_sum(dim_sum(vai_Bar))*1.0e-14 
      ext_mo_Bar(imo1:imo2) = dim_sum(dim_sum(ext_Bar))*1.0e-12
 
      vhi_mo_ArcOc(imo1:imo2) = dim_sum(dim_sum(vhi_ArcOc))*1.0e-13 
      vai_mo_ArcOc(imo1:imo2) = dim_sum(dim_sum(vai_ArcOc))*1.0e-14 
      ext_mo_ArcOc(imo1:imo2) = dim_sum(dim_sum(ext_ArcOc))*1.0e-12

      vhi_mo_Kara(imo1:imo2) = dim_sum(dim_sum(vhi_Kara))*1.0e-13
      vai_mo_Kara(imo1:imo2) = dim_sum(dim_sum(vai_Kara))*1.0e-14
      ext_mo_Kara(imo1:imo2) = dim_sum(dim_sum(ext_Kara))*1.0e-12

      vhi_mo_Laptev(imo1:imo2) = dim_sum(dim_sum(vhi_Laptev))*1.0e-13
      vai_mo_Laptev(imo1:imo2) = dim_sum(dim_sum(vai_Laptev))*1.0e-14
      ext_mo_Laptev(imo1:imo2) = dim_sum(dim_sum(ext_Laptev))*1.0e-12

      vhi_mo_Chukchi(imo1:imo2) = dim_sum(dim_sum(vhi_Chukchi))*1.0e-13
      vai_mo_Chukchi(imo1:imo2) = dim_sum(dim_sum(vai_Chukchi))*1.0e-14
      ext_mo_Chukchi(imo1:imo2) = dim_sum(dim_sum(ext_Chukchi))*1.0e-12

      vhi_mo_Sib(imo1:imo2) = dim_sum(dim_sum(vhi_Sib))*1.0e-13 
      vai_mo_Sib(imo1:imo2) = dim_sum(dim_sum(vai_Sib))*1.0e-14 
      ext_mo_Sib(imo1:imo2) = dim_sum(dim_sum(ext_Sib))*1.0e-12

      vhi_mo_Beau(imo1:imo2) = dim_sum(dim_sum(vhi_Beau))*1.0e-13 
      vai_mo_Beau(imo1:imo2) = dim_sum(dim_sum(vai_Beau))*1.0e-14 
      ext_mo_Beau(imo1:imo2) = dim_sum(dim_sum(ext_Beau))*1.0e-12
 
      vhi_mo_CArc(imo1:imo2) = dim_sum(dim_sum(vhi_CArc))*1.0e-13 
      vai_mo_CArc(imo1:imo2) = dim_sum(dim_sum(vai_CArc))*1.0e-14 
      ext_mo_CArc(imo1:imo2) = dim_sum(dim_sum(ext_CArc))*1.0e-12

      vhi_mo_Bering(imo1:imo2) = dim_sum(dim_sum(vhi_Bering))*1.0e-13 
      vai_mo_Bering(imo1:imo2) = dim_sum(dim_sum(vai_Bering))*1.0e-14 
      ext_mo_Bering(imo1:imo2) = dim_sum(dim_sum(ext_Bering))*1.0e-12

      vhi_mo_Okhotsk(imo1:imo2) = dim_sum(dim_sum(vhi_Okhotsk))*1.0e-13 
      vai_mo_Okhotsk(imo1:imo2) = dim_sum(dim_sum(vai_Okhotsk))*1.0e-14 
      ext_mo_Okhotsk(imo1:imo2) = dim_sum(dim_sum(ext_Okhotsk))*1.0e-12

      vhi_mo_Hudson(imo1:imo2) = dim_sum(dim_sum(vhi_Hudson))*1.0e-13 
      vai_mo_Hudson(imo1:imo2) = dim_sum(dim_sum(vai_Hudson))*1.0e-14 
      ext_mo_Hudson(imo1:imo2) = dim_sum(dim_sum(ext_Hudson))*1.0e-12

      vhi_mo_CAArch(imo1:imo2) = dim_sum(dim_sum(vhi_CAArch))*1.0e-13 
      vai_mo_CAArch(imo1:imo2) = dim_sum(dim_sum(vai_CAArch))*1.0e-14 
      ext_mo_CAArch(imo1:imo2) = dim_sum(dim_sum(ext_CAArch))*1.0e-12

      if (snow_exists) then
         vhs_mo_nh(imo1:imo2) = dim_sum(dim_sum(vhs_nh))*1.0e-13
         vhs_mo_Lab(imo1:imo2) = dim_sum(dim_sum(vhs_Lab))*1.0e-13
         vhs_mo_GIN(imo1:imo2) = dim_sum(dim_sum(vhs_GIN))*1.0e-13
         vhs_mo_Bar(imo1:imo2) = dim_sum(dim_sum(vhs_Bar))*1.0e-13
         vhs_mo_ArcOc(imo1:imo2) = dim_sum(dim_sum(vhs_ArcOc))*1.0e-13
         vhs_mo_Kara(imo1:imo2) = dim_sum(dim_sum(vhs_Kara))*1.0e-13
         vhs_mo_Laptev(imo1:imo2) = dim_sum(dim_sum(vhs_Laptev))*1.0e-13
         vhs_mo_Chukchi(imo1:imo2) = dim_sum(dim_sum(vhs_Chukchi))*1.0e-13
         vhs_mo_Sib(imo1:imo2) = dim_sum(dim_sum(vhs_Sib))*1.0e-13
         vhs_mo_Beau(imo1:imo2) = dim_sum(dim_sum(vhs_Beau))*1.0e-13
         vhs_mo_CArc(imo1:imo2) = dim_sum(dim_sum(vhs_CArc))*1.0e-13
         vhs_mo_Bering(imo1:imo2) = dim_sum(dim_sum(vhs_Bering))*1.0e-13
         vhs_mo_Okhotsk(imo1:imo2) = dim_sum(dim_sum(vhs_Okhotsk))*1.0e-13
         vhs_mo_Hudson(imo1:imo2) = dim_sum(dim_sum(vhs_Hudson))*1.0e-13
         vhs_mo_CAArch(imo1:imo2) = dim_sum(dim_sum(vhs_CAArch))*1.0e-13
      end if
;********************************
; southern hemisphere monthly ...
;********************************

      vhi_mo_sh(imo1:imo2) = dim_sum(dim_sum(vhi_sh))*1.0e-13 
      vai_mo_sh(imo1:imo2) = dim_sum(dim_sum(vai_sh))*1.0e-14 
      ext_mo_sh(imo1:imo2) = dim_sum(dim_sum(ext_sh))*1.0e-12

      vhi_mo_Wed(imo1:imo2) = dim_sum(dim_sum(vhi_Wed))*1.0e-13 
      vai_mo_Wed(imo1:imo2) = dim_sum(dim_sum(vai_Wed))*1.0e-14 
      ext_mo_Wed(imo1:imo2) = dim_sum(dim_sum(ext_Wed))*1.0e-12

      vhi_mo_Ross(imo1:imo2) = dim_sum(dim_sum(vhi_Ross))*1.0e-13 
      vai_mo_Ross(imo1:imo2) = dim_sum(dim_sum(vai_Ross))*1.0e-14 
      ext_mo_Ross(imo1:imo2) = dim_sum(dim_sum(ext_Ross))*1.0e-12

      vhi_mo_Ind(imo1:imo2) = dim_sum(dim_sum(vhi_Ind))*1.0e-13 
      vai_mo_Ind(imo1:imo2) = dim_sum(dim_sum(vai_Ind))*1.0e-14 
      ext_mo_Ind(imo1:imo2) = dim_sum(dim_sum(ext_Ind))*1.0e-12

      vhi_mo_Pac(imo1:imo2) = dim_sum(dim_sum(vhi_Pac))*1.0e-13 
      vai_mo_Pac(imo1:imo2) = dim_sum(dim_sum(vai_Pac))*1.0e-14 
      ext_mo_Pac(imo1:imo2) = dim_sum(dim_sum(ext_Pac))*1.0e-12

      vhi_mo_BAm(imo1:imo2) = dim_sum(dim_sum(vhi_BAm))*1.0e-13 
      vai_mo_BAm(imo1:imo2) = dim_sum(dim_sum(vai_BAm))*1.0e-14 
      ext_mo_BAm(imo1:imo2) = dim_sum(dim_sum(ext_BAm))*1.0e-12

      if (snow_exists) then
         vhs_mo_sh(imo1:imo2) = dim_sum(dim_sum(vhs_sh))*1.0e-13
         vhs_mo_Wed(imo1:imo2) = dim_sum(dim_sum(vhs_Wed))*1.0e-13
         vhs_mo_Ross(imo1:imo2) = dim_sum(dim_sum(vhs_Ross))*1.0e-13
         vhs_mo_Ind(imo1:imo2) = dim_sum(dim_sum(vhs_Ind))*1.0e-13
         vhs_mo_Pac(imo1:imo2) = dim_sum(dim_sum(vhs_Pac))*1.0e-13
         vhs_mo_BAm(imo1:imo2) = dim_sum(dim_sum(vhs_BAm))*1.0e-13
      end if

   end do ; nn

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

   vhi_mo_nh!0 = "time"
   vhi_mo_nh&time = time
   vhi_mo_nh@units = "1.E+13 m3"
   vhi_mo_nh@coordinates = "time"
   vhi_mo_nh@missing_value = 1.e+30
   vhi_mo_nh@_FillValue = 1.e+30
   vhi_mo_nh@cell_methods = "time:mean"

   vhs_mo_nh!0 = "time"
   vhs_mo_nh&time = time
   vhs_mo_nh@units = "1.E+13 m3"
   vhs_mo_nh@coordinates = "time"
   vhs_mo_nh@missing_value = 1.e+30
   vhs_mo_nh@_FillValue = 1.e+30
   vhs_mo_nh@cell_methods = "time:mean"

   copy_VarMeta(vai_mo_nh,vai_mo_Lab)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Lab)
   copy_VarMeta(ext_mo_nh,ext_mo_Lab)

   copy_VarMeta(vai_mo_nh,vai_mo_GIN)
   copy_VarMeta(ext_mo_nh,ext_mo_GIN)
   copy_VarMeta(vhi_mo_nh,vhi_mo_GIN)

   copy_VarMeta(vai_mo_nh,vai_mo_Bar)
   copy_VarMeta(ext_mo_nh,ext_mo_Bar)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Bar)

   copy_VarMeta(vai_mo_nh,vai_mo_ArcOc)
   copy_VarMeta(ext_mo_nh,ext_mo_ArcOc)
   copy_VarMeta(vhi_mo_nh,vhi_mo_ArcOc)

   copy_VarMeta(vai_mo_nh,vai_mo_Kara)
   copy_VarMeta(ext_mo_nh,ext_mo_Kara)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Kara)

   copy_VarMeta(vai_mo_nh,vai_mo_Laptev)
   copy_VarMeta(ext_mo_nh,ext_mo_Laptev)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Laptev)

   copy_VarMeta(vai_mo_nh,vai_mo_Chukchi)
   copy_VarMeta(ext_mo_nh,ext_mo_Chukchi)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Chukchi)

   copy_VarMeta(vai_mo_nh,vai_mo_Sib)
   copy_VarMeta(ext_mo_nh,ext_mo_Sib)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Sib)

   copy_VarMeta(vai_mo_nh,vai_mo_Beau)
   copy_VarMeta(ext_mo_nh,ext_mo_Beau)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Beau)

   copy_VarMeta(vai_mo_nh,vai_mo_CArc)
   copy_VarMeta(ext_mo_nh,ext_mo_CArc)
   copy_VarMeta(vhi_mo_nh,vhi_mo_CArc)

   copy_VarMeta(vai_mo_nh,vai_mo_Bering)
   copy_VarMeta(ext_mo_nh,ext_mo_Bering)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Bering)

   copy_VarMeta(vai_mo_nh,vai_mo_Okhotsk)
   copy_VarMeta(ext_mo_nh,ext_mo_Okhotsk)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Okhotsk)

   copy_VarMeta(vai_mo_nh,vai_mo_Hudson)
   copy_VarMeta(ext_mo_nh,ext_mo_Hudson)
   copy_VarMeta(vhi_mo_nh,vhi_mo_Hudson)

   copy_VarMeta(vai_mo_nh,vai_mo_CAArch)
   copy_VarMeta(ext_mo_nh,ext_mo_CAArch)
   copy_VarMeta(vhi_mo_nh,vhi_mo_CAArch)

   if (snow_exists) then
      vhs_mo_nh!0 = "time"
      vhs_mo_nh&time = time
      vhs_mo_nh@units = "1.E+13 m3"
      vhs_mo_nh@coordinates = "time"
      vhs_mo_nh@missing_value = 1.e+30
      vhs_mo_nh@_FillValue = 1.e+30
      vhs_mo_nh@cell_methods = "time:mean"

      copy_VarMeta(vhs_mo_nh,vhs_mo_Lab)
      copy_VarMeta(vhs_mo_nh,vhs_mo_GIN)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Bar)
      copy_VarMeta(vhs_mo_nh,vhs_mo_ArcOc)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Kara)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Laptev)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Chukchi)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Sib)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Beau)
      copy_VarMeta(vhs_mo_nh,vhs_mo_CArc)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Bering)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Okhotsk)
      copy_VarMeta(vhs_mo_nh,vhs_mo_Hudson)
      copy_VarMeta(vhs_mo_nh,vhs_mo_CAArch)
   end if

   vai_mo_sh!0 = "time"
   vai_mo_sh&time = time
   vai_mo_sh@units = "1.E+12 m2"
   vai_mo_sh@coordinates = "time"
   vai_mo_sh@missing_value = 1.e+30
   vai_mo_sh@_FillValue = 1.e+30
   vai_mo_sh@cell_methods = "time:mean"

   ext_mo_sh!0 = "time"
   ext_mo_sh&time = time
   ext_mo_sh@units = "1.E+12 m2"
   ext_mo_sh@coordinates = "time"
   ext_mo_sh@missing_value = 1.e+30
   ext_mo_sh@_FillValue = 1.e+30
   ext_mo_sh@cell_methods = "time:mean"

   vhi_mo_sh!0 = "time"
   vhi_mo_sh&time = time
   vhi_mo_sh@units = "1.E+13 m3"
   vhi_mo_sh@coordinates = "time"
   vhi_mo_sh@missing_value = 1.e+30
   vhi_mo_sh@_FillValue = 1.e+30
   vhi_mo_sh@cell_methods = "time:mean"

   copy_VarMeta(vai_mo_sh,vai_mo_Wed)
   copy_VarMeta(ext_mo_sh,ext_mo_Wed)
   copy_VarMeta(vhi_mo_sh,vhi_mo_Wed)

   copy_VarMeta(vai_mo_sh,vai_mo_Ross)
   copy_VarMeta(ext_mo_sh,ext_mo_Ross)
   copy_VarMeta(vhi_mo_sh,vhi_mo_Ross)

   copy_VarMeta(vai_mo_sh,vai_mo_Ind)
   copy_VarMeta(ext_mo_sh,ext_mo_Ind)
   copy_VarMeta(vhi_mo_sh,vhi_mo_Ind)

   copy_VarMeta(vai_mo_sh,vai_mo_Pac)
   copy_VarMeta(ext_mo_sh,ext_mo_Pac)
   copy_VarMeta(vhi_mo_sh,vhi_mo_Pac)

   copy_VarMeta(vai_mo_sh,vai_mo_BAm)
   copy_VarMeta(ext_mo_sh,ext_mo_BAm)
   copy_VarMeta(vhi_mo_sh,vhi_mo_BAm)

   if (snow_exists) then
      vhs_mo_sh!0 = "time"
      vhs_mo_sh&time = time
      vhs_mo_sh@units = "1.E+13 m3"
      vhs_mo_sh@coordinates = "time"
      vhs_mo_sh@missing_value = 1.e+30
      vhs_mo_sh@_FillValue = 1.e+30
      vhs_mo_sh@cell_methods = "time:mean"

      copy_VarMeta(vhs_mo_sh,vhs_mo_Wed)
      copy_VarMeta(vhs_mo_sh,vhs_mo_Ross)
      copy_VarMeta(vhs_mo_sh,vhs_mo_Ind)
      copy_VarMeta(vhs_mo_sh,vhs_mo_Pac)
      copy_VarMeta(vhs_mo_sh,vhs_mo_BAm)
   end if

   fout = addfile(out_dir+"/ice_vol_MHolland_regions_"+case_names(ii)+"_"+styrs_proc+".nc","c")
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

   filevardef(fout,"vai_mo_nh",typeof(vai_mo_nh),getvardims(vai_mo_nh))
   filevardef(fout,"ext_mo_nh",typeof(ext_mo_nh),getvardims(ext_mo_nh))
   filevardef(fout,"vhi_mo_nh",typeof(vhi_mo_nh),getvardims(vhi_mo_nh))
   filevardef(fout,"vai_mo_sh",typeof(vai_mo_sh),getvardims(vai_mo_sh))
   filevardef(fout,"ext_mo_sh",typeof(ext_mo_sh),getvardims(ext_mo_sh))
   filevardef(fout,"vhi_mo_sh",typeof(vhi_mo_sh),getvardims(vhi_mo_sh))

   filevardef(fout,"vai_mo_Lab",typeof(vai_mo_Lab),getvardims(vai_mo_Lab))
   filevardef(fout,"ext_mo_Lab",typeof(ext_mo_Lab),getvardims(ext_mo_Lab))
   filevardef(fout,"vhi_mo_Lab",typeof(vhi_mo_Lab),getvardims(vhi_mo_Lab))

   filevardef(fout,"vai_mo_GIN",typeof(vai_mo_GIN),getvardims(vai_mo_GIN))
   filevardef(fout,"ext_mo_GIN",typeof(ext_mo_GIN),getvardims(ext_mo_GIN))
   filevardef(fout,"vhi_mo_GIN",typeof(vhi_mo_GIN),getvardims(vhi_mo_GIN))

   filevardef(fout,"vai_mo_Bar",typeof(vai_mo_Bar),getvardims(vai_mo_Bar))
   filevardef(fout,"ext_mo_Bar",typeof(ext_mo_Bar),getvardims(ext_mo_Bar))
   filevardef(fout,"vhi_mo_Bar",typeof(vhi_mo_Bar),getvardims(vhi_mo_Bar))

   filevardef(fout,"vai_mo_ArcOc",typeof(vai_mo_ArcOc),getvardims(vai_mo_ArcOc))
   filevardef(fout,"ext_mo_ArcOc",typeof(ext_mo_ArcOc),getvardims(ext_mo_ArcOc))
   filevardef(fout,"vhi_mo_ArcOc",typeof(vhi_mo_ArcOc),getvardims(vhi_mo_ArcOc))

   filevardef(fout,"vai_mo_Kara",typeof(vai_mo_Kara),getvardims(vai_mo_Kara))
   filevardef(fout,"ext_mo_Kara",typeof(ext_mo_Kara),getvardims(ext_mo_Kara))
   filevardef(fout,"vhi_mo_Kara",typeof(vhi_mo_Kara),getvardims(vhi_mo_Kara))

   filevardef(fout,"vai_mo_Laptev",typeof(vai_mo_Laptev),getvardims(vai_mo_Laptev))
   filevardef(fout,"ext_mo_Laptev",typeof(ext_mo_Laptev),getvardims(ext_mo_Laptev))
   filevardef(fout,"vhi_mo_Laptev",typeof(vhi_mo_Laptev),getvardims(vhi_mo_Laptev))

   filevardef(fout,"vai_mo_Chukchi",typeof(vai_mo_Chukchi),getvardims(vai_mo_Chukchi))
   filevardef(fout,"ext_mo_Chukchi",typeof(ext_mo_Chukchi),getvardims(ext_mo_Chukchi))
   filevardef(fout,"vhi_mo_Chukchi",typeof(vhi_mo_Chukchi),getvardims(vhi_mo_Chukchi))

   filevardef(fout,"vai_mo_Sib",typeof(vai_mo_Sib),getvardims(vai_mo_Sib))
   filevardef(fout,"ext_mo_Sib",typeof(ext_mo_Sib),getvardims(ext_mo_Sib))
   filevardef(fout,"vhi_mo_Sib",typeof(vhi_mo_Sib),getvardims(vhi_mo_Sib))

   filevardef(fout,"vai_mo_Beau",typeof(vai_mo_Beau),getvardims(vai_mo_Beau))
   filevardef(fout,"ext_mo_Beau",typeof(ext_mo_Beau),getvardims(ext_mo_Beau))
   filevardef(fout,"vhi_mo_Beau",typeof(vhi_mo_Beau),getvardims(vhi_mo_Beau))

   filevardef(fout,"vai_mo_CArc",typeof(vai_mo_CArc),getvardims(vai_mo_CArc))
   filevardef(fout,"ext_mo_CArc",typeof(ext_mo_CArc),getvardims(ext_mo_CArc))
   filevardef(fout,"vhi_mo_CArc",typeof(vhi_mo_CArc),getvardims(vhi_mo_CArc))

   filevardef(fout,"vai_mo_Bering",typeof(vai_mo_Bering),getvardims(vai_mo_Bering))
   filevardef(fout,"ext_mo_Bering",typeof(ext_mo_Bering),getvardims(ext_mo_Bering))
   filevardef(fout,"vhi_mo_Bering",typeof(vhi_mo_Bering),getvardims(vhi_mo_Bering))

   filevardef(fout,"vai_mo_Okhotsk",typeof(vai_mo_Okhotsk),getvardims(vai_mo_Okhotsk))
   filevardef(fout,"ext_mo_Okhotsk",typeof(ext_mo_Okhotsk),getvardims(ext_mo_Okhotsk))
   filevardef(fout,"vhi_mo_Okhotsk",typeof(vhi_mo_Okhotsk),getvardims(vhi_mo_Okhotsk))

   filevardef(fout,"vai_mo_Hudson",typeof(vai_mo_Hudson),getvardims(vai_mo_Hudson))
   filevardef(fout,"ext_mo_Hudson",typeof(ext_mo_Hudson),getvardims(ext_mo_Hudson))
   filevardef(fout,"vhi_mo_Hudson",typeof(vhi_mo_Hudson),getvardims(vhi_mo_Hudson))

   filevardef(fout,"vai_mo_CAArch",typeof(vai_mo_CAArch),getvardims(vai_mo_CAArch))
   filevardef(fout,"ext_mo_CAArch",typeof(ext_mo_CAArch),getvardims(ext_mo_CAArch))
   filevardef(fout,"vhi_mo_CAArch",typeof(vhi_mo_CAArch),getvardims(vhi_mo_CAArch))

   if (snow_exists) then
      filevardef(fout,"vhs_mo_nh",typeof(vhs_mo_nh),getvardims(vhs_mo_nh))
      filevardef(fout,"vhs_mo_sh",typeof(vhs_mo_sh),getvardims(vhs_mo_sh))
      filevardef(fout,"vhs_mo_Lab",typeof(vhs_mo_Lab),getvardims(vhs_mo_Lab))
      filevardef(fout,"vhs_mo_GIN",typeof(vhs_mo_GIN),getvardims(vhs_mo_GIN))
      filevardef(fout,"vhs_mo_Bar",typeof(vhs_mo_Bar),getvardims(vhs_mo_Bar))
      filevardef(fout,"vhs_mo_ArcOc",typeof(vhs_mo_ArcOc),getvardims(vhs_mo_ArcOc))
      filevardef(fout,"vhs_mo_Kara",typeof(vhs_mo_Kara),getvardims(vhs_mo_Kara))
      filevardef(fout,"vhs_mo_Laptev",typeof(vhs_mo_Laptev),getvardims(vhs_mo_Laptev))
      filevardef(fout,"vhs_mo_Chukchi",typeof(vhs_mo_Chukchi),getvardims(vhs_mo_Chukchi))
      filevardef(fout,"vhs_mo_Sib",typeof(vhs_mo_Sib),getvardims(vhs_mo_Sib))
      filevardef(fout,"vhs_mo_Beau",typeof(vhs_mo_Beau),getvardims(vhs_mo_Beau))
      filevardef(fout,"vhs_mo_CArc",typeof(vhs_mo_CArc),getvardims(vhs_mo_CArc))
      filevardef(fout,"vhs_mo_Bering",typeof(vhs_mo_Bering),getvardims(vhs_mo_Bering))
      filevardef(fout,"vhs_mo_Okhotsk",typeof(vhs_mo_Okhotsk),getvardims(vhs_mo_Okhotsk))
      filevardef(fout,"vhs_mo_Hudson",typeof(vhs_mo_Hudson),getvardims(vhs_mo_Hudson))
      filevardef(fout,"vhs_mo_CAArch",typeof(vhs_mo_CAArch),getvardims(vhs_mo_CAArch))
   end if

   filevardef(fout,"vai_mo_Wed",typeof(vai_mo_Wed),getvardims(vai_mo_Wed))
   filevardef(fout,"ext_mo_Wed",typeof(ext_mo_Wed),getvardims(ext_mo_Wed))
   filevardef(fout,"vhi_mo_Wed",typeof(vhi_mo_Wed),getvardims(vhi_mo_Wed))

   filevardef(fout,"vai_mo_Ross",typeof(vai_mo_Ross),getvardims(vai_mo_Ross))
   filevardef(fout,"ext_mo_Ross",typeof(ext_mo_Ross),getvardims(ext_mo_Ross))
   filevardef(fout,"vhi_mo_Ross",typeof(vhi_mo_Ross),getvardims(vhi_mo_Ross))

   filevardef(fout,"vai_mo_Ind",typeof(vai_mo_Ind),getvardims(vai_mo_Ind))
   filevardef(fout,"ext_mo_Ind",typeof(ext_mo_Ind),getvardims(ext_mo_Ind))
   filevardef(fout,"vhi_mo_Ind",typeof(vhi_mo_Ind),getvardims(vhi_mo_Ind))

   filevardef(fout,"vai_mo_Pac",typeof(vai_mo_Pac),getvardims(vai_mo_Pac))
   filevardef(fout,"ext_mo_Pac",typeof(ext_mo_Pac),getvardims(ext_mo_Pac))
   filevardef(fout,"vhi_mo_Pac",typeof(vhi_mo_Pac),getvardims(vhi_mo_Pac))

   filevardef(fout,"vai_mo_BAm",typeof(vai_mo_BAm),getvardims(vai_mo_BAm))
   filevardef(fout,"ext_mo_BAm",typeof(ext_mo_BAm),getvardims(ext_mo_BAm))
   filevardef(fout,"vhi_mo_BAm",typeof(vhi_mo_BAm),getvardims(vhi_mo_BAm))

   if (snow_exists) then
      filevardef(fout,"vhs_mo_Wed",typeof(vhs_mo_Wed),getvardims(vhs_mo_Wed))
      filevardef(fout,"vhs_mo_Ross",typeof(vhs_mo_Ross),getvardims(vhs_mo_Ross))
      filevardef(fout,"vhs_mo_Ind",typeof(vhs_mo_Ind),getvardims(vhs_mo_Ind))
      filevardef(fout,"vhs_mo_Pac",typeof(vhs_mo_Pac),getvardims(vhs_mo_Pac))
      filevardef(fout,"vhs_mo_BAm",typeof(vhs_mo_BAm),getvardims(vhs_mo_BAm))
   end if

   filevarattdef(fout,"vhi_mo_nh",vhi_mo_nh)
   filevarattdef(fout,"vai_mo_nh",vai_mo_nh)
   filevarattdef(fout,"ext_mo_nh",ext_mo_nh)
   filevarattdef(fout,"vhi_mo_sh",vhi_mo_sh)
   filevarattdef(fout,"vai_mo_sh",vai_mo_sh)
   filevarattdef(fout,"ext_mo_sh",ext_mo_sh)

   filevarattdef(fout,"vhi_mo_Lab",vhi_mo_Lab)
   filevarattdef(fout,"vai_mo_Lab",vai_mo_Lab)
   filevarattdef(fout,"ext_mo_Lab",vai_mo_Lab)

   filevarattdef(fout,"vhi_mo_GIN",vhi_mo_GIN)
   filevarattdef(fout,"vai_mo_GIN",vai_mo_GIN)
   filevarattdef(fout,"ext_mo_GIN",ext_mo_GIN)

   filevarattdef(fout,"vhi_mo_Bar",vhi_mo_Bar)
   filevarattdef(fout,"vai_mo_Bar",vai_mo_Bar)
   filevarattdef(fout,"ext_mo_Bar",ext_mo_Bar)

   filevarattdef(fout,"vhi_mo_ArcOc",vhi_mo_ArcOc)
   filevarattdef(fout,"vai_mo_ArcOc",vai_mo_ArcOc)
   filevarattdef(fout,"ext_mo_ArcOc",ext_mo_ArcOc)

   filevarattdef(fout,"vhi_mo_Kara",vhi_mo_Kara)
   filevarattdef(fout,"vai_mo_Kara",vai_mo_Kara)
   filevarattdef(fout,"ext_mo_Kara",ext_mo_Kara)

   filevarattdef(fout,"vhi_mo_Laptev",vhi_mo_Laptev)
   filevarattdef(fout,"vai_mo_Laptev",vai_mo_Laptev)
   filevarattdef(fout,"ext_mo_Laptev",ext_mo_Laptev)

   filevarattdef(fout,"vhi_mo_Chukchi",vhi_mo_Chukchi)
   filevarattdef(fout,"vai_mo_Chukchi",vai_mo_Chukchi)
   filevarattdef(fout,"ext_mo_Chukchi",ext_mo_Chukchi)

   filevarattdef(fout,"vhi_mo_Sib",vhi_mo_Sib)
   filevarattdef(fout,"vai_mo_Sib",vai_mo_Sib)
   filevarattdef(fout,"ext_mo_Sib",ext_mo_Sib)

   filevarattdef(fout,"vhi_mo_Beau",vhi_mo_Beau)
   filevarattdef(fout,"vai_mo_Beau",vai_mo_Beau)
   filevarattdef(fout,"ext_mo_Beau",ext_mo_Beau)

   filevarattdef(fout,"vhi_mo_CArc",vhi_mo_CArc)
   filevarattdef(fout,"vai_mo_CArc",vai_mo_CArc)
   filevarattdef(fout,"ext_mo_CArc",ext_mo_CArc)

   filevarattdef(fout,"vhi_mo_Bering",vhi_mo_Bering)
   filevarattdef(fout,"vai_mo_Bering",vai_mo_Bering)
   filevarattdef(fout,"ext_mo_Bering",ext_mo_Bering)
   filevarattdef(fout,"vhi_mo_Okhotsk",vhi_mo_Okhotsk)
   filevarattdef(fout,"vai_mo_Okhotsk",vai_mo_Okhotsk)
   filevarattdef(fout,"ext_mo_Okhotsk",ext_mo_Okhotsk)
   filevarattdef(fout,"vhi_mo_Hudson",vhi_mo_Hudson)
   filevarattdef(fout,"vai_mo_Hudson",vai_mo_Hudson)
   filevarattdef(fout,"ext_mo_Hudson",vai_mo_Hudson)
   filevarattdef(fout,"vhi_mo_CAArch",vhi_mo_CAArch)
   filevarattdef(fout,"vai_mo_CAArch",vai_mo_CAArch)
   filevarattdef(fout,"ext_mo_CAArch",ext_mo_CAArch)

   if (snow_exists) then
      filevarattdef(fout,"vhs_mo_nh",vhs_mo_nh)
      filevarattdef(fout,"vhs_mo_sh",vhs_mo_sh)
      filevarattdef(fout,"vhs_mo_Lab",vhs_mo_Lab)
      filevarattdef(fout,"vhs_mo_GIN",vhs_mo_GIN)
      filevarattdef(fout,"vhs_mo_Bar",vhs_mo_Bar)
      filevarattdef(fout,"vhs_mo_ArcOc",vhs_mo_ArcOc)
      filevarattdef(fout,"vhs_mo_Kara",vhs_mo_Kara)
      filevarattdef(fout,"vhs_mo_Laptev",vhs_mo_Laptev)
      filevarattdef(fout,"vhs_mo_Chukchi",vhs_mo_Chukchi)
      filevarattdef(fout,"vhs_mo_Sib",vhs_mo_Sib)
      filevarattdef(fout,"vhs_mo_Beau",vhs_mo_Beau)
      filevarattdef(fout,"vhs_mo_CArc",vhs_mo_CArc)
      filevarattdef(fout,"vhs_mo_Bering",vhs_mo_Bering)
      filevarattdef(fout,"vhs_mo_Okhotsk",vhs_mo_Okhotsk)
      filevarattdef(fout,"vhs_mo_Hudson",vhs_mo_Hudson)
      filevarattdef(fout,"vhs_mo_CAArch",vhs_mo_CAArch)
   end if

   filevarattdef(fout,"vhi_mo_Wed",vhi_mo_Wed)
   filevarattdef(fout,"vai_mo_Wed",vai_mo_Wed)
   filevarattdef(fout,"ext_mo_Wed",ext_mo_Wed)

   filevarattdef(fout,"vhi_mo_Ross",vhi_mo_Ross)
   filevarattdef(fout,"vai_mo_Ross",vai_mo_Ross)
   filevarattdef(fout,"ext_mo_Ross",ext_mo_Ross)

   filevarattdef(fout,"vhi_mo_Ind",vhi_mo_Ind)
   filevarattdef(fout,"vai_mo_Ind",vai_mo_Ind)
   filevarattdef(fout,"ext_mo_Ind",ext_mo_Ind)

   filevarattdef(fout,"vhi_mo_Pac",vhi_mo_Pac)
   filevarattdef(fout,"vai_mo_Pac",vai_mo_Pac)
   filevarattdef(fout,"ext_mo_Pac",ext_mo_Pac)

   filevarattdef(fout,"vhi_mo_BAm",vhi_mo_BAm)
   filevarattdef(fout,"vai_mo_BAm",vai_mo_BAm)
   filevarattdef(fout,"ext_mo_BAm",ext_mo_BAm)

   if (snow_exists) then
      filevarattdef(fout,"vhs_mo_Wed",vhs_mo_Wed)
      filevarattdef(fout,"vhs_mo_Ross",vhs_mo_Ross)
      filevarattdef(fout,"vhs_mo_Ind",vhs_mo_Ind)
      filevarattdef(fout,"vhs_mo_Pac",vhs_mo_Pac)
      filevarattdef(fout,"vhs_mo_BAm",vhs_mo_BAm)
   end if

   fout->time = time
   fout->vhi_mo_nh = vhi_mo_nh
   fout->vai_mo_nh = vai_mo_nh
   fout->ext_mo_nh = ext_mo_nh
   fout->vhi_mo_sh = vhi_mo_sh
   fout->vai_mo_sh = vai_mo_sh
   fout->ext_mo_sh = ext_mo_sh

   fout->vhi_mo_Lab = vhi_mo_Lab
   fout->vai_mo_Lab = vai_mo_Lab
   fout->ext_mo_Lab = ext_mo_Lab
   fout->vhi_mo_GIN = vhi_mo_GIN
   fout->vai_mo_GIN = vai_mo_GIN
   fout->ext_mo_GIN = ext_mo_GIN
   fout->vhi_mo_Bar = vhi_mo_Bar
   fout->vai_mo_Bar = vai_mo_Bar
   fout->ext_mo_Bar = ext_mo_Bar

   fout->vhi_mo_ArcOc = vhi_mo_ArcOc
   fout->vai_mo_ArcOc = vai_mo_ArcOc
   fout->ext_mo_ArcOc = ext_mo_ArcOc

   fout->vhi_mo_Kara = vhi_mo_Kara
   fout->vai_mo_Kara = vai_mo_Kara
   fout->ext_mo_Kara = ext_mo_Kara

   fout->vhi_mo_Laptev = vhi_mo_Laptev
   fout->vai_mo_Laptev = vai_mo_Laptev
   fout->ext_mo_Laptev = ext_mo_Laptev

   fout->vhi_mo_Chukchi = vhi_mo_Chukchi
   fout->vai_mo_Chukchi = vai_mo_Chukchi
   fout->ext_mo_Chukchi = ext_mo_Chukchi

   fout->vhi_mo_Sib = vhi_mo_Sib
   fout->vai_mo_Sib = vai_mo_Sib
   fout->ext_mo_Sib = ext_mo_Sib
   fout->vhi_mo_Beau = vhi_mo_Beau
   fout->vai_mo_Beau = vai_mo_Beau
   fout->ext_mo_Beau = ext_mo_Beau
   fout->vhi_mo_CArc = vhi_mo_CArc
   fout->vai_mo_CArc = vai_mo_CArc
   fout->ext_mo_CArc = ext_mo_CArc

   fout->vhi_mo_Bering = vhi_mo_Bering
   fout->vai_mo_Bering = vai_mo_Bering
   fout->ext_mo_Bering = ext_mo_Bering
   fout->vhi_mo_Okhotsk = vhi_mo_Okhotsk
   fout->vai_mo_Okhotsk = vai_mo_Okhotsk
   fout->ext_mo_Okhotsk = ext_mo_Okhotsk
   fout->vhi_mo_Hudson = vhi_mo_Hudson
   fout->vai_mo_Hudson = vai_mo_Hudson
   fout->ext_mo_Hudson = ext_mo_Hudson
   fout->vhi_mo_CAArch = vhi_mo_CAArch
   fout->vai_mo_CAArch = vai_mo_CAArch
   fout->ext_mo_CAArch = ext_mo_CAArch

   if (snow_exists) then
      fout->vhs_mo_nh = vhs_mo_nh
      fout->vhs_mo_sh = vhs_mo_sh
      fout->vhs_mo_Lab = vhs_mo_Lab
      fout->vhs_mo_GIN = vhs_mo_GIN
      fout->vhs_mo_Bar = vhs_mo_Bar
      fout->vhs_mo_ArcOc = vhs_mo_ArcOc
      fout->vhs_mo_Kara = vhs_mo_Kara
      fout->vhs_mo_Laptev = vhs_mo_Laptev
      fout->vhs_mo_Chukchi = vhs_mo_Chukchi
      fout->vhs_mo_Sib = vhs_mo_Sib
      fout->vhs_mo_Beau = vhs_mo_Beau
      fout->vhs_mo_CArc = vhs_mo_CArc
      fout->vhs_mo_Bering = vhs_mo_Bering
      fout->vhs_mo_Okhotsk = vhs_mo_Okhotsk
      fout->vhs_mo_Hudson = vhs_mo_Hudson
      fout->vhs_mo_CAArch = vhs_mo_CAArch
   end if

   fout->vhi_mo_Wed = vhi_mo_Wed
   fout->vai_mo_Wed = vai_mo_Wed
   fout->ext_mo_Wed = ext_mo_Wed
   fout->vhi_mo_Ross = vhi_mo_Ross
   fout->vai_mo_Ross = vai_mo_Ross
   fout->ext_mo_Ross = ext_mo_Ross

   fout->vhi_mo_Ind = vhi_mo_Ind
   fout->vai_mo_Ind = vai_mo_Ind
   fout->ext_mo_Ind = ext_mo_Ind
   fout->vhi_mo_Pac = vhi_mo_Pac
   fout->vai_mo_Pac = vai_mo_Pac
   fout->ext_mo_Pac = ext_mo_Pac
   fout->vhi_mo_BAm = vhi_mo_BAm
   fout->vai_mo_BAm = vai_mo_BAm
   fout->ext_mo_BAm = ext_mo_BAm

   if (snow_exists) then
      fout->vhs_mo_Wed = vhs_mo_Wed
      fout->vhs_mo_Ross = vhs_mo_Ross
      fout->vhs_mo_Ind = vhs_mo_Ind
      fout->vhs_mo_Pac = vhs_mo_Pac
      fout->vhs_mo_BAm = vhs_mo_BAm
   end if

end do 						; end do ii
end
