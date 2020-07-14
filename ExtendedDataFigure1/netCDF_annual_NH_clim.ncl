; Reads in monthly ice volume timeseries 
; and creates netCDF file with 10-yr running annual daily/monthly 
; climatogy statistics for ice area, ice extent and ice volume for LE.
;
; LLandrum
; Feb, 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin
;  mname = "GFDL-CM3"
;  dirmname = "gfdl_cm3"
;  ncases = 20
;
;  mname = "MPI-ESM"
;  dirmname = "mpi"
;  ncases = 99
;  mname = "CSIRO-Mk3-6-0"
;  dirmname = "csiro_mk36"
;  ncases = 29
;  mname = "CanESM2"
;  dirmname = "canesm2"
;  ncases = 50
  mname = "GFDL-ESM2M"
  dirmname = "gfdl_esm2m"
  ncases = 30

  case_numbers = ispan(1,ncases,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  case1_names = mname+"_historical_rcp85_r"+case_numbers+"i1p1"
;  case1_names = mname+"_historical_r"+scnumbers+"i1850p3"
  case2_names = mname+"_rcp85_r"+scnumbers+"i2005p3"

  data_time_type = "monthly"			; choices are daily or monthly
  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"

  case1_yr1 = 1950
  case1_yr2 = 2100

  case2_yr1 = 1950
  case2_yr2 = 2100

  nyr_LE = case2_yr2 - case1_yr1 + 1
  year_LE = ispan(case1_yr2,case2_yr2,1)

  nyr_stats = 10
  yr1 = 1950
  yr2 = 2091

  nyrs = yr2 - yr1 + 1
  year = ispan(yr1,yr2,1)
  year!0 = "year"

  fice_vol1 = datadir + "ice_vol_MHolland_regions_"+case1_names+"_"+case1_yr1+"-"+case1_yr2+".nc"
  fice_vol2 = datadir + "ice_vol_MHolland_regions_"+case2_names+"_"+case2_yr1+"-"+case2_yr2+".nc"

  fclim = datadir+"monthly_clim."+nyr_stats+"_yrstats."+mname+"."+ncases+"_sim."+yr1+"-"+yr2+".nc"
  nt = 12
;---------------------------------------------------------------------------------
; read in data
;---------------------------------------------------------------------------------

  ext_LE = new((/nyr_LE*nt,ncases/),float)
  vhi_LE = ext_LE
  vai_LE = ext_LE

  if (data_time_type.eq."daily") then
     ext_var_name = "ext_d_nh"
     vhi_var_name = "vhi_d_nh"
     vai_var_name = "vai_d_nh"
  else
     ext_var_name = "ext_mo_nh"
     vhi_var_name = "vhi_mo_nh"
     vai_var_name = "vai_mo_nh"
  end if

  do i=0,ncases-1
      i2 = (case1_yr2 - yr1 + 1)*nt - 1
      f = addfile(fice_vol1(i),"r")
      ext_LE(:i2,i) = (/ f->$ext_var_name$ /)
      vhi_LE(:i2,i) = (/ f->$vhi_var_name$ /)
      vai_LE(:i2,i) = (/ f->$vai_var_name$ /)

      if (case2_yr1.gt.case1_yr2) then			; data span two files
         i3 = (case2_yr2 - yr1 + 1)*nt - 1
 
         f = addfile(fice_vol2(i),"r")
         ext_LE(i2+1:i3,i) = (/ f->$ext_var_name$ /)
         vhi_LE(i2+1:i3,i) = (/ f->$vhi_var_name$ /)
         vai_LE(i2+1:i3,i) = (/ f->$vai_var_name$ /)
      end if

   end do
   tmpext = f->$ext_var_name$
   tmpvai = f->$vai_var_name$
   tmpvhi = f->$vhi_var_name$

   ext_LE@units = tmpext@units
   vai_LE@units = tmpvai@units
   vhi_LE@units = tmpvhi@units

   delete(tmpext)
   delete(tmpvai)
   delete(tmpvhi)

;---------------------------------------------------------------------------------
; calculate climatologies and statistics
;---------------------------------------------------------------------------------

   ext_stats = new((/nyrs,nt,4/),float)
   ext_stats!0 = "year"
   ext_stats&year = year
   if (data_time_type.eq."daily") then
      ext_stats!1 = "day"
      dimNames = (/"year","day","moment"/)
   else
      ext_stats!1 = "month"
      dimNames = (/"year","month","moment"/)
   end if
   dimSizes = (/nyrs,nt,4/)
   dimUnlim = (/True,False,False/)
 
   ext_stats!2 = "moment"
   moment = ispan(1,4,1)
   ext_stats&moment = moment

   vai_stats = ext_stats
   vhi_stats = ext_stats

   ext_stats@units = ext_LE@units
   vai_stats@units = vai_LE@units
   vhi_stats@units = vhi_LE@units
   do iy=0,nyrs-1
      i1 = iy*nt
      i2 = i1 + (nyr_stats*nt) - 1
 
      tmpext = ext_LE(i1:i2,:)
      tmpvai = vai_LE(i1:i2,:)
      tmpvhi = vhi_LE(i1:i2,:)
 
      do j=0,nt-1
         ext1d = ndtooned(tmpext(j::nt,:))
         vai1d = ndtooned(tmpvai(j::nt,:))
         vhi1d = ndtooned(tmpvhi(j::nt,:))
 
         ext_stats(iy,j,:) = (/ dim_stat4(ext1d) /)
         vai_stats(iy,j,:) = (/ dim_stat4(vai1d) /)
         vhi_stats(iy,j,:) = (/ dim_stat4(vhi1d) /)
      end do
   end do
 
   statLabel = (/"Mean","Var","Skew","Kurtosis"/)
   statLabel!0 = "moment"
   statLabel&moment = moment

   ext_stats@long_name = "4 moments of "+data_time_type+" statistics for "+nyr_stats+" yrs of "+ncases+" LE simulations for NH sea ice extent"
   vai_stats@long_name = "4 moments of "+data_time_type+" statistics for "+nyr_stats+" yrs of "+ncases+" LE simulations for NH sea ice area"
   vhi_stats@long_name = "4 moments of "+data_time_type+" statistics for "+nyr_stats+" yrs of "+ncases+" LE simulations for NH ice volume"

   system("/bin/rm -f "+fclim)
   setfileoption("nc","Format","NetCDF4")
   fout = addfile(fclim,"c")
   setfileoption(fout,"DefineMode",True)

   fileAtt = True
   fileAtt@creation_date = systemfunc("date")
   fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
   fileattdef(fout,fileAtt)

   filedimdef(fout,dimNames,dimSizes,dimUnlim)

   filevardef(fout,"year",typeof(year),getvardims(year))
   filevardef(fout,"ext_stats",typeof(ext_stats),getvardims(ext_stats))
   filevardef(fout,"vai_stats",typeof(vai_stats),getvardims(vai_stats))
   filevardef(fout,"vhi_stats",typeof(vhi_stats),getvardims(vhi_stats))
   filevardef(fout,"statLabel",typeof(statLabel),getvardims(statLabel))

   filevarattdef(fout,"ext_stats",ext_stats)
   filevarattdef(fout,"vai_stats",vai_stats)
   filevarattdef(fout,"vhi_stats",vhi_stats)

   fout->year = (/ year /)
   fout->ext_stats = (/ ext_stats /)
   fout->vai_stats = (/ vai_stats /)
   fout->vhi_stats = (/ vhi_stats /)
   fout->statLabel = (/ statLabel /)

end

