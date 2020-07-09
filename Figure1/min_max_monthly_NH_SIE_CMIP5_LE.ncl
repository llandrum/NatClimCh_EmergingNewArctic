; Reads in monthly ice extent timeseries and creates
; netCDF file with annual minimum and maximum SIE, and
; months during which these occur, for nh from each of the LE simulations.
;
; LLandrum
; Jan, 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin
;  mname = "MPI-ESM"
;  dirmname = "mpi"
;  ncase = 99
;  mname = "CanESM2"
;  dirmname = "canesm2"
;  ncase = 50
  mname = "GFDL-ESM2M"
  dirmname = "gfdl_esm2m"
  ncase = 30
;
  case_numbers = ispan(1,ncase,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/ice/"+dirmname+"/"

  case_yr1 = 1950
  case_yr2 = 2100

 yr1 = case_yr1
 yr2 = case_yr2

  nyrs = yr2 - yr1 + 1
  year = ispan(yr1,yr2,1)
  year!0 = "year"

  if (mname.eq."MPI-ESM") then
     if (case_yr1.eq.1920) then
        case_names = mname+"_historical_r"+scnumbers+"i1850p3"
     else
        case_names = mname+"_rcp85_r"+scnumbers+"i2005p3"
     end if
  else
     case_names = mname+"_historical_rcp85_r"+case_numbers+"i1p1"
  end if

  fice_vol = datadir + "ice_vol_MHolland_regions_"+case_names+"_sitXsic_"+case_yr1+"-"+case_yr2+".nc"
  fmin_max_out = datadir+"monthly_NH_SIE_min_max_"+mname+"_LE."+ncase+"_sim.sitXsic."+yr1+"-"+yr2+".nc"

  ice_free_min = 1.			; definition of ice-free minimum (in units of 10^6 km^2)
;---------------------------------------------------------------------------------
; calcuate annual min, max SIE and months in which they occur
;---------------------------------------------------------------------------------
  nsim = ispan(0,ncase-1,1)
  case = case_numbers
  case!0= "nsim"
  case&nsim = nsim
  case@long_name = mname+"00+case"

  LE_min_ext = new((/nyrs,ncase/),float)
  LE_min_ext!0 = "year"
  LE_min_ext!1 = "nsim"
  LE_min_ext&year = year
  LE_min_ext&nsim = nsim

  LE_max_ext = LE_min_ext

  LE_min_month= new((/nyrs,ncase/),integer)
  LE_min_month!0 = "year"
  LE_min_month!1 = "nsim"
  LE_min_month&year = year
  LE_min_month&nsim = nsim

  LE_max_month = LE_min_month

  LE_min_nmonths = LE_min_month				; to be able to check if min/max occurs on more than 1 month/yr
  LE_max_nmonths = LE_min_month

  LE_icefree_month = LE_min_month                      ; also want to be able to compute statistics w.r.t. working definition of 
  LE_icefree_nmonths = LE_min_month                    ; "ice-free Arctic" (less than 1x10^6 km^2) so track this as well

  do i=0,ncase-1
     i1 = (yr1 - case_yr1)*12				; only grab data starting yr1
     f = addfile(fice_vol(i),"r")
     ext_mo_nh = f->ext_mo_nh(i1:)
     do iy=0,nyrs-1
        mo1 = iy*12
        mo2 = mo1+11
        tmpext = ext_mo_nh(mo1:mo2)

        if ((i.eq.0).and.(iy.eq.0)) then
           LE_max_ext@units = ext_mo_nh@units
           LE_min_ext@units = ext_mo_nh@units
        end if

        LE_max_ext(iy,i) = max(tmpext)
        ind_max = ind(tmpext.eq.max(tmpext))
        if (dimsizes(ind_max).eq.1) then
           LE_max_month(iy,i) = ind_max
           LE_max_nmonths(iy,i) = 1
        else
           i1 = ind_max(0)
           LE_max_month(iy,i) = i1
           LE_max_nmonths(iy,i) = dimsizes(ind_max)
        end if
        delete(ind_max)

        LE_min_ext(iy,i) = min(tmpext)
        ind_min = ind(tmpext.eq.min(tmpext))
        if (dimsizes(ind_min).eq.1) then
           LE_min_month(iy,i) = ind_min
           LE_min_nmonths(iy,i) = 1
        else
           i1 = ind_min(0)
           LE_min_month(iy,i) = i1
           LE_min_nmonths(iy,i) = dimsizes(ind_min)
        end if
        delete(ind_min)
        if (min(tmpext).lt.ice_free_min) then			; minimum falls below definition of "ice free arctic"        
           ind_min = ind(tmpext.lt.ice_free_min)
           if (dimsizes(ind_min).eq.1) then
              LE_icefree_month(iy,i) = ind_min
              LE_icefree_nmonths(iy,i) = 1
           else
              i1 = ind_min(0)
              LE_icefree_month(iy,i) = i1
              LE_icefree_nmonths(iy,i) = dimsizes(ind_min)
           end if
           delete(ind_min)
        end if
     end do
  end do

  LE_max_ext@units = ext_mo_nh@units
  LE_min_ext@units = ext_mo_nh@units

  system("/bin/rm -f "+fmin_max_out)
  setfileoption("nc","Format","NetCDF4")
  fout = addfile(fmin_max_out,"c")
  setfileoption(fout,"DefineMode",True)

  fileAtt = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
  fileattdef(fout,fileAtt)

  dimNames = (/"year","nsim"/)
  dimSizes = (/nyrs,ncase/)
   dimUnlim = (/False,False/)
   filedimdef(fout,dimNames,dimSizes,dimUnlim)

   filevardef(fout,"year",typeof(year),getvardims(year))
   filevardef(fout,"case",typeof(case),getvardims(case))
   filevardef(fout,"min_ext",typeof(LE_min_ext),getvardims(LE_min_ext))
   filevardef(fout,"min_month",typeof(LE_min_month),getvardims(LE_min_month))
   filevardef(fout,"min_nmonths",typeof(LE_min_nmonths),getvardims(LE_min_nmonths))
   filevardef(fout,"icefree_month",typeof(LE_icefree_month),getvardims(LE_icefree_month))
   filevardef(fout,"icefree_nmonths",typeof(LE_icefree_nmonths),getvardims(LE_icefree_nmonths))
   filevardef(fout,"max_ext",typeof(LE_max_ext),getvardims(LE_max_ext))
   filevardef(fout,"max_month",typeof(LE_max_month),getvardims(LE_max_month))
   filevardef(fout,"max_nmonths",typeof(LE_max_nmonths),getvardims(LE_max_nmonths))

   filevarattdef(fout,"min_ext",LE_min_ext)
   filevarattdef(fout,"min_month",LE_min_month)
   filevarattdef(fout,"min_nmonths",LE_min_nmonths)
   filevarattdef(fout,"icefree_month",LE_icefree_month)
   filevarattdef(fout,"icefree_nmonths",LE_icefree_nmonths)

   filevarattdef(fout,"max_ext",LE_max_ext)
   filevarattdef(fout,"max_month",LE_max_month)
   filevarattdef(fout,"max_nmonths",LE_max_nmonths)

   fout->year = (/ year /)
   fout->case = (/ case /)
   fout->min_ext = (/ LE_min_ext /)
   fout->min_month = (/ LE_min_month /)
   fout->min_nmonths = (/ LE_min_nmonths /)
   fout->icefree_month = (/ LE_icefree_month /)
   fout->icefree_nmonths = (/ LE_icefree_nmonths /)
   fout->max_ext = (/ LE_max_ext /)
   fout->max_month = (/ LE_max_month /)
   fout->max_nmonths = (/ LE_max_nmonths /)

end

