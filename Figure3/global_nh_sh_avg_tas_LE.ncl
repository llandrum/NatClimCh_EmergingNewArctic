;------------------------------------------------------
; Calculates and saves average tas timeseries files
; for all simulations of user-specified CMIP5 MMLE.
; Area averages include global, SH (55S-90S), and multiple
; NH averages (north of 45, 50, 55, 60, 65, 70, 75, 80)
;
; LLandrum
; Mar 2020
;------------------------------------------------------
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

  imodel = ind(model_names.eq.mname)
  dirmname = dir_names(imodel)
  ncase = model_nsim(imodel)
  case_yr1 = model_case_yr1(imodel)
  case_yr2 = model_case_yr2(imodel)

  case_numbers = ispan(1,ncase,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  LEdir = "/glade/collections/cdg/data/CLIVAR_LE/"+dirmname+"_lens/Amon/tas/"
  datadir = "/glade/p/cesm/pcwg/landrum/CMIP5/MMLE/atm/"+dirmname+"/"
;  datadir = "/glade/scratch/landrum/CMIP5/atm/"+dirmname+"/"

  yr1_proc = 1920                                                  ; start and end yrs for pre-processed time series
  if (case_yr1.gt.1920) then
     yr1_proc = case_yr1
  end if
  yr2_proc = case_yr2
  nyears = yr2_proc - yr1_proc + 1
  offset = (yr1_proc - case_yr1)*12

  plot_regions = (/"Global","NH","SH"/)
  nregions = dimsizes(plot_regions)

  styrs_proc = sprinti("%0.4i",yr1_proc) + "01-"+sprinti("%0.4i",yr2_proc)+"12"

  case_names = mname+"_historical_rcp85_r"+case_numbers+"i1p1"
  fdata = LEdir+"tas_Amon_"+case_names+"_"+case_yr1+"01-"+case_yr2+"12.nc"
  fts = datadir+"/"+mname+"."+ncase+"sim.tas_wgtavg."+styrs_proc+".nc"

  if (mname.eq."CESM1-CAM5") then
     fdata(0) = LEdir+"tas_Amon_"+case_names(0)+"_185001-"+case_yr2+"12.nc"
     offset0 = (yr1_proc - 1850)*12
     offset1 = offset
     offset = offset0
  end if
;--------------------------calculate and save wgt avg tas---------------------------
  do ii = 0,ncase-1
     print("calculating area totals for "+case_names(ii))
     f = addfile(fdata(ii),"r")
     time = f->time(offset:)
     tas = f->tas(offset:,:,:)
 
     if (ii.eq.0) then
        if (mname.eq."CESM1-CAM5") then
           offset = offset1
        end if
        lat = f->lat
        lon = f->lon
        pi=4.*atan(1.0)
        rad=(pi/180.)
        coswgt=cos(rad*lat)
        coswgt!0 = "lat"
        coswgt&lat= lat

        tas_mo_global = new((/nyears*12,ncase/),float)
        tas_mo_global!0 = "time"
        tas_mo_global&time = time
        tas_mo_global!1 = "case"
        tas_mo_global&case = case_numbers
        tas_mo_global@long_name = tas@long_name
        tas_mo_global@units = tas@units

        tas_mo_nh45 = tas_mo_global
        tas_mo_nh50 = tas_mo_global
        tas_mo_nh55 = tas_mo_global
        tas_mo_nh60 = tas_mo_global
        tas_mo_nh65 = tas_mo_global
        tas_mo_nh70 = tas_mo_global
        tas_mo_nh75 = tas_mo_global
        tas_mo_nh80 = tas_mo_global
        tas_mo_sh55 = tas_mo_global
     end if

     tas_mo_global(:,ii) = wgt_areaave(tas,coswgt,1.0,0)
     tas_mo_sh55(:,ii) = wgt_areaave(tas(time|:,{lat|-90.:-55.},lon|:),coswgt({-90.:-55.}),1.0,0)
     tas_mo_nh45(:,ii) = wgt_areaave(tas(time|:,{lat|45.:90.},lon|:),coswgt({45.:90.}),1.0,0)
     tas_mo_nh50(:,ii) = wgt_areaave(tas(time|:,{lat|50.:90.},lon|:),coswgt({50.:90.}),1.0,0)
     tas_mo_nh55(:,ii) = wgt_areaave(tas(time|:,{lat|55.:90.},lon|:),coswgt({55.:90.}),1.0,0)
     tas_mo_nh60(:,ii) = wgt_areaave(tas(time|:,{lat|60.:90.},lon|:),coswgt({60.:90.}),1.0,0)
     tas_mo_nh65(:,ii) = wgt_areaave(tas(time|:,{lat|65.:90.},lon|:),coswgt({65.:90.}),1.0,0)
     tas_mo_nh70(:,ii) = wgt_areaave(tas(time|:,{lat|70.:90.},lon|:),coswgt({70.:90.}),1.0,0)
     tas_mo_nh75(:,ii) = wgt_areaave(tas(time|:,{lat|75.:90.},lon|:),coswgt({75.:90.}),1.0,0)
     tas_mo_nh80(:,ii) = wgt_areaave(tas(time|:,{lat|80.:90.},lon|:),coswgt({80.:90.}),1.0,0)
  end do
;
;--------------------------save to netCDF---------------------------
;
  fout = addfile(fts,"c")
  setfileoption(fout,"DefineMode",True)

  fileAtt = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
  fileattdef(fout,fileAtt)

  dimNames = (/"time","case"/)
  dimSizes = (/-1,ncase/)
  dimUnlim = (/True,False/)
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout,"time",typeof(time),getvardims(time))

  filevardef(fout,"tas_mo_global",typeof(tas_mo_global),getvardims(tas_mo_global))
  filevardef(fout,"tas_mo_sh55",typeof(tas_mo_sh55),getvardims(tas_mo_sh55))
  filevardef(fout,"tas_mo_nh45",typeof(tas_mo_nh45),getvardims(tas_mo_nh45))
  filevardef(fout,"tas_mo_nh50",typeof(tas_mo_nh50),getvardims(tas_mo_nh50))
  filevardef(fout,"tas_mo_nh55",typeof(tas_mo_nh55),getvardims(tas_mo_nh55))
  filevardef(fout,"tas_mo_nh60",typeof(tas_mo_nh60),getvardims(tas_mo_nh60))
  filevardef(fout,"tas_mo_nh65",typeof(tas_mo_nh65),getvardims(tas_mo_nh65))
  filevardef(fout,"tas_mo_nh70",typeof(tas_mo_nh70),getvardims(tas_mo_nh70))
  filevardef(fout,"tas_mo_nh75",typeof(tas_mo_nh75),getvardims(tas_mo_nh75))
  filevardef(fout,"tas_mo_nh80",typeof(tas_mo_nh80),getvardims(tas_mo_nh80))

  filevarattdef(fout,"time",time)

  filevarattdef(fout,"tas_mo_global",tas_mo_global)
  filevarattdef(fout,"tas_mo_sh55",tas_mo_sh55)
  filevarattdef(fout,"tas_mo_nh45",tas_mo_nh45)
  filevarattdef(fout,"tas_mo_nh50",tas_mo_nh50)
  filevarattdef(fout,"tas_mo_nh55",tas_mo_nh55)
  filevarattdef(fout,"tas_mo_nh60",tas_mo_nh60)
  filevarattdef(fout,"tas_mo_nh65",tas_mo_nh65)
  filevarattdef(fout,"tas_mo_nh70",tas_mo_nh70)
  filevarattdef(fout,"tas_mo_nh75",tas_mo_nh75)
  filevarattdef(fout,"tas_mo_nh80",tas_mo_nh80)

  fout->time  = time
  fout->tas_mo_global = tas_mo_global
  fout->tas_mo_sh55 = tas_mo_sh55
  fout->tas_mo_nh45 = tas_mo_nh45
  fout->tas_mo_nh50 = tas_mo_nh50
  fout->tas_mo_nh55 = tas_mo_nh55
  fout->tas_mo_nh60 = tas_mo_nh60
  fout->tas_mo_nh65 = tas_mo_nh65
  fout->tas_mo_nh70 = tas_mo_nh70
  fout->tas_mo_nh75 = tas_mo_nh75
  fout->tas_mo_nh80 = tas_mo_nh80

end
