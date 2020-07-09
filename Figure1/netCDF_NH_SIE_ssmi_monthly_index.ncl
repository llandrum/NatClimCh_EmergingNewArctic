; Reads in monthly NH SIE from ssmi data (in excel csv form) and 
; saves to netCDF as well as calculates annual min/max SIE and saves
; that as well.
;
; LLandrum
; April, 2020
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin

  yr1_ssmi = 1978
  yr2_ssmi = 2020
  nyr_ssmi = yr2_ssmi - yr1_ssmi + 1
;  data_dir = "/glade/p/cesm/pcwg/mholland/ssmi_ice_index/NH/csv/"
  data_dir = "/glade/p/cesm/pcwg/landrum/ice_index/"

;---------------------------------------------------------------------------------
; read in data
;---------------------------------------------------------------------------------

  nh_ext = new(12*nyr_ssmi,float)
  nh_area = nh_ext

  index_year = new(12*nyr_ssmi,integer)
  index_month = index_year

  do i=0,11
     imonth = i+1
     if (i.lt.9) then
        fssmi = data_dir+"N_0"+imonth+"_extent_v3.0.csv"
     else
        fssmi = data_dir+"N_"+imonth+"_extent_v3.0.csv"
     end if
     nh_sie_data = asciiread(fssmi,-1,"string")
     year_data = stringtointeger(str_get_cols(nh_sie_data(1:),0,3))
     nyr_data = dimsizes(year_data)
     month = stringtointeger(str_get_cols(nh_sie_data(1:),6,7))
     if (num(month.ne.imonth).ge.1) then
        print("double check data for month "+imonth)
        return
     end if

     i1 = (year_data(0)-yr1_ssmi)*12 + i
     i2 = (year_data(nyr_data-1)-yr1_ssmi+1)*12-1
     nh_ext(i1:i2:12) = stringtofloat(str_get_cols(nh_sie_data(1:),33,37))
     nh_area(i1:i2:12) = stringtofloat(str_get_cols(nh_sie_data(1:),41,45))
     index_year(i1:i2:12) = (/ year_data /)
     index_month(i1:i2:12) = (/ month /)

     delete(nh_sie_data)
     delete(year_data)
     delete(month)
  end do

  nh_ext = where(nh_ext.lt.0.,nh_ext@_FillValue,nh_ext)
  nh_area = where(nh_area.lt.0.,nh_area@_FillValue,nh_area)

  yr1 = yr1_ssmi + 1			; first and last years are not full years of data - so do not use
  yr2 = yr2_ssmi - 1			; first and last years are not full years of data - so do not use
  nyr = yr2 - yr1 + 1
  year = ispan(yr1,yr2,1)

  fmin_max_out = "/glade/p/cesm/pcwg/landrum/SSMI/NH/NH_SIE_min_max_nsidc.from_monthly_sea_ice_index."+yr1+"-"+yr2+".nc"

  max_ext = new(nyr,float)
  max_ext!0 = "year"
  max_ext&year = year
  max_ext@units = "10^6 sq km"
  min_ext = max_ext

  min_month = new(nyr,integer)
  min_month!0 = "year"
  min_month&year = year
  max_month = min_month

  year!0 = "year"
  
  monthly_ext = new(12,float)

  do i=0,nyr-1
     indyr = ind(index_year.eq.year(i))
     tmpext = nh_ext(indyr)
     tmpmonth = index_month(indyr)
     do im=0,11
        indmnth = ind(tmpmonth.eq.(im+1))
        monthly_ext(im) = dim_avg(tmpext(indmnth))
        delete(indmnth)
     end do
     min_ext(i) = min(monthly_ext)
     max_ext(i) = max(monthly_ext)
     indmin = ind(monthly_ext.eq.min_ext(i))
     min_month(i) = indmin+1

     indmax = ind(monthly_ext.eq.max_ext(i))
     max_month(i) = indmax+1
     delete(indyr)
     delete(tmpext)
     delete(tmpmonth)
  end do
;---------------------------------------------------------------------------------
; save netCDF
;---------------------------------------------------------------------------------

  system("/bin/rm -f "+fmin_max_out)
  setfileoption("nc","Format","NetCDF4")
  fout = addfile(fmin_max_out,"c")
  setfileoption(fout,"DefineMode",True)

  fileAtt = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
  fileattdef(fout,fileAtt)

  file_title = "Yearly maximum, minimum NH SIE and month of min/max calculated from SSMI daily data from file: "+fssmi
  fileAtt@title = file_title
 
  dimNames = (/"year"/)
  dimSizes = (/nyr/)
  dimUnlim = (/False/)
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout,"year",typeof(year),getvardims(year))
  filevardef(fout,"min_ext",typeof(min_ext),getvardims(min_ext))
  filevardef(fout,"min_month",typeof(min_month),getvardims(min_month))
  filevardef(fout,"max_ext",typeof(max_ext),getvardims(max_ext))
  filevardef(fout,"max_month",typeof(max_month),getvardims(max_month))

  filevarattdef(fout,"min_ext",min_ext)
  filevarattdef(fout,"min_month",min_month)
  filevarattdef(fout,"max_ext",max_ext)
  filevarattdef(fout,"max_month",max_month)

  fout->year = (/ year /)
  fout->min_ext = (/ min_ext /)
  fout->min_month = (/ min_month /)
  fout->max_ext = (/ max_ext /)
  fout->max_month = (/ max_month /)

end

