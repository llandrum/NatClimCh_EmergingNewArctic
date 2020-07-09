; modified for use with CMIP5 MMLE (only CanESM2 has daily precip and snow)
; CMIP5				LENS
;  pr				PRECT
;  prsn				(PRECSC + PRECSL)
; CMIP5 output is a flux in kg/m2-s. LENS output is in m/s
; For comparison, will need to convert from kg/m2-s by multiplying by density of water
; This depends somewhat on temperature. For now simplify and do not take into account
; temperature and assume want liquid-water-equivalent for snow, so use 1000 kg/m3
;
; Feb 2020
; LLandrum
;
; Calculates and saves yearly day of freeze, thaw onset and earliest
; day of year for when relatively percent of total precip is rain. 
; This is NH!
;
; LLandrum
; April, 2018
;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/glade/u/home/landrum/ncl/asphilli_util.ncl"
; ================================================;
begin
  min_snow_percent = 65.     ; percent of total precip that must be snow for "freeze up onset" definition
  ndays_snow = 7             ; # of days in a row (or precipitation-days in a row) that must exceed min_snow_ratio for freeze onset declaration
  lat_min = 50.              ; looking at Arctic so retrieve only data north of this latitude
  precip_min = 0.2           ; minimum precipitation class (in mm/day) to count towards total days of precip phase

  mname = "CanESM2"
  dirmname = "canesm2"
  ncases = 50

  case_numbers = ispan(1,ncases,1)
  scnumbers = sprinti("%0.3i",case_numbers)

  case_yr1 = 1950
  case_yr2 = 2100

  yr1 = 1950
;  yr2 = 1951
  yr2 = 2100
  nyrs = yr2 - yr1 + 1
  ndays = nyrs*365
  year = ispan(yr1,yr2,1)

  CMIP5dir = "/glade/collections/cdg/data/CLIVAR_LE/"+dirmname+"_lens/day/"
  outdir = "/glade/scratch/landrum/CMIP5/atm/precip/"+mname+"/"

  if (.not.(isfilepresent(outdir))) then
     system("mkdir "+outdir)
  end if

  fdata_pr = CMIP5dir+"pr/pr_day_"+mname+"_historical_rcp85_r"+case_numbers+"i1p1_"+case_yr1+"0101-"+case_yr2+"1231.nc"
  fdata_prsn = CMIP5dir+"prsn/prsn_day_"+mname+"_historical_rcp85_r"+case_numbers+"i1p1_"+case_yr1+"0101-"+case_yr2+"1231.nc"
;  fdata_PRECSC = CMIP5dir+"PRECSC/"+case+".cam.h1.PRECSC."+syr1+"0101-"+syr2+"1231.nc"
;  fdata_PRECSL = CMIP5dir+"PRECSL/"+case+".cam.h1.PRECSL."+syr1+"0101-"+syr2+"1231.nc"
  fgrid = fdata_pr(0)

  ffreeze = outdir+mname+"."+ncases+"sim.NH.freeze_thaw_day.min"+min_snow_percent+"percent."+yr1+"-"+yr2+".nc"

  rho_water = 1000.					; kg/m3
  mps2mm = 60.*60.*24.*1000.				; to convert from m/s to mm for better plotting purposes
  kg2mm = mps2mm/rho_water
;------------------------------------------------------------------------------
; get indices then data
;------------------------------------------------------------------------------
  sim = ispan(0,ncases-1,1)
  case_numbers!0 = "sim"
  case_numbers&sim = sim

  f      = addfile (fgrid, "r")              ; get grid info and create variable arrays and get indices for grabbing only high lat NH
  tmplat = f->lat
  iNH = ind(tmplat.ge.lat_min)
  inh_lat = iNH(0)
  delete(tmplat)
  delete(iNH)

  lon = f->lon
  lat = f->lat(inh_lat:)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  first_rain = new((/nyrs,ncases,nlat,nlon/),integer)
  first_rain!0 = "year"
  first_rain&year = year
  first_rain!1 = "sim"
  first_rain&sim = sim
  first_rain!2 = "lat"
  first_rain&lat = lat
  first_rain!3 = "lon"
  first_rain&lon= lon

  first_snow = first_rain
  thaw_onset = first_rain
  freeze_onset = first_rain
  last_rain = first_rain

  first_rain@long_name = "First day of year with measureable precip. (>"+precip_min+"mm/day) when total rain exceeds "+min_snow_percent+"% of total precip"
  last_rain@long_name = "Last day of year with measureable precip. (>"+precip_min+"mm/day) when total rain exceeds "+min_snow_percent+"% of total precip"
  first_snow@long_name = "First day of year with measureable precip. (>"+precip_min+"mm/day) after first_rain when total snow exceeds "+min_snow_percent+"% of total precip"
  thaw_onset@long_name = "First day of year when total rain exceeds "+min_snow_percent+"% of total precip for "+ndays_snow+" consecutive days"
  freeze_onset@long_name = "First day of year after last "+ndays_snow+"-days rain event when total snow exceeds "+min_snow_percent+"% of total precip for "+ndays_snow+" consecutive days"

  dimNames = (/"year","sim","lat","lon"/)
  dimSizes = (/1,ncases,nlat,nlon/)
  dimUnlim = (/True,False,False,False/)

  total_snow = new((/ndays,nlat,nlon/),float)

  i1 = (yr1 - case_yr1)*365
  i2 = (yr2 - case_yr1 + 1)*365 - 1

  do ic=0,ncases-1

     print("computing rain/snow for "+mname+" "+case_numbers(ic))
     fPR = addfile(fdata_pr(ic),"r")
     fPSN = addfile(fdata_prsn(ic),"r")

     pr = kg2mm*(fPR->pr(i1:i2,inh_lat:,:))
     prsn = kg2mm*(fPSN->prsn(i1:i2,inh_lat:,:))
;     PRECSL = mps2mm*(fPSL->PRECSL(i1:i2,inh_lat:,:))

     total_snow = prsn

     tmp = total_snow
     tmp = (/ where(pr.ne.0.,pr,tmp@_FillValue) /)
     snow_percent = total_snow
     snow_percent = total_snow/tmp
     snow_percent = 100.*snow_percent

     total_rain = total_snow
     total_rain = (/ pr - total_snow /)
     rain_percent = total_rain
     rain_percent = total_rain/tmp
     rain_percent = 100.*rain_percent
     delete(tmp)

     do iyr=0,nyrs-1
        iday1 = 365*iyr
        iday2 = iday1 + 364
        tmp_snow = snow_percent(iday1:iday2,:,:)
        tmp_snow = where(tmp_snow.ge.min_snow_percent,1.,0.)

        tmp_precip = pr(iday1:iday2,:,:)
        tmp_meas_precip = tmp_precip
        tmp_meas_precip = where(tmp_precip.ge.precip_min,1.,0.)

        tmp_rain = rain_percent(iday1:iday2,:,:)
        tmp_rain = where(tmp_rain.ge.min_snow_percent,1.,0.)

        tmp_meas_rain = tmp_rain*tmp_meas_precip
        tmp_meas_snow = tmp_snow*tmp_meas_precip

        do i=0,nlat-1
           do j=0,nlon-1
              if (num(tmp_rain(:,i,j).eq.1.).ge.1) then
                 if (num(tmp_meas_rain(:,i,j).eq.1.).ge.1) then
                    tmpind = ind(tmp_meas_rain(:,i,j).eq.1.)
                    nrain = dimsizes(tmpind)
                    irain = tmpind(0)
                    first_rain(iyr,ic,i,j) = irain + 1				; day of first measureable rain of year
                    last_rain(iyr,ic,i,j) = tmpind(nrain-1) + 1
                    delete(tmpind)

                    if (num(tmp_meas_snow(irain:,i,j).eq.1.).ge.1) then
                       tmpind = ind(tmp_meas_snow(irain:,i,j).eq.1.)
                       first_snow(iyr,ic,i,j) = irain + tmpind(0) + 1                         ; day of first measureable snow (after first measureable rain) of year
                       delete(tmpind)
                    end if
                 end if

                 runavg_rain = runave(tmp_rain(:,i,j),ndays_snow,0)
                 if (num(runavg_rain.eq.1.).ge.1) then
                    tmpind = ind(runavg_rain.eq.1.)
                    ithaw = tmpind(0)
                    nrain = dimsizes(tmpind)
                    ilast_rain = tmpind(nrain-1)
                    thaw_onset(iyr,ic,i,j) = ithaw - (ndays_snow/2) + 1                       ; running average so beginning of ndays_snow of snow starts at

                    runavg_snow = runave(tmp_snow(:,i,j),ndays_snow,0)
                    numsnow = num(runavg_snow(ilast_rain:).eq.1.)
                    delete(tmpind)
                    if (numsnow.ge.1) then
                       tmpind = ind(runavg_snow(ilast_rain:).eq.1.)
                       freeze_onset(iyr,ic,i,j) = ilast_rain + tmpind(0) - (ndays_snow/2) + 1			; running average so beginning of ndays_snow of snow starts at 
                                                                ; ndays_snow/2 indices before running avg = ndays_snow
                       delete(tmpind)
                    end if
                 end if
              end if
           end do							; end nlon
        end do							; end nlat
     end do						; end nyrs
  end do						; end ncases
;------------------------------------------------------------------------------
; save netCDF
;------------------------------------------------------------------------------
  system("/bin/rm -f "+ffreeze)
  setfileoption("nc","Format","NetCDF4")
  fout = addfile(ffreeze,"c")
  setfileoption(fout,"DefineMode",True)

  file_title = "gridded LE first day of rain, day of onset of spring thaw and day of onset of fall freeze for CMIP5 "+mname+" ensemble"
  fileAtt = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@program = systemfunc("pwd")+"/"+get_script_name()
  fileAtt@title = file_title
  fileattdef(fout,fileAtt)

  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout,"lat",typeof(lat),getvardims(lat))
  filevardef(fout,"lon",typeof(lon),getvardims(lon))
  filevardef(fout,"first_rain",typeof(first_rain),getvardims(first_rain))
  filevardef(fout,"last_rain",typeof(last_rain),getvardims(last_rain))
  filevardef(fout,"first_snow",typeof(first_snow),getvardims(first_snow))
  filevardef(fout,"thaw_onset",typeof(thaw_onset),getvardims(thaw_onset))
  filevardef(fout,"freeze_onset",typeof(freeze_onset),getvardims(freeze_onset))
  filevardef(fout,"case_numbers",typeof(case_numbers),getvardims(case_numbers))

  filevarattdef(fout,"lat",lat)
  filevarattdef(fout,"lon",lon)
  filevarattdef(fout,"first_rain",first_rain)
  filevarattdef(fout,"last_rain",last_rain)
  filevarattdef(fout,"first_snow",first_snow)
  filevarattdef(fout,"thaw_onset",thaw_onset)
  filevarattdef(fout,"freeze_onset",freeze_onset)

  fout->lat = (/ lat /)
  fout->lon = (/ lon /)
  fout->first_rain = (/ first_rain /)
  fout->last_rain = (/ last_rain /)
  fout->first_snow = (/ first_snow /)
  fout->thaw_onset = (/ thaw_onset /)
  fout->freeze_onset = (/ freeze_onset /)
  fout->case_numbers = (/ case_numbers /)
  delete(fout)


end

