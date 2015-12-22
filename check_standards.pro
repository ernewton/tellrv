PRO check_standards, ccorr_fxn=ccorfxn, contf=contf, showplot=showplot

  ; read in atmospheric transmission spectrum (Lord, 1992)
  ; file is big, best to read it in only once
  ; you can find atrans in the spextool libary at Spextool/data/atrans.fits
  atrans=MRDFITS('/home/enewton/pro/Spextool/data/atrans.fits',0)

  ; list of standard stars and their RVs
  READCOL, 'standards.txt', stars, sptypes, rvs, skipline=1, format='(A,I,F)'
  spec_dir = 'spec/'
  
  ; RV standard: non-telluric corrected spectrum
  ; absolute wavelength calibration is set by the telluric features
  file = spec_dir+stars[0]+".fits" 
  std0 = MRDFITS(file,0,shdr)
  
  ; RV standard: telluric-corrected spectrum
  ; cross-correlation with standard star is done on final data product
; METHOD 1: use the standard spectrum as-is
;   file_tc = spec_dir+stars[0]+"_tc.fits"
;   std_tc0 = MRDFITS(file_tc, /silent)
;   std_tc0[*,0,*]=std0[*,0,*]	; want to use original wavelength array
;   stdrv = rvs[0]
;   print, "NOTE: RVs from Newton et al. (2014) include a -2.6 km/s systematic RV correction, not included here."

; METHOD 2 (preferred): use the at-rest zero-velocity standard created by makemeastandwich
; (seems to mitigate the random error that gives rise to the systematic correction)
  wlcal = 1
  stdrv = 0.
  file_tc = spec_dir+'J0727+0513_rest.fits'
  std_tc0 = MRDFITS(file_tc,0, shdr, /silent)
  std0 = std_tc0

  FOR i=0,N_ELEMENTS(stars)-1 DO BEGIN 
    std = std0
    std_tc = std_tc0

    ; read in data for this star
    star = stars[i]
    file = spec_dir+star+".fits" 
    data = MRDFITS(file,0,hdr, /silent)
    file_tc = spec_dir+star+"_tc.fits"
    data_tc = MRDFITS(file_tc, /silent)

    print, "==========="
    print, "Star is: ", star
    print, "File is: ", file

    data_tc[*,0,*]=data[*,0,*] ;  wavelength was modified during reduction

    temp = FLTARR(3)
    FOR order=0,2 DO BEGIN
      ; get pertinant variables for SpeX
      order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument="spex"
      ; calculate the RV
      NIR_RV, data_tc[*,*,order],hdr, data[*,*,order], $
	std_tc[*,*,order],shdr, std[*,*,order], $
	wlcal=wlcal, stdrv=stdrv, $ ; already wavelength calibrated?
	atrans=atrans, $
	pixscale=pixscale, polydegree=polydegree, $
	spixscale=pixscale, spolydegree=polydegree, $ ; standard is from same set-up 
	wrange=wrange, trange=trange, $
	mshft=mshft, smshft=smshft, rv=myrv, $
	showplot=showplot, torest=rv0, $
	ccorr_fxn=ccor_fxn, $
	quiet=1
	
      if order LT 3 then temp[order] = myrv

    ENDFOR
    
    ; Newton et al. (2014) adopts the median of the first three orders as the RV
    print, "RV is:   ", rvs[i]
    print, "RV NIR:  ", median(temp)
    print, "==========="
   
  ENDFOR
  
  
END


