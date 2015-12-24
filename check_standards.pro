;============================================
; Test code functionality and usage example

PRO check_standards, ccorr=ccorr, contf=contf, showplot=showplot

  quiet=1

  ; read in atmospheric transmission spectrum (Lord, 1992)
  ; file is big, best to read it in only once
  ; you can find atrans in the spextool libary at Spextool/data/atrans.fits
  atrans=MRDFITS('$SPEX_DIR/data/atrans.fits',0, silent=quiet)

  ; list of standard stars and their RVs
  READCOL, 'standards.txt', stars, sptypes, rvs, skipline=1, format='(A,I,F)', silent=quiet
  spec_dir = 'spec/'
  
  ; RV standard: non-telluric corrected spectrum
  ; absolute wavelength calibration is set by the telluric features
  
  ; RV standard: telluric-corrected spectrum
  ; cross-correlation with standard star is done on final data product
  
; ; METHOD 1: use the standard spectrum as-is
;   file = spec_dir+stars[0]+".fits" 
;   std0 = MRDFITS(file,0,shdr)
;   file_tc = spec_dir+stars[0]+"_tc.fits"
;   std_tc0 = MRDFITS(file_tc, /silent)
;   std_tc0[*,0,*]=std0[*,0,*]	; want to use original wavelength array
;   stdrv = rvs[0]

; ; METHOD 2: use the wavelength-calibrated standard star
;   wlcal = 1
;   atrest = 0
;   stdrv = rvs[0]
;   file_tc = spec_dir+'J0727+0513_wlcal.fits'
;   std_tc0 = MRDFITS(file_tc,0, shdr, /silent)
;   std0 = std_tc0

; ; METHOD 3: use the at-rest standard star created by makemeastandwich
  wlcal = 1
  atrest = 1
  stdrv = rvs[0]
  file_tc = spec_dir+'J0727+0513_rest.fits'
  std_tc0 = MRDFITS(file_tc,0, shdr, silent=quiet)
  std0 = std_tc0

  
  FOR i=0,N_ELEMENTS(stars)-1 DO BEGIN 
    std = std0
    std_tc = std_tc0

    ; read in data for this star
    star = stars[i]
    file = spec_dir+star+".fits" 
    data = MRDFITS(file,0,hdr, silent=quiet)
    file_tc = spec_dir+star+"_tc.fits"
    data_tc = MRDFITS(file_tc, silent=quiet)

    print, "==========="
    print, "Star is: ", star
    print, "File is: ", file

    data_tc[*,0,*]=data[*,0,*] ;  wavelength was modified during reduction

    ; J, H, and K bands independently
    temp = FLTARR(3)
    FOR order=0,2 DO BEGIN
      ; get pertinant variables for SpeX
      order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument="spex"
      ; calculate the RV
      NIR_RV, data_tc[*,*,order],hdr, data[*,*,order], $
	std_tc[*,*,order],shdr, std[*,*,order], $
	wlcal=wlcal, atrest=atrest, stdrv=stdrv, $ ; already wavelength calibrated?
	atrans=atrans, $
	pixscale=pixscale, polydegree=polydegree, $
	spixscale=pixscale, spolydegree=polydegree, $ ; standard is from same set-up 
	wrange=wrange, trange=trange, $
	ccorr=ccorr, contf=contf, $
	showplot=showplot, quiet=quiet, $
	rv = myrv
      if order LT 3 then temp[order] = myrv

    ENDFOR
    
    ; Newton et al. (2014) adopts the median of the first three orders as the RV
    print, "RV is:   ", rvs[i]
    print, "RV NIR:  ", median(temp)
    print, "==========="
   
  ENDFOR

  print, "NOTE: RVs from Newton et al. (2014) include a -2.6 km/s systematic RV correction, not included here."
  
  
END


