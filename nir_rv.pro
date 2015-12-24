;+
; NAME:
;	NIR_RV
; PURPOSE:
;	Determine a star's radial velocity
; EXPLANATION:
;	Performs wavelength calibration and cross-correlation with a standard
; CALLING SEQUENCE:
;      NIR_RV, data_tc, data_hdr, data, std_tc, std_hdr, 
;		[std, /WVCAL, STDRV=, ATRANS=,
;		POLYDEGREE=, PIXSCALE=, SPOLYDEGREE=, SPIXSCALE,
;		TRANGE=, WRANGE=, 
;		CORR_RANGE=, MAXSHIFT=, CCORR=,
;		FRAC=, BIN=, /CONTF,
;		SMSHFT=, MSHFT=, SHFTARR=,
;		RV=, RELRV=, TOREST=,
;		/QUIET, /SHOWPLOT]
;
; INPUTS:
;	data_tc = array containing telluric-corrected science target data, [[wave],[flux]]
;	data_hdr = FITS header for data
;	data = array containing non telluric-corrected science target data
;	std_tc = array containinng telluric-corrected standard star data
;	std_hdr = FITS header for std
;	
; OPTIONAL INPUTS:
;	std = array containing non telluric-corrected standard star data 
;		(required only if wavelength calibration for standard has not been done, see /WVCAL)
;
; OPTIONAL KEYWORD INPUTS:
;	wvcal = flag indicating wavelength calibration has been done for stanadrd [0]
;	atrest = flag indicating RV standard is already at rest velocity
;	rvstd = known RV of standard star [0]
;	atrans = atmospheric transmission spectrum or file name [read from $SPEX_DIR/Spextool/data/atrans.fits]
;	polydegree = Legendre polynomial degree for model on non-telluric corrected science star data [4]
;	pixscale = science pixel scale in microns/pix [estimated from spectrum]
;	spolydegree = same as above but for standard star
;	spixscale = same as above but for standard star
;	trange = wavelength range used for telluric modeling [data[0,0],data[-1,0]]
;	wrange = wavelength range used for cross-correlation [data[0,0],data[-1,0]]
;	corr_range = number of oversampled pixels for cross-correlation [about 20]
;	maxshift = maximum allowable shift in telluric modeling [0.008 microns]
;	ccorr = choice of cross-correlation routine [c_correlate]
;	frac, bin = options for continuum removal
;	contf = flag indicating the contf routine can be used [0]	
;	
; OUTPUTS:
;	smshft = shift in microns required to put standard star to absolute wavelength frame
;	mshft = same as above but for science star
;	shftarr = shift required to put spectrum at rest wavelengths (including wave cal and RV) 
;	rv = the absolute RV of the science star
;	relrv = the RV of the science star relative to the standard
;	torest = the RV required to shift to rest frame (doesn't include wavelength calibration!)
;
; EXAMPLE:
;	ALSO SEE CHECK_STANDARDS.PRO FOR MORE INFORMATION
;
; 	; read in data
; 	data = MRDFITS('spec/J0727+0513.fits',0,hdr)
; 	data_tc  = MRDFITS('spec/J0727+0513_tc.fits')
; 	data_tc[*,0,*]=data[*,0,*] ; wavelength was modified during reduction
; 
; 	std = MRDFITS('spec/J0727+0513.fits',0,shdr)
; 	std_tc = MRDFITS('spec/J0727+0513_tc.fits',0,shdr)
; 	std_tc[*,0,*]=std[*,0,*] ; wavelength was modified during reduction
; 
; 	; get pertinent info
; 	order = 1
; 	order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument="spex"
; 	
; 	; run program
; 	NIR_RV, data_tc[*,*,order], hdr, data[*,*,order], $
; 	  std_tc[*,*,order], shdr, std[*,*,order], stdrv=18.2, $
; 	  pixscale=pixscale, polydegree=polydegree, $
; 	  spixscale=pixscale, spolydegree=polydegree, $
; 	  wrange=wrange, trange=trange, $
; 	  shftarr=shftarr, rv=myrv
; 	print, "RV: ", myrv ; answer should be 18.2(+-0.1)
; 	data_rest = data_tc[*,*,order]
; 	data_rest[*,0] = data_rest + shftarr ; ready to use
;
; METHOD:
;	Uses telluric features to determine absolute wavelength calibration for the science
;	target (and optionally the standard star). Cross-correlates to determine the 
; 	radial velocity relative to the standard star. Calculates absolute radial velocity.
;
; PROCEDURES USED:
;	tellrv
;
;-



PRO NIR_RV, mydata_tc, hdr, mydata, $
	mystd_tc, shdr, mystd, $
	wlcal=wlcal, atrest=atrest, stdrv=stdrv, $
	atrans=atrans, $
	polydegree=plorder, pixscale=pixscale, $
	spolydegree=s_plorder, spixscale=s_pixscale, $
	trange=trange, wrange=wrange, $
	corr_range=corr_range, maxshift=maxshift, ccorr=ccorr, $
	frac=frac, sbin=sbin, contf=contf, $
	smshft=s_mshft, mshft=mshft, $ ; linear shift in microns
	shftarr=shft, $
	rv=rv, relrv=rv0, torest=torest, $
	quiet=quiet, showplot=showplot

	IF ~KEYWORD_SET(quiet) THEN quiet=0
	IF ~KEYWORD_SET(pixscale) THEN pixscale = abs(median(mydata_tc[0:-2,0]-mydata_tc[1:-1,0]))
	IF ~KEYWORD_SET(spixscale) THEN spixscale = abs(median(mystd_tc[0:-2,0]-mystd_tc[1:-1,0]))
	IF ~KEYWORD_SET(plorder) THEN plorder = 4
	IF ~KEYWORD_SET(s_plorder) THEN s_plorder = 4
	
	data = mydata
	data_tc = mydata_tc
	IF KEYWORD_SET(mystd) THEN std = mystd
	std_tc = mystd_tc

	IF ~KEYWORD_SET(atrans) THEN BEGIN
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: Reading atrans from spex directory."
		atrans=MRDFITS('$SPEX_DIR/data/atrans.fits',0, silent=quiet) 
	ENDIF ELSE IF size(atrans,/type) EQ 7 THEN BEGIN; string! 
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: Reading atrans from file name supplied."		
		atrans=MRDFITS(atrans,0, silent=quiet)
	ENDIF ELSE $
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: atrans array supplied."

		
	; parameter for oversampling the spectrum; value set by minimum needed for spex echelle
	oversamp=6.d
	
	; maximum shift to consider; in spex echelle, no shifts beyond 0.0006 (and very few >0.0004)
	if ~KEYWORD_SET(maxshift) THEN $
	        maxshift=0.0008
	if ~KEYWORD_SET(quiet) THEN print, "NIR_RV: max shift in microns is ", maxshift
	
	; 
	; SCIENCE TARGET
	;
	
	; data: get wavelength calibration by modeling the telluric features
	maxshft=maxshift/pixscale*oversamp
	TELL_MODEL, atrans, data, hdr, data_new, atrans_new=atrans_new, $
	  plorder=plorder, trange=trange, oversamp=oversamp, $
	  pixscale=pixscale, maxshft=maxshft, showplot=showplot, $
	  res=res, shft=mshft, origcont=origcont, quiet=quiet
	; shift telluric corrected data to absolute wavelength
	data_tc_new = data_tc
	data_tc_new[*,0] = data_tc[*,0]+mshft
	
	IF KEYWORD_SET(showplot) THEN BEGIN
		erase & multiplot, /default
		plot, data[*,0], data[*,1]/origcont, xrange=trange, /xsty, yrange=[0,1.5], title='Resulting absolute wavelength calibration (data)'
		oplot, data_new[*,0], data_new[*,1], co=7
		oplot, atrans_new[*,0], atrans_new[*,1], co=2, linestyle=0
		oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
		al_legend, ['input data','shifted, normalized, interpolated data', 'scaled, interpolated atrans'], color=[1,7,2], linestyle=0, /right, /top
		wait, 2
	ENDIF

	;
	; STANDARD STAR
	;
	
	; standard: get wavelength calibration by modeling the telluric features
	IF ~KEYWORD_SET(wlcal) THEN BEGIN 
		maxshft=maxshift/s_pixscale*oversamp
		TELL_MODEL,atrans, std, shdr, std_new, atrans_new=atrans_new, $
		  plorder=plorder, trange=trange, oversamp=oversamp, $
		  pixscale=s_pixscale, maxshft=maxshft, showplot=showplot, $
		  res=s_res, shft=s_mshft, origcont=s_origcont, quiet=quiet
		; shift telluric corrected standard to absolute wavelength
		std_tc_new = std_tc
		std_tc_new[*,0] = std_tc[*,0]+s_mshft	
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: Wavecal done for standard."
	ENDIF ELSE BEGIN ; user says this has already been done
		std_tc_new = std_tc
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: wlcal keyword set. No calibration done."
	ENDELSE
	IF KEYWORD_SET(showplot) THEN BEGIN
		erase & multiplot, /default
		plot, data[*,0], data[*,1]/origcont, xrange=trange, /xsty, yrange=[0,1.5], title='Resulting absolute wavelength calibration (standard)'
		oplot, data_new[*,0], data_new[*,1], co=7
		oplot, atrans_new[*,0], atrans_new[*,1], co=2, linestyle=0
		oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
		al_legend, ['input  data','shifted, normalized, interpolated data', 'scaled, interpolated atrans'], color=[1,7,2], linestyle=0, /right, /top
		wait, 2
	ENDIF

	;
	; GET RADIAL VELOCITY
	;

	; now cross-correlate data and standard to get relative radial velocity
	ERN_RV, data_tc_new, std_tc_new, wrange=wrange, pixscale=pixscale, rv0=rv0, showplot=showplot, corr_range=corr_range, ccorr=ccorr, contf=contf, frac=frac, sbin=sbin, quiet=quiet

	IF KEYWORD_SET(atrest) THEN BEGIN
		bc = GET_HELIO(hdr, quiet=quiet) 
		rv = rv0 + bc ; actual RV
		torest = rv0 ; the required correction
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: atrest keyword set. No RV for standard star."
	ENDIF ELSE BEGIN
		rv=DO_RVSHIFTS(rv0, hdr, shdr, bc=bc, rv_std=stdrv, quiet=quiet)
		torest=rv-bc
		IF ~KEYWORD_SET(quiet) THEN print, "NIR_RV: Adjusting for RV of standard star using velocity provided and barycentric velocity."
	ENDELSE
	shft = mshft - (mshft+data_tc[*,0])*torest/(3.e5) ; this is what you add to get to rest wavelengths! Fixed to treat linear shift correctly.
	
END




