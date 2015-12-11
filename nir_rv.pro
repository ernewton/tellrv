


PRO NIR_RV, order, mydata, hdr, $
	mydata_tc, mystd, mystd_tc, shdr, stdrv=stdrv, $
	atrans=atrans, $
	polydegree=plorder, pixscale=pixscale, $
	spolydegree=s_plorder, spixscale=s_pixscale, $
	trange=trange, wrange=wrange, $
	smshft=s_mshft, mshft=mshft, $ ; linear shift in microns
	shftarr=shft, $
	rv=rv, torest=torest, $
	showplot=showplot, chi=chi, corr_range=corr_range, maxshift=maxshift, ccorr_fxn=ccorr_fxn, $
	contf=contf, frac=frac, sbin=sbin

	data = mydata
	data_tc = mydata_tc
	std = mystd
	std_tc = mystd_tc

	IF ~KEYWORD_SET(atrans) THEN $
		atrans=XMRDFITS('/home/enewton/pro/Spextool/data/atrans.fits',0)

	; parameter for oversampling the spectrum; value set by minimum needed for spex echelle
	oversamp=6.d
	
	; maximum shift to consider; in spex echelle, no shifts beyond 0.0006 (and very few >0.0004)
	if ~KEYWORD_SET(maxshift) THEN $
	        maxshift=0.0008
	maxshft=maxshift/pixscale*oversamp

	; first, get wavelength calibration by modeling the data
	TELL_MODEL, order, atrans, data, hdr, data_new, atrans_new=atrans_new, $
		plorder=plorder, trange=trange, oversamp=oversamp, $
		pixscale=pixscale, maxshft=maxshft, showplot=showplot, $
		res=res, shft=mshft, origcont=origcont

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

	; now do RV

	; first, run wavelength calibration

	IF NOT KEYWORD_SET(std) THEN BEGIN
		std=MRDFITS('/home/enewton/Mdwarf-metallicities/irtf/data/14feb2012/proc/J0727+0513.fits',0,shdr1)
		std_tc=MRDFITS('/home/enewton/Mdwarf-metallicities/irtf/data/14feb2012/proc/J0727+0513_tc.fits',0,shdr)
		std_tc[*,0,*]=std[*,0,*]	; want to use original wavelength array
	ENDIF

	maxshft=maxshift/s_pixscale*oversamp

	; get wavelength calibration for the standaerd by modeling the data

	TELL_MODEL, order, atrans, std, shdr, std_new, atrans_new=atrans_new, $
	plorder=plorder, trange=trange, oversamp=oversamp, $
	pixscale=s_pixscale, maxshft=maxshft, showplot=showplot, $
	res=s_res, shft=s_mshft, origcont=s_origcont

	IF KEYWORD_SET(showplot) THEN BEGIN

		erase & multiplot, /default
		plot, data[*,0], data[*,1]/origcont, xrange=trange, /xsty, yrange=[0,1.5], title='Resulting absolute wavelength calibration (standard)'
		oplot, data_new[*,0], data_new[*,1], co=7
		oplot, atrans_new[*,0], atrans_new[*,1], co=2, linestyle=0
		oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
		al_legend, ['input data','shifted, normalized, interpolated data', 'scaled, interpolated atrans'], color=[1,7,2], linestyle=0, /right, /top
		
		wait, 2

	ENDIF

	; shift telluric corrected data to absolute wavelength
	data_tc_new = data_tc
	data_tc_new[*,0] = data_tc[*,0]+mshft

	; shift telluric corrected standard to absolute wavelength
	std_tc_new = std_tc
	std_tc_new[*,0] = std_tc[*,0]+s_mshft		

	ERN_RV, data_tc_new, std_tc_new, wrange=wrange, pixscale=pixscale, rv0=rv0, showplot=showplot, chi=chi, corr_range=corr_range, ccorr_fxn=ccorr_fxn, contf=contf, frac=frac, sbin=sbin

	rv=DO_RVSHIFTS(rv0, hdr, shdr, bc=bc, rv_std=stdrv)
	torest=rv-bc
; 		print, 'shifts', shft, s_shft
	; shift is now the exact vector you need

	shft = mshft - data_tc[*,0]*rv0/(3.e5)
; 	print, 'RVs', rv, bc, torest

END




