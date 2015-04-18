PRO NIR_RV, order, atrans, data, hdr, $
	data_tc, std, std_tc, shdr, stdrv=stdrv, $
	plorder=plorder, trange=trange, wrange=wrange, $
	shft=shft, s_shft=s_shft, rv=rv, torest=torest, $
	showplot=showplot, chi=chi, corr_range=corr_range

	IF ~KEYWORD_SET(atrans) THEN $
		atrans=XMRDFITS('/home/enewton/pro/Spextool/data/atrans.fits',0)


	; default settings
	IF NOT KEYWORD_SET(order) THEN order=0

	; get pixel scale
	CASE order OF
		0: pixscale = FXPAR(hdr, 'DISPO03')
		1: pixscale = FXPAR(hdr, 'DISPO04')
		2: pixscale = FXPAR(hdr, 'DISPO05')
		3: pixscale = FXPAR(hdr, 'DISPO06')
		4: pixscale = FXPAR(hdr, 'DISPO07')
		ELSE: message, 'Order is', order, 'This is an invalid order'
	ENDCASE
	oversamp=6.d
	maxshft=0.0008/pixscale*oversamp ; in all of the full sample, nothing beyond 0.0006 (and very few >0.0004)

	; default polynomial orders
	IF NOT KEYWORD_SET(plorder) THEN BEGIN
	CASE order OF
		0:	plorder = 5
		1:	plorder = 4
		2:	plorder = 5
		3:	plorder = 5
		4:	plorder = 4
	ENDCASE
	ENDIF

	; default range for doing telluric matching
	IF NOT KEYWORD_SET(trange) THEN BEGIN
	CASE order OF
		0:	trange = [1.995, 2.4]
		1:	trange = [1.43, 1.8]
		2:	trange = [1.142,1.35]
		3:	trange = [0.94, 1.17]
		4:	trange = [0.89, 0.99]
	ENDCASE
	ENDIF

	; default range for radial velocities
	IF NOT KEYWORD_SET(wrange) THEN BEGIN			
	CASE order OF
		0:	wrange = [2.18, 2.41]
		1:	wrange = [1.49, 1.73]
		2:	wrange = [1.15,1.32]
		3:	wrange = [1.0, 1.1]
; 		4:	wrange = [0.82, 1.01]
; 		4:	wrange = [0.95, 1.01]
		4:	wrange = [0.82, 0.93]
	ENDCASE
	ENDIF

	; some regions should be masked out
; 	IF NOT KEYWORD_SET(roi) THEN BEGIN
; 	CASE order OF
; 		0:	
; 		1:	
; 		2:	mask=[[1.2675,1.270]]
; 		3:	
; 		4:	
; 	ENDCASE
; 	ENDIF

	; first, get wavelength calibration by modeling the data
	TELL_MODEL, order, atrans, data, hdr, data_new, atrans_new=atrans_new, roi=roi, $
		plorder=plorder, trange=trange, wrange=wrange, oversamp=oversamp, $
		pixscale=pixscale, maxshft=maxshft, showplot=showplot, $
		res=res, shft=shft, origcont=origcont

	IF KEYWORD_SET(showplot) THEN BEGIN

; 		erase & multiplot, [1,2]
; 		plot, data[*,0], data[*,1]/origcont, yrange=[0,1.2], xrange=trange
; 		oplot, atrans[*,0], atrans[*,1], co=2, linestyle=0
; 		al_legend, ['unshifted, normalized data', 'original atrans'], color=[1,2], linestyle=0, /right, /bottom
; 		multiplot

		erase & multiplot, /default
		plot, data[*,0], data[*,1]/origcont,yrange=[0,1.2], xrange=trange
		oplot, data_new[*,0], data_new[*,1], co=7
		oplot, data[*,0]+shft, data[*,1]/origcont, co=5, linestyle=1
		oplot, atrans_new[*,0], atrans_new[*,1], co=2, linestyle=0
		oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
		al_legend, ['shifted, normalized, interpolated data', 'scaled, interpolated atrans'], color=[1,2], linestyle=0, /right, /bottom

		wait, 2

	ENDIF

	; now do RV
	IF NOT KEYWORD_SET(norv) THEN BEGIN

		; first, run wavelength calibration
; 		IF NOT KEYWORD_SET(s_shft) THEN BEGIN

		IF NOT KEYWORD_SET(std) THEN BEGIN
			std=MRDFITS('/home/enewton/Mdwarf-metallicities/irtf/data/14feb2012/proc/J0727+0513.fits',0,shdr1)
			std_tc=MRDFITS('/home/enewton/Mdwarf-metallicities/irtf/data/14feb2012/proc/J0727+0513_tc.fits',0,shdr)
			std_tc[*,0,*]=std[*,0,*]	; want to use original wavelength array
		ENDIF

			; get pixel scale
			CASE order OF
				0: s_pixscale = FXPAR(shdr, 'DISPO03')
				1: s_pixscale = FXPAR(shdr, 'DISPO04')
				2: s_pixscale = FXPAR(shdr, 'DISPO05')
				3: s_pixscale = FXPAR(shdr, 'DISPO06')
				4: s_pixscale = FXPAR(shdr, 'DISPO07')
				ELSE: message, 'Order is', order, 'This is an invalid order'
			ENDCASE
			maxshft=0.0008/s_pixscale*oversamp
	
			; get wavelength calibration for the standaerd by modeling the data

			TELL_MODEL, order, atrans, std, shdr, std_new, atrans_new=atrans_new, roi=roi, $
			plorder=plorder, trange=trange, wrange=wrange, oversamp=oversamp, $
			pixscale=s_pixscale, maxshft=maxshft, showplot=showplot, $
			res=s_res, shft=s_shft, origcont=s_origcont
	
			IF KEYWORD_SET(showplot) THEN BEGIN
		
				erase & multiplot, /default
				plot, std[*,0], std[*,1]/s_origcont,yrange=[0,1.2], xrange=trange
				oplot, std_new[*,0], std_new[*,1], co=7
				oplot, std[*,0]+s_shft, std[*,1]/s_origcont, co=5, linestyle=1
				oplot, atrans_new[*,0], atrans_new[*,1], co=2, linestyle=0
				oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
				oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
				al_legend, ['shifted, normalized, interpolated data', 'scaled, interpolated atrans'], color=[1,2], linestyle=0, /right, /bottom
				
				wait, 2
		
			ENDIF

; 		ENDIF

		; shift telluric corrected data
		data_tc_new = data_tc
		data_tc_new[*,0] = data_tc[*,0]+shft

		; shift telluric corrected standard
		std_tc_new = std_tc
		std_tc_new[*,0] = std_tc[*,0]+s_shft		

		RV, data_tc_new, std_tc_new, wrange=wrange, pixscale=pixscale, rv0=rv0, showplot=showplot, chi=chi, corr_range=corr_range

		rv=DO_BARY(rv0, hdr, shdr, bc=bc, rv_std=stdrv)
		torest=rv-bc
; 		print, 'shifts', shft, s_shft
		; shift is now the exact vector you need
		print, shft
		print, mean(data_tc[*,0]*rv0/(3.e5))
		shft = shft - data_tc[*,0]*rv0/(3.e5)
		print, 'RVs', rv, bc, torest

	ENDIF
END



;============================================
; PRO TELL_MODEL
; Modify the atmospheric transmission spectrum until it matches the observation to find the necessary wavelength shift


PRO TELL_MODEL, order, atrans, data, hdr, roi=roi, $
	data_new, atrans_new=atrans_new, $
	plorder=plorder, trange=trange, wrange=wrange, maxshft=maxshft, $
	oversamp=oversamp, pixscale=pixscale, $
	res=res, shft=shft, origcont=origcont, $
	showplot=showplot

	; interpolate data and atrans onto new wavelength grid
	; data_interp, atrans_interp are interpolated fluxes
	; wl_vector is interpolated wavelengths
	TELLSPEC_INTERP, data, atrans, lambda_interp, data_interp, atrans_interp, pixscale, oversamp
	roi=WHERE(lambda_interp GT trange[0] AND lambda_interp LT trange[1])

	; initialize MPFIT
	fa = {LAMBDA:lambda_interp, DATA:data_interp, ATRANS:atrans_interp, ROI:roi, $
		PIXSCALE:pixscale, OVERSAMP:oversamp} 
	base={VALUE:1.d, FIXED:0., LIMITED:[0.,0.], LIMITS:[0.,0.]}
	parinfo=REPLICATE(base,plorder+2.)

	; mpfit can get caught in local minima since telluric features are regularly spaced. Real shifts won't be far enough for this to matter, but for testing I need to start at a reasonable distance from the true answer. This is realistic.
	IF KEYWORD_SET(testoffset) THEN $
		 parinfo[0].value=testoffset/pixscale*oversamp+0.0005*RANDOMN(seed)
	parinfo[0].value=0.d
	parinfo[1].value=2.d			; 2 is typical for all but the K band.
	parinfo[0].limited=[1.,1.]		; limit the shift in pixels to...
	parinfo[0].limits=[-maxshft,maxshft]	; ... 0.0015 microns
; 	parinfo[1].limited[0]=0			; lower limit on the scaling
; 	parinfo[1].limits[0]=0.5

	; run MPFIT
	res = MPFIT('tellfuncdos',parinfo=parinfo,functargs=fa, dof=dof, bestnorm=chi2,covar=covar, quiet=1)

; 	parinfo.fixed=1.
; 	parinfo.value=res
; 	parinfo[0].fixed=0.
; 	res = MPFIT('tellfuncdos',parinfo=parinfo,functargs=fa, dof=dof, bestnorm=chi2,covar=covar);, quiet=0)

	; get result
	; res[0] is shift in pixels
	; res[1] is atrans flux scaling
	; remaining are Legendre polynomial coefficients
	diff = TELLFUNCDOS(res, lambda=lambda_interp, atrans=atrans_interp, data=data_interp, roi=roi, model=model, shft=shft, cont=cont, pixscale=pixscale, oversamp=oversamp)

	data_new = [[lambda_interp+shft],[data_interp/cont]]
	atrans_new = [[lambda_interp],[atrans_interp^res[1]]]

	; want to save continuum on original grid as well
	n=N_ELEMENTS(data[*,0])
	pippo=SCALE_VECTOR(FINDGEN(N_ELEMENTS(cont)),0,n-1)
	origcont = INTERPOL(cont, pippo, FINDGEN(n))
; 	plot, lambda_interp, data_interp/cont, yrange=[0,1.2]
; 	oplot, data[*,0], data[*,1]/origcont, co=2
; 	print, parinfo[0].value
; 	wait, 1

	IF KEYWORD_SET(showplot) THEN BEGIN

		print, 'max shift', maxshft
		print, 'shift in pixels', res[0], size(res[0])
		print, 'shift in microns', shft, size(shft)

		erase & multiplot, [1,2]
		plot, lambda_interp, data_interp
		oplot, [trange[0], trange[0]], [-20,5000], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [-20,5000], co=4, linestyle=2
		IF order EQ 4 THEN adj=1. ELSE adj=2.
		oplot, lambda_interp, atrans_interp*cont+adj, co=3, linestyle=2
		oplot, lambda_interp, model, co=2
		al_legend, ['original interpolated data', 'unshifted atrans model', 'shifted atrans model'], color=[1,3,2], linestyle=[0,2,2], /right
		multiplot

		plot, lambda_interp, data_interp/cont, yrange=[0,1.1], xrange=trange
		oplot, [0,3],[1,1], co=4, linestyle=2
		oplot, lambda_interp+shft, data_interp/cont, co=7
		oplot, lambda_interp, (atrans_interp)^res[1], co=2, linestyle=2
		oplot, atrans[*,0], (atrans[*,1])^res[1], co=2, linestyle=1
		oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
		al_legend, ['unshifted, normalized data', 'shifted, normalized data', 'original atrans','original, interpolated atrans'], color=[1,7,2,2], linestyle=[0,0,1,2], /right, /bottom

; 		wait, 2
	ENDIF
	multiplot,/default


END




;============================================
; FUNCTION TELLFUNCDOS
; Function for MPFIT



FUNCTION TELLFUNCDOS, p, lambda=lambda, atrans=atrans, data=data, roi=roi, model=model, cont=cont, pixscale=pixscale, oversamp=oversamp, shft=shft

	IF NOT KEYWORD_SET(roi) THEN roi=FINDGEN(N_ELEMENTS(data))

	; scale atrans by a constant to account for precipital water vapor and airmass differnces
	atrans_new=atrans^(p[1])

	shft = p[0]*pixscale/oversamp
	wl_shift = lambda + shft
	atrans_new=INTERPOL(atrans_new, lambda, wl_shift)

	; p is polynomial coefficients
	; x is vector of pixel positions
	x = SCALE_VECTOR(FINDGEN(N_ELEMENTS(data)), -1, 1) 
	poly=LEGPOLY(x,p[2:*])

	atrans_curved = atrans_new*poly

	IF KEYWORD_SET(mask) THEN $
		diff=((data-atrans_curved)*mask)[roi] $
	ELSE $
		diff=(data-atrans_curved)[roi] 

; 	plot, lambda, data
; 	oplot, lambda, atrans_curved, co=2
; 	wait, .1

	model=atrans_curved
	cont=poly		
	RETURN, diff

END




;============================================
; PRO TELLSPEC_INTERP
; Interpolate data and atrans onto supersampled, uniformly spaced grids



PRO TELLSPEC_INTERP, data, atrans, wl_vector, data_interp, atrans_interp, pixscale, oversamp

	; wavelength range for new data
	start_wl = MIN(data[*,0])
	end_wl = MAX(data[*,0])

	; new oversampled wavelength vector on which to interpolate all data 
	wl_vector = SCALE_VECTOR(FIX(FINDGEN((end_wl-start_wl)*oversamp/pixscale)), start_wl, end_wl)
	
	; interpolate atrans and object flux onto wl_vector
	atrans_interp = INTERPOL(atrans[*,1],atrans[*,0],wl_vector, /spline) 
	data_interp = INTERPOL(data[*,1],data[*,0],wl_vector, /spline)

; plot, data[*,0], data[*,1]
; oplot, wl_vector, data_interp, co=2
; wait,1
; plot, wl_vector, atrans_interp, /nodata
; oplot, atrans[*,0], atrans[*,1]
; oplot, wl_vector, atrans_interp, co=2
; wait,1

END




;============================================
; FUNCTION LEGPOLY
; Legendre polynomial of arbitary order



FUNCTION LEGPOLY, x, p

	poly=0.
	FOR i=0,N_ELEMENTS(p)-1 DO poly=poly + p[i]*legendre(x,i)
	RETURN, poly

END



;============================================
; PRO RV
; Find the offset between data and a standard star



PRO RV, data, std, pixscale=pixscale, wrange=wrange, showplot=showplot, rv0=rv0, chi=chi, corr_range=corr_range

	IF KEYWORD_SET(showplot) THEN showplot=2

	; select wavelength range (log lambda)
	IF KEYWORD_SET(wrange) THEN BEGIN
		start_wl = ALOG(wrange[0])
		end_wl = ALOG(wrange[1])
		roi = WHERE(data[*,0] GT wrange[0] AND data[*,0] LT wrange[1])
	ENDIF ELSE BEGIN
		start_wl = ALOG(MIN(data[*,0]))	
		end_wl = ALOG(MAX(data[*,0]))	
		roi = FINDGEN(N_ELEMENTS(data[*,0]))
	ENDELSE

	; create oversampled grid uniformly spaced in log lambda	
	oversamp=6. 
	wl_vector = SCALE_VECTOR(FINDGEN((end_wl-start_wl)*oversamp*mean(data[roi,0])/pixscale), start_wl, end_wl) 

	; interpolate object and standard onto new grid
	int_obj = INTERPOL(data[*,1],ALOG(data[*,0]),wl_vector, /spline) 
	int_std = INTERPOL(std[*,1],ALOG(std[*,0]),wl_vector, /spline)	

	; remove continuum
	flat_obj=FLATTEN(int_obj, showplot=showplot)
	flat_std=FLATTEN(int_std, showplot=showplot)

	IF KEYWORD_SET(showplot) THEN BEGIN
	
		print, wl_vector[1]-wl_vector[0]
		print, wl_vector[0], wl_vector[1]
		print, start_wl, end_wl
		
		; this plot checks the interpolations
		plot, wl_vector, int_obj
		oplot, wl_vector, int_std-0.4
		oplot, ALOG(data[*,0]), data[*,1], co=3
		oplot, ALOG(data[*,0]), data[*,1]-0.4, co=4
		oplot, [ALOG(2.2063),ALOG(2.2063)],[0,5],linestyle=2
		oplot, [ALOG(2.26267),ALOG(2.26267)],[0,5],linestyle=2
		wait, 2
	
	ENDIF

	IF NOT KEYWORD_SET(corr_range) THEN corr_range=20
	xcorl, flat_std, flat_obj, corr_range, shft, chi, minchi, plot=showplot, print=showplot
	IF KEYWORD_SET(showplot) THEN wait, 2

	; these are the pixel arrays
	pix_fiducial = MAKE_ARRAY(N_ELEMENTS(wl_vector),/index)
	pix_shifted = pix_fiducial + shft
	wl_shifted = INTERPOL(wl_vector, pix_fiducial, pix_shifted)

	; offset and relative RV
	; if shift is positive, shifting object to redder wavelengths = it was blue shifted = negative RV

	offset=-shft*pixscale/(oversamp*mean(data[roi,0]))  	; delta lambda/lambda
	RV0=3.0e5*offset					; velocity
print, offset, rv0, rv0/(3e5)

	IF KEYWORD_SET(showplot) THEN BEGIN
	
		erase & multiplot, [1,3]
		
		plot, pix_fiducial, int_obj/MEAN(int_obj), /xsty
		oplot, pix_fiducial, int_std/MEAN(int_std)-0.4
		oplot, pix_shifted, int_obj/MEAN(int_obj)-0.2, co=3
		oplot, [N_ELEMENTS(wl_vector)/5,N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		oplot, [2*N_ELEMENTS(wl_vector)/5,2*N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		oplot, [3*N_ELEMENTS(wl_vector)/5,3*N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		oplot, [4*N_ELEMENTS(wl_vector)/5,4*N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		
		multiplot
		
		plot, wl_vector, flat_obj, yrange=[0.3,1.2], /xsty
		oplot, wl_vector, flat_std-0.4
		oplot, wl_shifted, flat_obj-0.2, co=3
		oplot, [ALOG(2.2063),ALOG(2.2063)],[0,5],linestyle=2
		oplot, [ALOG(2.26267),ALOG(2.26267)],[0,5],linestyle=2

		multiplot
		
		plot, data[*,0], data[*,1]/MEAN(data[roi,1]), xrange=wrange, /xsty
		oplot, std[*,0], std[*,1]/MEAN(std[WHERE(FINITE(std[*,1])),1])-0.4
		oplot, data[*,0]-offset*data[*,0], data[*,1]/MEAN(data[roi,1])-0.2, co=3
; 		oplot, [2.2063,2.2063],[0,5],linestyle=2
; 		oplot, [2.26267,2.26267],[0,5],linestyle=2
; 		oplot, [2.29,2.29],[0,5],linestyle=2
; 	oplot, [1.1383850,1.1383850], [0,100],color=!col.magenta
; 	oplot, [1.1408517,1.1408517], [0,100],color=!col.magenta
; 	oplot, [1.5028185,1.5028185], [0,100],color=!col.magenta
; 	oplot, [1.1692427, 1.1692427], [0,100],color=!col.magenta ; K+Fe
; 	oplot, [1.1778406, 1.1778406], [0,100],color=!col.magenta ; K
; 	oplot, [1.1801096, 1.1801096], [0,100],color=!col.magenta ; Fe
; 	oplot, [1.1833522,1.1833522], [0,100],color=!col.magenta
; 	oplot, [1.1891897, 1.1891897], [0,100],color=!col.magenta
; 	oplot, [1.1974243, 1.1974243], [0,100],color=!col.magenta
; 	oplot, [1.4878380,1.4878380], [0,100],color=!col.magenta ; Mg
; 	oplot, [1.5167645, 1.5167645], [0,100],color=!col.magenta
; 	oplot, [1.5028185,1.5028185], [0,100],color=!col.magenta ; Mg
; 	oplot, [1.5046375, 1.5046375], [0,100],color=!col.magenta ; Mg
; 	oplot, [1.5771145, 1.5771145], [0,100],color=!col.magenta ; Mg
; 	oplot, [1.5894133, 1.5894133], [0,100],color=!col.magenta ; Si
; 	oplot, [1.6723865, 1.6723865], [0,100],color=!col.magenta ; Al
; 	oplot, [1.6755865, 1.6755865], [0,100],color=!col.magenta ; Al
; 	oplot, [1.7110102, 1.7110102], [0,100],color=!col.magenta
; ; 	oplot, [1.9314969,1.9314969], [0,100],color=!col.magenta ; Ca
; 	oplot, [1.9458160,1.9458160], [0,100],color=!col.magenta ; Ca
; 	oplot, [1.9509361,1.9509361], [0,100],color=!col.magenta ; Ca
; 	oplot, [1.9781980, 1.9781980], [0,100],color=!col.magenta ; Ca
; 	oplot, [1.9817054, 1.9817054], [0,100],color=!col.magenta ; Ca
; ; 	oplot, [1.9865941, 1.9865941], [0,100],color=!col.magenta
; 	oplot, [2.1661178, 2.1661178], [0,100],color=!col.magenta ; Brackett Gamma
; 	oplot, [2.2062996, 2.2062996], [0,100],color=!col.magenta ; Na
; 	oplot, [2.2089866, 2.2089866], [0,100],color=!col.magenta ; Na
; 	oplot, [2.2626701, 2.2626701], [0,100],color=!col.magenta ; Ca
; 	oplot, [2.2656055, 2.2656055], [0,100],color=!col.magenta ; Ca

		print, 'The radial velocity is ', RV0
		
; 		wait, 2
	stop
	ENDIF

END



;============================================
; FUNCTION FLATTEN
; Remove continuum from a spectrum



FUNCTION FLATTEN, int, showplot=showplot


	IF KEYWORD_SET(showplot) THEN showplot=2

	contf, int, c1, nord=4, sbin=10,frac=0.5, plot=showplot, mask=fin, xspline=xspline, yspline=yspline
	t1=spline(xspline,yspline,FINDGEN(N_ELEMENTS(int)))
	
	contf, int/t1, c2, nord=4, sbin=6,frac=0.2, plot=showplot, mask=fin, xspline=xspline, yspline=yspline
	t2=spline(xspline,yspline,FINDGEN(N_ELEMENTS(int)))
	
	flat=int/t1/t2
	IF KEYWORD_SET(showplot) THEN wait,2

	RETURN, flat

END



;============================================
; FUNCTION GET_BARY
; Use hfits fx commands to get info out of header and make barycentric correction



FUNCTION DO_BARY, rv0, hdr, shdr, rv_std=rv_std, bc=hcorr


	; barycentric correction
	hcorr=GET_BARY(hdr)
	
	; offset from standard star
	hcorr_std=GET_BARY(shdr)
	IF NOT KEYWORD_SET(rv_std) THEN rv_std=18.2
	RVoff=rv_std-hcorr_std
	
	; barycentric corrected absolute RV
	RV=RV0+hcorr+RVoff
	
	RETURN, rv

END



;============================================
; FUNCTION GET_BARY
; Use hfits fx commands to get info out of header and make barycentric correction



FUNCTION GET_BARY, hdr, degrees=degrees

	; from date and UT time, get Julian day
	date = FXPAR(hdr, 'DATE_OBS')
	date_ex = STRSPLIT(date,'-',/EXTRACT) ; 0 = year, 1 = month, 2 = day
	time = FXPAR(hdr, 'TIME_OBS')
	time_ex = STRSPLIT(time,':',/EXTRACT) ; 0 = hour, 1 = min, 2 = sec
	jdate = julday(date_ex[1],date_ex[2],date_ex[0],time_ex[0],time_ex[1],time_ex[2])

	; ra and dec in degrees
	ra = FXPAR(hdr, 'RA')
	dec = FXPAR(hdr, 'DEC')
	IF NOT KEYWORD_SET(degrees) THEN BEGIN
		ra_ex = STRSPLIT(ra,'[+:]',/EXTRACT)
		ra_deg = 15 * (ra_ex[0] + ra_ex[1]/60. + ra_ex[2]/3600.)
		dec_ex = STRSPLIT(dec,'[+-:]',/EXTRACT)
		dec_deg = dec[0] + dec_ex[1]/60. + dec_ex[2]/3600.
		; if first character in dec string was negative, make negative
		char = STRMID(dec,0,1)
		IF STRMATCH(char,'\-') EQ 1 THEN dec_deg = -dec_deg
	ENDIF

	; telescope longitude and latitude 
	tel = FXPAR(hdr,'TELESCOP')
	IF STRMATCH(tel,'NASA IRTF') THEN BEGIN
		lat = 19. + 49./60. + 34.38594/3600. ; north
		lon = 360. - (155. + 28./60. + 19.19564/3600.) ; west
		alt = 13674.76*0.3048 ; feet to meters
	ENDIF ELSE message, 'Telescope invalid'
	hcorr = heliocentric(ra_deg,dec_deg,jd=jdate,longitude=lon,latitude=lat, altitude=alt)

	baryvel, jdate, 0, vh, vb
	ra_rad=ra_deg/!RADEG
	dec_rad=dec_deg/!RADEG
	vel = vb[0]*cos(dec_rad)*cos(ra_rad) + $   ;Project velocity toward star
               vb[1]*cos(dec_rad)*sin(ra_rad) + vb[2]*sin(dec_rad) 
; 	print, 'CHECKING', -hcorr, vel

RETURN, -hcorr

END
