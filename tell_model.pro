;============================================
; FUNCTION LEGPOLY
; Legendre polynomial of arbitary order



FUNCTION LEGPOLY, x, p

	poly=0.
	FOR i=0,N_ELEMENTS(p)-1 DO poly=poly + p[i]*legendre(x,i)
	RETURN, poly

END



;============================================
; FUNCTION TELL_FUNC
; Function for MPFIT


FUNCTION TELL_FUNC, p, lambda=lambda, atrans=atrans, data=data, model=model, cont=cont, pixscale=pixscale, oversamp=oversamp, shft=shft, mask=mask

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

	diff=(data-atrans_curved)
	IF KEYWORD_SET(mask) THEN $
		diff=(diff*mask)		

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



PRO TELLSPEC_INTERP, data, atrans, wl_vector, data_interp, atrans_interp, pixscale, oversamp, trange=trange

	; wavelength range for new data
	; NOTE! This differs slightly from what was used in Newton et al.
	;   (2014,2015), for which the region of interest was selected in 
	;   'tell_func'
	IF N_ELEMENTS(trange) NE 2 THEN BEGIN
		start_wl = MIN(data[*,0])
		end_wl = MAX(data[*,0])
	ENDIF ELSE BEGIN
		start_wl = MIN(trange)
		end_wl = MAX(trange)
	ENDELSE

	; new oversampled wavelength vector on which to interpolate all data 
	wl_vector = SCALE_VECTOR(FIX(FINDGEN((end_wl-start_wl)*oversamp/pixscale)), start_wl, end_wl)
	
	; interpolate atrans and object flux onto wl_vector
	roi = WHERE(atrans[*,0] GE start_wl-0.1 AND atrans[*,0] LT end_wl+0.1 AND FINITE(atrans[*,1]))
	atrans_interp = INTERPOL(atrans[roi,1],atrans[roi,0],wl_vector, /spline) 

	roi = WHERE(data[*,0] GE start_wl-0.1 AND data[*,0] LT end_wl+0.1 AND FINITE(data[*,1]) AND data[*,1] GT 0)
	data_interp = INTERPOL(data[roi,1],data[roi,0],wl_vector, /spline)

; plot, data[*,0], data[*,1]
; oplot, wl_vector, data_interp, co=2
; wait,1
; plot, wl_vector, atrans_interp, /nodata
; oplot, atrans[*,0], atrans[*,1]
; oplot, wl_vector, atrans_interp, co=2
; wait,1

END








;============================================
; PRO TELL_MODEL
; Modify the atmospheric transmission spectrum until it matches the observation to find the necessary wavelength shift


PRO TELL_MODEL, atrans, data, $
	data_new, atrans_new=atrans_new, $
	plorder=plorder, trange=trange, maxshft=maxshft, $
	oversamp=oversamp, pixscale=pixscale, $
	res=res, shft=shft, origcont=origcont, $
	showplot=showplot, quiet=quiet
	
	IF ~KEYWORD_SET(quiet) THEN quiet=0
	IF ~KEYWORD_SET(trange) THEN trange=[data[0,0],data[-1,0]]
	IF ~KEYWORD_SET(pixscale) THEN pixscale = MEAN(data[1:-1,0]-data[0:-2,0])
	IF ~KEYWORD_SET(oversamp) THEN oversamp = 1.
	IF ~KEYWORD_SET(plorder) THEN plorder = 5
	IF ~KEYWORD_SET(maxshft) THEN maxshft = pixscale*5.

	; interpolate data and atrans onto new wavelength grid
	; data_interp, atrans_interp are interpolated fluxes
	; wl_vector is interpolated wavelengths
	TELLSPEC_INTERP, data, atrans, lambda_interp, data_interp, atrans_interp, pixscale, oversamp, trange=trange

	; initialize MPFIT
	fa = {LAMBDA:lambda_interp, DATA:data_interp, ATRANS:atrans_interp, $
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
; not sure whether to include
	parinfo[1].limited=[1.,0.]			; lower limit on the scaling
	parinfo[1].limits=[0.5,10.]

	; run MPFIT
	res = MPFIT('tell_func',parinfo=parinfo,functargs=fa, dof=dof, bestnorm=chi2,covar=covar, quiet=1)

	; get result
	; res[0] is shift in pixels
	; res[1] is atrans flux scaling
	; remaining are Legendre polynomial coefficients
	diff = TELL_FUNC(res, lambda=lambda_interp, atrans=atrans_interp, data=data_interp, model=model, shft=shft, cont=cont, pixscale=pixscale, oversamp=oversamp)

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
		plot, lambda_interp, data_interp, /xsty, xrange=trange, /ynozero
		oplot, [trange[0], trange[0]], [-20,5000], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [-20,5000], co=4, linestyle=2
; 		IF order EQ 4 THEN adj=1. ELSE adj=2.
		oplot, lambda_interp, atrans_interp*cont+2., co=3, linestyle=2
		oplot, lambda_interp, model, co=2
		al_legend, ['original interpolated data', 'unshifted atrans model', 'shifted atrans model'], color=[1,3,2], linestyle=[0,2,2], /right
		multiplot

		plot, lambda_interp, data_interp/cont, xrange=trange, /xsty, /ynozero
		oplot, [0,3],[1,1], co=4, linestyle=2
		oplot, lambda_interp+shft, data_interp/cont, co=7
		oplot, lambda_interp, (atrans_interp)^res[1], co=2, linestyle=2
		oplot, atrans[*,0], (atrans[*,1])^res[1], co=2, linestyle=1
		oplot, [trange[0], trange[0]], [0,2], co=4, linestyle=2
		oplot, [trange[1], trange[1]], [0,2], co=4, linestyle=2
		al_legend, ['unshifted, normalized data', 'shifted, normalized data', 'original atrans','original, interpolated atrans'], color=[1,7,2,2], linestyle=[0,0,1,2], /right, /bottom

		multiplot,/default
		wait, 2

	ENDIF


END




