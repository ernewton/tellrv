FUNCTION SPLINE_CONT, int, frac=frac, sbin=sbin, showplot=showplot

	IF ~KEYWORD_SET(frac) THEN frac1=0.2
	IF ~KEYWORD_SET(sbin) THEN sbin1=10 ELSE sbin1=sbin
	ordrd = SORT(int)
	
	xall = FINDGEN(N_ELEMENTS(int))
	x = FLTARR(sbin1) ; spline arrays
	y = FLTARR(sbin1) ; spline arrays
	
	binsize = long(N_ELEMENTS(int)/sbin1)
	FOR b=0,sbin1-1 DO BEGIN
		i1 = binsize*b
		i2 = binsize*(b+1)
		IF i2 GE N_ELEMENTS(int) THEN i2=N_ELEMENTS(int)-1
		bin = xall[i1:i2]
		fbin = int[bin]

		; highest frac fraction of points in bin
		lim = (fbin[SORT(fbin)])[fix((1.-frac)*N_ELEMENTS(fbin))]
		use = WHERE(fbin GT lim, count)
		; five points in each bin at least...
		WHILE count LT 5. DO BEGIN	
			lim = (fbin[SORT(fbin)])[fix((1.-frac+0.1)*N_ELEMENTS(fbin))]
			use = WHERE(int GT lim, count)
			stop
		ENDWHILE
		
		x[b] = MEDIAN(bin[use])
		y[b] = MEDIAN(fbin[use])
	ENDFOR
	
	c = SPLINE(x,y,FINDGEN(N_ELEMENTS(int)))

	IF KEYWORD_SET(showplot) THEN BEGIN
		plot, int
		oplot, x, y
		oplot, c
	ENDIF
	
	RETURN, c
	
END

;============================================
; FUNCTION FLATTEN
; Remove continuum from a spectrum



FUNCTION FLATTEN, int, showplot=showplot, contf=contf, frac=frac, sbin=sbin
	
	IF ~KEYWORD_SET(sbin) THEN BEGIN
	  sbin1 = 10
	  sbin2 = 6
	ENDIF ELSE BEGIN
	  sbin1=sbin
	  sbin2=sbin
	ENDELSE
	
	IF ~KEYWORD_SET(frac) THEN BEGIN
	  frac1 = 0.4
	  frac2 = 0.2
	ENDIF ELSE BEGIN
	  frac1 = frac
	  frac2 = frac
	ENDELSE
	
; 	IF KEYWORD_SET(showplot) THEN showplot=2

	IF KEYWORD_SET(contf) THEN BEGIN
	
		CONTF, int, c1, nord=4, sbin=sbin1,frac=frac1, plot=0, mask=fin, xspline=xspline, yspline=yspline
		t1=SPLINE(xspline,yspline,FINDGEN(N_ELEMENTS(int)))
		
		CONTF, int/t1, c2, nord=4, sbin=sbin2,frac=frac2, plot=0, mask=fin, xspline=xspline, yspline=yspline
		t2=SPLINE(xspline,yspline,FINDGEN(N_ELEMENTS(int)))
		
		flat=int/t1/t2
		
	ENDIF ELSE BEGIN
	
		t1 = SPLINE_CONT(int, sbin=sbin1, frac=frac1, showplot=0)
		t2 = SPLINE_CONT(int/t1, sbin=sbin2, frac=frac2, showplot=0)
		flat=int/t1/t2
	
	ENDELSE

	; sometimes goes bad...
	flat[WHERE(flat LT 0)] = 0
	flat[WHERE(flat GT 2)] = 0
	IF KEYWORD_SET(showplot) THEN BEGIN
		plot, flat
		wait, 2
	ENDIF
	
	RETURN, flat

END




;============================================
; PRO ERN_RV
; Find the offset between data and a standard star



PRO ERN_RV, data, std, rv0=rv0, chi0=chi0, $
	showplot=showplot, $ ; show plots in output?
	pixscale=pixscale, $ ; pixel scale
	wrange=wrange, $ ; wavelenghts to flatten over,
	oversamp=oversamp, $ ; oversampling multiple
	zero=zero, nan=nan, $ ; bad data flags
	ccorr=ccorr, $ ; choosing cross-correlation routine
	corr_range=corr_range, $ ; cross-corr range (assumes microns)
	contf=contf, frac=frac, sbin=sbin, $ ; for flattening spectrum
	quiet=quiet, $ ;mask=mymask, $
	flatten=flatten, $  ; flatten spec?
	svd=svd ; flatten using SVD?

	IF KEYWORD_SET(showplot) THEN showplot=1
	IF ~KEYWORD_SET(pixscale) THEN pixscale = abs(median(data[0:-2,0]-data[1:-1,0]))
	IF ~KEYWORD_SET(wrange) THEN wrange = [data[0,0],data[-1,0]]
	IF ~KEYWORD_SET(oversamp) THEN oversamp=6. 
	IF N_ELEMENTS(flatten) EQ 0 THEN flatten=1
	IF N_ELEMENTS(mymask) EQ 0 THEN $
	  mask = FLTARR(N_ELEMENTS(data))+1 $
	ELSE IF N_ELEMENTS(mask) EQ N_ELEMENTS(data[*,0]) THEN $
	  mask = mymask $
	ELSE message, "ERN_RV: Mask doesn't match data"

	; select wavelength range (log lambda)
	start_wl = ALOG(wrange[0])
	end_wl = ALOG(wrange[1])
	roi = WHERE(data[*,0] GE wrange[0] AND data[*,0] LE wrange[1])

	; create oversampled grid uniformly spaced in log lambda	
	wl_vector = SCALE_VECTOR(FINDGEN((end_wl-start_wl)*oversamp*mean(data[roi,0])/pixscale), start_wl, end_wl) 

	IF KEYWORD_SET(zero) THEN BEGIN
	  roi = roi[WHERE(data[roi,1] GT 0)]
	  IF ~KEYWORD_SET(quiet) THEN print, "ERN_RV: using non-zero fluxes."
	ENDIF ELSE IF KEYWORD_SET(nan) THEN BEGIN
	  roi = roi[WHERE(FINITE(data[roi,1]))]
	  IF ~KEYWORD_SET(quiet) THEN print, "ERN_RV: using finite fluxes."
	ENDIF ELSE BEGIN
	  IF ~KEYWORD_SET(quiet) THEN print, "ERN_RV: using all provided data."
	ENDELSE

	; interpolate object and standard onto new grid
	int_obj = INTERPOL(data[*,1],ALOG(data[*,0]),wl_vector, /spline) 
	sroi = WHERE(std[*,0] GE wrange[0] AND std[*,0] LE wrange[1])
	int_std = INTERPOL(std[*,1],ALOG(std[*,0]),wl_vector, /spline)	

	; remove continuum
	IF ~KEYWORD_SET(quiet) THEN BEGIN
	  IF KEYWORD_SET(contf) THEN print, "ERN_RV: using contf routine for flattening" $
	  ELSE print, "ERN_RV: using spline flattening routine"
	  IF KEYWORD_SET(frac) THEN print, "ERN_RV: using user-supplied value for frac in FLATTEN. Are you sure?"
	  IF KEYWORD_SET(sbin) THEN print, "ERN_RV: using user-supplied value for sbin in FLATTEN. Are you sure?"
	ENDIF
	
	IF KEYWORD_SET(flatten) THEN BEGIN
	  IF KEYWORD_SET(svd) THEN BEGIN
	    res = SVDFIT(wl_vector, int_obj, 3, /legendre, YFIT=fit)
	    flat_obj = int_obj/fit
	    plot, wl_vector, flat_obj 
 	    res = SVDFIT(wl_vector, int_std, 3, /legendre, YFIT=fit)
 	    flat_std = int_std/fit
	    oplot, wl_vector, flat_std, color=2
	    wait, 1
	  ENDIF ELSE BEGIN
	    flat_obj=FLATTEN(int_obj, showplot=showplot, contf=contf, frac=frac, sbin=sbin)
	    flat_std=FLATTEN(int_std, showplot=showplot, contf=contf, frac=frac, sbin=sbin)
	  ENDELSE
	ENDIF ELSE BEGIN
	  flat_obj = int_obj/median(int_obj)
	  flat_std = int_std/median(int_std)
	ENDELSE
; 	IF KEYWORD_SET(showplot) THEN BEGIN
; 			
; 		; this plot checks the interpolations
; 		plot, wl_vector, int_obj
; 		oplot, wl_vector, int_std-0.4
; 		oplot, ALOG(data[*,0]), data[*,1], co=3
; 		oplot, ALOG(data[*,0]), data[*,1]-0.4, co=4
; 		oplot, [ALOG(2.2063),ALOG(2.2063)],[0,5],linestyle=2
; 		oplot, [ALOG(2.26267),ALOG(2.26267)],[0,5],linestyle=2
; 		wait, 2
; 	
; 	ENDIF

	IF NOT KEYWORD_SET(corr_range) THEN corr_range=fix(20*.0005/pixscale)
	IF NOT KEYWORD_SET(ccorr) THEN ccorr='c_correlate'
	IF ccorr EQ 'xcorl' THEN BEGIN
		IF ~KEYWORD_SET(quiet) THEN print, "ERN_RV: Using xcorl"
		xcorl, flat_std, flat_obj, corr_range, shft, chi, minchi, plot=showplot, print=showplot
		chi0 = 1.-minchi/(total(flat_obj*flat_obj))
	ENDIF ELSE IF ccorr EQ 'c_correlate' THEN BEGIN
		IF ~KEYWORD_SET(quiet) THEN print, "ERN_RV: Using c_correlate + modifications"
		lag = findgen(corr_range*2)-float(corr_range)
		lag[where(lag EQ 0)] = 0.01 ; hack because it spikes at 0 for no reason
		result = c_correlate(flat_obj, flat_std, lag) ; opposite order
		pk = MAX(result,p)

		IF KEYWORD_SET(showplot) THEN BEGIN
			plot, lag, result, /ynozero, xrange=[-20,20]
			wait, 1
		ENDIF
		; if not at ends, quadratic offset to get better peak
		IF (p GT 0) AND (p LT (N_ELEMENTS(lag)-1)) THEN BEGIN
			aa = result[p]
			bb = 0.5*(result[p+1] - result[p-1])
			cc = 0.5*(result[p+1] + result[p-1] - 2.0*aa)
			offset = -0.5*bb/cc
		ENDIF ELSE $
			offset = 0.
		shft = lag[p] + offset
		chi0 = 1*pk
		
	ENDIF ELSE IF ccorr EQ 'cross_correlate' THEN BEGIN
		IF ~KEYWORD_SET(quiet) THEN print, "ERN_RV: Using cross_correlate"
		cross_correlate, flat_std, flat_obj, shft, result, width=corr_range*2
		chi0 = 1*result
	ENDIF ELSE message, 'Cross-correlation routine not implemented: ', ccorr

	; these are the pixel arrays
	pix_fiducial = MAKE_ARRAY(N_ELEMENTS(wl_vector),/index)
	pix_shifted = pix_fiducial + shft
	wl_shifted = INTERPOL(wl_vector, pix_fiducial, pix_shifted)

	; offset and relative RV
	; if shift is positive, shifting object to redder wavelengths = it was blue shifted = negative RV

	offset=-shft*pixscale/(oversamp*mean(data[roi,0]))  	; delta lambda/lambda
	RV0=3.0e5*offset					; velocity

	IF KEYWORD_SET(showplot) THEN BEGIN
	wait, 2
		erase & multiplot, [1,3]
		
		plot, pix_fiducial, int_obj/MEAN(int_obj), /xsty, title='Shifting observed to match standard'
		oplot, pix_shifted, int_obj/MEAN(int_obj)-0.2, co=3
		oplot, pix_fiducial, int_std/MEAN(int_std)-0.4, co=2
		oplot, [N_ELEMENTS(wl_vector)/5,N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		oplot, [2*N_ELEMENTS(wl_vector)/5,2*N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		oplot, [3*N_ELEMENTS(wl_vector)/5,3*N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		oplot, [4*N_ELEMENTS(wl_vector)/5,4*N_ELEMENTS(wl_vector)/5], [0,2],linestyle=2
		
		multiplot
		
		plot, wl_vector, flat_obj, yrange=[0.3,1.2], /xsty
		oplot, wl_shifted, flat_obj-0.2, co=3
		oplot, wl_vector, flat_std-0.4, co=2
		oplot, [ALOG(2.2063),ALOG(2.2063)],[0,5],linestyle=2
		oplot, [ALOG(2.26267),ALOG(2.26267)],[0,5],linestyle=2

		multiplot
		
		plot, data[*,0], data[*,1]/MEAN(data[roi,1]), xrange=wrange, /xsty
		oplot, data[*,0]-offset*data[*,0], data[*,1]/MEAN(data[roi,1])-0.2, co=3
		oplot, std[*,0], std[*,1]/MEAN(std[WHERE(FINITE(std[roi,1])),1])-0.4, co=2

; 		print, 'The radial velocity is ', RV0
		al_legend, ['Object','Shifted object', 'Standard'], color=[1,3,2], linestyle=[0,0,0], /right, /bottom
		multiplot, /default
		
		wait, 2
	ENDIF

END



