FUNCTION SPLINE_CONT, int, frac=frac, sbin=sbin, showplot=showplot

	IF ~KEYWORD_SET(frac) THEN frac=0.2
	IF ~KEYWORD_SET(sbin) THEN sbin=10
	ordrd = SORT(int)
	
	xall = FINDGEN(N_ELEMENTS(int))
	x = FLTARR(sbin) ; spline arrays
	y = FLTARR(sbin) ; spline arrays
	
	binsize = fix(N_ELEMENTS(int)/sbin)
	FOR b=0,sbin-1 DO BEGIN
		i1 = binsize*b
		i2 = binsize*(b+1)
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


	IF KEYWORD_SET(showplot) THEN showplot=2

	IF KEYWORD_SET(contf) THEN BEGIN
	
		CONTF, int, c1, nord=4, sbin=10,frac=0.5, plot=showplot, mask=fin, xspline=xspline, yspline=yspline
		t1=SPLINE(xspline,yspline,FINDGEN(N_ELEMENTS(int)))
		
		CONTF, int/t1, c2, nord=4, sbin=6,frac=0.2, plot=showplot, mask=fin, xspline=xspline, yspline=yspline
		t2=SPLINE(xspline,yspline,FINDGEN(N_ELEMENTS(int)))
		
		flat=int/t1/t2
		
	ENDIF ELSE BEGIN
	
		t1 = SPLINE_CONT(int, sbin=10, frac=0.4, showplot=showplot)
		t2 = SPLINE_CONT(int/t1, sbin=10, frac=0.2, showplot=showplot)
		flat=int/t1/t2
	
	ENDELSE
	IF KEYWORD_SET(showplot) THEN wait,2

	RETURN, flat

END




;============================================
; PRO RV
; Find the offset between data and a standard star



PRO ERN_RV, data, std, rv0=rv0, chi=chi,  $
	showplot=showplot, $ ; show plots in output?
	pixscale=pixscale, $ ; pixel scale
	wrange=wrange, $ ; wavelenghts to flatten over,
	oversamp=oversamp, $ ; oversampling multiple
	xcorl=xcorl, ccorr=ccorr, $ ; choosing cross-correlation routine
	corr_range=corr_range, $ ; cross-corr range
	contf=contf, frac=frac, sbin=sbin ; for flattening spectrum

	IF KEYWORD_SET(showplot) THEN showplot=1

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
	IF ~KEYWORD_SET(oversamp) THEN $
		oversamp=6. 
	wl_vector = SCALE_VECTOR(FINDGEN((end_wl-start_wl)*oversamp*mean(data[roi,0])/pixscale), start_wl, end_wl) 

	; interpolate object and standard onto new grid
	int_obj = INTERPOL(data[*,1],ALOG(data[*,0]),wl_vector, /spline) 
	int_std = INTERPOL(std[*,1],ALOG(std[*,0]),wl_vector, /spline)	

	; remove continuum
	flat_obj=FLATTEN(int_obj, showplot=showplot, contf=contf, frac=frac, sbin=sbin)
	flat_std=FLATTEN(int_std, showplot=showplot, contf=contf, frac=frac, sbin=sbin)

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
	IF KEYWORD_SET(xcorl) THEN BEGIN
		xcorl, flat_std, flat_obj, corr_range, shft, chi, minchi, plot=showplot, print=showplot
	ENDIF ELSE IF KEYWORD_SET(ccorr) THEN BEGIN
		lag = indgen(corr_range*2)-corr_range
		result = c_correlate(flat_obj, flat_std, lag) ; opposite order
		pk = MAX(result,p)
		; if not at ends, quadratic offset to get better peak
		IF (p GT 0) AND (p LT (N_ELEMENTS(lag)-1)) THEN BEGIN
			aa = result[p]
			bb = 0.5*(result[p+1] - result[p-1])
			cc = 0.5*(result[p+1] + result[p-1] - 2.0*aa)
			offset = -0.5*bb/cc
		ENDIF ELSE $
			offset = 0.
		shft = lag[p] + offset
	ENDIF ELSE BEGIN
		cross_correlate, flat_std, flat_obj, shft, result, width=corr_range*2
	ENDELSE

	; these are the pixel arrays
	pix_fiducial = MAKE_ARRAY(N_ELEMENTS(wl_vector),/index)
	pix_shifted = pix_fiducial + shft
	wl_shifted = INTERPOL(wl_vector, pix_fiducial, pix_shifted)

	; offset and relative RV
	; if shift is positive, shifting object to redder wavelengths = it was blue shifted = negative RV

	offset=-shft*pixscale/(oversamp*mean(data[roi,0]))  	; delta lambda/lambda
	RV0=3.0e5*offset					; velocity

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

		print, 'The radial velocity is ', RV0
		
; 		wait, 2
		stop
	ENDIF

END



