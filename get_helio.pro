

;============================================
; FUNCTION GET_HELIO
; Use hfits fx commands to get info out of header and make helio(bary)centric correction



FUNCTION GET_HELIO, hdr, degrees=degrees, barycentric=barycentric, force=force, quiet=quiet

	IF N_ELEMENTS(barycentric) EQ 0 THEN barycentric = 1

	; telescope longitude and latitude 
	tel = FXPAR(hdr,'TELESCOP')
	CASE tel OF
	  'NASA IRTF': BEGIN
		lat = 19. + 49./60. + 34.38594/3600. ; north
		lon = 360. - (155. + 28./60. + 19.19564/3600.) ; west
		alt = 13674.76*0.3048 ; feet to meters
		date = FXPAR(hdr, 'DATE_OBS')
		time = FXPAR(hdr, 'TIME_OBS')
		degrees = 0
		END	
	  'Baade': BEGIN
		; helio already in header
		IF ~KEYWORD_SET(force) THEN RETURN, FXPAR(hdr, 'HELIO')

		print, "GET_HELIO: calculating from header info for FIRE."
		lat = FXPAR(hdr, 'SITELAT')
		lon = FXPAR(hdr, 'SITELONG')
		alt = FXPAR(hdr, 'SITEALT')
		date = FXPAR(hdr, 'DATE-OBS')
		time = FXPAR(hdr, 'UT-TIME')
		degrees = 1				
		END
	   'TILLINGHAST': BEGIN
		;
		IF ~KEYWORD_SET(force) THEN RETURN, FXPAR(hdr, 'BCV')
		message, 'GET_HELIO: Not implemented.'
		END
	   ELSE: message, 'Telescope invalid'
	ENDCASE


	; from date and UT time, get Julian day
	date_ex = STRSPLIT(date,'-',/EXTRACT) ; 0 = year, 1 = month, 2 = day
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
	ENDIF ELSE BEGIN
		ra_deg = ra
		dec_deg = dec
	ENDELSE

	IF KEYWORD_SET(barycentric) THEN BEGIN
		baryvel, jdate, 0, vh, vb
		ra_rad=ra_deg/!RADEG
		dec_rad=dec_deg/!RADEG
		vel = vb[0]*cos(dec_rad)*cos(ra_rad) + $   ;Project velocity toward star
		      vb[1]*cos(dec_rad)*sin(ra_rad) + vb[2]*sin(dec_rad) 
		IF ~KEYWORD_SET(quiet) THEN print, "GET_HELIO: returning barycentric."
		RETURN, vel
	ENDIF ELSE BEGIN
		hcorr = heliocentric(ra_deg,dec_deg,jd=jdate,longitude=lon,latitude=lat, altitude=alt)
		IF ~KEYWORD_SET(quiet) THEN print, "GET_HELIO: returning heliocentric."
		RETURN, -hcorr
	ENDELSE

END
