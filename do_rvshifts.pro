;============================================
; FUNCTION DO_RVSHIFTS
; get and apply heliocentric correction



FUNCTION DO_RVSHIFTS, rv0, hdr, shdr, rv_std=rv_std, bc=hcorr


	; barycentric correction
	hcorr=GET_HELIO(hdr)
	
	; offset from standard star
	hcorr_std=GET_HELIO(shdr)
	IF NOT KEYWORD_SET(rv_std) THEN rv_std=18.2
	RVoff=rv_std-hcorr_std
	
	; heliocentric corrected absolute RV
	RV=RV0+hcorr+RVoff
	
	RETURN, rv

END

