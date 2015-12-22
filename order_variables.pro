;============================================
; get order variables for RV fitting

PRO order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument=instrument

  IF KEYWORD_SET(instrument) THEN BEGIN
    IF instrument EQ "spex" THEN BEGIN
      spex_order_variables, hdr,order, wrange, trange, pixscale, polydegree 
    ENDIF ELSE IF instrument EQ "fire" THEN BEGIN
      fire_order_variables, hdr,order, wrange, trange, pixscale, polydegree 
    ENDIF ELSE BEGIN
      message, "Invalid telescope"
    ENDELSE
  ENDIF ELSE BEGIN
      spex_order_variables, hdr,order, wrange, trange, pixscale, polydegree
  ENDELSE

END