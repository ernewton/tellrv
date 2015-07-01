; ##################
; get order variables for RV fitting
; ##################

PRO order_variables, hdr, order, wrange, trange, pixscale, polydegree, telescope=telescope

  IF KEYWORD_SET(telescope) THEN BEGIN
    IF telescope EQ "irtf" THEN BEGIN
      irtf_order_variables, hdr,order, wrange, trange, pixscale, polydegree 
    ENDIF ELSE BEGIN
      message, "Invalid telescope"
    ENDELSE
  ENDIF ELSE BEGIN
      irtf_order_variables, hdr,order, wrange, trange, pixscale, polydegree
  ENDELSE

END