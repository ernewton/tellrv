; ##################
; get order variables for RV fitting
; ##################

PRO order_variables, order, wrange, trange, telescope=telescope

  IF KEYWORD_SET(telescope) THEN BEGIN
    IF telescope EQ "irtf" THEN BEGIN
      print, "IRTF"
      irtf_order_variables, order, wrange, trange 
    ENDIF ELSE BEGIN
      print, "Invalid telescope"
    ENDELSE
  ENDIF ELSE BEGIN
      print, "IRTF"
      irtf_order_variables, order, wrange, trange 
  ENDELSE

END