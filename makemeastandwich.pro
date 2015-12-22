; Makes the RV standard on which this code has been tested.
; note that a -2.6 km/s offset correction was applied in Newton et al. (2014)
; if using the standard as created here, this is not necessary.
; the difference arises from slightly different treatment of the way
; the "zero RV" standard is applied in this code, which mitigates the
; random errors that are present.

pro makemeastandwich

  file = "spec/J0727+0513.fits" 
  std0 = MRDFITS(file,0,shdr)

  file_tc = "spec/J0727+0513_tc.fits"
  std_tc0 = MRDFITS(file_tc, /silent)
  std_tc0[*,0,*] = std0[*,0,*]
  
  hdr=shdr

  for i=0,5 do begin

    order = i
    order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument="spex"
    
    std = std0
    std_tc = std_tc0
    data = std0
    data_tc = std_tc0
    
    NIR_RV, data[*,*,order], hdr, data_tc[*,*,order],  $
      std[*,*,order], shdr, std_tc[*,*,order],  $
      atrans=atrans, stdrv=0., $
      pixscale=pixscale, polydegree=polydegree, $
      spixscale=pixscale, spolydegree=polydegree, $ ; standard is same
      wrange=wrange, trange=trange, $
      ccorr_fxn='xcorl', $
      rv=rv, $ ; true rv of star, will be 18.2
      relrv=relrv, $ ; rv relative to standard, will be 0.
      torest=torest, $ ; rv as observed from Earth
      mshft=mshft, $  ; shift required to get to absolute wavelength calibration
      shftarr=shftarr ; total shift to zero RV
    
    std_tc0[*,0,order] = std_tc0[*,0,order]+mshft ; CORRECTS FOR WAVELENGTH CALIBRATION
    print, relrv
    print, rv
    print, torest

  endfor
  
  std_tc0[*,0,*] = std_tc0[*,0,*]*(1.-(18.2-GET_HELIO(hdr))/(3.e5)) ; CORRECTS FOR KNOWN RADIAL VELOCITY
  mwrfits, std_tc0, "spec/J0727+0513_rest.fits", shdr, /create

end