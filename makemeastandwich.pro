pro makemeastandwich

  file = "spec/J0727+0513.fits" 
  std0 = MRDFITS(file,0,shdr)

  file_tc = "spec/J0727+0513_tc.fits"
  std_tc0 = MRDFITS(file_tc, /silent)
  std_tc0[*,0,*] = std0[*,0,*]

  rv = 18.2
  helio = GET_HELIO(shdr)
  
  hdr=shdr

  for i=0,5 do begin

    order = i
    order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument="spex"
    
    std = std0
    std_tc = std_tc0
    data = std0
    data_tc = std_tc0
    
    NIR_RV, order, data[*,*,order], hdr, $
      data_tc[*,*,order], std[*,*,order], std_tc[*,*,order], shdr, $
      atrans=atrans, $
      pixscale=pixscale, polydegree=polydegree, $
      spixscale=pixscale, spolydegree=polydegree, $ ; standard is from same set-up 
      wrange=wrange, trange=trange, $
      shft=myshft, s_shft=s_shft, rv=myrv, $
      showplot=1, torest=rv0, $
      ccorr=ccorr, xcorl=xcorl, $
      contf=contf, frac=frac, sbin=sbin

    stop

    lambda = std0[*,0,i]*(1-(rv-helio)/(3.e5))
    plot, std0[*,0,i], std_tc0[*,1,i]
    std_tc0[*,0,i] = lambda
    oplot, std_tc0[*,0,i], std_tc0[*,1,i], co=2
  
  endfor
  
  mwrfits, std_tc0, "spec/J0727+0513_rest.fits", shdr, /create

end