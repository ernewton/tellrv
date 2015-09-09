PRO check_standards, xcorl=xcorl, ccorr=ccorr, contf=contf, showplot=showplot

  ;file is big, best to read it in only once
  atrans=XMRDFITS('/home/enewton/pro/Spextool/data/atrans.fits',0)

  ; list of standard stars and their RVs
  READCOL, 'standards.txt', stars, sptypes, rvs, skipline=1, format='(A,I,F)'
  spec_dir = 'spec/'
  
  ; non-telluric corrected spectrum
  file = spec_dir+stars[0]+".fits" 
  std0 = MRDFITS(file,0,shdr)
  
  ; telluric-corrected spectrum
  file_tc = spec_dir+stars[0]+"_tc.fits"
  std_tc0 = MRDFITS(file_tc, /silent)
  std_tc0[*,0,*]=std0[*,0,*]	; want to use original wavelength array
  stdrv = rvs[0]

  i = 0
  FOREACH star, stars DO BEGIN 
    std = std0
    std_tc = std_tc0

    file = spec_dir+star+".fits" 
    data = MRDFITS(file,0,hdr, /silent)
    file_tc = spec_dir+star+"_tc.fits"
    data_tc = MRDFITS(file_tc, /silent)

    print, "Star is: ", star
    print, "File is: ", file
    print, "RV is:   ", rvs[i]

    order = 2
    data_tc[*,0,*]=data[*,0,*]

    order_variables, hdr, order, wrange, trange, pixscale, polydegree, telescope="irtf"
    NIR_RV, order, data[*,*,order], hdr, $
      data_tc[*,*,order], std[*,*,order], std_tc[*,*,order], shdr, $
      atrans=atrans, $
      pixscale=pixscale, polydegree=polydegree, $
      spixscale=pixscale, spolydegree=polydegree, $ ; standard is from same set-up 
      wrange=wrange, trange=trange, $
      shft=myshft, s_shft=s_shft, rv=myrv, $
      showplot=showplot, torest=rv0, $
      ccorr=ccorr, xcorl=xcorl, $
      contf=contf, frac=frac, sbin=sbin

    print, "Measured absolute RV:", myrv
    print, "RV offset to rest:", rv0
    i++
    
  ENDFOREACH
  
END


