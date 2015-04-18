PRO check_standards

  ;file is big, best to read it in only once
  atrans=XMRDFITS('/home/enewton/pro/Spextool/data/atrans.fits',0)

  READCOL, 'standards.txt', stars, sptypes, rvs, skipline=1, format='(A,I,F)'
  spec_dir = 'spec/'
  
  file = spec_dir+stars[0]+".fits" 
  std0 = MRDFITS(file,0,shdr)
  file_tc = spec_dir+stars[0]+"_tc.fits"
  std_tc0 = MRDFITS(file_tc)
  std_tc0[*,0,*]=std0[*,0,*]	; want to use original wavelength array
  stdrv = rvs[0]
  print, file
  print, stdrv

  FOREACH star, stars DO BEGIN 
    std = std0
    std_tc = std_tc0

    file = spec_dir+star+".fits" 
    data = MRDFITS(file,0,hdr)
    file_tc = spec_dir+star+"_tc.fits"
    data_tc = MRDFITS(file_tc)

    order = 2
    data_tc[*,0,*]=data[*,0,*]
    print, file_tc
    print, data[0,*,order]
    print, data_tc[0,*,order]
    print, std[0,*,order]
    print, std_tc[0,*,order]

    order_variables, order, wrange, trange, telescope="irtf"
    NIR_RV, order, atrans, data[*,*,order], hdr, $
      data_tc[*,*,order], std[*,*,order], std_tc[*,*,order], shdr, $
      wrange=wrange, trange=trange, $
      shft=myshft, s_shft=s_shft, rv=myrv, $
      showplot=0, torest=rv0

    print, star, myrv, rv0
  
  ENDFOREACH
  
END


