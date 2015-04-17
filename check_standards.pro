PRO check_standards

  ;file is big, best to read it in only once
  atrans=XMRDFITS('/home/enewton/pro/Spextool/data/atrans.fits',0)

  READCOL, 'standards.txt', stars, sptypes, rvs, skipline=1, format='(A,I,F)'
  spec_dir = 'spec/'

  file = spec_dir+stars[4]+".fits" 
  std = MRDFITS(file,0,shdr)
  file_tc = spec_dir+stars[4]+"_tc.fits"
  std_tc = MRDFITS(file)
  stdrv = rvs[4]

  FOREACH star, stars DO BEGIN 

    file = spec_dir+star+".fits" 
    data = MRDFITS(file,0,hdr)
    file_tc = spec_dir+star+"_tc.fits"
    data_tc = MRDFITS(file)

    order = 2
    order_variables, order, wrange, trange, telescope="irtf"
    TELL_RV, order, atrans, data[*,*,order], hdr, $
      data_tc[*,*,order], std[*,*,order], std_tc[*,*,order], shdr, $
      wrange=wrange, trange=trange, $
      shft=myshft, s_shft=s_shft, rv=myrv, $
      showplot=0

    print, star, myrv
  
  ENDFOREACH
  
END


