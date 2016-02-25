PRO fire_order_variables, hdr, order, wrange, trange, pixscale, polydegree

CASE order OF
	12: BEGIN
		wrange = [2.19, 2.22]
		trange = [2.16, 2.20]
	END
	15: BEGIN
		wrange = [1.68,1.72] ; extends over order boundary
		trange = [1.72,1.76]
	END
	16: BEGIN
		wrange = [1.62, 1.66]
		trange = [1.63, 1.67]
	END
	17: BEGIN
		wrange = [1.53, 1.565]
		trange = [1.515, 1.55]
	END
	20: BEGIN
		wrange = [1.29, 1.32]
		trange = [1.29, 1.32]
	END
	21: BEGIN
		wrange = [1.23, 1.255]
		trange = [1.23, 1.265]
	END
	22: BEGIN
		wrange = [1.175, 1.20]
		trange = [1.175, 1.21]
	END
	24: BEGIN
		wrange = [1.07, 1.1]
		trange = [1.08, 1.11]
	END
	27: BEGIN ; wavecal may be bad
		wrange = [0.98, 0.995] ; extends over order boundary
		trange = [0.97, 0.982] 
	END
	28: BEGIN
		wrange = [0.915,0.93]
		trange = [0.915,0.93]
	END
	29: BEGIN
		wrange = [0.875, 0.895] ; extends over order boundary
		trange = [0.885,0.91]
	END

	ELSE: BEGIN
	      print, 'Order is', order
	      message, 'This is an invalid order'
	END
ENDCASE

pixscale = fxpar(hdr, 'CDELT1')
polydegree = 5

END