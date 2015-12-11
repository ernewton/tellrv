PRO fire_order_variables, hdr, order, wrange, trange, pixscale, polydegree

CASE order OF
	12: BEGIN
		wrange = [2.19, 2.22]
; 		trange = [2.24, 2.28]
		trange = [2.16, 2.20]
	END
	15: BEGIN
		wrange = [1.715,1.745] 
		trange = [1.72,1.8]
	END
	16: BEGIN
		wrange = [1.62,  1.66]
		trange = [1.63, 1.67]
	END
	17: BEGIN
		wrange = [1.49, 1.53]
		trange = [1.53, 1.56]
	END
	20: BEGIN
		wrange = [1.28, 1.3]
		trange = [1.3, 1.31]
	END
	21: BEGIN
		wrange = [1.235, 1.255]
		trange = [1.255, 1.265]
	END
	22: BEGIN
		wrange = [1.165, 1.2]
		trange = [1.175, 1.21]
	END
	24: BEGIN
		wrange = [1.07, 1.1]
		trange = [1.09, 1.12]
	END
	27: BEGIN
		wrange = [0.98, 1.0]
		trange = [0.96, 0.975] 
	END

	ELSE: message, 'Order is', order, 'This is an invalid order'
ENDCASE

pixscale = fxpar(hdr, 'CDELT1')
polydegree = 5

END