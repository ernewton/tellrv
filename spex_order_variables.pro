PRO spex_order_variables, hdr, order, wrange, trange, pixscale, polydegree

CASE order OF
	0: BEGIN
		wrange = [2.18, 2.41]
		trange = [1.995, 2.4]
		pixscale = FXPAR(hdr, 'DISPO03')
		polydegree = 5
	END
	1: BEGIN
		wrange = [1.49, 1.73]
		trange = [1.43, 1.8]
		pixscale = FXPAR(hdr, 'DISPO04')
		polydegree = 4
	END
	2: BEGIN
		wrange = [1.15,1.32]
		trange = [1.142,1.35]
		pixscale = FXPAR(hdr, 'DISPO05')
		polydegree = 5
	END
	3: BEGIN
		wrange = [1.0, 1.1]
		trange = [0.94, 1.17]
		pixscale = FXPAR(hdr, 'DISPO06')
		polydegree = 5
	END
	4: BEGIN
		wrange = [0.82, 0.92]
		trange = [0.89, 0.99]
		pixscale = FXPAR(hdr, 'DISPO07')
		polydegree = 4
	END
	5: BEGIN
		wrange = [0.81, 0.9]
		trange = [0.81, 0.9]
		pixscale = FXPAR(hdr, 'DISPO08')
		polydegree = 5
	END
	ELSE: message, 'Order is', order, 'This is an invalid order'
ENDCASE

END