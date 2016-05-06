PRO fire_order_variables, hdr, order, wrange, trange, pixscale, polydegree

CASE order OF
	; 0 -- wavecal usually fails
	1: BEGIN
		wrange = [2.19, 2.22]
		trange = [2.16, 2.20]
	END
	2: BEGIN
		trange = [2.03, 2.12]
		wrange = [1.95, 2.10]
	END
	; 3 -- too much atm; no lines
	4: BEGIN
		wrange = [1.68,1.72]
		trange = [1.72,1.76]
	END
	5: BEGIN
		wrange = [1.62, 1.66]
		trange = [1.63, 1.67]
	END
	6: BEGIN
		wrange = [1.53, 1.565]
		trange = [1.515, 1.55]
	END
	7: BEGIN ; 7 -- too much atm; but Mg line
		wrange = [1.48, 1.51]
		trange = [1.49, 1.51]
	END
	; 8 -- too much atm., better portion overlaps with 9
	9: BEGIN
		wrange = [1.29, 1.32]
		trange = [1.29, 1.32]
	END
	10: BEGIN
		wrange = [1.23, 1.255]
		trange = [1.23, 1.265]
	END
	11: BEGIN
		wrange = [1.175, 1.20]
		trange = [1.175, 1.21]
	END
	12: BEGIN
		wrange = [1.1355, 1.17]
		trange = [1.135, 1.16]
	END
	13: BEGIN
		wrange = [1.07, 1.1]
		trange = [1.08, 1.11]
	END
	; 14 -- no atm. absorption!
	15: BEGIN ; wavecal may be bad
		wrange = [0.995, 1.03]
		trange = [0.975, 0.985] 
	END
	; 16 -- too much atm.
	17: BEGIN
		wrange = [0.905,0.93]
		trange = [0.905,0.93]
	END
	; 18 -- overlaps with 17
	; 19 -- no atm.
	; 20 -- no atm.

	ELSE: BEGIN
	      print, 'Order is', order
	      message, 'This is an invalid order'
	END
ENDCASE

pixscale = fxpar(hdr, 'CDELT1')
polydegree = 5

END