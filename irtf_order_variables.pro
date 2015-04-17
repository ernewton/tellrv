PRO irtf_order_variables, order, wrange, trange

CASE order OF
	0: BEGIN
		wrange = [2.18, 2.41]
		trange = [1.995, 2.4]
	END
	1: BEGIN
		wrange = [1.49, 1.73]
		trange = [1.43, 1.8]
	END
	2: BEGIN
		wrange = [1.15,1.32]
		trange = [1.142,1.35]
	END
	3: BEGIN
		wrange = [1.0, 1.1]
		trange = [0.94, 1.17]
	END
	4: BEGIN
		wrange = [0.82, 0.92]
		trange = [0.89, 0.99]
	END
ENDCASE

END