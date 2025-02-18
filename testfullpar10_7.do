
mata:
version 13
mata clear
mata set matastrict on

//void myfun_combined(string scalar varlist, string scalar touse)
void testfullpar(		real matrix Z			)
{
	real matrix  X, Y, Xdnp, Xcnp, gridnp, gridid, grid
	real scalar dv, N, quant, ymin, ymax, xdnum, t_xdnum, xcnum
	
	
	
	quant = 0.5
	xdnum = 3
	t_xdnum = 3
	xcnum = 1
	
	
	
	
	Y = Z[.,1]
	X = Z[|1,2\.,.|]

	N = rows(X)
	X = (J(N,1,1) , X)
	dv = cols(X)
	
	ymin = min(Y)
	ymax = max(Y)
	
	
	st_numscalar("e(ymin)",ymin)
	st_numscalar("e(ymax)",ymax)
	st_numscalar("e(xdnum)",xdnum)
	st_numscalar("e(quant)",quant)
	//gridid = ((grid[.,1] + grid[.,2]) :< 2)
	//grid = select(grid,gridid)
	grid = grid'
	gridnp = grid[|1,1\xdnum+xcnum,.|]
	
	//rseed(13579)
	
	//--------------------------------------------------------------------------
	//defining parameters
	real scalar G, JJ, B, boots, mm, b0, d0, lower, upper, step, dd, dbb, i, j, b, gg, R
	real matrix l, ciboot1, phi, nb_mse, disd, disbd, disbb, disdbias, disbbbias, disbdbias
				
	G = 40
	l = (0.65,0.85,1.15,1.45)
	JJ = cols(l)
	B = 500
	boots = 550
	R = 200
	mm = 1.2
	b0 = 0
	d0 = 1
	
	ciboot1 = J(2*dv-2,2,0)
	
	lower = min((80/boots, 0.1))
	//lower = 0.05
	upper = 0.3
	step = (upper - lower)/G
	
	phi = (1,1,0,1,1)
	dbb = sum(phi)
	dd = dv - 1 - dbb

	pointer(pointer(real colvector)) colvector par_boot, par_bootb
	pointer(real matrix) rowvector  V_homb, V_mse
	//pointer(real matrix) rowvector  V_homd
	par_boot = J(B,1,NULL)
	par_bootb = J(B,1,NULL)
	V_homb = J(1,G,NULL)
	//V_homd = J(1,G,NULL)
	V_mse = J(1,G,NULL)
	
	
	real scalar nbtemp, tau, tempdis_d
	real matrix tempmsef, tempVd, tempVb, tempmse
	real matrix par_mse, par_hom, Nb_hom, chid, tempV
	//real matrix  chibb
	real scalar median_d1, median_d2, median_b1, tempdis, tempdis_b, tempnb
	real matrix tempchi_bootd, tempchi_bootbd, tempchi_bootbb
	tempchi_bootd = J(B,G,0)
	tempchi_bootbd = J(B,G,0)
	tempchi_bootbb = J(B,G,0)
	
	real matrix Zz, tempboot, tempZ

	for(b=1;b<=B;b++){
		timer_on(2)
		//idz = 1::N
		//tempidz = jumble(idz)
		//seltempidz = tempidz[|1,1\boots,1|]
		//Zz = Z[seltempidz,.]
		tempZ = jumble(Z)
		Zz = tempZ[|1,1\boots,.|]
		X = Zz[|1,2\.,.|]
		X = (J(boots,1,1), X)
		Y = Zz[|1,1\.,1|]
		par_boot[b] = &(J(G,1,NULL))
		par_bootb[b] = &(J(G,1,NULL))
		
		tempmsef = (0)
		tempmse = (0)
		
		for(gg=1;gg<=G;gg++){	
			tau = lower + step * gg
			myfun_hetero(tau,mm,b0,d0,X,Y,l,tempmsef,tempdis)
			tempchi_bootd[b,gg] = tempdis
			(*par_boot[b])[gg] = &(tempmsef[|1,1\dv-1,1|])
		    myfun_hom(tau,mm,phi,X,Y,l,tempmse,tempdis_d,tempdis_b)
			(*par_bootb[b])[gg] = &(tempmse)
			tempchi_bootbd[b,gg] = tempdis_d
			tempchi_bootbb[b,gg] = tempdis_b
		}
		timer_off(2)

	}
	}
end
