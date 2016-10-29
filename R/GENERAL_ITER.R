################################################################
###### GENERAL_ITER: function to execute iterations          ###
###### for MATLDA and MATLDA_GRID; gets passed initial       ###
###### iterations. For speed-ups, can make modifications     ###
###### to MATLDA functions and use GENERAL_ITERS             ###
###### in parallel                                           ### 
################################################################
###### Inputs: see MATLDA documentation                      ###
################################################################
##### Outputs: passed back to MATLDA functions               ###
################################################################

GENERAL_ITER = function(X, class, lambda1, lambda2,  quiet, UInit, VInit, U, V, 
	meanstacked,mean.array, M.update, D, weightmat, cov.tol, m.tol, full.tol, init.objfncval, k.iter){

	## priliminaries
	nc = count(class)[,2]

	C = length(nc)
	K = (C*(C-1))/2
	N = sum(nc)
	p = dim(X)[2]
	r = dim(X)[1]
	
	## normalize covariance matrices
	out1 = sum(abs(U))
	U.inner = (r/out1)*U
	V.inner = (out1/r)*V
	Uinv.inner = solve(U.inner)
	Vinv.inner = solve(V.inner)

	## set step size for AMA
	out1 = min(eigen(U.inner)$val)
	out2 = min(eigen(V.inner)$val)
	rho = (4/C)*(min(nc)/N)*out1*out2
	rho = rho/10

	## initialize
	k.outer = 1
	outer.converge = FALSE

	## compute obj func values
	old.objfncval = EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.inner,  weightmat, lambda1, lambda2)$out
	
	## print first objective function value 
	if (!quiet) {
		cat("k.outer=", k.outer,"objfncval = ", old.objfncval,"\n")
	}
	
	old.objfncval = init.objfncval

	## iterating 
	while (!outer.converge) {

			if (lambda1!=0){

				out = M_Update_Rcpp(Mean=meanstacked, 
				 	D=D, Uinv=Uinv.inner, Vinv=Vinv.inner, nc = nc, rho=rho, weightmat = weightmat, 
				 	lambda=lambda1, C=C, r=r, mtol=m.tol)
			
				D = out$D

				M = THRESH_MEAN(out$M, out$G, C)
				
				M.update = array(0, dim=c(r, p, C))
				for(c in 1:C){
					M.update[,,c] = M[((c-1)*r +1):(c*r),]
				}
				
			} else {
					M.update = mean.array
			}

			if (!quiet) {
	 			cat("obj.func after M update=", EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.inner, weightmat, lambda1, lambda2)$out, "\n")
	 		}

	 		## center Y and update V and U
			Xc = CENTER_X(X, class, M.update)
       		S.v = V_SAMPLE(Xc, U=U.inner)

       		if (lambda2!=0) {

  				V.t = glasso(S.v, rho=lambda2, trace=0,  penalize.diagonal=TRUE)$wi
  				
  				if (!quiet) {
  					cat("obj.func after V update=", EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.t, weightmat, lambda1, lambda2, S.v = S.v)$out, "\n")
  				}

				S.u = U_SAMPLE(Xc, V=V.t)

				lambda.u = (lambda2/p)*sum(abs(V.t))
  				U.t = glasso(S.u, rho=lambda.u, trace=0, penalize.diagonal=TRUE)$wi
  				if (!quiet) {
	 				cat("obj.func after U update=", (new.objfncval <- EVAL_OBJ_FUNC(X, class, M.update, U.t, V.t, weightmat, lambda1, lambda2, S.u = S.u)$out), "\n")
  				} else {
  					new.objfncval = EVAL_OBJ_FUNC(X, class, M.update, U.t, V.t, weightmat, lambda1, lambda2, S.u = S.u)$out
  				}

 	 		} else {

 	 				V.t = solve(S.v)
 	 				if (!quiet) {
						cat("obj.func after V update=", EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.t, weightmat, lambda1, lambda2, S.v = S.v)$out, "\n")
					}

					S.u = U_SAMPLE(Xc, V=V.t)
	 				U.t = solve(S.u)
	 				if (!quiet) {
	 					cat("obj.func after U update=", (new.objfncval <- EVAL_OBJ_FUNC(X, class, M.update, U.t, V.t, weightmat, lambda1, lambda2, S.u = S.u)$out), "\n")
	 				} else {
	 					new.objfncval = EVAL_OBJ_FUNC(X, class, M.update, U.t, V.t, weightmat, lambda1, lambda2, S.u = S.u)$out
	 				}

	 		}

		
		k.outer = k.outer + 1
		## check 
		resid = (old.objfncval - new.objfncval)/abs(init.objfncval)
		## print if not quiet

		if (!quiet) {
			cat("k.outer=", k.outer, "resid=", resid,"\n")
		}
		 if (resid < full.tol) {
		 		outer.converge=TRUE
		}

		if (k.outer > k.iter) {
		 	outer.converge=TRUE
		} 

		old.objfncval = new.objfncval

		## normalize
		out1 = sum(abs(U.t))

	  	V.inner = (V.t*out1)/r
	  	U.inner = (U.t*r)/out1

	  	Uinv.inner = solve(U.inner)
		Vinv.inner = solve(V.inner)
	
	}


return(list("Mean" = M.update, "U" = U.inner, "V" = V.inner))

}

