#########################################################################
#### MATLDA_GRID: function fit the penalized matrix normal        #######
#########################################################################
###		X: an r x p x N array of training data
###		class: N-vector of trainig class labels {1,...,C}
###		lambda1: vector of nonneg candidate tps for lam1
###		lambda2: vector of nonneg candidate tps for lam2
###		quiet: should print iterations and obj fnc vals?
###		Xval: r x p x Ntest array of test data
###		classval: Ntest-vector of test class data labels
###		k.iter: max iterations for full algorithm
###		cov.tol: tolerance for glasso algorithms
###		m.tol: tolerance for AMA algorithms
###		full.tol: tolerance for full algorithm
####################################################################
###		Val: misclassification rates on test set
###		Mean: r x p x C array of estiamted class means
###		U:  estimated U
###		V:  estimated V
###		pi.list:  marginal class probabilities from training data
####################################################################

MatLDA = function(X, class, lambda1, lambda2, quiet=TRUE,
	Xval=NULL, classval= NULL,
	k.iter = 100, cov.tol=1e-5, m.tol=1e-5, full.tol=1e-6){

	## sort classes so n1 > .... > nC
	nc.orig = count(class)[,2]
	N = sum(nc.orig)

	class.orig = class
	class.val.orig =  classval

	out = rank(-nc.orig, ties.method = "first")
	class = rep(0, N)
	for (gg in 1:N) {
		class[gg] = out[class.orig[gg]]
	}

	nc = count(class)[,2]

	C = length(nc)
	K = (C*(C-1))/2
	N = sum(nc)
	p = dim(X)[2]
	r = dim(X)[1]
	pihat = nc/N

	if (r < p) {
		warning("Row dimension is greater than column dimension; consider transposing data.")
	}
	if(!is.null(Xval)){
		classval = rep(0, dim(Xval)[3])
		for (gg in 1:dim(Xval)[3]) {
			classval[gg] = out[class.val.orig[gg]]
		}
	}


	## initialize U and V at diagonal
	out = MN_MLE(X, class)
	V.unscaled = diag(diag(out$V))
	U.unscaled =  diag(diag(out$U))
	out1 = sum(abs(U.unscaled))
	UInit = (r/out1)*U.unscaled
	VInit = (out1/r)*V.unscaled

	Vinv = solve(VInit)
	Uinv = solve(UInit)

	## rho.max is maximum value for rho to ensure convergence
	## rho < rho.max
	rho.max = (4*min(diag(UInit))*min(diag(VInit))*min(nc))/(C*N)
	rho = rho.max/10

	## get marginal statistics
	mean.array = array(0, dim=c(r,p,C))
	for (k in 1:C) {
		mean.array[,,k] = apply(X[,,class==k], c(1,2), mean)
	}

	## compute weight matrix
	weightmat = matrix(0, nrow=K*r, ncol=p)
	for (c in 1:(C-1)) {
		for (cc in (c+1):C) {
		k = (C - 1)*(c-1) + (cc-1) - (c*(c-1)/2)
		weightmat[((k-1)*r + 1):(k*r), ] = 1/abs(mean.array[,,c] - mean.array[,,cc])
		}
	}

	init.objfncval = EVAL_OBJ_FUNC(X, class, mean.array, UInit, VInit,  weightmat, lambda1, lambda2)$out

	if (!is.null(Xval)) {
				Val = matrix(0, nrow=length(lambda1), ncol=length(lambda2))
			}

	outer.converge = FALSE
	meanstacked = STACKED(mean.array)
	M.update = mean.array

	for (j in 1:length(lambda1)) {

		if (lambda1[j]!= 0) {
		## update means  ######
				D = matrix(0, nrow=K*r, ncol=p)
				out = M_Update_Rcpp(Mean=meanstacked,
					 D=D, Uinv=Uinv, Vinv=Vinv, nc = nc, rho=rho, weightmat = weightmat,
				  lambda=lambda1[j], C=C, r=r, mtol = m.tol)

				D = out$D

				M = THRESH_MEAN(out$M, out$G, C)

				M.update = array(0, dim=c(r, p, C))
				for (c in 1:C) {
					M.update[,,c] = M[((c-1)*r +1):(c*r),]
				}

		} else {
				M.update = mean.array
		}

		## center and compute path of covariances
		Xc = CENTER_X(X, class, M.update)
       	S.v = V_SAMPLE(X.c=Xc, U = UInit)
  		V.temp = glassopath(S.v, rholist=lambda2, trace=0, thr=cov.tol, penalize.diagonal=TRUE)$wi

  		for (i in 1:length(lambda2)) {
  			if (lambda2[i]!=0) {
  				if (!quiet) {
  					cat("obj.func after V update=", EVAL_OBJ_FUNC(X, class, M.update, UInit, V.temp[,,i], weightmat, lambda1[j], lambda2[i], S.v = S.v)$out, "\n")
  				}

  				lambda.u = (lambda2[i]*sum(abs(V.temp[,,i]))/p)
  				S.u = U_SAMPLE(X.c = Xc, V = V.temp[,,i])
  				U.temp = glasso(S.u, rho = lambda.u, thr=cov.tol, trace=0, penalize.diagonal=TRUE)$wi

  				if (!quiet) {
					cat("obj.func after U update=", EVAL_OBJ_FUNC(X, class, M.update, U.temp, V.temp[,,i], weightmat, lambda1[j], lambda2[i], S.u = S.u)$out, "\n")
				}
				out.inner = GENERAL_ITER(X, class, lambda1 = lambda1[j], lambda2 = lambda2[i],
				quiet=quiet, UInit = UInit, VInit = VInit, U=U.temp, V=V.temp[,,i],
						  meanstacked = meanstacked, mean.array= mean.array, M.update = M.update, D=D, weightmat=weightmat,
						  full.tol = full.tol, cov.tol=cov.tol, m.tol=m.tol, init.objfncval=init.objfncval, k.iter=k.iter)

			} else {
				V.temp = solve(S.v)
				S.u = U_SAMPLE(X.c = Xc, V = V.temp)
				U.temp = solve(S.u)

				out.inner = GENERAL_ITER(X, class, lambda1 = lambda1[j], lambda2 = lambda2[i],
				quiet=quiet, UInit = UInit, VInit = VInit, U=U.temp, V=V.temp,
						 meanstacked = meanstacked, mean.array= mean.array, M.update = M.update, D=D, weightmat=weightmat,
						  full.tol = full.tol, cov.tol=cov.tol, m.tol=m.tol, init.objfncval=init.objfncval, k.iter=k.iter)
			}

  			## If Val set, compute misclassification rate to pass
  			if (!is.null(Xval)) {
				Val[j,i] = MN_CLASSIFICATION(X=Xval, class=classval, out.inner$Mean, out.inner$U, out.inner$V, pi.list=pihat, C=C)$misclass
			} else {
				Val = NULL
			}

			MeanTemp = out.inner$Mean
			for (kk in 1:C) {
				ind = unique(class.orig[which(class==kk)])
				MeanTemp[,,ind] = out.inner$Mean[,,kk]
			}
		}
	}

		out = list("Mean"= MeanTemp, "U"= out.inner$U, "V" = out.inner$V, "Val" = Val, "pi.list" = nc.orig/N)
		class(out) = "MN"

return(out)

}

