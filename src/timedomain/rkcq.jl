

"""
    RungeKuttaConvolutionQuadrature{T,N,NN}

T: the value type of the basis function.
N: the number of stages.
NN: N*N.

Performs a convolution quadrature on a laplaceKernel to represent an operator
in time domain using an implicit Runge-Kutta method.

laplaceKernel: function of the Laplace variable s that returns an IntegralOperator.
A, b: Coefficient matrix and vectors from the Butcher tableau.
Δt: time step.
zTransformedTermCount: Number of terms in the inverse Z-transform.
contourRadius: radius of circle used as integration contour for the inverse Z-transform.
"""
struct RungeKuttaConvolutionQuadrature{}
	timedomainKernel :: AbstractSpaceTimeOperator # function of s that returns an IntegralOperator
end
#scalartype(rkcq::RungeKuttaConvolutionQuadrature{LK}) where {LK} = Complex


# M = H*diagm(D)*invH
struct DiagonalizedMatrix{T,N,NN}
	H    :: SArray{Tuple{N,N},Complex{T},2,NN}
	invH :: SArray{Tuple{N,N},Complex{T},2,NN}
	D    :: SVector{N,Complex{T}}
end

# M = H*diagm(D)*invH
function diagonalizedmatrix(M :: SArray{Tuple{N,N},Complex{T},2,NN}) where {T,N,NN}
	ef = eigen(Array{Complex{T},2}(M));

	efValues = SVector{N,Complex{T}}(ef.values)  :: SVector{N,Complex{T}};
	efVectors = SArray{Tuple{N,N},Complex{T},2,NN}(ef.vectors) :: SArray{Tuple{N,N},Complex{T},2,NN};
	return DiagonalizedMatrix(efVectors, inv(efVectors), efValues);
end

function assemble(rkcq :: RungeKuttaConvolutionQuadrature,
                  testfns :: SpaceTimeBasis,
                  trialfns :: SpaceTimeBasis)

	@warn "staged assemble of the left-hand side"
	sol = rkcq.timedomainKernel.speed_of_light
	numdiffweak = rkcq.timedomainKernel.ws_diffs
	numdiffhyper = rkcq.timedomainKernel.hs_diffs

	@show numdiffhyper
	@show numdiffweak

	@info "converting time-domain kernel to Laplace-domain kernel"
	HsEFIO(s::T) where {T} = MWSingleLayer3D(s/sol, T(0), T(sol))
	WsEFIO(s::T) where {T} = MWSingleLayer3D(s/sol, 1/T(sol), T(0))
	sWsEFIO(s::T) where {T} = MWSingleLayer3D(s/sol, s/T(sol), T(0))
	ssWsEFIO(s::T) where {T} = MWSingleLayer3D(s/sol, s*s/T(sol), T(0))
	WsKernel = WsEFIO
	HsKernel = HsEFIO
	sWsKernel = sWsEFIO
	ssWsKernel = ssWsEFIO
	
	A = testfns.time.A
	b = testfns.time.b
	Δt = testfns.time.Δt
	Q = testfns.time.zTransformedTermCount
	rho = testfns.time.contourRadius
	p = length(b) # stage count

	test_spatial_basis  = testfns.space
	trial_spatial_basis = trialfns.space

	#build quasi Helmholtz Projectors
	Γ = trial_spatial_basis.geo

    ∂Γ = boundary(Γ)

    setminus(A,B) = submesh(!in(B), A)

    edges = setminus(skeleton(Γ,1), ∂Γ)
    verts = setminus(skeleton(Γ,0), skeleton(∂Γ,0))

    Σ = Matrix(connectivity(Γ, edges, sign))
    Λ = Matrix(connectivity(verts, edges, sign))

	#adapt Projector size to stage count "Large Time Step and DC Stable TD-EFIE Discretized With Implicit Runge–Kutta Methods"
	I = LinearAlgebra.I
    PΣ = Σ * pinv(Σ'*Σ) * Σ'
    PΛH = I - PΣ

	Ip = diagm(@SVector ones(p))
	@show Ip
	@show size(PΣ)

    ℙΣ = kron(PΣ,Ip)
    ℙΛH = kron(I - PΣ,Ip)


    #MR = γ * PΣ + PΛH
    #ML = PΣ + 1/γ * PΛH

	# Compute the Z transformed sequence.
	# Assume that the operator applied on the conjugate of s is the same as the
	# conjugate of the operator applied on s,
	# so that only half of the values are computed
	Qmax = Q>>1+1
	M = numfunctions(test_spatial_basis)
	N = numfunctions(trial_spatial_basis)
	Tz = ComplexF64
	#Tz = promote_type(scalartype(rkcq), scalartype(testfns), scalartype(trialfns))
	#Tz = promote_type(scalartype(testfns), scalartype(trialfns))

	Zz = Vector{Array{Tz,2}}(undef,Qmax)
	MTsM = Vector{Array{Tz,2}}(undef,Qmax)
	MThM = Vector{Array{Tz,2}}(undef,Qmax)
	MDTsM = Vector{Array{Tz,2}}(undef,Qmax)
	MDTsM2 = Vector{Array{Tz,2}}(undef,Qmax)
	MD2TsM = Vector{Array{Tz,2}}(undef,Qmax)
	smat = Vector{Array{Tz,2}}(undef,Qmax)

	blocksEigenvalues_weakly = Vector{Array{Tz,2}}(undef,p)
	blocksEigenvalues_sweakly = Vector{Array{Tz,2}}(undef,p)
	blocksEigenvalues_ssweakly = Vector{Array{Tz,2}}(undef,p)
	blocksEigenvalues_hyper = Vector{Array{Tz,2}}(undef,p)
	tmpDiag_weakly = Vector{Tz}(undef,p)
	tmpDiag_sweakly = Vector{Tz}(undef,p)
	tmpDiag_ssweakly = Vector{Tz}(undef,p)
	tmpDiag_hyper = Vector{Tz}(undef,p)
	for q = 0:Qmax-1
		# Build a temporary matrix for each eigenvalue
		s = laplace_to_z(rho, q, Q, Δt, A, b)
		sFactorized = diagonalizedmatrix(s)
		for (i,sD) in enumerate(sFactorized.D)
			blocksEigenvalues_weakly[i] = assemble((WsKernel(sD)), test_spatial_basis, trial_spatial_basis)
			blocksEigenvalues_sweakly[i] = assemble((sWsKernel(sD)), test_spatial_basis, trial_spatial_basis)
			blocksEigenvalues_ssweakly[i] = assemble((ssWsKernel(sD)), test_spatial_basis, trial_spatial_basis)
			blocksEigenvalues_hyper[i] = assemble((HsKernel(sD)), test_spatial_basis, trial_spatial_basis)
		end

		# Compute the Z transformed matrix by block
		D = diagm(sFactorized.D)

		Zz[q+1] = zeros(Tz, M*p, N*p)
		MTsM[q+1] = zeros(Tz, M*p, N*p)
		MThM[q+1] = zeros(Tz, M*p, N*p)
		MDTsM[q+1] = zeros(Tz, M*p, N*p)
		MDTsM2[q+1] = zeros(Tz, M*p, N*p)
		MD2TsM[q+1] = zeros(Tz, M*p, N*p)
		smat[q+1] = zeros(Tz, M*p, N*p)

		for m = 1:M
			for n = 1:N
				for i = 1:p
					tmpDiag_weakly[i] = blocksEigenvalues_weakly[i][m,n]
					tmpDiag_sweakly[i] = blocksEigenvalues_sweakly[i][m,n]
					tmpDiag_ssweakly[i] = blocksEigenvalues_ssweakly[i][m,n]
					tmpDiag_hyper[i] = blocksEigenvalues_hyper[i][m,n]
				end
				D_weakly = SVector{p,Tz}(tmpDiag_weakly)
				D_sweakly = SVector{p,Tz}(tmpDiag_sweakly)
				D_ssweakly = SVector{p,Tz}(tmpDiag_ssweakly)
				D_hyper = SVector{p,Tz}(tmpDiag_hyper)

				MTsM[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D_weakly) * sFactorized.invH
				MThM[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D_hyper) * sFactorized.invH
				MDTsM[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D_sweakly) * sFactorized.invH
				MDTsM2[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D_sweakly) * sFactorized.invH
				MD2TsM[q+1][(m-1)*p.+(1:p),(n-1)*p.+(1:p)] = sFactorized.H * diagm(D_ssweakly) * sFactorized.invH
			end
		end
		a = 1.0
		smat[q+1] = D
		Zz[q+1] = 1/a*ℙΛH*MTsM[q+1]*ℙΛH + a*ℙΣ*MThM[q+1]*ℙΣ + a/sol^2*ℙΣ*MD2TsM[q+1]*ℙΣ + 1/sol*ℙΣ*MDTsM[q+1]*ℙΛH + 1/sol*ℙΛH*MDTsM[q+1]*ℙΣ
	end

	# return the inverse Z transform
	kmax = Q
	T = real(Tz)
	Z = zeros(T, M*p, N*p, kmax)
	for q = 0:kmax-1
		Z[:,:,q+1] = real_inverse_z_transform(q, rho, Q, Zz)
	end
	ZC = ConvolutionOperators.DenseConvOp(Z)
	return ZC

end
