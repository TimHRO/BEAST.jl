


"""
    StagedTimeStep{T,N,NN}

T: the value type of the basis function.
N: the number of stages.
NN: the number of stages squared
Each time step has intermediary stages given by the vertor c in a Butcher tableau (A,b,c)
"""
struct StagedTimeStep{T, N, NN}
	Δt :: T
	Nt :: Int
	c  :: SVector{N,T}
	A  :: SArray{Tuple{N,N},T,2,NN}
	b  :: SVector{N,T}
	zTransformedTermCount :: Int
	contourRadius         :: T
end

scalartype(sts :: StagedTimeStep{T, N, NN}) where {T, N, NN} = T
temporalbasis(sts :: StagedTimeStep{T, N, NN}) where {T, N, NN} = timebasisdelta(sts.Δt, sts.Nt)

numfunctions(s::StagedTimeStep) = s.Nt

numstages(s) = 1
numstages(s::StagedTimeStep) = size(s.c,1)
