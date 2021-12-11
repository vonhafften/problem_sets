# Problem Set 3
# Computational Economics
# Taught by JF Houde

# Alex von Hafften


# IV regression function
function ivreg(mY::Array{Float64}, mVariables::Array{Float64}, mInstruments::Array{Float64}, mWeight::Array{Float64})  
    mInvXZX=inv((mVariables'* mInstruments)*mWeight*(mInstruments'* mVariables))
    if mInvXZX==0
        vIV=fill!(zeros(size(mVariables)[2], 1), NaN)
    else 
        vIV=mInvXZX*(mVariables'mInstruments)*mWeight*(mInstruments'mY)
    end
    return vIV
end

function idio_value(vParam, t, aProductID)
    rowid = aProductID[t]
    mMu = zeros(length(rowid), Sim)

    for i = 1:length(vParam)
        mMu += vParam[i]*aZ[t]
    end

    return exp.(mMu)
end

function demand(mMu, vDelta, t, vParam)
    
    eV = exp.(vDelta[rowid]) .* mMu
    mS = eV ./ (1 .+ sum(eV, dims = 1))
    vShat = mean(mS, dims = 2)

    return vShat 
    
end


function inverse(aDelta, vParam, eps1, eps, vDelta0, aProductID, vShare)
    vDelta = vDelta0
    mJacobian = 1
    maxit = 10
    vIT = zeros(T)
    maxT = T

    # Parallelize the inversion across markets. When possible use Nb processors = T (31)
    # need to parallelize
    for t = 1
        # Pre-compute the heterogeneity component of utility (indepedent of delta)
        mMu = idio_value(vParam, t, aProductID)
        
        rowid = aProductID[t]
        vIT[t] = 0
        f = 1000

        # Evalute the demand without the jacobian matrix if the norm is larger than 1 
        while true
            if norm(f) > eps1
                vShat = demand(mMu, vDelta, t, vParam)
                f = log.(vShare[rowid]).-log.(vShat) # Zero function: Log difference in shares
                vDelta[rowid] = vDelta[rowid] + f # Contraction mapping step
            else
                # need to change to newton step                 
                break
            end
            vIT[t]+=1
            if t==1
                println("t = ",t," it ",vIT[t]," norm : ",norm(f));
            end
            if vIT[t]>maxit
                break
            end
        end
    end
end







