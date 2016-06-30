
# Linear algebra routines and wrappers


function findortho(A)
  m,n = size(A)
  return sparse(randn(n))
end

"""
Finds (an estimate for) the (small) nullspace of a large sparse matrix inductively,
by solving linear equations only
"""
function sparse_nullspace(A)
  max_iter = 100 # the maximum size of the nullspace
  m,n = size(A)
  rhs = [zeros(m);1]
  ret =0
  for idx = 1:max_iter
    ortho = findortho(A) # we guess that ortho is not in nullspace
    #println(size(ret))
    B= idx==1 ? [A;ortho'] : [A;ortho';ret']
    vv = B\rhs

    #println(norm(B*vv-rhs))
    if norm(B*vv-rhs) > 1e-8
        break
    end
    ret = idx==1? vv : [ret vv]
    rhs = [rhs;0]
  end
  return ret
end


function sparse_nullspace_dim4(A)
  B = copy(A)
  n = size(A,1)
  B=[B;sprandn(100,size(B,2), 0.2)]
  for rkdef = 0:100
    F = lufact(B[1:n+rkdef,:])
    if all(diag(F[:U]) .!= 0)
      return rkdef
    end
  end
  return 99999
  error("max iterations reached")
end

function sparse_nullspace_dim(A)
  max_iter = 100 # the maximum size of the nullspace
  m,n = size(A)
  println("Solver called on A of size $(size(A))...")
  rhs = [zeros(m);1]
  ret =0
  retn = 0
  for idx = 1:max_iter
    println("    $idx")
    #ortho = sparse(randn(n)) # we guess that ortho is not in nullspace
    ortho = zeros(n)
    ind1 = 1+round(rand()*(n-1))
    ortho[ind1]=1
    ortho=sparse(ortho)
    #println(size(ret))
    B= idx==1 ? [A;ortho'] : [A;ortho';ret']
    vv = 0
    try
      vv = B\rhs
    catch
      println("Singular")
      return retn
    end

    #println(norm(B*vv-rhs))
    if norm(B*vv-rhs) > 1e-8
        break
    end
    ret = idx==1? vv : [ret vv]
    rhs = [rhs;0]
    retn += 1
  end
  return retn
end

function first_nzidx(v)
  rs = rowvals(v)
  for (i,x) in enumerate(nonzeros(v))
    if abs(x) >0
      return rs[i]
    end
  end
  return -1
end

function myelim(A)
  m,n = size(A)
  mm = min(m,n)
  # remove empty columns
  l = [!isempty(nzrange(A,i)) for i in 1:n]
  AA = A[:,l]
  m,n = size(AA)
  mm = min(m,n)

  # create subatrices by start row
  rs = rowvals(AA)
  start_rows = [first_nzidx(AA[:,i]) for i in 1:n]

  srd = [ s => find(start_rows.==s)  for s in start_rows]

  rank_count = 0
  println(start_rows)

  for sidx=1:m
    #println(sidx)
    # eliminate all but possibly one column
    if rem(sidx,100)==0
      println(sidx)
    end
    if !haskey(srd,sidx)
      continue
    else
      rank_count +=1
      #println("ha $(size(srd[sidx]))")
    end
    b = srd[sidx]
    nb = length(b)
    # find pivot -> sparsest colum
    #println(typeof(B))
    nnzs = [nnz(AA[:,b[j]]) for j in 1:nb]
    pivotidx = indmin(nnzs)
    #println("..$(B[sidx,1])")
    #B[:,[1,j]] = B[:,[j,1]]
    # eliminate
    w = AA[:,b[pivotidx]]
    for j=1:nb
      if j!= pivotidx
        #println(B[sidx,1])
        v = AA[:,b[j]]
        v -= w*v[sidx]/w[sidx]
        AA[:,b[j]] = v
        new_sidx = first_nzidx(v)

        if new_sidx>0
          # find non-zeros
          if haskey(srd,new_sidx)
            push!(srd[new_sidx], b[j])
          else
            srd[new_sidx] = [b[j]]
          end
        end
      end
    end
  end

    rank_count
    #return srd
end
