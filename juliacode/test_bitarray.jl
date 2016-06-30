
function ATest()
A=BitMatrix(100,100)

p=false
for i=1:100
  for j=1:100
    A[i,j] = true
    p = A[j,i]
  end
end

end


function BTest()
A = Matrix{Bool}(100,100)
p=false
for i=1:100
  for j=1:100
    A[i,j] = true
    p = A[j,i]
  end
end

end

A= wheel_graph(12)

function listtest()
  count=0
  for (u,v) in edges(A)
    count +=1
  end
end

function cortest()
  count=0
  for (u,v) in edges22(A)
    count +=1
  end
end



function run1000k(f)
  for k=1:100000
    f()
  end
end
