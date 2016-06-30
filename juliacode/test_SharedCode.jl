using Base.Test


function testPermSign()
    p = [1,2,3]
    @test permSign(p) == 1
    p = [1,3,2]
    @test permSign(p) == -1
    p = [3,2,1]
    @test permSign(p) == -1
    p = [3,2,1,5,6,4]
    @test permSign(p) == -1
    p = [3,5,1,2,6,4]
    @test permSign(p) == 1
    println( "permSign : success" )
end

function testGetFString()
    @test getFString(8,10,7)=="aaaaaaaaaaaaaaaaaa"
    println( "getFString : success")
end

function testInducedPerm()
    perm=[2, 9, 10, 3, 6, 4, 5,8,7 ]
    lst = [4, 5, 6]
    #print inducedPerm(perm, lst)
    @test inducedPerm(perm, lst) == [1,3,2] #TODO
    perm=[2, 9, 10, 3, 6, 4, 5,8,7 ]
    lst = [1, 2, 5]
    #print inducedPerm(perm, lst)
    @test inducedPerm(perm, lst) == [1,3,2]
    perm=[2, 9, 10, 3, 6, 4, 5,8,7 ]
    lst = [2,3,4]
    #print inducedPerm(perm, lst)
    @test inducedPerm(perm, lst) == [3,1,2]
    println( "inducedPerm : success")
end

function testG6Code()
    G = small_graph(3)
    add_edge!(G,1,2)
    add_edge!(G,2,3)
    add_edge!(G,1,3)
    @test generate_graph6(G) == "Bw"

    s = "Bw"
    G = parse_graph6(s)
    @test num_vertices(G) == 3
    s2 = generate_graph6(G)
    @test s == s2

    s = "IiIOOC@gG"
    G = parse_graph6(s)
    s2 = generate_graph6(G)
    #println(G)
    @test s == s2

    println( "G6code : success")
end

testPermSign()

testInducedPerm()

testG6Code()
