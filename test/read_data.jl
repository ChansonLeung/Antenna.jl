using antenna


(vecx, vecy, vecz, point_p) = include("D:/ChansonDocument/project/2112_7_机翼工况计算/yang/read_data.jl")
point_p = map(x->(x=x[1], y=x[2], z=x[3]), point_p)
I = antenna.Iₛ.(point_p, 0,0, k)
point = map(p -> [p.x,p.y,p.z] , point_p )
return (point,  angle.(I))