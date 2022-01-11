using PlotlyJS


vec_p = [[10,0,0], [0,0,10]]
vec_coord = [[1,0,0], [0,0,1]]

rotate45 = [1 0 0 
            0 1/sqrt(2) -1/sqrt(2)
             0 1/sqrt(2) 1/sqrt(2)]
rotate0 = Matrix(1.0I, (3,3))

# vec_p_rotate = [rotate0] .* vec_p
# vec_coord_rotate = [rotate0] .* vec_coord
vec_p_rotate = [rotate45] .* vec_p
vec_coord_rotate = [rotate45] .* vec_coord

x = [i[1] for i in vec_p_rotate] 
y = [i[2] for i in vec_p_rotate] 
z = [i[3] for i in vec_p_rotate] 

u = [i[1] for i in vec_coord_rotate] 
v = [i[2] for i in vec_coord_rotate] 
w = [i[3] for i in vec_coord_rotate] 

plot(
    cone(
    x=x,
    y=y,
    z=z,
    u=u,
    v=v,
    w=w,
    sizemode="absolute",
    sizeref=1
    )
)
