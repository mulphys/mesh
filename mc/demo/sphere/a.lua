-- Lua script.
p=tetview:new()
p:load_mesh("sphere.366")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
