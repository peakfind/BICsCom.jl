@testset "Transparent boundary condition" begin
   n = 9
   k = 0.7
   
   # wave vector for the plane wave 
   kx = 0.0
   ky = sqrt(k^2 - kx^2)
   
   sq = Square([0.0, 0.0], π)
   sp = get_samplingpoints(sq, n)
   
   xt = sp.cc[1, 2n+1:3n]
   yt = sp.cc[2, 2n+1:3n]
   xb = sp.cc[1, 1:n]
   yb = sp.cc[2, 1:n]
   
   # upward plane wave
   pwu = exp.(im * kx * xt + im * ky * yt)
   pwun = im * ky * pwu

   # downward plane wave
   pwd = exp.(im * kx * xb - im * ky * yb)
   pwdn = im * ky * pwd
   
   # for odd n
   modes = collect(-(n - 1)÷2 : (n - 1)÷2)
   β = [beta_m(m, k) for m in modes]
   
   Tt = assemble_tbc(xt, β, modes)
   Tb = assemble_tbc(xb, β, modes)

   # numerical solution
   pwuN = Tt * pwu
   pwdN = Tb * pwd
   
   @test pwuN ≈ pwun
   @test pwdN ≈ pwdn
end 