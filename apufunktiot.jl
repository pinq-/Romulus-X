# Funktio joka muuttaa asteet radiaaneiksi
function func_radians(kulma::Float64)
    kulma * (pi/180.0)
end

# Funktio joka muuttaa radiaanit asteiksi
function func_asteet(kulma::Float64)
    kulma * (180.0/pi)
end

# Funktio joka muuttaa polaarikoordinaatin carteesiseen koordinaatistoon
function func_cartesis(R::Float64, Θ::Float64)
    x::Float64 = R * cos(Θ)
    y::Float64 = R * sin(Θ)
    return x,y
end

# Funktio joka muuttaa karteesiset koordinaatit polaariseen koordinaatistoon
function func_polar(X::Float64, Y::Float64)
    R = sqrt(X^2 + Y^2)
    Θ = atan(Y/X)
    return R,Θ
end

# Laskee rotaatiomatriisin eulerin parametreille
function func_rotation_matrix(e1::Float64,e2::Float64,e3::Float64,e4::Float64)
    Rotation_matrix_1 = Array(Float64, 3,3)

#=    # Wikipedia
    Rotation_matrix_1[1,1]=e1^2 + e2^2 - e3^2 - e4^2
    Rotation_matrix_1[1,2]=2*(e2*e3 - e1*e4)
    Rotation_matrix_1[1,3]=2*(e2*e4 + e1*e3)
    Rotation_matrix_1[2,1]=2*(e2*e3 + e1*e4)
    Rotation_matrix_1[2,2]=e1^2 - e2^2 + e3^2 - e4^2
    Rotation_matrix_1[2,3]=2*(e3*e4 - e1*e2)
    Rotation_matrix_1[3,1]=2*(e2*e4 - e1*e3)
    Rotation_matrix_1[3,2]=2*(e3*e4 + e1*e2)
    Rotation_matrix_1[3,3]=e1^2 - e2^2 - e3^2 + e4^2
=#
    # Mathworld.Wolfram
    Rotation_matrix_1[1,1]=e1^2 + e2^2 - e3^2 - e4^2
    Rotation_matrix_1[1,2]=2*(e2*e3 + e1*e4)
    Rotation_matrix_1[1,3]=2*(e2*e4 - e1*e3)
    Rotation_matrix_1[2,1]=2*(e2*e3 - e1*e4)
    Rotation_matrix_1[2,2]=e1^2 - e2^2 + e3^2 - e4^2
    Rotation_matrix_1[2,3]=2*(e3*e4 + e1*e2)
    Rotation_matrix_1[3,1]=2*(e2*e4 + e1*e3)
    Rotation_matrix_1[3,2]=2*(e3*e4 - e1*e2)
    Rotation_matrix_1[3,3]=e1^2 - e2^2 - e3^2 + e4^2

#=    # arizona.edu
    Rotation_matrix_1[1,1]=2*(e1^2 + e2^2 - 1/2)
    Rotation_matrix_1[1,2]=2*(e2*e3 - e1*e4)
    Rotation_matrix_1[1,3]=2*(e2*e4 + e1*e3)

    Rotation_matrix_1[2,1]=2*(e2*e3 + e1*e4)
    Rotation_matrix_1[2,2]=2*(e1^2 + e3^2 - 1/2)
    Rotation_matrix_1[2,3]=2*(e3*e4 - e1*e2)

    Rotation_matrix_1[3,1]=2*(e2*e4 - e1*e3)
    Rotation_matrix_1[3,2]=2*(e3*e4 + e1*e2)
    Rotation_matrix_1[3,3]=2*(e1^2 + e4^2 - 1/2)
=#
    return Rotation_matrix_1
end

function func_rotation_matrix!(Rotation_matrix_1::Array,e1::Float64,e2::Float64,e3::Float64,e4::Float64)
    Rotation_matrix_1[1,1]=e1^2 + e2^2 - e3^2 - e4^2
    Rotation_matrix_1[1,2]=2*(e2*e3 - e1*e4)
    Rotation_matrix_1[1,3]=2*(e2*e4 + e1*e3)
    Rotation_matrix_1[2,1]=2*(e2*e3 + e1*e4)
    Rotation_matrix_1[2,2]=e1^2 - e2^2 + e3^2 - e4^2
    Rotation_matrix_1[2,3]=2*(e3*e4 - e1*e2)
    Rotation_matrix_1[3,1]=2*(e2*e4 - e1*e3)
    Rotation_matrix_1[3,2]=2*(e3*e4 + e1*e2)
    Rotation_matrix_1[3,3]=e1^2 - e2^2 - e3^2 + e4^2

    return nothing
end
