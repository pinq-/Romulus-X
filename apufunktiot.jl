# Funktio joka muuttaa asteet radiaaneiksi
function func_radians(kulma::Float64)
    kulma * (pi/180.0)
end

# Funktio joka muuttaa radiaanit asteiksi
function func_asteet(kulma::Float64)
    kulma * (180.0/pi)
end

# Funktio joka muuttaa polaarikoordinaatin carteesiseen koordinaatistoon
function func_cartesis(R::Float64,Θ::Float64)
    x::Float64 = R * cos(Θ)
    y::Float64 = R * sin(Θ)
    return x,y
end

#=function func_cartesis_X(R::Float64,Θ::Float64)
    x::Float64 = R * cos(Θ)
    return x
end=#
