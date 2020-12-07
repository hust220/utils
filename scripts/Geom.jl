# Geometry operation

using LinearAlgebra
using Statistics

function distance(c1::Vector{<:Real}, c2::Vector{<:Real})
    return sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2 + (c1[3] - c2[3])^2)
end

function angle(c1::Vector{<:Real}, c2::Vector{<:Real}, c3::Vector{<:Real})
    return angle(c1 - c2, c3 - c2)
end

function angle(vec_a::Vector{<:Real}, vec_b::Vector{<:Real})
    return acos(dot(vec_a, vec_b) / (norm(vec_a) * norm(vec_b)))
end

function dihedral(vec_a, vec_b, vec_c)
    return atan(dot(cross(cross(vec_a, vec_b), cross(vec_b, vec_c)), vec_b / norm(vec_b)),
    dot(cross(vec_a, vec_b), cross(vec_b, vec_c)))
end                

function dihedral(c1::Vector{<:Real}, c2::Vector{<:Real}, c3::Vector{<:Real}, c4::Vector{<:Real})
    return dihedral(c2 - c1, c3 - c2, c4 - c3)
end

abstract type Rotation end

struct Suppos <: Rotation
    o1::Vector{Float64}
    rot::Matrix{Float64}
    o2::Vector{Float64}
    rmsd::Float64
end

struct RotateAlong <: Rotation
    o1::Vector{Float64}
    rot::Matrix{Float64}
    o2::Vector{Float64}
end

# Align 3 points
struct Align3P <: Rotation
    o1::Vector{Float64}
    rot::Matrix{Float64}
    o2::Vector{Float64}
end

function rotate(c, rot)
    x = c[1] * rot[1,1] + c[2] * rot[2,1] + c[3] * rot[3,1]
    y = c[1] * rot[1,2] + c[2] * rot[2,2] + c[3] * rot[3,2]
    z = c[1] * rot[1,3] + c[2] * rot[2,3] + c[3] * rot[3,3]
    c[1], c[2], c[3] = x, y, z
end

function apply(r::Rotation, c)
    for i = 1:3
        c[i] -= r.o1[i]
    end
    rotate(c, r.rot)
    for i = 1:3
        c[i] += r.o2[i]
    end
end

function suppos(x_, y_)
    len = size(x_, 1)

    x = copy(x_)
    y = copy(y_)

    c1 = mean(x, dims=1)
    c2 = mean(y, dims=1)
    for i = 1:len
        for j = 1:3
            x[i,j] -= c1[j]
            y[i,j] -= c2[j]
        end
    end
    
    g = transpose(x) * y;
    u, s, v = svd(g);
    
    identity = Matrix{Float64}(I, 3, 3)
    if det(g) < 0 identity[3,3] = -1 end
    rot = u * identity * transpose(v)

    d = x * rot - y;
    rmsd = 0;
    for i = 1:len
        for j = 1:3
            rmsd += d[i,j] * d[i,j]
        end    
    end    
    rmsd = sqrt(rmsd / len)

    return Suppos([c1...], rot, [c2...], rmsd)
end

function rmsd(x, y)
    return suppos(x, y).rmsd
end

function xrot(c::Real, s::Real) return Float64[1 0 0; 0 c s; 0 -s c] end
function yrot(c::Real, s::Real) return Float64[c 0 -s; 0 1 0; s 0 c] end
function zrot(c::Real, s::Real) return Float64[c s 0; -s c 0; 0 0 1] end

function rotatealong(beg_, end_, ang_::Real)
    l = [end_[i] - beg_[i] for i in 1:3]

    r1 = sqrt(l[1] * l[1] + l[2] * l[2]);
    c1 = l[2] / r1;
    s1 = l[1] / r1;

    r2 = sqrt(l[1] * l[1] + l[2] * l[2] + l[3] * l[3]);
    c2 = l[3] / r2;
    s2 = sqrt(l[1] * l[1] + l[2] * l[2]) / r2;

    rot = Matrix{Float64}(I, 3, 3)
    if r1 != 0 rot = rot * zrot(c1, s1) end
    if r2 != 0 rot = rot * xrot(c2, s2) end
    rot = rot * zrot(cos(ang_), sin(ang_))
    if r2 != 0 rot = rot * xrot(c2, -s2) end
    if r1 != 0 rot = rot * zrot(c1, -s1) end

    return RotateAlong([beg_...], rot, [beg_...])
end

function align3p(x, y)
    x2 = x[2,:] - x[1,:]
    y2 = y[2,:] - y[1,:]

    ra1 = rotatealong(zeros(3), cross(x2, y2), angle(x2, y2))

    x3 = x[3,:] - x[1,:]
    apply(ra1, x3)
    y3 = y[3,:] - y[1,:]

    ra2 = rotatealong(zeros(3), y2, dihedral(x3, zeros(3), y2, y3))

    return Align3P(x[1,:], ra1.rot * ra2.rot, y[1,:])
end
