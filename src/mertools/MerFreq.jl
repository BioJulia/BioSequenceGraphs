
"""
    Kmer frequency type
"""
struct MerFreq{M<:AbstractMer}
    mer::M
    count::UInt8
end

const DNAMerFreq{K} = MerFreq{DNAMer{K}}
const RNAMerFreq{K} = MerFreq{RNAMer{K}}

Base.isless(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer} = x.mer < y.mer
Base.:(>)(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer} = x.mer > y.mer
Base.:(==)(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer} = x.mer == y.mer

function merge(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer}
    newcount = min(UInt16(x.count) + UInt16(y.count), typemax(x.count))
    return MerFreq{M}(x.mer, newcount)
end

#=
v = DNAMer{4}[]
push!(v, mer"AAAA")
push!(v, mer"AAAA")
push!(v, mer"AAAA")
push!(v, mer"ATAG")
push!(v, mer"GGGG")
push!(v, mer"GGGG")
push!(v, mer"GGGT")
push!(v, mer"AGGT")
push!(v, mer"AGGT")
push!(v, mer"AGGT")
push!(v, mer"AGGT")
sort!(v)
a = collapse_into_freqs!(v, DNAMerFreq{4}[])
a′ = copy(a)

v2 = DNAMer{4}[]
push!(v2, mer"AAAA")
push!(v2, mer"AAAA")
push!(v2, mer"AAAA")
push!(v2, mer"ATAG")
push!(v2, mer"GGGG")
push!(v2, mer"GGGG")
push!(v2, mer"GGGT")
push!(v2, mer"AGGT")
push!(v2, mer"ACGT")
push!(v2, mer"AAGT")
push!(v2, mer"AAGT")
sort!(v2)
b = collapse_into_freqs!(v2, DNAMerFreq{4}[])
b′ = copy(b)

merge_into!(a′, b′)

z = MerFreq{4}[]

merge_into!(z, a′)
=#

"""
    merge_into!(a::Vector{MerFreq{M}}, b::Vector{MerFreq{M}}) where {M<:AbstractMer}

Collapse the 
"""
function merge_into!(a::Vector{MerFreq{M}}, b::Vector{MerFreq{M}}) where {M<:AbstractMer}
    a_i = firstindex(a)
    a_end = lastindex(a) + 1
    b_i = b_i2 = firstindex(b)
    b_end = lastindex(b) + 1
    
    # Merge, accumulating counts on `a`.
    while b_i < b_end
        while a_i < a_end && a[a_i] < b[b_i]
            a_i = a_i + 1
        end
        if a_i < a_end && a[a_i] == b[b_i]
            # Combine entries
            a[a_i] = merge(a[a_i], b[b_i])
            a_i = a_i + 1
            b_i = b_i + 1
        end
        while b_i < b_end && (a_i == a_end || b[b_i] < a[a_i])
            b[b_i2] = b[b_i]
            b_i2 = b_i2 + 1
            b_i = b_i + 1
        end
        
    end
    # Shrink `b` to the size of the remaining contents.
    resize!(b, b_i2 - 1)
    
    # Expand `a` to allow the insertion of unique values in `b`.
    oldsize = length(a)
    resize!(a, oldsize + length(b))
    r_a = oldsize
    
    # Merge-sort from the bottom into `a`.
    wr_a = lastindex(a)
    rend_a = firstindex(a)
    r_b = lastindex(b)
    r_end_b = firstindex(b)
    while wr_a >= rend_a
        if r_b >= r_end_b && (r_a < rend_a || b[r_b] > a[r_a])
            a[wr_a] = b[r_b]
            r_b = r_b - 1
        else
            a[wr_a] = a[r_a]
            r_a = r_a - 1
        end
        wr_a = wr_a - 1
    end
    empty!(b)
    return a
end

function collapse_into_freqs!(mers::Vector{M}, result::Vector{MerFreq{M}}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    stop = lastindex(mers) + 1
    empty!(result)
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mers[ri] == mers[wi]
            ci = ci + one(UInt16)
        end
        push!(result, MerFreq{M}(mers[wi], ci))
        wi = ri
    end
    return result
end

function collapse_into_freqs(mers::Vector{M}) where {M<:AbstractMer}
    return collapse_into_freqs!(mers, Vector{MerFreq{M}}())
end
