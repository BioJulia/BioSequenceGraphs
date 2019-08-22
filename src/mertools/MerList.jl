struct MerFreq{K}
    mer::DNAMer{K}
    count::UInt8
end

Base.isless(x::MerFreq{K}, y::MerFreq{K}) where K = x.mer < y.mer
Base.:(>)(x::MerFreq{K}, y::MerFreq{K}) where K = x.mer > y.mer
Base.:(==)(x::MerFreq{K}, y::MerFreq{K}) where K = x.mer == y.mer

function merge(x::MerFreq{K}, y::MerFreq{K}) where K
    newcount = min(UInt16(x.count) + UInt16(y.count), typemax(x.count))
    return MerFreq{K}(x.mer, newcount)
end

struct MerFreqList{K}
    kmers::Vector{MerFreq{K}}
end


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
a = _collapse_into_freqs!(v, MerFreq{4}[])
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
b = _collapse_into_freqs!(v2, MerFreq{4}[])
b′ = copy(b)

merge_into!(a′, b′)

"""
    merge_into!(a::Vector{MerFreq{K}}, b::Vector{MerFreq{K}}) where {K}

Collapse the 
"""
function merge_into!(a::Vector{MerFreq{K}}, b::Vector{MerFreq{K}}) where {K}
    a_i = firstindex(a)
    a_end = lastindex(a) + 1
    b_i = b_i2 = firstindex(b)
    b_end = lastindex(b) + 1
    
    # Merge, accumulating counts on `a`.
    @info "Merging, accumulating counts on a"
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
    
    @debug "Shrink b to size of remaining contents"
    # Shrink `b` to the size of the remaining contents.
    
    resize!(b, b_i2 - 1)
    
    # Expand `a` to allow the insertion of unique values in `b`.
    @debug "Expand a to allow insertion of unique values in b"
    
    oldsize = length(a)
    resize!(a, oldsize + length(b))
    r_a = oldsize
    
    @debug "Expanded a" a b
    
    # Merge-sort from the bottom into `a`.
    wr_a = lastindex(a)
    rend_a = firstindex(a)
    
    r_b = lastindex(b)
    r_end_b = firstindex(b)
    @debug "Merge sort from the bottom into `a`"
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
end

function build_freq_list(::Type{DNAMer{K}}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}, batch_size::Int = 1000000) where {K}
    chunk_mers = Vector{DNAMer{K}}()
    chunk_freqs = Vector{MerFreq{K}}()
    sizehint!(batch_mers, batch_size)
    
    read = first(range)
    lastread = last(range)
    read_sequence = LongDNASeq()
    
    while read <= lastread
        while read <= lastread && length(batch_mers) < batch_size
            # Collect mers into a batch
            for mer in each(DNAMer{K}, load_sequence!(bufds, read, read_sequence))
                push!(batch_mers, canonical(mer))
            end
        end
        @info string("Collected a chunk of ", batch_size, " mers")
        # Sort and collapse the batch into a MerList
        sort!(batch_mers)
        collapse_into_freqs!(chunk_mers, chunk_freqs)
        
    end
end

function collapse_into_freqs!(mers::Vector{DNAMer{K}}, result::Vector{MerFreq{K}}) where {K}
    wi = 1
    ri = 1
    stop = lastindex(mers) + 1
    empty!(result)
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mers[ri] == mers[wi]
            ci = ci + one(UInt16)
        end
        push!(result, MerFreq{K}(mers[wi], ci))
        wi = ri
    end
    return result
end

