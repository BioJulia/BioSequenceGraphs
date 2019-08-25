
function build_freq_list(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}, chunk_size::Int = 1000000) where {M<:AbstractMer}
    chunk_mers = Vector{M}()
    chunk_freqs = Vector{MerFreq{M}}()
    sizehint!(chunk_mers, chunk_size)
    output = Vector{MerFreq{M}}()
    
    read = first(range)
    lastread = last(range)
    read_sequence = LongDNASeq()
    
    while read <= lastread
        while read <= lastread && length(chunk_mers) < chunk_size
            # Collect mers into a batch
            for mer in each(M, load_sequence!(sbuf, read, read_sequence))
                push!(chunk_mers, canonical(mer))
            end
            read = read + 1
        end
        # Sort and collapse the batch into a MerList
        collapse_into_freqs!(chunk_mers, chunk_freqs)
        merge_into_sorted!(output, chunk_freqs)
        empty!(chunk_mers)
    end
    return output
end

function build_freq_list2(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}, chunk_size::Int = 1000000) where {M<:AbstractMer}
    chunk_freqs = Vector{MerFreq{M}}()
    sizehint!(chunk_freqs, chunk_size)
    output = Vector{MerFreq{M}}()
    
    read = first(range)
    lastread = last(range)
    read_sequence = LongDNASeq()
    
    while read <= lastread
        while read <= lastread && length(chunk_freqs) < chunk_size
            # Collect mers into a batch
            for mer in each(M, load_sequence!(sbuf, read, read_sequence))
                push!(chunk_freqs, MerFreq{M}(canonical(mer), 0x01))
            end
            read = read + 1
        end
        # Sort and collapse the batch of mers.
        collapse!(chunk_freqs)
        merge_into_sorted!(output, chunk_freqs)
        empty!(chunk_freqs)
    end
    return output
end

function build_freq_list3(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}, chunk_size::Int) where {M<:AbstractMer}
    @info chunk_size
    
    chunk_mers = Vector{M}(undef, chunk_size)
    chunk_freqs = Vector{MerFreq{M}}()
    sizehint!(chunk_freqs, chunk_size)
    #sizehint!(chunk_mers, chunk_size)
    output = Vector{MerFreq{M}}()
    
    #@info "mers_per_read" 
    #@info "chunk_mers" chunk_mers
    
    read = first(range)
    lastread = last(range)
    read_sequence = LongDNASeq()
    
    mergen = each(M, load_sequence!(sbuf, read, read_sequence))
    mernext = iterate(mergen)
    chunkfill = 1
    
    @inbounds while read <= lastread
        while mernext !== nothing && chunkfill <= chunk_size
            chunk_mers[chunkfill] = canonical(mernext[1])
            chunkfill = chunkfill + 1
            mernext = iterate(mergen, mernext[2])
        end
        @info "finished inner while" chunkfill mernext
        if chunkfill > chunk_size
            #@info "Finished a chunk"
            collapse_into_freqs!(chunk_mers, chunk_freqs)
            merge_into_sorted!(output, chunk_freqs)
            chunkfill = 1
        else
            read = read + 1
            if read <= lastread
                mergen = each(M, load_sequence!(sbuf, read, read_sequence))
                mernext = iterate(mergen)
            end
        end
    end
    return output
end

function build_freq_list3(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}) where {M<:AbstractMer}
    return build_freq_list3(M, sbuf, range, length(range) * (300 - 31 + 1))
end

function build_freq_list4(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}) where {M<:AbstractMer}
    chunk_mers = Vector{M}(undef, length(range) * (300 - 31 + 1))
    wi = firstindex(chunk_mers)
    read_sequence = LongDNASeq()
    @inbounds for i in range
        for mer in each(M, load_sequence!(sbuf, i, read_sequence))
            chunk_mers[wi] = canonical(mer)
            wi = wi + 1
        end
    end
    resize!(chunk_mers, wi)
    return collapse_into_freqs(chunk_mers)
end

function build_freq_list5(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}) where {M<:AbstractMer}
    chunk_mers = Vector{MerFreq{M}}(undef, length(range) * (300 - 31 + 1))
    wi = firstindex(chunk_mers)
    read_sequence = LongDNASeq()
    @inbounds for i in range
        for mer in each(M, load_sequence!(sbuf, i, read_sequence))
            chunk_mers[wi] = MerFreq{M}(canonical(mer), 0x01)
            wi = wi + 1
        end
    end
    resize!(chunk_mers, wi)
    return collapse!(chunk_mers)
end

function build_freq_list6(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}, chunk_size::Int) where {M<:AbstractMer}
    
    chunk_mers = Vector{M}(undef, chunk_size)
    chunk_freqs = Vector{MerFreq{M}}()
    sizehint!(chunk_freqs, chunk_size)
    output = Vector{MerFreq{M}}()
    
    read = first(range)
    lastread = last(range)
    read_sequence = LongDNASeq()
    
    mergen = each(M, load_sequence!(sbuf, read, read_sequence))
    mernext = iterate(mergen)
    chunkfill = 1
    
    @inbounds while read <= lastread
        while mernext !== nothing && chunkfill <= chunk_size
            chunk_mers[chunkfill] = canonical(mernext[1])
            chunkfill = chunkfill + 1
            mernext = iterate(mergen, mernext[2])
        end
        if isnothing(mernext)
            read = read + 1
            if read <= lastread
                mergen = each(M, load_sequence!(sbuf, read, read_sequence))
                mernext = iterate(mergen)
            end
        end
        
        if chunkfill > chunk_size || read > lastread
            collapse_into_freqs!(chunk_mers, chunk_freqs)
            merge_into_sorted!(output, chunk_freqs)
            chunkfill = 1
        end
    end
    return output
end

function build_freq_list6(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}) where {M<:AbstractMer}
    return build_freq_list6(M, sbuf, range, length(range) * (300 - 31 + 1))
end