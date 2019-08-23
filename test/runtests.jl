module TestBioSequenceGraphs

using BioSequenceGraphs, BioSequences, Test

import BioSequenceGraphs.MerFreq

# write your own tests here
@testset "MerFreq" begin
    v = DNAMer{4}[
        mer"AAAA", mer"AAAA", mer"AAAA", mer"ATAG",
        mer"GGGG", mer"GGGG", mer"GGGT", mer"AGGT",
        mer"AGGT", mer"AGGT", mer"AGGT"
    ]
    
    v2 = DNAMer{4}[
        mer"AAAA", mer"AAAA", mer"AAAA", mer"ATAG",
        mer"GGGG", mer"GGGG", mer"GGGT", mer"AGGT",
        mer"ACGT", mer"AAGT", mer"AAGT"
    ]
    
    @test collapse_into_freqs(v) == [
        MerFreq{DNAMer{4}}(mer"AAAA", 3),
        MerFreq{DNAMer{4}}(mer"AGGT", 4),
        MerFreq{DNAMer{4}}(mer"ATAG", 1),
        MerFreq{DNAMer{4}}(mer"GGGG", 2),
        MerFreq{DNAMer{4}}(mer"GGGT", 1)
    ]
    
    @test collapse_into_freqs(v2) == [
        MerFreq{DNAMer{4}}(mer"AAAA", 3),
        MerFreq{DNAMer{4}}(mer"AAGT", 2),
        MerFreq{DNAMer{4}}(mer"ACGT", 1),
        MerFreq{DNAMer{4}}(mer"AGGT", 1),
        MerFreq{DNAMer{4}}(mer"ATAG", 1),
        MerFreq{DNAMer{4}}(mer"GGGG", 2),
        MerFreq{DNAMer{4}}(mer"GGGT", 1)
    ]
end

end # module