
@testset "Hamming Alignment" begin
    # Test cases for hamming_align function
    # Correct signature: hamming_align(query, ref, max_error_rate, ref_search_range, max_start_pos, min_end_pos, trim_side)

    # 1. Exact match
    q = "AAAA"
    r = "TTAAAAgg"
    # Search range effectively covering everything
    res = BioDemuX.hamming_align(q, r, 0.2, 1:length(r), length(r), 1, nothing)
    @test res == (0.0, 3, 6)

    # 2. One mismatch (allowed)
    q = "AAAA"
    r = "TTAATAgg" # "AATA" vs "AAAA" -> 1 mismatch / 4 = 0.25. If max_error=0.3, should pass
    res = BioDemuX.hamming_align(q, r, 0.3, 1:length(r), length(r), 1, nothing)
    @test res == (0.25, 3, 6)

    # 3. Mismatch exceeding error rate
    res = BioDemuX.hamming_align(q, r, 0.2, 1:length(r), length(r), 1, nothing)
    @test res == (Inf, -1, -1)

    # 4. 'N' in query (wildcard)
    q = "ANNA"
    r = "TTAATAgg" # "AATA" vs "ANNA" -> match, match, match, match (N matches anything)
    res = BioDemuX.hamming_align(q, r, 0.0, 1:length(r), length(r), 1, nothing)
    @test res == (0.0, 3, 6)

    # 5. 'N' in reference (mismatch unless query is 'N')
    q = "AAAA"
    r = "TTANAAgg" # "ANAA" vs "AAAA" -> N in ref is not wildcard for A
    # If N in ref does not match A in query, it's a mismatch. 
    # Current implementation: if q_char != r_char && q_char != 'N', mismatch.
    # So 'A' != 'N' -> mismatch.
    res = BioDemuX.hamming_align(q, r, 0.0, 1:length(r), length(r), 1, nothing)
    @test res == (Inf, -1, -1)

    # 6. Boundary conditions
    q = "AAAA"
    r = "AAAA"
    res = BioDemuX.hamming_align(q, r, 0.0, 1:4, 4, 1, nothing)
    @test res == (0.0, 1, 4)

    # 7. Trimming side 3 (keep better match to the right? No, standard tie breaking)
    # hamming_align logic: if score == best_score && trim_side == 3 && j > best_start -> update
    # Case where same best score appears twice
    q = "AA"
    r = "AATAA" # Matches at 1 (AA) and 4 (AA)
    # If trim_side=3, prefer rightmost (later start)
    res = BioDemuX.hamming_align(q, r, 0.0, 1:5, 5, 1, 3)
    @test res == (0.0, 4, 5)

    # If trim_side=nothing or 5, prefer leftmost (first found)
    res_def = BioDemuX.hamming_align(q, r, 0.0, 1:5, 5, 1, nothing)
    @test res_def == (0.0, 1, 2)
end
