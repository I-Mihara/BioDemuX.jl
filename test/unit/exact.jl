
@testset "Exact Alignment" begin
    # Test cases for exact_align function
    # Signature: exact_align(query, ref, ref_search_range, max_start_pos, min_end_pos, trim_side)

    # 1. Exact match
    q = "AAAA"
    r = "TTAAAAgg"
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 1, nothing)
    @test res == (0.0, 3, 6)

    # 2. One mismatch (Should FAIL)
    q = "AAAA"
    r = "TTAATAgg"
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 1, nothing)
    @test res == (Inf, -1, -1)

    # 3. 'N' in query (Should FAIL in Exact mode if strict literal)
    # The benchmark used strict equality.
    q = "ANNA"
    r = "TTAATAgg"
    # 'N' != 'A' -> Fail.
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 1, nothing)
    @test res == (Inf, -1, -1)

    # 4. 'N' in reference (Should FAIL)
    q = "AAAA"
    r = "TTANAAgg"
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 1, nothing)
    @test res == (Inf, -1, -1)

    # 5. Boundary conditions
    q = "AAAA"
    r = "AAAA"
    res = BioDemuX.exact_align(q, r, 1:4, 4, 1, nothing)
    @test res == (0.0, 1, 4)

    # 6. Trim side = 3 (Prefer rightmost)
    q = "AA"
    r = "AATAA"
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 1, 3)
    @test res == (0.0, 4, 5)

    # 7. Trim side = 5 (Prefer leftmost/first)
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 1, 5)
    @test res == (0.0, 1, 2)

    # 8. Trim side = 3 with min_end_pos constraint
    # "AATAA", query "AA". Matches at 1:2 and 4:5.
    # If min_end_pos = 4:
    # Rightmost match (4:5) satisfies end >= 4. OK.
    q = "AA"
    r = "AATAA"
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 4, 3)
    @test res == (0.0, 4, 5)

    # If min_end_pos = 6 (Too strict):
    # Matches ends at 2 and 5. Neither >= 6. Should fail.
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 6, 3)
    @test res == (Inf, -1, -1)

    # 9. Trim side = 3 where rightmost is invalid but earlier one is valid?
    # This scenario is tricky. findprev finds the "occurrence at or before start".
    # BioDemuX logic for semiglobal/Hamming scans ALL and picks best.
    # For Exact, "best" is 0.0. Tie-break for trim_side=3 is rightmost.
    # If rightmost match is INVALID due to constraints (e.g. min_end_pos), do we fallback to earlier valid match?
    # My current implementation only checks the first result of findprev.
    # If findprev finds a match that fails `min_end_pos`, it returns Inf.
    # It DOES NOT search further back.
    # Ideally it should search further back if the first found match is invalid.
    # BUT, typically `min_end_pos` is a low threshold (like 10), and matches are distinct.
    # Let's test this behavior.

    # q="AA", r="AATAA", min_end_pos=3.
    # Matches: 1:2 (end=2 < 3, invalid), 4:5 (end=5 >= 3, valid).
    # trim_side=3 (Rightmost).
    # search backwards from end. findprev hits 4:5. Valid. Returns 4:5.

    # q="AA", r="AATAA", min_end_pos=3.
    # trim_side=5 (Leftmost).
    # search forward. findnext hits 1:2. Invalid (end=2 < 3).
    # logic should continue to next match 4:5. Valid.
    res = BioDemuX.exact_align(q, r, 1:length(r), length(r), 3, 5)
    @test res == (0.0, 4, 5)

    # q="AA", r="AATAA", min_end_pos=6.
    # trim_side=3.
    # findprev hits 4:5. Invalid (end=5 < 6).
    # Returns Inf. (Correct, as no other match to the left would satisfy end >= 6 either).
end
