"""up-set at j, exactly as typeD.tex def:upset has it: there is a degree D with
X_j[d] = 0 for d < D, X_j[D] an ARBITRARY subspace, and X_j[d] full for d > D.
Reading degrees from the top: a run of full slots, then at most one partial,
then zeros."""
def is_upset(col, P):
    seen_partial = False
    for d in sorted(P, reverse=True):
        c = col.get(d, 0)
        if c == P[d] and not seen_partial:
            continue                      # still in the full run
        if c == 0:
            seen_partial = True           # from here down everything must be 0
            continue
        if 0 < c < P[d] and not seen_partial:
            seen_partial = True           # the one allowed partial slot
            continue
        return False                      # nonzero below the partial slot
    return True
