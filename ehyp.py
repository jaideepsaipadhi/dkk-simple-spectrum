"""
HYPOTHESIS (uniform ADE):  let L(Gamma,k) = the number of leading 1s in the
colour-t multiplicity profile of Pi e_k (t = trivalent node).  Then

        T acts with simple spectrum   <=>   v_t <= L(Gamma,k).

In type D, (Pi e_k)_t = (1,2^{k-1},1) so L = 1 for k >= 2 ... no: L counts
leading 1s, which is 1 for k>=2 and the whole thing for k=1.  The type-D
theorem says the threshold is 2, so the correct statement uses the
BOTTOM-TRUNCATION formulation: V_t is the greedy bottom-truncation of
(Pi e_k)_t to total v_t, and simple <=> that truncation is all-ones, i.e.

        v_t <= (sum of the leading run of 1s)  =  L(Gamma,k).

For type D k>=2 the leading run is a single 1 of length 1, giving threshold 1
-- which contradicts the proved threshold 2.  So the profile alone is not
enough; the *degrees* matter.  This script measures the truth directly:
for each (type,k) it reports the largest v_t among simple elements and the
smallest v_t among non-simple ones, i.e. whether v_t is a complete invariant
in type E as it is in type D.
"""
import json
from collections import defaultdict

TRIV = {"D": lambda n: n - 2, "E": lambda n: 4}


def main():
    data = json.load(open("data.json"))
    try:
        data += [json.loads(l) for l in open("extra.jsonl")]
    except OSError:
        pass
    buckets = defaultdict(lambda: (set(), set()))
    for d in data:
        typ = d["typ"]
        if typ == "A":
            continue
        t = TRIV[typ](d["n"])
        vt = d["v"][t - 1]
        s, ns = buckets[(typ, d["n"], d["k"])]
        (s if d["simple"] else ns).add(vt)
    print(f"{'case':12} {'v_t simple':22} {'v_t non-simple':22} separated?")
    allsep = True
    for key in sorted(buckets):
        s, ns = buckets[key]
        sep = (not s or not ns) or (max(s) < min(ns))
        allsep &= sep
        name = f"{key[0]}{key[1]} om{key[2]}"
        print(f"{name:12} {str(sorted(s)):22} {str(sorted(ns)):22} "
              f"{'YES thr=' + str(max(s) if s else '-') if sep else 'NO -- overlap'}")
    print()
    print("v_t is a complete invariant in every case" if allsep
          else "v_t alone does NOT decide simplicity")


if __name__ == "__main__":
    main()
