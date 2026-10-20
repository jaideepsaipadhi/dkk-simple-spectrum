"""Ground truth + candidate features for a type-E simple-spectrum criterion.

Type D's answer has the shape: write mu in a natural basis, count two things, require
r + 2m <= 2.  We need the analogue for E.  This dumps, for every extremal weight, the
verdict plus a battery of candidate statistics of mu and of the c-trace, so a separating
rule can be searched for rather than guessed.
usage: python3 efeat.py [TYPE=E6] [OUT=e6.tsv]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace

from socle import soc_sigma

def main():
    typ = sys.argv[1] if len(sys.argv) > 1 else "E6"
    out = sys.argv[2] if len(sys.argv) > 2 else "e6.tsv"
    T, N = typ[0], int(typ[1:])
    rs = RootSystem(T, N)
    rows = []
    for k in rs.I:
        lam = rs.fundamental(k)
        for mu, w in orbit_with_words(rs, lam).items():
            if not w: continue
            tr = trace(rs, w, lam)
            cs = [c for _, _, c, _ in tr]
            cols = [i for _, i, c, _ in tr if c >= 2]
            X, rowsd, dv = soc_sigma(rs, k, w)
            simple = all(c <= 1 for x in rowsd.values() for c in x.values())
            neg = [x for x in mu if x < 0]
            rows.append(dict(
                ka=k, L=len(w), simple=int(simple),
                maxc=max(cs), n2=sum(1 for c in cs if c == 2),
                n3=sum(1 for c in cs if c == 3), n4p=sum(1 for c in cs if c >= 4),
                exc=sum(c-1 for c in cs),                  # total excess over minuscule
                nneg=len(neg), sneg=-sum(neg), minmu=min(mu),

                bigcols="|".join(map(str, sorted(set(cols)))),
                mu="|".join(map(str, mu)), v="|".join(map(str, dv)),
            ))
        print(f"  om{k} done ({len(rows)} rows)", flush=True)
    keys = list(rows[0].keys())
    with open(out, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows: f.write("\t".join(str(r[k]) for k in keys) + "\n")
    ns = sum(r['simple'] for r in rows)
    print(f"{typ}: {len(rows)} extremal weights, {ns} with simple spectrum -> {out}")

if __name__ == "__main__":
    main()
