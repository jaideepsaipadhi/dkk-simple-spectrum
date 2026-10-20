"""Elements whose c-trace is minuscule except for a single entry 3.

DFS from lambda, extending only by i with c = <alpha_i^vee, mu> > 0 (which is exactly the
condition that the length of the coset representative goes up), and keeping only branches
whose c-values all lie in {1,3} with at most `maxthree` threes. This reaches E7 and E8,
where enumerating the whole Weyl orbit does not.
usage: python3 e3dfs.py [E6|E7|E8] [maxthree=1] [depthcap=40]"""
import sys
from rootsys import RootSystem
from dfs3 import solve

def main():
    typ = sys.argv[1] if len(sys.argv) > 1 else "E6"
    T, N = typ[0], int(typ[1:])
    maxthree = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    rs = RootSystem(T, N)
    print(f"{typ}: c-traces in {{1,3}} with <= {maxthree} threes, depth <= {cap}")
    for k in rs.I:
        lam = rs.fundamental(k)
        hits = []
        seen = set()
        def dfs(mu, cs, word):
            key = (mu, cs.count(3))
            if key in seen: return
            seen.add(key)
            if cs.count(3) >= 1:
                hits.append((tuple(mu), tuple(cs), tuple(word)))
            if len(cs) >= cap: return
            for i in rs.I:
                c = mu[i-1]
                if c not in (1, 3): continue
                if c == 3 and cs.count(3) >= maxthree: continue
                dfs(rs.act(i, mu), cs + [c], word + [i])
        dfs(tuple(lam), [], [])
        if not hits: continue
        # verdict depends only on v
        out = []
        vcache = {}
        for mu, cs, word in hits:
            v = rs.root_coords(tuple(lam[j] - mu[j] for j in range(N)))
            v = tuple(int(x) for x in v)
            if v not in vcache:
                r = solve(rs, k, v)
                if not r: vcache[v] = None
                else:
                    rows = list(r.values())[0]
                    vcache[v] = all(c <= 1 for x in rows.values() for c in x.values())
            if vcache[v] is None: continue
            out.append((vcache[v], len(cs), mu, v, cs, word))
        ns = sum(1 for o in out if o[0])
        print(f"  lam=om{k}: {len(out)} elements, {ns} with simple spectrum")
        for o in sorted(out, key=lambda z: (not z[0], z[1]))[:6]:
            w = "".join("s"+str(i) for i in reversed(o[5]))
            print(f"      simple={str(o[0]):5} L={o[1]:2} mu={o[2]} v={o[3]} c={list(o[4])}")
            print(f"           w = {w}")

if __name__ == "__main__":
    main()
