"""Systematic hypothesis search on the labelled dataset."""
import json
from rootsys import RootSystem

DATA = json.load(open("data.json"))
RS = {}
def rs_of(d):
    key = (d['typ'], d['n'])
    if key not in RS: RS[key] = RootSystem(*key)
    return RS[key]

def arms(rs, t):
    out = []
    for nb in sorted(rs.adj[t]):
        seen = {t}; cur = nb; path = [nb]
        while True:
            nxt = [x for x in rs.adj[cur] if x not in seen]
            seen.add(cur)
            if len(nxt) != 1: break
            cur = nxt[0]; path.append(cur)
        out.append(path)
    return out

def tri(rs):
    c = [i for i in rs.I if len(rs.adj[i]) >= 3]
    return c[0] if c else None

CAND = {}

def reg(name):
    def deco(f):
        CAND[name] = f; return f
    return deco

@reg("cmax==1 or v_t<=2")
def p1(d, rs, t, ar, v):
    return d['cmax'] == 1 or (t is None) or v[t-1] <= 2

@reg("cmax==1 or exists arm with exact staircase")
def p2(d, rs, t, ar, v):
    if d['cmax'] == 1 or t is None: return True
    vt = v[t-1]
    for A in ar:
        if len(A) < vt - 1: continue
        if all(v[A[j]-1] == vt - 1 - j for j in range(vt - 1)): return True
    return False

@reg("cmax==1 or exists arm with >= staircase")
def p3(d, rs, t, ar, v):
    if d['cmax'] == 1 or t is None: return True
    vt = v[t-1]
    for A in ar:
        if len(A) < vt - 1: continue
        if all(v[A[j]-1] >= vt - 1 - j for j in range(vt - 1)): return True
    return False

@reg("cmax==1 or TWO arms with >= staircase")
def p4(d, rs, t, ar, v):
    if d['cmax'] == 1 or t is None: return True
    vt = v[t-1]
    cnt = 0
    for A in ar:
        if len(A) < vt - 1: continue
        if all(v[A[j]-1] >= vt - 1 - j for j in range(vt - 1)): cnt += 1
    return cnt >= 2

@reg("cmax==1 or ALL arms >= staircase")
def p5(d, rs, t, ar, v):
    if d['cmax'] == 1 or t is None: return True
    vt = v[t-1]
    for A in ar:
        if len(A) < vt - 1: return False
        if not all(v[A[j]-1] >= vt - 1 - j for j in range(vt - 1)): return False
    return True

@reg("cmax==1 or v_i<=2 for ALL i")
def p6(d, rs, t, ar, v):
    return d['cmax'] == 1 or all(x <= 2 for x in v)

@reg("cmax==1 or v_t <= 1+min arm len")
def p7(d, rs, t, ar, v):
    if d['cmax'] == 1 or t is None: return True
    return v[t-1] <= 1 + min(len(A) for A in ar)

@reg("cmax==1 or v_t <= 1+2nd smallest arm len")
def p8(d, rs, t, ar, v):
    if d['cmax'] == 1 or t is None: return True
    L = sorted(len(A) for A in ar)
    return v[t-1] <= 1 + L[1]

for name, f in CAND.items():
    exc = 0; byt = {}
    first = []
    for d in DATA:
        rs = rs_of(d); t = tri(rs); ar = arms(rs, t) if t else []
        pred = f(d, rs, t, ar, d['v'])
        if pred != d['simple']:
            exc += 1
            byt[d['typ']] = byt.get(d['typ'], 0) + 1
            if len(first) < 2:
                first.append((d['typ'], d['n'], d['k'], d['v'], d['cmax'], d['simple'], pred))
    print(f"{name:<42} exceptions {exc:4d}/{len(DATA)}  by type {byt}")
    for x in first: print(f"      {x}")
