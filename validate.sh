#!/usr/bin/env bash
# validate.sh -- check the proof on your own machine.
#
#   ./validate.sh            the audit and the cited scripts   (~5 minutes)
#   ./validate.sh quick      the audit alone                   (~2 minutes)
#   ./validate.sh full       ... and the exploration log too
#
# Only the first two matter: they are the audit and the scripts the manuscript
# cites.  "full" additionally runs the exploration log -- scratch files from the
# search phase that nothing in the proof depends on -- and several of those are
# slow, so each script is capped (TIMEOUT seconds, default 120) and anything
# over the cap is reported as "slow" rather than allowed to block.
#
# Needs python3 and nothing else: no packages, no network, no compiler.
# Exit status is 0 only if every check passes.

set -u
cd "$(dirname "$0")"

PY=${PYTHON:-python3}
MODE=${1:-default}
FAILED=""
T0=$SECONDS

# every script the manuscript cites, in dependency order
CITED="proof.py verify.py sweepall.py inwardlemma.py master.py upsets.py starchar.py localstar.py reduce.py canon.py strands.py thresh.py epsilon.py rppmod.py eproj.py epscrit.py filltop.py headrule.py homological.py inwardproof.py socdim.py squeeze.py staircase.py trunctest.py violators.py bigrade.py claimG.py"

hdr () { printf '\n\033[1m%s\033[0m\n\n' "$1"; }

SLOW=""
CAP=0                          # 0 = no cap; set per section below
EXPCAP=${TIMEOUT:-120}         # cap for the exploration log only
HAVE_TIMEOUT=0
command -v timeout >/dev/null 2>&1 && HAVE_TIMEOUT=1

run () {                       # run <script>
    printf '  %-18s ' "$1"
    if [ "$CAP" != 0 ] && [ "$HAVE_TIMEOUT" = 1 ]; then
        OUT=$(timeout "$CAP" "$PY" "$1" 2>&1); RC=$?
    else
        OUT=$("$PY" "$1" 2>&1); RC=$?
    fi
    if [ "$RC" = 124 ]; then
        printf 'slow  (over %ss, skipped)\n' "$CAP"
        SLOW="$SLOW $1"
        return
    fi
    if [ "$RC" = 0 ]; then
        printf 'ok    %s\n' "$(printf '%s' "$OUT" | tail -1 | cut -c1-52)"
    else
        printf 'FAIL\n'
        printf '%s\n' "$OUT" | tail -12 | sed 's/^/        /'
        FAILED="$FAILED $1"
    fi
}

"$PY" - <<'EOF' || { echo "python3 too old: need 3.8+"; exit 2; }
import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)
EOF

hdr "THE AUDIT  --  every link in the chain, plus negative controls"
run audit.py

if [ "$MODE" != "quick" ]; then
    hdr "THE CITED SCRIPTS  --  what the manuscript points at"
    for s in $CITED; do
        run "$s"
    done
fi

if [ "$MODE" = "full" ]; then
    hdr "EVERYTHING ELSE  --  the exploration log (not load-bearing)"
    CAP=$EXPCAP
    for s in *.py; do
        case " audit.py $CITED " in
            *" $s "*) continue ;;
        esac
        run "$s"
    done
fi

hdr "RESULT"
if [ -n "$SLOW" ]; then
    printf '  slow, skipped:%s\n' "$SLOW"
fi
if [ -z "$FAILED" ]; then
    printf '  everything passed  (%ds)\n\n' "$((SECONDS - T0))"
    exit 0
else
    printf '  FAILED:%s   (%ds)\n\n' "$FAILED" "$((SECONDS - T0))"
    exit 1
fi
