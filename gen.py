"""Keyed generator for Cross-Polarity Acquisition Grouping (creator side).

A batch is a pile of acquisitions, not a list of compounds. Each compound in a batch is acquired zero to
three times in positive mode and zero to three times in negative mode, so the correspondence between the
two polarities is many-to-many and the answer is a partition of the batch rather than a matching.

Spectra come from a withheld fragmentation grammar: a hidden vocabulary of substructure units, each with
two different fragment signatures (one per polarity), mode-specific neutral losses, sequence ions in both
directions and neighbour-context effects. A compound is a chain of units.

Three degradation modalities are applied per acquisition and never labelled:
  co-isolation    the isolation window admitted a neighbouring compound of the same batch; its fragments
                  appear at reduced intensity, so the spectrum is a mixture
  saturation      the most intense peaks are clipped flat at a detector ceiling, destroying the intensity
                  profile while leaving m/z intact
  low injection   only the strongest few peaks rise above noise
"""
import hashlib, hmac, math, os, random
from pathlib import Path

N_UNITS = 64
CHAIN = (3, 6)
POS_LOSSES_N, NEG_LOSSES_N = 8, 8
MZ_LO, MZ_HI = 40.0, 1250.0
N_BATCHES = 900

COMPOUNDS_PER_BATCH = (10, 15)
REPS_PER_FAMILY = 3
ACQ_PER_MODE = (0, 3)
P_COISOLATION = 0.15
P_SATURATION = 0.15
P_LOWINJECTION = 0.15


def load_secret():
    p = os.environ.get("XPOL3_SECRET_FILE") or str(Path(__file__).resolve().parent / "SECRET_DO_NOT_SHARE.txt")
    return bytes.fromhex(Path(p).read_text().strip())


def rng_for(key, label):
    return random.Random(int.from_bytes(hmac.new(key, label.encode(), hashlib.sha256).digest()[:8], "big"))


def opaque(key, label, n=12):
    return hmac.new(key, label.encode(), hashlib.sha256).hexdigest()[:n]


def grammar(key):
    """The hidden chemistry: units, per-mode signatures, losses and pair-context effects."""
    r = rng_for(key, "grammar")
    units = []
    for u in range(N_UNITS):
        mass = r.uniform(44.0, 176.0)
        sig = {}
        for mode in ("POS", "NEG"):
            silent = r.random() < 0.15
            n_frag = 0 if silent else r.randint(2, 5)
            frags = []
            for _ in range(n_frag):
                frac = r.uniform(0.35, 0.98)
                shift = r.choice([1.007, -1.007, 15.011, 17.003, -17.003, 28.031, -18.011, 0.0])
                frags.append((mass * frac + shift, r.uniform(0.15, 1.0)))
            sig[mode] = {"silent": silent, "frags": frags,
                         "loss_prop": [r.random() for _ in range(POS_LOSSES_N if mode == "POS" else NEG_LOSSES_N)],
                         "end_boost": r.uniform(1.0, 2.2)}
        units.append({"mass": mass, "sig": sig})
    losses = {"POS": [r.uniform(14.0, 60.0) for _ in range(POS_LOSSES_N)],
              "NEG": [r.uniform(14.0, 60.0) for _ in range(NEG_LOSSES_N)]}
    context = {"POS": {}, "NEG": {}}
    for mode in context:
        for _ in range(int(N_UNITS * N_UNITS * 0.12)):
            a, b = r.randrange(N_UNITS), r.randrange(N_UNITS)
            context[mode][(a, b)] = r.choice([0.25, 0.4, 0.6, 1.6, 2.2])
    similar = {}
    for u in range(N_UNITS):
        similar[u] = sorted((v for v in range(N_UNITS) if v != u and abs(units[v]["mass"] - units[u]["mass"]) < 6.0),
                            key=lambda v: abs(units[v]["mass"] - units[u]["mass"]))
    return {"units": units, "losses": losses, "context": context, "similar": similar}


def compound_mass(G, chain):
    return sum(G["units"][u]["mass"] for u in chain) + 1.0


def clean_peaks(G, key, chain, mode, label):
    """The undegraded fragment list of one acquisition, as (mz, linear intensity) before normalisation."""
    r = rng_for(key, label)
    units = G["units"]
    M = compound_mass(G, chain)
    precursor = M + (1.007 if mode == "POS" else -1.007)
    peaks = []
    K = len(chain)
    for i, u in enumerate(chain):
        sg = units[u]["sig"][mode]
        if sg["silent"]:
            continue
        pos_factor = sg["end_boost"] if i in (0, K - 1) else 1.0
        ctx = 1.0
        if i > 0:
            ctx *= G["context"][mode].get((chain[i - 1], u), 1.0)
        if i < K - 1:
            ctx *= G["context"][mode].get((u, chain[i + 1]), 1.0)
        for mz, base in sg["frags"]:
            inten = base * pos_factor * ctx * math.exp(r.gauss(0, 0.35))
            if r.random() < 0.85:
                peaks.append((mz, inten))
            for k, lm in enumerate(G["losses"][mode]):
                if r.random() < 0.35 * sg["loss_prop"][k]:
                    peaks.append((mz - lm, inten * r.uniform(0.1, 0.6)))
    adduct = 1.007 if mode == "POS" else -1.007
    cum = 0.0
    seq = chain if mode == "POS" else chain[::-1]
    for u in seq[:-1]:
        cum += units[u]["mass"]
        if r.random() < 0.90:
            peaks.append((cum + adduct, r.uniform(0.05, 0.5)))
    cum2 = 0.0
    for u in seq[::-1][:-1]:
        cum2 += units[u]["mass"]
        if r.random() < 0.60:
            peaks.append((cum2 + adduct, r.uniform(0.05, 0.5)))
    if r.random() < 0.6:
        peaks.append((precursor, r.uniform(0.05, 0.8)))
    for lm in G["losses"][mode][:3]:
        if r.random() < 0.3:
            peaks.append((precursor - lm, r.uniform(0.05, 0.4)))
    for _ in range(r.randint(3, 8)):
        peaks.append((r.uniform(MZ_LO, min(MZ_HI, precursor + 5)), r.uniform(0.01, 0.08)))
    return [(mz + r.gauss(0, 0.003), it) for mz, it in peaks if MZ_LO <= mz <= MZ_HI and it > 0]


def acquisition(G, key, chain, mode, label, contaminant_chain=None):
    """One acquisition, with its degradations applied. Returns (peaks, flags) - flags are creator-side."""
    r = rng_for(key, "deg:" + label)
    peaks = clean_peaks(G, key, chain, mode, label)
    flags = []
    if contaminant_chain is not None and r.random() < P_COISOLATION:
        share = r.uniform(0.10, 0.40)
        cpk = clean_peaks(G, key, contaminant_chain, mode, label + ":contam")
        if cpk:
            top_own = max(it for _, it in peaks) if peaks else 1.0
            top_c = max(it for _, it in cpk)
            peaks = peaks + [(mz, it / top_c * top_own * share) for mz, it in cpk]
            flags.append("coisolated")
    if not peaks:
        peaks = [(compound_mass(G, chain), 1.0)]
    top = max(it for _, it in peaks)
    if r.random() < P_SATURATION:
        ceiling = top * r.uniform(0.25, 0.55)
        peaks = [(mz, min(it, ceiling)) for mz, it in peaks]
        flags.append("saturated")
        top = max(it for _, it in peaks)
    out = {}
    for mz, it in peaks:
        k = round(mz, 4)
        out[k] = max(out.get(k, 0.0), it / top * 999.0)
    rows = sorted(((mz, it) for mz, it in out.items() if it >= 1.0), key=lambda x: x[0])[:80]
    if r.random() < P_LOWINJECTION and len(rows) > 6:
        keep = r.randint(5, 12)
        rows = sorted(sorted(rows, key=lambda x: -x[1])[:keep], key=lambda x: x[0])
        flags.append("low_injection")
    if not rows:
        rows = [(compound_mass(G, chain), 999.0)]
    return [(mz, math.log1p(it)) for mz, it in rows], flags


def family(G, key, n):
    """Mass-matched, spectrally confusable compounds: a seed chain plus unit substitutions."""
    r = rng_for(key, f"family:{n}")
    K = r.randint(*CHAIN)
    seed = [r.randrange(N_UNITS) for _ in range(K)]
    members, tries = [tuple(seed)], 0
    while len(members) < 18 and tries < 300:
        tries += 1
        c = list(r.choice(members))
        for _ in range(r.randint(1, 2)):
            i = r.randrange(K)
            sim = G["similar"][c[i]]
            if sim:
                c[i] = r.choice(sim[:6])
        t = tuple(c)
        if t not in members:
            members.append(t)
    return members


def batch(G, key, n, rep=0):
    """One batch: compounds with several acquisitions each, degraded, anonymised into slots.
    Batches sharing a family index share their compounds and differ only in the acquisitions drawn."""
    r = rng_for(key, f"batch:{n}:{rep}")
    chains = family(G, key, n)
    n_comp = min(len(chains), r.randint(*COMPOUNDS_PER_BATCH))
    chosen = chains[:n_comp]
    plan = []
    for j, chain in enumerate(chosen):
        np_, nn = r.randint(*ACQ_PER_MODE), r.randint(*ACQ_PER_MODE)
        if np_ + nn == 0:
            np_ = 1
        plan.append([chain, np_, nn])
    # the batch must be able to score both terms: at least two compounds present in both polarities,
    # and at least one compound acquired twice in some polarity
    both = [p for p in plan if p[1] > 0 and p[2] > 0]
    while len(both) < 2:
        p = r.choice(plan)
        p[1] = max(1, p[1]); p[2] = max(1, p[2])
        both = [q for q in plan if q[1] > 0 and q[2] > 0]
    if not any(p[1] > 1 or p[2] > 1 for p in plan):
        p = r.choice(both)
        p[r.choice([1, 2])] += 1
    acqs = []
    for j, (chain, np_, nn) in enumerate(plan):
        others = [c for c in chosen if c != chain]
        for mode, count in (("POS", np_), ("NEG", nn)):
            for a in range(count):
                contam = r.choice(others) if others else None
                pk, flags = acquisition(G, key, chain, mode, f"acq:{n}:{rep}:{j}:{mode}:{a}", contam)
                acqs.append({"compound": j, "mode": mode, "peaks": pk, "flags": flags,
                             "precursor": compound_mass(G, chain) + (1.007 if mode == "POS" else -1.007)})
    pos = [a for a in acqs if a["mode"] == "POS"]
    neg = [a for a in acqs if a["mode"] == "NEG"]
    r.shuffle(pos); r.shuffle(neg)
    for i, a in enumerate(pos):
        a["slot"] = "P%02d" % i
    for i, a in enumerate(neg):
        a["slot"] = "N%02d" % i
    return {"case_id": "batch_" + opaque(key, f"case:{n}:{rep}"), "block_id": "fam_" + opaque(key, f"fam:{n}", 8),
            "acquisitions": pos + neg,
            "group_of": {a["slot"]: "c%d" % a["compound"] for a in pos + neg}}
