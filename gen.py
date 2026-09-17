"""Keyed synthetic fragmentation grammar for Cross-Polarity Spectral Correspondence (creator side).

Nothing here is a record of any public archive. A hidden vocabulary of substructure units is drawn once
from a withheld secret. Every unit has a mass and two DIFFERENT fragment signatures, one per ionization
mode, plus mode-specific neutral losses; some units are silent in one mode. A compound is a short chain
of units. Its positive-mode and negative-mode spectra are generated from the same chain through the two
signature tables, with sequence ions, neutral losses, neighbour-context intensity effects, noise peaks and
intensity jitter. The two spectra therefore share almost no peaks, yet correspond through the hidden unit
grammar, which can only be learned from paired training examples.
"""
import hashlib, hmac, math, os, random
from pathlib import Path

N_UNITS = 64
CHAIN = (3, 6)
POS_LOSSES_N, NEG_LOSSES_N = 8, 8
MZ_LO, MZ_HI = 40.0, 1250.0


def load_secret():
    p = os.environ.get("XPOL2_SECRET_FILE") or str(Path(__file__).resolve().parent / "SECRET_DO_NOT_SHARE.txt")
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
                # a fragment is a fixed sub-mass of the unit plus a mode-specific adduct shift
                frac = r.uniform(0.35, 0.98)
                shift = r.choice([1.007, -1.007, 15.011, 17.003, -17.003, 28.031, -18.011, 0.0])
                frags.append((mass * frac + shift, r.uniform(0.15, 1.0)))
            sig[mode] = {"silent": silent, "frags": frags,
                         "loss_prop": [r.random() for _ in range(POS_LOSSES_N if mode == "POS" else NEG_LOSSES_N)],
                         "end_boost": r.uniform(1.0, 2.2)}
        units.append({"mass": mass, "sig": sig})
    losses = {"POS": [r.uniform(14.0, 60.0) for _ in range(POS_LOSSES_N)],
              "NEG": [r.uniform(14.0, 60.0) for _ in range(NEG_LOSSES_N)]}
    # neighbour context: some ordered unit pairs damp or boost each other's fragments, per mode
    context = {"POS": {}, "NEG": {}}
    for mode in context:
        for _ in range(int(N_UNITS * N_UNITS * 0.12)):
            a, b = r.randrange(N_UNITS), r.randrange(N_UNITS)
            context[mode][(a, b)] = r.choice([0.25, 0.4, 0.6, 1.6, 2.2])
    # units of similar mass, for building confusable batches by substitution
    similar = {}
    for u in range(N_UNITS):
        similar[u] = sorted((v for v in range(N_UNITS) if v != u and abs(units[v]["mass"] - units[u]["mass"]) < 6.0),
                            key=lambda v: abs(units[v]["mass"] - units[u]["mass"]))
    return {"units": units, "losses": losses, "context": context, "similar": similar}


def compound_mass(G, chain):
    return sum(G["units"][u]["mass"] for u in chain) + 18.011 * (len(chain) - 1) * 0.0 + 1.0


def spectrum(G, key, chain, mode, label):
    """One acquisition of a compound in one mode. Deterministic given the key and label."""
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
            # neutral losses from this fragment
            for k, lm in enumerate(G["losses"][mode]):
                if r.random() < 0.35 * sg["loss_prop"][k]:
                    peaks.append((mz - lm, inten * r.uniform(0.1, 0.6)))
    # sequence ions: the mode's own direction (prefix in positive mode, suffix in negative mode) is
    # the primary signal; a lighter pass in the opposite direction is also observed (fragmentation does
    # not stop cleanly at the chain boundary), which spreads structurally meaningful mass more evenly
    # across the range instead of leaving it concentrated near each unit's own fragment mass
    adduct = 1.007 if mode == "POS" else -1.007
    cum = 0.0
    seq = chain if mode == "POS" else chain[::-1]
    for i, u in enumerate(seq[:-1]):
        cum += units[u]["mass"]
        if r.random() < 0.90:
            peaks.append((cum + adduct, r.uniform(0.05, 0.5)))
    cum2 = 0.0
    for i, u in enumerate(seq[::-1][:-1]):
        cum2 += units[u]["mass"]
        if r.random() < 0.60:
            peaks.append((cum2 + adduct, r.uniform(0.05, 0.5)))
    # precursor and precursor losses
    if r.random() < 0.6:
        peaks.append((precursor, r.uniform(0.05, 0.8)))
    for lm in G["losses"][mode][:3]:
        if r.random() < 0.3:
            peaks.append((precursor - lm, r.uniform(0.05, 0.4)))
    # noise peaks
    for _ in range(r.randint(3, 8)):
        peaks.append((r.uniform(MZ_LO, min(MZ_HI, precursor + 5)), r.uniform(0.01, 0.08)))
    peaks = [(mz + r.gauss(0, 0.003), it) for mz, it in peaks if MZ_LO <= mz <= MZ_HI and it > 0]
    if not peaks:
        peaks = [(precursor, 1.0)]
    top = max(it for _, it in peaks)
    out = {}
    for mz, it in peaks:
        k = round(mz, 4)
        out[k] = max(out.get(k, 0.0), it / top * 999.0)
    rows = sorted(out.items())
    rows = [(mz, it) for mz, it in rows if it >= 1.0][:80]
    return precursor, [(mz, math.log1p(it)) for mz, it in rows]


def family(G, key, n):
    """A batch family: a seed chain and variants made by substituting units of similar mass, so the
    members are mass-matched and spectrally confusable. Returns list of chains."""
    r = rng_for(key, f"family:{n}")
    K = r.randint(*CHAIN)
    seed = [r.randrange(N_UNITS) for _ in range(K)]
    members = [tuple(seed)]
    tries = 0
    while len(members) < 16 and tries < 200:
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
