"""
Convert a Z_37 lift solution (connection sets D, from Modal LIFT_SOLUTION.json) into the
668x668 Hadamard matrix and verify it.

Chain:  D -> A (333x333 srg adjacency)
        -> S = J - I - 2A          (Seidel matrix; S^2 = 333 I - J,  S*1 = 0)
        -> C = [[0, j^T],[j, S]]   (334x334 SYMMETRIC CONFERENCE MATRIX: C^2 = 333 I_334)
        -> H = [[C+I, C-I],[C-I, -C-I]]   (668x668, entries +-1,  H H^T = 668 I)

Validated on P(5)->H(12) and P(13)->H(28) before use.
"""
import numpy as np, json, sys
import conf_core as cc

def srg_to_conference(A):
    v = A.shape[0]
    S = np.ones((v, v), dtype=np.int64) - np.eye(v, dtype=np.int64) - 2 * A   # J - I - 2A
    C = np.zeros((v + 1, v + 1), dtype=np.int64)
    C[0, 1:] = 1; C[1:, 0] = 1; C[1:, 1:] = S
    return C

def conference_to_hadamard(C):
    n = C.shape[0]
    I = np.eye(n, dtype=np.int64)
    top = np.hstack([C + I, C - I])
    bot = np.hstack([C - I, -C - I])
    return np.vstack([top, bot])

def verify_hadamard(H):
    n = H.shape[0]
    return bool(np.array_equal(np.abs(H), np.ones((n, n), dtype=np.int64))) and \
           bool(np.array_equal(H @ H.T, n * np.eye(n, dtype=np.int64)))

def paley_circulant_adj(p):
    qr = set((x * x) % p for x in range(1, p))
    A = np.zeros((p, p), dtype=np.int64)
    for a in range(p):
        for d in qr:
            A[a, (a + d) % p] = 1
    return A

def validate():
    for p in (5, 13, 17):  # P(p) conference graphs -> C(p+1) -> H(2p+2)
        A = paley_circulant_adj(p)
        C = srg_to_conference(A)
        cc_ok = bool(np.array_equal(C @ C, (p) * np.eye(p + 1, dtype=np.int64))) and bool(np.array_equal(C, C.T))
        H = conference_to_hadamard(C)
        print(f"  P({p}) -> C({p+1}) conference(C^2={p}I,sym)={cc_ok} -> H({2*p+2}) verified={verify_hadamard(H)}")

def emit_from_solution(path="LIFT_SOLUTION.json", out="hadamard_668.csv"):
    sol = json.load(open(path))
    D = [[set(b) for b in row] for row in sol["solution"]]
    A = cc.build_adjacency(D, 37, 9)
    assert cc.srg_identity_direct(A, 166, 82, 83), "solution is NOT a valid srg(333,166,82,83)!"
    C = srg_to_conference(A)
    assert np.array_equal(C @ C, 333 * np.eye(334, dtype=np.int64)), "C(334) not a conference matrix!"
    H = conference_to_hadamard(C)
    assert verify_hadamard(H), "H(668) failed HH^T=668I!"
    np.savetxt(out, H, fmt="%d", delimiter=",")
    print(f"*** H(668) CONSTRUCTED AND VERIFIED -> {out}  (HH^T = 668 I, all entries +-1) ***")
    return H

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "emit":
        emit_from_solution()
    else:
        print("=== VALIDATE conference->Hadamard chain ===")
        validate()
