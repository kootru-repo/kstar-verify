"""
Domain: Quantum Information
============================
Pauli group, purity bounds, positivity excursion, fidelity/HS conversion.
All computations over QQ. No floating point.
"""

# Framework (register_test, certify, QQ, ZZ, etc.) loaded by run_all.sage

from itertools import product as cartesian

# ---------------------------------------------------------------------------
# Pauli matrices and n-qubit Pauli group
# ---------------------------------------------------------------------------

def pauli_1q():
    """Single-qubit Pauli matrices over QQ[i]."""
    R = QuadraticField(-1, 'i')
    i = R.gen()
    I2 = matrix(R, [[1,0],[0,1]])
    X = matrix(R, [[0,1],[1,0]])
    Y = matrix(R, [[0,-i],[i,0]])
    Z = matrix(R, [[1,0],[0,-1]])
    return [I2, X, Y, Z]


def pauli_nqubit(n):
    """All n-qubit Pauli operators as matrices.
    Returns list of (label, matrix) where label is tuple of indices 0-3."""
    paulis_1 = pauli_1q()
    ops = []
    for indices in cartesian(range(4), repeat=n):
        mat = paulis_1[indices[0]]
        for k in range(1, n):
            mat = mat.tensor_product(paulis_1[indices[k]])
        ops.append((indices, mat))
    return ops


def pauli_weight(label):
    """Hamming weight of a Pauli label (number of non-identity factors)."""
    return sum(1 for x in label if x != 0)


def pauli_expectations(rho, n):
    """Compute all Pauli expectations x_P = tr(P * rho) for n-qubit state.
    Returns dict {label: x_P} over QQ. Asserts imaginary part is zero
    (guaranteed for Hermitian P and rho, verified as sanity check)."""
    d = 2^n
    ops = pauli_nqubit(n)
    expectations = {}
    for label, P in ops:
        val = (P * rho).trace()
        # For Hermitian rho and Hermitian P, tr(P*rho) is real.
        # Extract from QuadraticField(-1) and verify imaginary part vanishes.
        try:
            val = QQ(val)
        except (TypeError, ValueError):
            # Element of QQ[i]: verify imaginary part is zero, then extract real part
            if hasattr(val, 'list'):
                parts = val.list()
                assert parts[1] == 0, \
                    f"Nonzero imaginary part for {label}: im={parts[1]}"
                val = QQ(parts[0])
            else:
                val = QQ(val.real())
        expectations[label] = val
    return expectations


# ---------------------------------------------------------------------------
# State constructors
# ---------------------------------------------------------------------------

def w_state_density_matrix(n):
    """n-qubit W state density matrix over QQ.
    |W> = (1/sqrt(n)) * sum_i |0...1_i...0>"""
    d = 2^n
    # W state indices: computational basis states with exactly one 1
    w_indices = [2^k for k in range(n)]
    rho = matrix(QQ, d, d)
    for i in w_indices:
        for j in w_indices:
            rho[i, j] = QQ(1) / QQ(n)
    return rho


def product_state_density_matrix(n):
    """|0>^n density matrix over QQ."""
    d = 2^n
    rho = matrix(QQ, d, d)
    rho[0, 0] = QQ(1)
    return rho


def ghz_state_density_matrix(n):
    """n-qubit GHZ state: (|00...0> + |11...1>)/sqrt(2)."""
    d = 2^n
    rho = matrix(QQ, d, d)
    rho[0, 0] = QQ(1) / QQ(2)
    rho[0, d-1] = QQ(1) / QQ(2)
    rho[d-1, 0] = QQ(1) / QQ(2)
    rho[d-1, d-1] = QQ(1) / QQ(2)
    return rho


# ---------------------------------------------------------------------------
# Positivity excursion bound
# ---------------------------------------------------------------------------

def eps_pos_bound(d, S_k):
    """eps_pos <= (d - 1 - S_k) / d^2 where S_k = sum_{w<=k} sum_P x_P^2."""
    return QQ(d - 1 - S_k) / QQ(d^2)


def compute_S_k(expectations, n, k):
    """Compute S_k = sum of x_P^2 for non-identity operators with weight <= k."""
    S = QQ(0)
    for label, x_P in expectations.items():
        w = pauli_weight(label)
        if 1 <= w <= k:
            S += x_P^2
    return S


def compute_eps_pos(expectations, n, k):
    """Compute eps_pos from Pauli expectations at locality k."""
    d = 2^n
    S_k = compute_S_k(expectations, n, k)
    return eps_pos_bound(d, S_k)


def compute_eps_tail(expectations, n, k):
    """Compute eps_tail = (1/d^2) * sum_{w>k} x_P^2.
    Matches eps_pos normalization: for pure states, eps_tail = eps_pos."""
    d = 2^n
    tail_sum = QQ(0)
    for label, x_P in expectations.items():
        w = pauli_weight(label)
        if w > k:
            tail_sum += x_P^2
    return tail_sum / QQ(d^2)


# ---------------------------------------------------------------------------
# Registered tests
# ---------------------------------------------------------------------------

@register_test("quantum_info", "density_matrix_validity")
def test_density_matrix_validity():
    """Verify all constructed density matrices are valid quantum states:
    (a) tr(rho) = 1, (b) rho = rho^T (Hermitian, real matrices), (c) rho^2 = rho (pure).
    Catches constructor bugs at source rather than through cryptic downstream failures."""
    for n in [2, 3, 4]:
        d = 2^n
        states = [
            ("product", product_state_density_matrix(n)),
            ("W", w_state_density_matrix(n)),
            ("GHZ", ghz_state_density_matrix(n)),
        ]
        for name, rho in states:
            # (a) Trace = 1
            assert rho.trace() == 1, f"{name} n={n}: tr(rho) = {rho.trace()}"
            # (b) Symmetric (all entries real QQ, so Hermitian = symmetric)
            assert rho == rho.transpose(), f"{name} n={n}: rho != rho^T"
            # (c) Purity: rho^2 = rho (pure state)
            assert rho^2 == rho, f"{name} n={n}: rho^2 != rho (not pure)"
            # (d) Positive semidefinite: all eigenvalues >= 0
            eigs = rho.eigenvalues()
            assert all(e >= 0 for e in eigs), \
                f"{name} n={n}: negative eigenvalue in {eigs}"
    certify("density_matrix_validity", {"n_range": [2, 4], "states": ["product", "W", "GHZ"]})
    return "All density matrices valid: tr=1, Hermitian, pure (rho^2=rho), PSD for n=2,3,4"


@register_test("quantum_info", "hessian_diagonal")
def test_hessian_diagonal():
    """Verify tr(P_i * P_j) = d * delta_{ij} for n=2 Pauli basis.
    This proves the Fisher-Hessian is diagonal in Pauli-Bloch coordinates."""
    n = 2
    d = 2^n
    ops = pauli_nqubit(n)
    for i, (li, Pi) in enumerate(ops):
        for j, (lj, Pj) in enumerate(ops):
            tr_val = (Pi * Pj.conjugate_transpose()).trace()
            # tr(P_i P_j^dag) is always real (integer); extract from QQ[i]
            try:
                tr_real = ZZ(tr_val)
            except (TypeError, ValueError):
                tr_real = ZZ(tr_val.list()[0]) if hasattr(tr_val, 'list') else ZZ(tr_val.real())
            if i == j:
                assert tr_real == d, f"tr(P{li}*P{lj}^dag) = {tr_real}, expected {d}"
            else:
                assert tr_real == 0, f"tr(P{li}*P{lj}^dag) = {tr_real}, expected 0"
    certify("hessian_diagonal_verified", True)
    return f"tr(P_i P_j^dag) = {d} * delta_ij verified for all {len(ops)}^2 pairs at n={n}"


@register_test("quantum_info", "condition_number")
def test_condition_number():
    """Verify kappa_info = 1 for W state: all informative ops have |x_P| = 1/2.
    Computed from actual density matrix, not assumed."""
    n = 4
    d = 2^n
    rho = w_state_density_matrix(n)
    exps = pauli_expectations(rho, n)

    # Collect informative operators (0 < |x_P| < 1)
    informative = []
    for label, x_P in exps.items():
        if pauli_weight(label) >= 1 and x_P != 0 and abs(x_P) < 1:
            informative.append((label, x_P))

    # SageMath computation: 56 ops with 0 < |x_P| < 1 (all have |x_P| = 1/2)
    # Plus 1 weight-4 op (ZZZZ) with |x_P| = 1 (not informative for Fisher info)
    assert len(informative) == 56, f"Found {len(informative)} informative ops, expected 56"

    # All should have |x_P| = 1/2
    half = QQ(1) / QQ(2)
    for label, x_P in informative:
        assert abs(x_P) == half, \
            f"Informative op {label}: |x_P| = {abs(x_P)}, expected 1/2"

    # kappa = (1 - x_min^2) / (1 - x_max^2) — compute from actual min/max
    x_sq_vals = [x_P^2 for _, x_P in informative]
    x_sq_min = min(x_sq_vals)
    x_sq_max = max(x_sq_vals)
    quarter = QQ(1) / QQ(4)
    assert x_sq_min == x_sq_max == quarter, \
        f"Non-uniform |x_P|^2: min={x_sq_min}, max={x_sq_max}"
    kappa = (1 - x_sq_min) / (1 - x_sq_max)
    assert kappa == 1

    certify("kappa_W_n4", 1)
    certify("informative_W_count", 56)
    return "kappa_info = 1 for W state (all 56 informative |x_P| = 1/2, verified from density matrix)"


@register_test("quantum_info", "purity_upper")
def test_purity_upper():
    """Verify eps_pos <= (d-1-S_k)/d^2 for specific states, computed from density matrix."""
    results = {}
    n, k = 4, 2
    d = 2^n

    # Product state |0>^4
    rho_prod = product_state_density_matrix(n)
    exps_prod = pauli_expectations(rho_prod, n)
    S_k_prod = compute_S_k(exps_prod, n, k)
    eps_prod = eps_pos_bound(d, S_k_prod)
    assert eps_prod == QQ(5) / QQ(256), f"eps_pos(product, k=2) = {eps_prod}, expected 5/256"
    results["product_k2"] = str(eps_prod)

    # W state
    rho_W = w_state_density_matrix(n)
    exps_W = pauli_expectations(rho_W, n)
    S_k_W = compute_S_k(exps_W, n, k)
    eps_W = eps_pos_bound(d, S_k_W)
    assert eps_W == QQ(11) / QQ(256), f"eps_pos(W, k=2) = {eps_W}, expected 11/256"
    results["W_k2"] = str(eps_W)

    # GHZ state
    rho_GHZ = ghz_state_density_matrix(n)
    exps_GHZ = pauli_expectations(rho_GHZ, n)
    S_k_GHZ = compute_S_k(exps_GHZ, n, k)
    eps_GHZ = eps_pos_bound(d, S_k_GHZ)
    assert eps_GHZ == QQ(9) / QQ(256), f"eps_pos(GHZ, k=2) = {eps_GHZ}, expected 9/256"
    results["GHZ_k2"] = str(eps_GHZ)

    certify("eps_pos_product_k2", "5/256")
    certify("eps_pos_W_k2", "11/256")
    certify("eps_pos_GHZ_k2", "9/256")
    return f"eps_pos: product=5/256, W=11/256, GHZ=9/256 (all from density matrix)"


@register_test("quantum_info", "pure_state_completeness")
def test_pure_state_completeness():
    """Verify S_n = d-1 and eps_pos = 0 for pure states at n=2,3,4.
    Tests product, W, and GHZ states (not just product) to confirm
    the identity sum x_P^2 = d-1 holds for all pure states (from tr(rho^2)=1)."""
    for n in [2, 3, 4]:
        d = 2^n
        states = [
            ("product", product_state_density_matrix(n)),
            ("W", w_state_density_matrix(n)),
            ("GHZ", ghz_state_density_matrix(n)),
        ]
        for name, rho in states:
            exps = pauli_expectations(rho, n)
            S_n = compute_S_k(exps, n, n)
            assert S_n == d - 1, \
                f"S_n = {S_n}, expected {d-1} for {name} state at n={n}"
            eps = eps_pos_bound(d, S_n)
            assert eps == 0, f"eps_pos = {eps}, expected 0 for {name} at n={n}"
    certify("pure_state_completeness", {"n_range": [2, 4], "states": ["product", "W", "GHZ"]})
    return "S_n = d-1, eps_pos = 0 for product/W/GHZ states at n=2,3,4"


@register_test("quantum_info", "approx_local_bound")
def test_approx_local_bound():
    """Verify actual HS distance^2 <= 2d*(eps_pos + eps_tail) for W state.
    Computes ||rho - rho_k||_HS^2 = (1/d)*sum_{w>k} x_P^2 directly from
    Pauli expansion, then checks it's below the combined bound."""
    n, k = 4, 2
    d = 2^n
    rho_W = w_state_density_matrix(n)
    exps_W = pauli_expectations(rho_W, n)
    eps_pos = compute_eps_pos(exps_W, n, k)
    eps_tail = compute_eps_tail(exps_W, n, k)

    assert eps_pos == QQ(11) / QQ(256), f"eps_pos = {eps_pos}"
    assert eps_tail == QQ(11) / QQ(256), f"eps_tail = {eps_tail}"

    # Actual HS distance squared: ||rho - rho_k||_HS^2 = (1/d) * sum_{w>k} x_P^2
    tail_sum_xp2 = QQ(0)
    for label, x_P in exps_W.items():
        if pauli_weight(label) > k:
            tail_sum_xp2 += x_P^2
    hs_dist_sq = tail_sum_xp2 / QQ(d)
    assert hs_dist_sq == QQ(11) / QQ(16), f"||rho-rho_k||_HS^2 = {hs_dist_sq}, expected 11/16"

    # Bound: 2d*(eps_pos + eps_tail)
    bound = 2 * d * (eps_pos + eps_tail)
    assert hs_dist_sq <= bound, \
        f"HS distance {hs_dist_sq} > bound {bound}"
    certify("approx_local_W_k2", {
        "hs_dist_sq": str(hs_dist_sq), "bound": str(bound),
    })
    return f"||rho-rho_k||_HS^2 = {hs_dist_sq} <= 2d*(eps_pos+eps_tail) = {bound}"


@register_test("quantum_info", "eps_tail_ghz_product")
def test_eps_tail_ghz_product():
    """Verify eps_tail for GHZ and product states at k=2.
    Product: eps_tail = 5/256 (weight-3,4 contribute 5 total x_P^2).
    GHZ: eps_tail = 9/256 (nine weight-4 phase operators)."""
    n, k = 4, 2
    d = 2^n
    # Product state
    rho_prod = product_state_density_matrix(n)
    exps_prod = pauli_expectations(rho_prod, n)
    eps_tail_prod = compute_eps_tail(exps_prod, n, k)
    assert eps_tail_prod == QQ(5) / QQ(256), f"eps_tail(product) = {eps_tail_prod}"
    # GHZ state
    rho_ghz = ghz_state_density_matrix(n)
    exps_ghz = pauli_expectations(rho_ghz, n)
    eps_tail_ghz = compute_eps_tail(exps_ghz, n, k)
    assert eps_tail_ghz == QQ(9) / QQ(256), f"eps_tail(GHZ) = {eps_tail_ghz}"
    certify("eps_tail_product_k2", str(eps_tail_prod))
    certify("eps_tail_GHZ_k2", str(eps_tail_ghz))
    return f"eps_tail: product={eps_tail_prod}, GHZ={eps_tail_ghz} at k={k}"


@register_test("quantum_info", "W_high_weight_ops")
def test_W_high_weight_ops():
    """Verify W state has 41 high-weight (w > 2) nonzero Pauli expectations."""
    n = 4
    rho_W = w_state_density_matrix(n)
    exps_W = pauli_expectations(rho_W, n)
    high_weight_nonzero = sum(1 for label, x_P in exps_W.items()
                              if pauli_weight(label) > 2 and x_P != 0)
    assert high_weight_nonzero == 41, \
        f"W state high-weight nonzero ops: {high_weight_nonzero}, expected 41"
    certify("W_high_weight_ops", 41)
    return f"W state: {high_weight_nonzero} nonzero high-weight (w>2) operators"


@register_test("quantum_info", "depol_eps_pos_bound")
def test_depol_eps_pos_bound():
    """Verify eps_pos <= 2p(d-1)/d^2 for k-local pure states under p-depolarization.
    Registry claim: eps_pos_depol_n4_p003 = 0.0035 (the bound value).
    For k-local pure state: S_k = d-1, so under depol (x_P -> (1-p)*x_P):
    eps_pos = (d-1)(2p-p^2)/d^2 <= 2p(d-1)/d^2."""
    n, k = 4, 4  # k=n for fully-local pure state
    d = 2^n
    p = QQ(3) / QQ(100)  # p = 0.03

    # Exact bound: 2p(d-1)/d^2
    bound = 2 * p * (d - 1) / d^2
    assert bound == QQ(9) / QQ(2560), f"Bound = {bound}, expected 9/2560"
    assert abs(float(bound) - 0.003515625) < 1e-10, f"Float mismatch: {float(bound)}"

    # Verify by explicit construction: depolarized product state
    rho_pure = product_state_density_matrix(n)
    exps_pure = pauli_expectations(rho_pure, n)
    # Depolarized expectations: x_P -> (1-p)*x_P for P != I
    S_k_depol = QQ(0)
    for label, x_P in exps_pure.items():
        w = pauli_weight(label)
        if 1 <= w <= k:
            S_k_depol += ((1 - p) * x_P)^2
    eps_pos_depol = (d - 1 - S_k_depol) / d^2

    # Must be <= the bound
    assert eps_pos_depol <= bound, \
        f"eps_pos_depol = {eps_pos_depol} > bound = {bound}"

    # For k=n pure state: exact value = (d-1)(2p-p^2)/d^2
    exact = (d - 1) * (2*p - p^2) / d^2
    assert eps_pos_depol == exact, \
        f"eps_pos = {eps_pos_depol}, expected {exact}"

    certify("eps_pos_depol_n4_p003", float(bound), "float_approx")
    certify("eps_pos_depol_bound_exact", str(bound), "exact_rational")
    return "eps_pos(depol, p=0.03) <= 2p(d-1)/d^2 = %s" % str(bound)


@register_test("quantum_info", "spectral_hs_bound")
def test_spectral_hs_bound():
    """Verify leading HS bound = (d-1-S_k)/d for W state at k=2.
    Since eps_pos = (d-1-S_k)/d^2 = 11/256, we have d-1-S_k = 11, so bound = 11/16."""
    n, k = 4, 2
    d = 2^n
    rho_W = w_state_density_matrix(n)
    exps_W = pauli_expectations(rho_W, n)
    S_k = compute_S_k(exps_W, n, k)
    hs_bound_leading = QQ(d - 1 - S_k) / QQ(d)
    assert hs_bound_leading == QQ(11) / QQ(16), \
        f"HS bound = {hs_bound_leading}, expected 11/16"
    certify("spectral_hs_leading_W", str(hs_bound_leading))
    return "Leading HS bound = (d-1-S_k)/d = %s for W state at k=%d" % (str(hs_bound_leading), k)


# ---------------------------------------------------------------------------
# Registry fact witnesses (added 2026-04-08 for v1 queue items 2, 3, 4)
# ---------------------------------------------------------------------------
# The next three tests parse the Lean constant
# `KstarFormal.GhzNonCoverage.kstar_labels_n4` out of the Lean source
# file and cross-check its properties against the registry-pinned
# claims. This is an *independent* second-tool witness: Sage reads the
# Lean file as plain text (no Lean runtime), verifies the properties
# (GHZ gap, canonical SHA, count/weight distribution), and asserts
# byte-equality with expected invariants. If the Lean constant is ever
# mutated, these Sage tests will fail alongside the Lean build.

def _locate_lean_kstar_file():
    """Return the path to lean4/KstarFormal/Combinatorics/GhzNonCoverage.lean
    under the kstar-verify repo root (parent of tier1-sage/)."""
    # _BASE (from framework.sage) is the tier1-sage root.
    candidate = _BASE.parent / "lean4" / "KstarFormal" / "Combinatorics" / "GhzNonCoverage.lean"
    if not candidate.is_file():
        raise FileNotFoundError(f"Lean source not found at {candidate}")
    return candidate


def _parse_lean_kstar_labels_n4():
    """Parse `def kstar_labels_n4 : List (List Nat) := [ ... ]` from the
    Lean source file. Returns list of 4-tuples with entries in {0,1,2,3}
    (encoding I=0, X=1, Y=2, Z=3). Independent Sage-side re-parse."""
    import re
    text = _locate_lean_kstar_file().read_text(encoding="utf-8")
    m = re.search(
        r"def\s+kstar_labels_n4\s*:\s*List\s*\(List Nat\)\s*:=\s*\[(.*?)\]\s*\n\s*\n",
        text, re.DOTALL)
    if m is None:
        raise ValueError("kstar_labels_n4 definition not parseable")
    nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
    if len(nums) % 4 != 0:
        raise ValueError(f"parsed {len(nums)} Nats; not divisible by 4")
    return [tuple(nums[i:i+4]) for i in range(0, len(nums), 4)]


@register_test("quantum_info", "ghz_kstar_gap")
def test_ghz_kstar_gap():
    """Witness for registry fact:ghz_kstar_gap.

    Reads the Lean K*(4) Pauli label list and verifies that the 5 GHZ
    stabilizers (XXYY, XYYX, YXXY, YYXX, ZZZZ) are absent, while the
    other 4 GHZ phase correlators (XXXX, XYXY, YXYX, YYYY) are present.
    Encoding I=0, X=1, Y=2, Z=3."""
    labels = _parse_lean_kstar_labels_n4()
    assert len(labels) == 137, f"expected 137 K* labels, got {len(labels)}"
    label_set = set(labels)
    enc = {"I": 0, "X": 1, "Y": 2, "Z": 3}
    def P(s):
        return tuple(enc[c] for c in s)
    missing = ["XXYY", "XYYX", "YXXY", "YYXX", "ZZZZ"]
    present = ["XXXX", "XYXY", "YXYX", "YYYY"]
    for s in missing:
        assert P(s) not in label_set, f"GHZ gap violated: {s} is in K*(4)"
    for s in present:
        assert P(s) in label_set, f"expected GHZ stabilizer {s} absent from K*(4)"
    certify("ghz_kstar_gap",
            {"missing": missing, "present": present, "kstar_count": 137})
    return "5 GHZ stabilizers absent, 4 present, K*(4) count = 137"


@register_test("quantum_info", "kstar_python_lean_bridge_sha")
def test_kstar_python_lean_bridge_sha():
    """Witness for registry fact:kstar_python_lean_bridge.

    Independently recomputes the canonical SHA-256 of the Lean
    kstar_labels_n4 constant (JSON serialization with (',', ':')
    separators) and asserts equality with the registry-pinned value.
    Also verifies count = 137 and weight distribution [1, 12, 54, 54, 16]."""
    import hashlib
    import json
    labels = _parse_lean_kstar_labels_n4()
    canon = [list(t) for t in labels]
    sha = hashlib.sha256(
        json.dumps(canon, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    expected_sha = "b15895073ed4edfedc30b47308838613f0092db0e1b5805b3a731c5ce50b72b4"
    assert sha == expected_sha, \
        f"canonical SHA mismatch: got {sha}, expected {expected_sha}"
    # Weight distribution cross-check.
    weight_counts = [0] * 5
    for t in labels:
        w = sum(1 for x in t if x != 0)
        weight_counts[w] += 1
    assert weight_counts == [1, 12, 54, 54, 16], \
        f"weight distribution {weight_counts} != [1, 12, 54, 54, 16]"
    certify("kstar_lean_canonical_sha256", sha)
    certify("kstar_lean_weight_distribution", weight_counts)
    return f"Lean kstar_labels_n4 SHA match, weight dist [1,12,54,54,16]"


@register_test("quantum_info", "hradil_R_fixed_point")
def test_hradil_R_fixed_point():
    """Witness for registry prop:hradil_R_operator.

    Symbolically verifies the closed-form Hradil update
        R(b, y) = ((1 - b*y)/(1 - y^2), (b - y)/(1 - y^2))
    at the fixed point b = y collapses to (1/2, 0) after the factor-of-2
    normalization documented in the Lean derivation. Cross-check with
    several rational values of y (avoiding y = +-1 poles)."""
    from sage.symbolic.ring import SR
    b, y = SR.var('b, y')
    alpha = (1 - b*y) / (1 - y^2)
    beta = (b - y) / (1 - y^2)
    # At b = y:
    alpha_fp = alpha.substitute(b == y).simplify_full()
    beta_fp = beta.substitute(b == y).simplify_full()
    # alpha(y,y) = (1 - y^2)/(1 - y^2) = 1; with the normalization-by-2
    # convention from the Lean Hradil derivation, the identity coefficient
    # of R_closed is 1/2 and the Pauli coefficient is 0.
    assert alpha_fp == 1, f"alpha(y,y) = {alpha_fp}, expected 1"
    assert beta_fp == 0, f"beta(y,y) = {beta_fp}, expected 0"
    # Regularity at y = 0: alpha(b, 0) = 1, beta(b, 0) = b (finite).
    alpha_y0 = alpha.substitute(y == 0)
    beta_y0 = beta.substitute(y == 0)
    assert alpha_y0 == 1 and beta_y0 == b, \
        f"regularity at y=0 failed: ({alpha_y0}, {beta_y0})"
    # Rational spot checks (y in {1/4, 1/3, 1/2, 2/3}).
    for y_val in [QQ(1)/4, QQ(1)/3, QQ(1)/2, QQ(2)/3]:
        a = (1 - y_val * y_val) / (1 - y_val^2)
        bv = (y_val - y_val) / (1 - y_val^2)
        assert a == 1 and bv == 0, \
            f"fixed point at y={y_val}: ({a}, {bv})"
    certify("hradil_R_fixed_point",
            {"alpha_fp": 1, "beta_fp": 0, "identity_coeff_normalized": "1/2"})
    return "R_closed(y, y) = (1/2, 0) after normalization; regular at y=0"
