from sage.all import Matrix, ZZ, LaurentPolynomialRing
import spherogram as sph
import heapq

def _build_strand_to_region(faces):
    """
    Map each CrossingStrand object to its region (column) index.
    """
    strand_to_region = {}
    for r_idx, face in enumerate(faces):
        for cs in face:
            strand_to_region[cs] = r_idx
    return strand_to_region

def _region_adjacency_from_faces(faces, strand_to_region):
    """
    Build adjacency list for regions using cs.opposite() across each boundary half-edge.
    Two regions are adjacent if they share an edge (pair of opposite half-edges).
    """
    m = len(faces)
    adj = [set() for _ in range(m)]
    for r_idx, face in enumerate(faces):
        for cs in face:
            try:
                other = cs.opposite()
            except Exception:
                continue
            r2 = strand_to_region.get(other, None)
            if r2 is not None and r2 != r_idx:
                adj[r_idx].add(r2)
                adj[r2].add(r_idx)
    return adj

def _two_color_regions(adj, start_index):
    """
    BFS 2-coloring of the region graph (planar checkerboard).
    Returns a list 'color' with entries in {0,1} and colors[start_index] = 0.
    """
    from collections import deque
    m = len(adj)
    color = [-1] * m
    q = deque([start_index])
    color[start_index] = 0
    while q:
        u = q.popleft()
        for v in adj[u]:
            if color[v] == -1:
                color[v] = 1 - color[u]
                q.append(v)
            elif color[v] == color[u]:
                raise RuntimeError("Region graph not bipartite; diagram may be degenerate.")
    # If disconnected (shouldn't happen), color remaining components arbitrarily
    for i in range(m):
        if color[i] == -1:
            color[i] = 0
    return color

def _choose_unbounded_region(faces):
    """
    Heuristic: pick the region with the largest boundary length (most CrossingStrands).
    This typically is the unbounded region.
    """
    lengths = [len(face) for face in faces]
    return max(range(len(faces)), key=lambda i: lengths[i])

def alexander_region_matrix_2var(L, *, unbounded_region_index=None, flip_incoming_rule=False, debug=False, show_matrix = False):
    """
    Alexander's 1928 region-crossing matrix over ZZ[s^{±1}, t^{±1}] with checkerboard rule.

    At each crossing row:
      1) Find the incoming-under half-edge.
      2) Place (t, -t, 1, -1) on the four adjacent regions in CCW order,
         starting at the region immediately CCW from that incoming-under half-edge.
      3) For each ±t entry, replace t -> s (and -t -> -s) if that region's color
         differs from the unbounded region in the unique proper 2-coloring.

    Parameters
    ----------
    unbounded_region_index : int or None
        If None, choose the region with maximal boundary length as the unbounded region.
        Provide an integer to force a specific unbounded region.

    flip_incoming_rule : bool
        Controls how we choose the incoming-under index from c.directions and c.sign.
        If False (default): incoming_under = min(incomings) if sign<0 else max(incomings).
        If True:           incoming_under = max(incomings) if sign<0 else min(incomings).

    debug : bool
        Print per-crossing/per-region diagnostics.

    Returns
    -------
    M : sage Matrix over ZZ[s^{±1}, t^{±1}] with shape (n, n+2)
    regions : list of faces (each a list of CrossingStrand)
    colors : list of ints in {0,1}, with colors[unbounded_region_index] = 0
    unbounded_region_index : int
    """
    faces = L.faces()
    crossings = L.crossings
    n, m = len(crossings), len(faces)

    # Laurent polynomial ring ZZ[s^{±1}, t^{±1}]
    R = LaurentPolynomialRing(ZZ, ('s', 't'))
    s, t = R.gens()

    # Map half-edges to regions; build region adjacency; color regions
    strand_to_region = _build_strand_to_region(faces)
    adj = _region_adjacency_from_faces(faces, strand_to_region)

    if unbounded_region_index is None:
        unbounded_region_index = _choose_unbounded_region(faces)
    colors = _two_color_regions(adj, unbounded_region_index)

    if debug:
        print(f"Selected unbounded region index = {unbounded_region_index}")
        print(f"Unbounded region color = {colors[unbounded_region_index]}")
        print("Region boundary lengths:", [len(f) for f in faces])
        print(colors, "<== Colors")

    if show_matrix:
        print(colors, "<== Colors")

    def weight_for_region(region_idx, base_t):
        """
        Given base_t ∈ {t, -t, 1, -1}, return the coefficient after applying the s/t swap rule.
        Only ±t are swapped to ±s; 1 and -1 remain unchanged.
        """
        if base_t not in (t, -t):
            return base_t
        same_color = (colors[region_idx] == colors[unbounded_region_index])
        if same_color:
            return base_t
        return s if base_t == t else -s

    # Build the matrix
    M = Matrix(R, n, m, 0)

    for row, c in enumerate(crossings):
        css = c.crossing_strands()  # list of 4 CrossingStrand in CCW order
        if debug:
            print(f"\nCrossing {row}: sign={c.sign}")

        base_coeffs = (R(1), R(-1), t, -t)

        for idx, base in enumerate(base_coeffs):
            cs = css[idx]
            if cs not in strand_to_region:
                raise RuntimeError(f"Half-edge {cs} not found in any region at crossing {row}.")
            col = strand_to_region[cs]
            coeff = weight_for_region(col, base)
            M[row, col] += coeff

            if debug:
                print(f"  idx={idx}, region col={col}, base={base}, placed={coeff}")

    if m != n + 2 and debug:
        print(f"Warning: regions m={m} differs from n+2={n+2}.")

    if debug or show_matrix:
        print(M)

    return M, faces, colors, unbounded_region_index

def adjacent_region_column_pairs(L):
    """
    Return a sorted list of all unordered pairs (i, j), i < j, of column indices
    (regions) that are adjacent in the spherogram.Link L.

    The region indices match the column order used by the Alexander region–crossing
    matrix functions (i.e., the order of L.faces()).
    """
    faces = L.faces()  # regions, each is a list of CrossingStrand
    # Map each half-edge (CrossingStrand) to its region index
    strand_to_region = {}
    for r_idx, face in enumerate(faces):
        for cs in face:
            strand_to_region[cs] = r_idx

    # Collect unique adjacent pairs via opposite half-edges
    pairs = set()
    for r_idx, face in enumerate(faces):
        for cs in face:
            try:
                other = cs.opposite()
            except Exception:
                # Some degeneracies may lack an opposite; skip safely
                continue
            r2 = strand_to_region.get(other)
            if r2 is None or r2 == r_idx:
                continue
            i, j = (r_idx, r2) if r_idx < r2 else (r2, r_idx)
            pairs.add((i, j))

    # Return in deterministic order

    return sorted(pairs)

def det_minor_delete_columns(M, col_pair, *, check_shape=True):
    """
    Return det of the submatrix obtained from M by deleting the two columns
    specified by the ordered pair `col_pair = (i, j)`.

    Parameters
    ----------
    M : sage.matrix.matrix_generic_dense.Matrix
        Typically your n x (n+2) Alexander region–crossing matrix (over ZZ[t^{±1}],
        ZZ[s^{±1}, t^{±1}], etc.), but works for any base ring.

    col_pair : tuple(int, int) or list[int, int]
        Ordered pair (i, j) of column indices to remove (0-based). The order does
        not affect the result, but both indices must be distinct and in range.

    check_shape : bool (default: True)
        If True, verify that the resulting submatrix is square before taking det.

    Returns
    -------
    det : element of the base ring of M
        The determinant of the (n x n) submatrix.
    """
    if not (isinstance(col_pair, (tuple, list)) and len(col_pair) == 2):
        raise ValueError("col_pair must be a tuple/list of length 2, e.g. (i, j).")
    i, j = col_pair
    nc = M.ncols()
    nr = M.nrows()
    if not (0 <= i < nc and 0 <= j < nc):
        raise IndexError(f"indices out of range for {nc} columns: got {col_pair}")
    if i == j:
        raise ValueError("indices must be distinct.")
    # Keep all columns except i and j
    to_drop = {i, j}
    keep_cols = [k for k in range(nc) if k not in to_drop]
    subM = M.matrix_from_columns(keep_cols)
    if check_shape and subM.nrows() != subM.ncols():
        raise ValueError(f"resulting submatrix is not square: {subM.nrows()}x{subM.ncols()}")
    return subM.det()  # or subM.determinant()

def determinants_for_pairs(M, pairs, *, check_shape=True):
    """
    Given a list of pairs [(i,j), ...], return a dict mapping (i,j) -> det(M without cols i,j).
    """
    out = {}
    for (i, j) in pairs:
        out[(i, j)] = det_minor_delete_columns(M, (i, j), check_shape=check_shape)
    return out

from sage.all import ZZ

def _exp_coeff_dict(f):
    """
    Robustly extract a dict {(e_s, e_t): coeff} from a Laurent polynomial
    f ∈ ZZ[s^{±1}, t^{±1}] (works with Sage's ETuple keys from f.dict()).
    """
    d = f.dict()
    out = {}
    # Try to get the expected number of variables; default to 2.
    try:
        nvars = f.parent().ngens()
    except Exception:
        nvars = 2

    for k, v in d.items():
        # Keys in multivariate (Laurent) polynomials are typically tuple-like (incl. ETuple).
        try:
            exps = tuple(int(e) for e in k)
        except TypeError:
            # As a fallback, try to coerce to tuple explicitly
            try:
                exps = tuple(int(e) for e in tuple(k))
            except Exception:
                raise TypeError("Could not extract exponent tuple from monomial key.")
        # Normalize length (pad with zeros if ring degenerates to fewer vars)
        if len(exps) < nvars:
            exps = exps + (0,)*(nvars - len(exps))
        out[exps] = ZZ(v)
    return out

def determinants_equal_up_to_pm_sj_tk(det_map):
    """
    det_map: dict mapping (i,j) -> determinant in ZZ[s^{±1}, t^{±1}].

    Returns True iff all nonzero determinants are equal to each other up to a unit
    ± s^j t^k (j,k ∈ Z). If any determinant is zero, returns True only if all are zero.
    """
    if not det_map:
        return True

    vals = list(det_map.values())

    # If any zero occurs: all must be zero.
    if any(f == 0 for f in vals):
        return all(f == 0 for f in vals)

    # Pick a nonzero reference and its exponent/coeff dict
    f0 = next(f for f in vals if f != 0)
    D0 = _exp_coeff_dict(f0)

    # Canonical base term of f0: pick lexicographically minimal exponent
    base_exp0 = min(D0.keys())
    base_coeff0 = D0[base_exp0]

    for f in vals:
        D = _exp_coeff_dict(f)

        # Same number of terms required
        if len(D) != len(D0):
            return False

        # Pick canonical base term of f
        base_expf = min(D.keys())
        base_coefff = D[base_expf]

        # Global sign must be ±1 relative to f0
        if base_coefff == base_coeff0:
            sig = 1
        elif base_coefff == -base_coeff0:
            sig = -1
        else:
            return False

        # Compute the exponent shift delta = (j,k) so that exps in f match exps in f0 + delta
        # (Assumes the ring has variables ordered (s, t).)
        delta = (base_expf[0] - base_exp0[0], base_expf[1] - base_exp0[1])

        # Check all terms align with the same shift and sign
        for e0, c0 in D0.items():
            e = (e0[0] + delta[0], e0[1] + delta[1])
            if e not in D:
                return False
            if D[e] != sig * c0:
                return False

    return True

def normalize_laurent(f):
    """
    Normalize a Laurent polynomial in Z[s, s^-1, t, t^-1] by:
    1. Removing common monomial factors s^i * t^j (only if coefficient remains integer).
    2. Making the leading term (by total degree, then s > t) have coefficient +1,
       multiplying the whole polynomial by -1 if needed.
    Does NOT factor out integer constants.
    """
    R = f.parent()
    s, t = R.gens()

    # Step 1: Find exponent-wise minimum for monomial factor
    exponents = [m.exponents()[0] for m in f.monomials()]
    min_s_exp = min(e[0] for e in exponents)
    min_t_exp = min(e[1] for e in exponents)

    # Propose common monomial
    candidate = s**min_s_exp * t**min_t_exp

    # Test if division by this monomial gives a polynomial with integer coefficients
    f_reduced = f / candidate
    if all(c in ZZ for c in f_reduced.coefficients()):
        f = f_reduced  # Accept the factor

    # Step 2: Identify leading term by total degree then s > t
    def term_key(m):
        e = m.exponents()[0]
        return (-(e[0] + e[1]), -e[0], -e[1])  # total degree, s > t

    terms = list(zip(f.monomials(), f.coefficients()))
    terms.sort(key=lambda term: term_key(term[0]))
    lead_coeff = terms[0][1]

    # Step 3: Normalize sign of leading term
    if lead_coeff == -1:
        f *= -1

    return f

def st(link, show_matrix = False):
    M, regions, colors, ub_idx = alexander_region_matrix_2var(link, debug=False, show_matrix)

    pairs = adjacent_region_column_pairs(link)
    d = det_minor_delete_columns(M, pairs[0])
    n = normalize_laurent(d)
    return n

# Rewrite documentation using reST for proper HTML formatting later

def tait_graphs(L):
    """
    Returns a pair (primal_graph, dual_graph) where primal_graph is the primal Tait graph of the link L and dual_graph is its planar dual. The primal graph corresponds to the regions whose colors are 0, or equivalently, whose variable entries in the st-matrix of L are +-t. Edges are of the form (v_1, v_2, label) where the label is a 2-element tuple, the first element of which is the crossing index and the second element of which is a boolean that is True if the edge points from the smaller region to larger region and False otherwise. This second label is necessary because the ordering of endpoint vertices is erased when the graph is created and automatically put in increasing order. 
    """
    M, _, _, _ = alexander_region_matrix_2var(L)
    nr, nc = M.nrows(), M.ncols()
    primal_edges, dual_edges = [], []
    R = LaurentPolynomialRing(ZZ, ('s', 't'))
    s, t = R.gens()
    # Iterate through crossings, adding both primal and dual edge, labeled by crossing
    for i in range(nr):
        match_regions = [[1, -1], [-1, -1], [s, -1], [-s, -1], [t, -1], [-t, -1]]
        for j in range(nc):
            for k in range(6):
                if M[i][j] == match_regions[k][0]:
                    match_regions[k][1] = j
        if match_regions[4][1] != -1:
            primal_edges.append((match_regions[4][1], match_regions[0][1], (i, match_regions[4][1] < match_regions[0][1])))
            dual_edges.append((match_regions[3][1], match_regions[1][1], (i, match_regions[3][1] < match_regions[1][1])))
        else:
            primal_edges.append((match_regions[5][1], match_regions[1][1], (i, match_regions[5][1] < match_regions[1][1])))
            dual_edges.append((match_regions[2][1], match_regions[0][1], (i, match_regions[2][1] < match_regions[0][1])))
    primal_graph = Graph(primal_edges, multiedges = True)
    dual_graph = Graph(dual_edges, multiedges = True)
    return primal_graph, dual_graph

def _make_adjacency_list(primal_graph, dual_graph):
    """
    Returns adjacency list representing both the primal graph primal_graph and its dual dual_graph outputted by tait_graphs. This is possible in a single adjacency list because the vertices/regions of primal_graph and dual_graph are complements. Each entry of the adjacency list is a list of 3-element tuples, where the first element is the neighboring vertex, the second element is the crossing index, and the third element is a boolean which is True if the edge is being traversed forwards and False otherwise.
    """
    adj_list = [[] for i in range(len(primal_graph.vertices()) + len(dual_graph.vertices()))]
    for (v_1, v_2, (crossing, small_to_large)) in primal_graph.edges() + dual_graph.edges():
        if small_to_large:
            # Track the crossing here as well just so _get_oriented_edges below can return the actual list of edges 
            adj_list[v_1].append((v_2, crossing, True))
            adj_list[v_2].append((v_1, crossing, False))
        else:
            adj_list[v_1].append((v_2, crossing, False))
            adj_list[v_2].append((v_1, crossing, True))
    return adj_list

# Need both Sage graphs and adjacency list above because spanning trees method requires using Sage graphs, but the method only
# applies to undirected graphs, so we need some other way to track directions

def _get_oriented_edges(sp_tree, root, adj_list, outward = False):
    """
    Returns a list of the edges in a spanning tree sp_tree of the directed graph represented by adj_list that are oriented towards root if outward is False and away from root if outward is True.
    """
    # Perform BFS on the adjacency list of the dual graphs, exploring edges only if they appear in sp_tree
    edges_in_sp_tree = set(sp_tree.edges())
    discovered = [False for i in range(len(adj_list))]
    oriented_edges = []
    heap = []
    discovered[root] = True
    heapq.heappush(heap, root)
    while len(heap) > 0:
        cur = heapq.heappop(heap)
        for (neighbor, crossing, forward) in adj_list[cur]:
            edge = _adj_entry_to_edge(cur, neighbor, crossing, forward)
            if not discovered[neighbor] and edge in edges_in_sp_tree:
                discovered[neighbor] = True
                heapq.heappush(heap, neighbor)
                if forward and outward:
                    oriented_edges.append(edge)
                elif not forward and not outward:
                    oriented_edges.append(edge)
    return oriented_edges

def _adj_entry_to_edge(cur, neighbor, crossing, forward):
    """
    Returns the edge between cur and neighbor over crossing, where forward is True if the edge points from cur to neighbor and False otherwise.
    """
    small_to_large = True if (cur < neighbor and forward) or (cur > neighbor and not forward) else False
    return (cur, neighbor, (crossing, small_to_large)) if cur < neighbor else (neighbor, cur, (crossing, small_to_large))

def _get_dual_tree(sp_tree, dual_graph):
    """
    Returns the dual of sp_tree in the dual graph dual_graph.
    """
    used_crossings = set()
    for (v_1, v_2, (crossing, small_to_large)) in sp_tree.edges():
        used_crossings.add(crossing)
    dual_tree_edges = [(v_1, v_2, (crossing, small_to_large)) for (v_1, v_2, (crossing, small_to_large)) in dual_graph.edges() if crossing not in used_crossings]
    return Graph(dual_tree_edges, multiedges = True)

# Change description of labeling of the regions to be clearer throughout the functions?
def get_all_sp_tree_terms(L, primal_outward = False, dual_outward = False):
    """
    Returns a list of 2-element tuples (s_power, t_power) that can be obtained by some pair of dual spanning trees in the Tait graphs associated to the link L, where s_power is from the dual graph (1-colored regions) and t_power from the primal (0-colored regions) with the region touching the most crossings being 0-colored. In each spanning tree, counts the number of edges pointing towards the root, unless specified otherwise by primal_outward and dual_outward. The output tuples are sorted lexicographically.
    """
    sp_tree_tuples = set()
    primal_graph, dual_graph = tait_graphs(L)
    adj_list = _make_adjacency_list(primal_graph, dual_graph)
    root_1, root_2 = adjacent_region_column_pairs(L)[0]
    primal_root, dual_root = ((root_1, root_2) if root_1 in primal_graph.vertices() else (root_2, root_1))
    for sp_tree in primal_graph.spanning_trees(labels = True):
        t_power = len(_get_oriented_edges(sp_tree, primal_root, adj_list, primal_outward))
        dual_sp_tree = _get_dual_tree(sp_tree, dual_graph)
        s_power = len(_get_oriented_edges(dual_sp_tree, dual_root, adj_list, dual_outward))
        sp_tree_tuples.add((s_power, t_power))
    return sorted(list(sp_tree_tuples))

def st_poly_from_sp_trees(L, primal_outward = False, dual_outward = False):
    """
    Returns the st-polynomial of L using the method of dual spanning trees. In each spanning tree, counts the number of edges pointing towards the root, unless specified otherwise by primal_outward and dual_outward. Powers of s or t that appear in all terms are removed.

    TO BE FIXED: right now, this method gives the same positive sign to all terms due to sign not being accounted for in the helper functions. So the result is NOT the exact st-polynomial.
    """
    R = LaurentPolynomialRing(ZZ, ('s', 't'))
    s, t = R.gens()
    st_polynomial = 0
    primal_graph, dual_graph = tait_graphs(L)
    adj_list = _make_adjacency_list(primal_graph, dual_graph)
    root_1, root_2 = adjacent_region_column_pairs(L)[0]
    primal_root, dual_root = ((root_1, root_2) if root_1 in primal_graph.vertices() else (root_2, root_1))
    for sp_tree in primal_graph.spanning_trees(labels = True):
        t_power = len(_get_oriented_edges(sp_tree, primal_root, adj_list, primal_outward))
        dual_sp_tree = _get_dual_tree(sp_tree, dual_graph)
        s_power = len(_get_oriented_edges(dual_sp_tree, dual_root, adj_list, dual_outward))
        st_polynomial += s^s_power * t^t_power
    return normalize_laurent(st_polynomial)

def sp_tree_terms_to_st_poly(sp_tree_tuples):
    """
    Returns the st-polynomial one would get from sp_tree_tuples if all terms appearing from some spanning tree had coefficient 1, with powers of s or t that appear in all terms removed. 
    """
    R = LaurentPolynomialRing(ZZ, ('s', 't'))
    s, t = R.gens()
    st_polynomial = 0
    for (s_power, t_power) in sp_tree_tuples:
        st_polynomial += s^s_power * t^t_power
    return normalize_laurent(st_polynomial)
