from sage.all import text3d

def newton_polytope(f):
    """
    Given a Laurent polynomial f in ZZ[s, s^-1, t, t^-1],
    return its Newton polytope as a 3D polytope lying in the plane
    x + y + z = deg(f), where each monomial s^a t^b is mapped to
    (a, b, deg(f) - a - b).
    """
    if f.is_zero():
        raise ValueError("Zero polynomial has no Newton polytope.")

    # Extract monomials with nonzero coefficients
    exponents = [tuple(m) for m, c in f.dict().items() if c != 0]

    # Compute degree as the max total degree (a + b) among terms
    degrees = [sum(m) for m in exponents]
    d = max(degrees)

    # Embed each exponent (a, b) as (a, b, d - a - b)
    lifted_points = [(a, b, d - a - b) for (a, b) in exponents]

    # Return the convex hull
    return Polyhedron(vertices=lifted_points)

def plot_conv_hull(points, labels):
    """
    Given a list of points in 3-space and their corresponding labels, creates a plot with labeled axes of the polytope that is the convex hull of the points. Points are labeled by their label.
    """
    # Create a 3D plot of the points
    point_plot = point3d(points, size=25, color='blue')
    
    # Create the Newton polytope as the convex hull of the lifted points
    p = Polyhedron(vertices=points)
    poly_plot = p.plot(opacity=0.4, color='green')

    # Create text labels for each point
    label_offset = (0.2, 0.2, 0.2) 
    labels_plot = sum(
        text3d(
            str(c), 
            (pt[0] + label_offset[0], pt[1] + label_offset[1], pt[2] + label_offset[2]), 
            color='black'
        ) 
        for c, pt in zip(labels, points)
    )
    
    combined_plot = point_plot + poly_plot + labels_plot
    
    # Determine axis ranges to create prominent tick marks
    if points:
        x_coords = [p[0] for p in points]
        y_coords = [p[1] for p in points]
        z_coords = [p[2] for p in points]
        
        # Define the range with a small buffer around the points
        min_x, max_x = min(x_coords) - 1, max(x_coords) + 1
        min_y, max_y = min(y_coords) - 1, max(y_coords) + 1
        min_z, max_z = min(z_coords) - 1, max(z_coords) + 1
        
        # Generate ticks at integer intervals across the range
        x_ticks = list(range(floor(min_x), ceil(max_x) + 1))
        y_ticks = list(range(floor(min_y), ceil(max_y) + 1))
        z_ticks = list(range(floor(min_z), ceil(max_z) + 1))
    else:
        # Default ticks for an empty plot
        x_ticks, y_ticks, z_ticks = [], [], []

    # Manually create prominent axes, ticks, and labels
    axes_plot = None
    if points:
        x_coords = [p[0] for p in points]
        y_coords = [p[1] for p in points]
        z_coords = [p[2] for p in points]
        
        # Define the range with a buffer around the points
        min_x, max_x = min(x_coords) - 1, max(x_coords) + 1
        min_y, max_y = min(y_coords) - 1, max(y_coords) + 1
        min_z, max_z = min(z_coords) - 1, max(z_coords) + 1
        
        # Create axis lines
        x_axis = line3d([(min_x, 0, 0), (max_x, 0, 0)], color='red', thickness=2)
        y_axis = line3d([(0, min_y, 0), (0, max_y, 0)], color='green', thickness=2)
        z_axis = line3d([(0, 0, min_z), (0, 0, max_z)], color='blue', thickness=2)
        
        # Create axis labels
        x_label = text3d('x', (max_x + 0.5, 0, 0), color='red')
        y_label = text3d('y', (0, max_y + 0.5, 0), color='green')
        z_label = text3d('z', (0, 0, max_z + 0.5), color='blue')
        
        axes_plot = x_axis + y_axis + z_axis + x_label + y_label + z_label
        
        # Generate ticks and tick labels at integer intervals
        tick_size = 0.1
        # X-axis ticks
        for t in range(floor(min_x), ceil(max_x) + 1):
            if t == 0: continue
            axes_plot += line3d([(t, -tick_size, 0), (t, tick_size, 0)], color='red', thickness=1)
            axes_plot += text3d(str(t), (t, -tick_size*2.5, 0), color='red')
        # Y-axis ticks
        for t in range(floor(min_y), ceil(max_y) + 1):
            if t == 0: continue
            axes_plot += line3d([(-tick_size, t, 0), (tick_size, t, 0)], color='green', thickness=1)
            axes_plot += text3d(str(t), (-tick_size*2.5, t, 0), color='green')
        # Z-axis ticks
        for t in range(floor(min_z), ceil(max_z) + 1):
            if t == 0: continue
            axes_plot += line3d([(0, -tick_size, t), (0, tick_size, t)], color='blue', thickness=1)
            axes_plot += text3d(str(t), (0, -tick_size*2.5, t), color='blue')

    if axes_plot:
        combined_plot += axes_plot
    
    combined_plot.show(frame = False)

# Add monomials onto labels?
def plot_newton_polytope(f):
    """
    Given a Laurent polynomial f in ZZ[s, s^-1, t, t^-1], plot its
    Newton polytope as a 3D polytope lying in the plane x + y + z = deg(f)
    along with the lattice points inside the polytope labeled with the 
    polynomial's coefficients. Each monomial s^a t^b is mapped to 
    (a, b, deg(f) - a - b).

    The plot includes labeled axes (x, y, z) with tick marks for clarity.
    """
    if f.is_zero():
        raise ValueError("The zero polynomial does not have a Newton polytope.")

    # Extract monomials and their corresponding non-zero coefficients
    poly_dict = f.dict()
    exponents = [tuple(m) for m, c in poly_dict.items() if c != 0]
    coefficients = [c for m, c in poly_dict.items() if c != 0]

    # Handle the case of a constant polynomial
    if not exponents:
        if f.is_constant():
            d = 0
            exponents = [(0, 0)]
            coefficients = [f.constant_coefficient()]
        else:
            print("Warning: No terms with non-zero coefficients found.")
            return

    # Compute degree as the max total degree (a + b) among terms
    degrees = [sum(m) for m in exponents]
    d = max(degrees) if degrees else 0

    # Embed each 2D exponent (a, b) as a 3D point (a, b, d - a - b)
    lifted_points = [(a, b, d - a - b) for (a, b) in exponents]
    
    plot_conv_hull(lifted_points, coefficients)

def plot_diffs_polytope(f, g):
    """
    Given Laurent polynomials f, g in ZZ[s, s^-1, t, t^-1], plot the common 
    Newton polytope as a 3D polytope lying in the plane x + y + z = deg(f)
    along with the lattice points inside the polytope labeled with the net change in coefficient 
    going from f to g. Each monomial s^a t^b is mapped to (a, b, deg(f) - a - b).

    The plot includes labeled axes (x, y, z) with tick marks for clarity.

    Requires: f and g have the same Newton polytope.
    """
    # Extract monomials, non-zero coefficients, and differences in coefficients
    f_exponents = [tuple(m) for m, c in f.dict().items() if c != 0]
    g_exponents = [tuple(m) for m, c in g.dict().items() if c != 0]
    f_coeffs = [c for m, c in f.dict().items() if c != 0]
    g_coeffs = [c for m, c in g.dict().items() if c != 0]
    if len(f_exponents) != len(g_exponents):
        return False
    coeffs_diffs = [g_coeffs[i] - f_coeffs[i] for i in range(len(f_coeffs))]

    # Compute degree as the max total degree (a + b) among terms
    degrees = [sum(m) for m in f_exponents]
    d = max(degrees) if degrees else 0

    # Embed each 2D exponent (a, b) as a 3D point (a, b, d - a - b)
    lifted_points = [(a, b, d - a - b) for (a, b) in f_exponents]
    
    plot_conv_hull(lifted_points, coeffs_diffs)