# ==============================================================================
# Example: Tetrahedral Mesh Generation from SDF Data
# ==============================================================================

using Implicit2TetMesh

# ------------------------------------------------------------------------------
# Step 1: Define Boundary Planes (optional)
# ------------------------------------------------------------------------------
# Boundary planes are used to align mesh nodes to specific geometric features,
# which are important for applying boundary conditions in FEM analysis.
#
# PlaneDefinition(normal, point, shape):
#   - normal: Vector defining plane orientation [x, y, z]
#   - point:  Reference point on the plane [x, y, z]
#   - shape:  Geometric shape defining the plane's bounded area

plane_definitions = [
    # Left boundary plane (x = 0):
    # Normal [-1,0,0] points in -X direction (into the domain)
    # Square(30.0) creates a 30×30 bounded region for constraint application
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),

    # Right boundary plane (x = 60):
    # Normal [1,0,0] points in +X direction (into the domain)
    # Square(5.0) creates a 5×5 bounded region (e.g., for load application)
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0)),
]

# ------------------------------------------------------------------------------
# Step 2: Generate Tetrahedral Mesh
# ------------------------------------------------------------------------------
# Main function that performs the complete mesh generation workflow:
# 1. Load SDF and grid data from JLD2 files
# 2. Generate initial mesh using selected discretization scheme
# 3. Warp nodes to isosurface (SDF = 0)
# 4. Slice ambiguous elements crossing the boundary
# 5. Apply boundary plane constraints
# 6. Correct mesh volume to match reference geometry
# 7. Export mesh to VTU format for visualization

const DATA_DIR = "data/beam"
grid_file = joinpath(DATA_DIR, "Z_beam_HEX8_FineGrid_B-1.0_smooth-1.jld2")
sdf_file = joinpath(DATA_DIR, "Z_beam_HEX8_FineSDF_B-1.0_smooth-1.jld2")

mesh = generate_tetrahedral_mesh(
    # Input files:
    grid_file,  # Grid coordinates (fine_grid array)
    sdf_file,   # SDF values (fine_sdf array)

    # Output file prefix (without extension):
    # Will generate: cantilever_beam_TriMesh-A15.vtu
    #                cantilever_beam_TriMesh-A15_cut.vtu
    "cantilever_beam",

    # Mesh generation options:
    options = MeshGenerationOptions(
        # ------------------------------------------------------------------
        # scheme: Tetrahedral discretization scheme
        # ------------------------------------------------------------------
        # "A15"      - A15 lattice structure (27 nodes per unit cell, 46 tets)
        #              Better surface representation, recommended for most cases
        # "Schlafli" - Schlafli orthoscheme (8 nodes per unit cell, 6 tets)
        #              Simpler structure, fewer elements
        scheme = "A15",

        # ------------------------------------------------------------------
        # warp_param: Surface node warping threshold (0.0 - 1.0+)
        # ------------------------------------------------------------------
        # Determines how far from boundary planes nodes will be warped:
        #   threshold_distance = warp_param × mesh.grid_step
        #
        # 0.0  - No warping, planes are ignored
        # 0.3  - Nodes within 30% of grid step from planes are warped (recommended)
        # 0.5  - More aggressive warping, larger influence zone
        # 1.0+ - Very aggressive, may affect mesh quality
        warp_param = 0.3,

        # ------------------------------------------------------------------
        # plane_definitions: Boundary plane constraints
        # ------------------------------------------------------------------
        # When provided, surface nodes near these planes will be warped to
        # align exactly with the plane surfaces. Essential for:
        #   - Applying Dirichlet boundary conditions (fixed displacements)
        #   - Defining load application surfaces
        #   - Enforcing symmetry planes
        plane_definitions = plane_definitions,

        # ------------------------------------------------------------------
        # correct_volume: Enable volume correction (true/false)
        # ------------------------------------------------------------------
        # true  - Adjusts surface node positions to match the reference volume
        #         from SDF data using a bisection method. Improves geometric
        #         accuracy at the cost of computation time (~20 iterations).
        # false - Skip volume correction (faster, less accurate)
        correct_volume = true,

        # ------------------------------------------------------------------
        # experimental_nzzz: Enable experimental NZZZ case warping
        # ------------------------------------------------------------------
        # Controls handling of elements with 1 negative + 3 zero SDF nodes:
        # true  - Uses relaxed safety parameters for NZZZ case warping
        #         (may improve surface quality but less conservative)
        # false - Conservative mode, discards NZZZ elements (safer, default)
        experimental_nzzz = true,
    ),
)

# ------------------------------------------------------------------------------
# Output Description
# ------------------------------------------------------------------------------
# The function generates two VTU files for ParaView visualization:
#
# 1. cantilever_beam_TriMesh-A15.vtu
#    - Initial mesh after isosurface extraction
#    - Contains SDF values at nodes
#
# 2. cantilever_beam_TriMesh-A15_cut.vtu (if plane_definitions are provided)
#    - Mesh with nodes warped to boundary planes
#    - Volume-corrected if correct_volume = true
#
# The returned 'mesh' object (BlockMesh type) contains:
#   mesh.X        - Node coordinates (Vector{SVector{3,Float64}})
#   mesh.IEN      - Element connectivity (Vector{Vector{Int64}})
#   mesh.node_sdf - SDF values at nodes (Vector{Float64})
#   mesh.INE      - Inverse connectivity (node → elements)

# ------------------------------------------------------------------------------
# Quality Metrics (printed during execution)
# ------------------------------------------------------------------------------
# The generation process reports:
# - Number of nodes and elements at each stage
# - Volume statistics and problematic elements
# - Surface node warping results
# - Volume correction convergence (if enabled)
# - Inverted element detection and fixing

# ------------------------------------------------------------------------------
# Typical Workflow After Mesh Generation
# ------------------------------------------------------------------------------
# 1. Visualize in ParaView:
#    - Load .vtu file
#    - Apply 'Surface' or 'Surface with edges' to see the boundary
#    - Color by 'sdf' field to verify isosurface
#
# 2. Export for FEM solver
#
# 3. Quality assessment:
#    using Implicit2TetMesh.Utils
#    assess_mesh_quality(mesh, "quality_report")
#    # Generates dihedral angle histogram and quality metrics

println("\n✓ Mesh generation complete!")
println("  Total nodes:    $(length(mesh.X))")
println("  Total elements: $(length(mesh.IEN))")
println("\nVisualize results in ParaView:")
println("  paraview cantilever_beam_TriMesh-A15_cut.vtu")
