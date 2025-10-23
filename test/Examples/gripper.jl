# ==============================================================================
# Example: Tetrahedral Mesh Generation for Gripper Geometry
# ==============================================================================

using Implicit2TetMesh

# ------------------------------------------------------------------------------
# Step 1: Define Boundary Planes
# ------------------------------------------------------------------------------
# Boundary planes are used to align mesh nodes to specific geometric features,
# which are important for applying boundary conditions in FEM analysis.
#
# PlaneDefinition(normal, point, shape):
#   - normal: Vector defining plane orientation [x, y, z]
#   - point:  Reference point on the plane [x, y, z]
#   - shape:  Geometric shape defining the plane's bounded area

plane_definitions = [
    # YZ Symmetry Plane (x = 0):
    # Normal [1,0,0] points in +X direction
    # Square(400.0) creates a large bounded region to cover entire YZ cross-section
    PlaneDefinition([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], Square(400.0)),

    # XZ Constrained Region (y ≈ 75):
    # Normal [0,1,0] points in +Y direction
    # Ellipse(8.0, 18.0) covers the contact/constraint zone
    PlaneDefinition([0.0, 1.0, 0.0], [8.0, 75.0, 114.5], Ellipse(8.0, 18.0)),

    # XY Force Application Plane (z = -90):
    # Normal [0,0,1] points in +Z direction
    # Square(200.0) covers the entire force application area
    PlaneDefinition([0.0, 0.0, 1.0], [0.0, 0.0, -90.0], Square(200.0)),
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

const DATA_DIR = "data/gripper"
grid_file = joinpath(DATA_DIR, "Z_robot_gripper_HEX8_FineGrid_B-2.0085_smooth-1.jld2")
sdf_file = joinpath(DATA_DIR, "Z_robot_gripper_HEX8_FineSDF_B-2.0085_smooth-1.jld2")

mesh = generate_tetrahedral_mesh(
    # Input files:
    grid_file,  # Grid coordinates (fine_grid array)
    sdf_file,   # SDF values (fine_sdf array)

    # Output file prefix (without extension):
    # Will generate: gripper_TriMesh-A15.vtu
    #                gripper_TriMesh-A15_cut.vtu
    "gripper",

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
        # 0.1  - Nodes within 10% of grid step from planes are warped (precise)
        # 0.5  - More aggressive warping, larger influence zone
        # 1.0+ - Very aggressive, may affect mesh quality
        warp_param = 0.1,

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
# 1. tet_chapadlo_B-2.0_TriMesh-A15.vtu
#    - Initial mesh after isosurface extraction
#    - Contains SDF values at nodes
#
# 2. tet_chapadlo_B-2.0_TriMesh-A15_cut.vtu (if plane_definitions are provided)
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
println("  paraview tet_chapadlo_B-2.0_TriMesh-A15_cut.vtu")
