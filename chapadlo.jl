#!/usr/bin/env julia

"""
Tetrahedral Mesh Generation Script for Chapadlo (Gripper) Geometry

This script generates a tetrahedral mesh from SDF (Signed Distance Function) data
for a gripper mechanism with specific boundary plane constraints.

"""

using Implicit2TetMesh

"""
    define_boundary_planes() -> Vector{PlaneDefinition}

Define boundary planes that constrain the mesh generation process.
These planes represent physical constraints and symmetry conditions for the gripper geometry.

Returns:
    Vector{PlaneDefinition}: Array of plane definitions for mesh constraints
"""
function define_boundary_planes()
    planes = PlaneDefinition[]
    
    # Plane 1: YZ Symmetry Plane (x = 0)
    # This plane enforces symmetry boundary condition along the X-axis
    # Normal vector points in +X direction, ensuring proper constraint orientation
    push!(planes, PlaneDefinition(
        [1.0, 0.0, 0.0],    # normal vector pointing in +X direction
        [0.0, 0.0, 0.0],    # reference point at origin
        Square(400.0)        # large square to cover entire YZ cross-section
    ))
    
    # Plane 2: XZ Constrained Region (y ≈ 75.024)
    # This plane defines the contact/constraint region for the gripper mechanism
    # Elliptical area covers the expected contact zone: z ∈ (92, 137), x ∈ (0, 16)
    push!(planes, PlaneDefinition(
        [0.0, 1.0, 0.0],     # normal vector pointing in +Y direction
        [8.0, 75.0, 114.5], # center point of constrained elliptical region
        Ellipse(8.0, 18.)   # ellipse: a=8 (x-semi-axis), b=22.5 (z-semi-axis)
    ))
    
    # Plane 3: XY Force Application Plane (z = 90)
    # This plane represents the surface where external forces are applied
    # Large square ensures coverage of the entire force application area
    push!(planes, PlaneDefinition(
        [0.0, 0.0, 1.0],    # normal vector pointing in +Z direction
        [0.0, 0.0, -90.0],   # reference point on the force application plane
        Square(200.0)        # large square to cover XY cross-section
    ))
    
    return planes
end

"""
    main()

Main function that orchestrates the tetrahedral mesh generation process.
Configures mesh generation options and executes the discretization workflow.
"""
function main()
    # Print header information
    println("="^60)
    println("Tetrahedral Mesh Generation for Chapadlo Geometry")
    println("="^60)
    
    # Define input data file paths
    # These files contain the fine grid coordinates and SDF values for the gripper geometry
    grid_file = "data/chapadlo/Z_chapadlo_FineGrid_B-2.0085_smooth-1_Interpolation.jld2"
    sdf_file = "data/chapadlo/Z_chapadlo_FineSDF_B-2.0085_smooth-1_Interpolation.jld2"
    
    # Define output file prefix (without extension - .vtu will be added automatically)
    output_prefix = "tet_chapadlo_B-2.0"
    
    # Get boundary plane definitions for mesh constraints
    plane_definitions = define_boundary_planes()
    
    println("Configuration:")
    println("  Grid file: $grid_file")
    println("  SDF file:  $sdf_file")
    println("  Output:    $output_prefix")
    println("  Boundary planes: $(length(plane_definitions)) defined")
    println()
    
    # Configure mesh generation options with specified parameters
    options = MeshGenerationOptions(
        scheme = "A15",              # Use A15 discretization scheme for better surface representation
        warp_param = 0.1,            # Small warping parameter for precise boundary alignment
        plane_definitions = plane_definitions,  # Apply boundary plane constraints
        quality_export = true,       # Export detailed quality metrics for analysis
        optimize = false,             # Enable mesh optimization for better element quality
        split_elements = true,       # Use element splitting for accurate isosurface representation
        correct_volume = true        # Apply volume correction to match reference geometry
    )
    
    # Display mesh generation settings
    println("Mesh Generation Settings:")
    println("  Discretization scheme: $(options.scheme)")
    println("  Warp parameter: $(options.warp_param)")
    println("  Quality export: $(options.quality_export)")
    println("  Mesh optimization: $(options.optimize)")
    println("  Element splitting: $(options.split_elements)")
    println("  Volume correction: $(options.correct_volume)")
    println()
    
    # Execute tetrahedral mesh generation
    println("Starting mesh generation process...")
    try
        mesh = generate_tetrahedral_mesh(
            grid_file,
            sdf_file,
            output_prefix;
            options = options
        )
        
        # Success message with mesh statistics
        println()
        println("✓ Mesh generation completed successfully!")
        println("  Final mesh statistics:")
        println("    Nodes: $(length(mesh.X))")
        println("    Elements: $(length(mesh.IEN))")
        println("  Output files generated:")
        println("    - $(output_prefix)_TriMesh-$(options.scheme).vtu")
        if options.plane_definitions !== nothing && options.warp_param > 0.0
            println("    - $(output_prefix)_TriMesh-$(options.scheme)_cut.vtu")
        end
        
    catch e
        # Error handling with detailed message
        println("✗ Mesh generation failed with error:")
        println("  Error type: $(typeof(e))")
        println("  Error message: $e")
        
        # Print stack trace for debugging
        println("\nStack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
    
    println("="^60)
end

# Execute main function when script is run directly
# if abspath(PROGRAM_FILE) == @__FILE__
    main()
# end
