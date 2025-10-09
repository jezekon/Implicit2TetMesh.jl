using Test
using JLD2
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using Implicit2TetMesh.Optimization
using Implicit2TetMesh.Modification
using Implicit2TetMesh.Utils


@testset "Implicit2TetMesh.jl" begin

    RUN_beam = false
    RUN_main = false
    RUN_main_param = false
    RUN_chapadlo = true

    if RUN_chapadlo
        function define_boundary_planes()
            planes = PlaneDefinition[]

            # Plane 1: YZ Symmetry Plane (x = 0)
            # This plane enforces symmetry boundary condition along the X-axis
            # Normal vector points in +X direction, ensuring proper constraint orientation
            push!(planes, PlaneDefinition(
                [1.0, 0.0, 0.0],    # normal vector pointing in +X direction
                [0.0, 0.0, 0.0],    # reference point at origin
                Square(400.0),        # large square to cover entire YZ cross-section
            ))

            # Plane 2: XZ Constrained Region (y ≈ 75.024)
            # This plane defines the contact/constraint region for the gripper mechanism
            # Elliptical area covers the expected contact zone: z ∈ (92, 137), x ∈ (0, 16)
            push!(
                planes,
                PlaneDefinition(
                    [0.0, 1.0, 0.0],     # normal vector pointing in +Y direction
                    [0.0, 75.0, 114.5], # center point of constrained elliptical region
                    Ellipse(18.0, 16.0),   # ellipse: a=8 (x-semi-axis), b=22.5 (z-semi-axis)
                ),
            )

            # Plane 3: XY Force Application Plane (z = 90)
            # This plane represents the surface where external forces are applied
            # Large square ensures coverage of the entire force application area
            push!(
                planes,
                PlaneDefinition(
                    [0.0, 0.0, 1.0],    # normal vector pointing in +Z direction
                    [0.0, 0.0, -90.0],   # reference point on the force application plane
                    Square(200.0),        # large square to cover XY cross-section
                ),
            )

            return planes
        end

        # Print header information
        println("="^60)
        println("Tetrahedral Mesh Generation for Chapadlo Geometry")
        println("="^60)

        # Define input data file paths
        # These files contain the fine grid coordinates and SDF values for the gripper geometry
        grid_file = "../data/chapadlo/Z_chapadlo_FineGrid_B-2.0085_smooth-1_Interpolation.jld2"
        sdf_file = "../data/chapadlo/Z_chapadlo_FineSDF_B-2.0085_smooth-1_Interpolation.jld2"

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
            warp_param = 0.3,            # Small warping parameter for precise boundary alignment
            plane_definitions = plane_definitions,  # Apply boundary plane constraints
            quality_export = true,       # Export detailed quality metrics for analysis
            optimize = false,             # Enable mesh optimization for better element quality
            split_elements = true,       # Use element splitting for accurate isosurface representation
            correct_volume = true,        # Apply volume correction to match reference geometry
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
                options = options,
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

    if RUN_beam

        taskName = "beam_vfrac_04_B-1.0_smooth-1_Interp"

        @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2" fine_grid
        @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2" fine_sdf

        scheme = "A15"
        # scheme = "Schlafli"

        # Definice rovin
        plane_definitions = [
            PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
            PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0)),
        ]
        warp_param = 0.3

        mesh = BlockMesh(fine_sdf, fine_grid)

        # Choose scheme: "A15" or "Schlafli"
        generate_mesh!(mesh, scheme)

        warp!(mesh, scheme) # Warp nearest nodes to the isocontour
        # export_mesh_vtu_quality(mesh, "$(taskName)_TriMesh-step_warp$(scheme).vtu")

        slice_ambiguous_tetrahedra!(mesh) # Remove elements outside the body
        # adjust_nodes_to_isosurface!(mesh) # Simple cut of elements to follow the isocontour

        update_connectivity!(mesh)

        export_mesh_vtu_quality(mesh, "$(taskName)_TriMesh-notOPT_$(scheme).vtu")

        TetMesh_volumes(mesh)

        remove_inverted_elements!(mesh)
        update_connectivity!(mesh)
        # optimize_mesh!(mesh, scheme)

        # Korekce objemu
        success = correct_mesh_volume!(
            mesh,
            fine_sdf,
            fine_grid,
            scheme,
            plane_definitions = plane_definitions,
        )
        remove_inverted_elements!(mesh)

        # Vyhodnocení přesnosti
        # assess_volume_accuracy(mesh, fine_sdf, fine_grid)

        export_mesh_vtu(mesh, "$(taskName)_TriMesh-$(scheme).vtu")
        # export_mesh_vtu_quality(mesh, "$(taskName)_TriMesh-volume_modif_$(scheme).vtu")
    end

    #     # Apply cutting planes only if they are defined
    #     if !isempty(plane_definitions)
    #       warp_mesh_by_planes_sdf!(mesh, plane_definitions, warp_param)
    #       update_connectivity!(mesh)
    #       export_mesh_vtu(mesh, "$(taskName)_TriMesh-$(scheme)_cut.vtu")
    #     end
    #
    #     TetMesh_volumes(mesh)
    #
    #     slice_mesh_with_plane!(mesh, "x", 0.977, -1)
    #     export_mesh_vtu(mesh, "$(taskName)_TriMesh-$(scheme)_plane.vtu")
    #     # assess_mesh_quality(mesh, "mesh_quality")
    # end

    if RUN_main
        mesh = generate_tetrahedral_mesh(
            "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
            "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
            "cantilever_beam_interp",
        )
    end

    if RUN_main_param

        plane_definitions = [
            PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
            PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0)),
        ]

        options=MeshGenerationOptions(
            scheme = "A15",
            optimize = true,
            split_elements = false,
            quality_export = true,
            warp_param = 0.0,
            plane_definitions = plane_definitions,
        )

        mesh = generate_tetrahedral_mesh(
            "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
            "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
            "cantilever_beam_interp_cut",
            options = options,
        )

        # assess_mesh_quality(mesh, "mesh_quality")
    end
end
