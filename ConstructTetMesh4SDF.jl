using Implicit2TetMesh

plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0))
]

mesh = generate_tetrahedral_mesh(
    "data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
    "data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
    "cantilever_beam_interp_cut",
    options=MeshGenerationOptions(
        scheme = "A15",
        warp_param = 0.3,
        plane_definitions = plane_definitions
    )
)

