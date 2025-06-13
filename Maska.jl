using Implicit2TetMesh

plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0))
]

mesh = generate_tetrahedral_mesh(
    "Maska_grid.jld2",
    "Maska_sdf.jld2",
    "Maska",
    options=MeshGenerationOptions(
        scheme = "A15",
    )
)

