# ----------------------------
# Definition of Schlafli orthoscheme connectivity
# ----------------------------
const schlafli_tet_connectivity = [
  [1, 2, 3, 7],  # Path 1: x, y, z
  [1, 6, 2, 7],  # Path 2: x, z, y
  [1, 3, 4, 7],  # Path 3: y, x, z
  [1, 4, 8, 7],  # Path 4: y, z, x
  [1, 5, 6, 7],  # Path 5: z, x, y
  [1, 8, 5, 7]   # Path 6: z, y, x
]
