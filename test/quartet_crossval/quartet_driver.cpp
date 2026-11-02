// Cross-validation driver: runs quartet's stuffing steps (unchanged quartet code, included
// as a source file so its static functions are reachable) on an SDF grid read from a file,
// and dumps the mesh after each step in full float precision.
//   input  (binary): int ni, nj, nk; float ox, oy, oz, dx; float phi[ni*nj*nk] (i fastest)
//   output (binary, one file per step <prefix>_<step>.bin):
//          int nv, nt; float X[3*nv]; float vphi[nv]; int T[4*nt]
// optimize = false, no feature matching (same as make_tet_mesh(mesh, sdf, false)).
#include "make_tet_mesh.cpp"
#include <cstdio>
#include <string>

static void dump(const std::string& file, const TetMesh& mesh, const std::vector<float>& vphi)
{
    FILE* f = fopen(file.c_str(), "wb");
    int nv = (int)mesh.vSize(), nt = (int)mesh.tSize();
    fwrite(&nv, sizeof(int), 1, f);
    fwrite(&nt, sizeof(int), 1, f);
    for (int i = 0; i < nv; ++i) fwrite(&mesh.V(i)[0], sizeof(float), 3, f);
    fwrite(vphi.data(), sizeof(float), nv, f);
    for (int t = 0; t < nt; ++t) fwrite(&mesh.T(t)[0], sizeof(int), 4, f);
    fclose(f);
    printf("  %-28s %9d vertices %9d tets\n", file.c_str(), nv, nt);
}

int main(int argc, char** argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s <sdf.bin> <out prefix>\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb");
    int n[3]; float o[4];
    if (fread(n, sizeof(int), 3, f) != 3 || fread(o, sizeof(float), 4, f) != 4) return 2;
    SDF sdf(Vec3f(o[0], o[1], o[2]), o[3], n[0], n[1], n[2]);
    for (int k = 0; k < n[2]; ++k) for (int j = 0; j < n[1]; ++j) for (int i = 0; i < n[0]; ++i)
        if (fread(&sdf.phi(i, j, k), sizeof(float), 1, f) != 1) return 3;
    fclose(f);

    std::string pre(argv[2]);
    initialize_exact();
    TetMesh mesh;
    std::vector<float> vphi;
    cut_lattice(sdf, mesh, vphi);            dump(pre + "_1_lattice.bin", mesh, vphi);
    warp_vertices(0.3f, mesh, vphi);         dump(pre + "_2_warp.bin", mesh, vphi);
    trim_spikes(mesh, vphi);                 dump(pre + "_3_trim.bin", mesh, vphi);
    remove_exterior_tets(mesh, vphi, sdf);   dump(pre + "_4_exterior.bin", mesh, vphi);
    return 0;
}
