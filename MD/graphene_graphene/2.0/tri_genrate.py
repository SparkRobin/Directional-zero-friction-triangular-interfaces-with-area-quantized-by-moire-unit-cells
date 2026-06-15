"""
Generate a 30deg right-angled triangular graphene-related structure and export to a LAMMPS data file.

1. Build a graphene-like lattice from two sublattices.
2. Rotate one layer by twist angle theta determined by (m, n).
3. Cut the rotated layer using a triangular / 30-degree right-angle boundary.
4. Plot bottom layer and cut top layer.
5. Assign z coordinates to bottom layer, top flake, and staged top flake.
6. Export structure into a LAMMPS .in data file.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path


def main():
    # =========================
    # Basic geometric parameters
    # =========================
    aCC = 1.42                  # C-C bond length
    dy = aCC * 1.5
    dx = aCC * np.sqrt(3)

    gra_b_list = []

    xi = 2

    # adjust for boundary
    imax = xi*50+30
    jmax = int(imax * 1.5)
    
    # Chiral/twist indices (m, n)
    m = 22
    n = m + 1

    # Twist angle formula
    theta = np.arccos(
        0.5 * (m**2 + n**2 + 4 * m * n) / (m**2 + n**2 + m * n)
    )
    theta_deg = theta / np.pi * 180.0

    print(f"xi = {xi}")
    print(f"m = {m}, n = {n}")
    print(f"theta_deg = {theta_deg:.10f}")

    # =========================================
    # Construct the first sublattice-like points
    # =========================================
    for i in range(-20, imax + 2):        # MATLAB: -20:imax+1
        for j in range(-20, jmax + 1):    # MATLAB: -20:jmax
            if j % 2 == 1:
                x = (i - 1) * dx
                y = (j - 1) * dy
            else:
                # Construct triangle lattice
                x = (i - 1) * dx + dx / 2
                y = (j - 1) * dy
            gra_b_list.append([x, y])

    gra_b = np.array(gra_b_list, dtype=float)

    # =========================================
    # Construct graphene lattice (second basis)
    # =========================================
    nB_half = gra_b.shape[0]
    gra_b_shifted = gra_b.copy()
    gra_b_shifted[:, 1] += aCC

    # Final bottom graphene lattice
    gra_b = np.vstack([gra_b, gra_b_shifted])

    # =========================================
    # Rotate to create twisted top layer
    # =========================================
    gra_t = (rotate_z(theta) @ gra_b.T).T

    # =========================================
    # Calculate supercell / cutting boundary
    # =========================================
    a1 = np.array([dx, 0.0])
    a2 = np.array([0.5 * dx, np.sqrt(3) * dx / 2.0])

    # L has 3 vertices in MATLAB; the last point remains [0, 0]
    L = np.zeros((3, 2), dtype=float)

    # 30-degree right-angle boundary
    L[0, :] = m * a1 + n * a2
    lambda_val = np.linalg.norm(L[0, :])

    # Scale by sqrt(3)/2 * xi
    L[0, :] = (np.sqrt(3) / 2.0) * xi * L[0, :]

    # Rotate by -30 degrees
    rotation_deg1 = np.deg2rad(-30.0)
    R1 = rotate_z(rotation_deg1)
    L[0, :] = R1 @ L[0, :]

    # Rotate by +90 degrees to get another edge vector
    rotation_deg2 = np.deg2rad(90.0)
    R2 = rotate_z(rotation_deg2)
    L[1, :] = R2 @ L[0, :]
    L[1, :] = L[1, :] * np.sqrt(3)

    # L[2, :] remains [0, 0], as in the MATLAB script
    # Polygon vertices: [L0, L1, origin]
    polygon = L.copy()

    # =========================================
    # Cut rotated top layer using polygon
    # =========================================
    gra_t_in_mask = points_in_polygon(gra_t[:, :2], polygon)
    gra_t_cut = gra_t[gra_t_in_mask, :]

    # =========================================
    # Plot
    # =========================================
    plt.figure(figsize=(8, 8))
    plt.scatter(gra_b[:, 0], gra_b[:, 1], s=4, c='blue', label='GraB')
    plt.scatter(gra_t_cut[:, 0], gra_t_cut[:, 1], s=4, c='red', label='GraTcut')
    plt.axis('equal')
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Bottom lattice and cut rotated top flake")
    plt.tight_layout()
    plt.show()

    # =========================================
    # Add z coordinates
    # =========================================
    gra_b = np.column_stack([gra_b, np.zeros(gra_b.shape[0])])  # z = 0
    gra_t_cut = np.column_stack([gra_t_cut, np.full(gra_t_cut.shape[0], 3.35)])  # z = 3.35

    # =========================================
    # Create GraT stage
    # =========================================
    distance = 0.0  # vertical distance between stage and top flake
    gra_t_stage = gra_t_cut.copy()
    gra_t_stage[:, 2] += distance

    # =========================================
    # Output / bookkeeping
    # =========================================
    nB = len(gra_b)
    nT = len(gra_t_cut)

    # MATLAB code had:
    # iG=nB+nT;
    # iG=nB+2*nT;
    # Final valid value is:
    iG = nB + 2 * nT

    nBonds = nT

    # Boundary of the graphene
    XGra = [np.min(gra_b[:, 0]), np.max(gra_b[:, 0]) + dx / 2.0]
    YGra = [np.min(gra_b[:, 1]), np.max(gra_b[:, 1]) + aCC / 2.0]
    ZGra = [-20.0, 20.0]

    filename = f"m{m}n{n}_tri30_{xi:.1f}.in"

    # =========================================
    # Write LAMMPS data file
    # =========================================
    with open(filename, "w", encoding="utf-8") as fid:
        fid.write("#Lammps data file\n\n")
        fid.write(f"{iG:5d} atoms\n")
        fid.write(f"{nBonds:5d} bonds\n")
        fid.write(f"{3:5d} atom types\n")
        fid.write(f"{1:5d} bond types\n\n")

        # boundaries
        fid.write(f"{XGra[0]:8.4f} {XGra[1]:8.4f} xlo xhi\n")
        fid.write(f"{YGra[0]:8.4f} {YGra[1]:8.4f} ylo yhi\n")
        fid.write(f"{ZGra[0]-10:8.4f} {ZGra[1]+10:8.4f} zlo zhi\n\n")

        fid.write("Masses\n\n")
        fid.write(" 1 12\n")
        fid.write(" 2 12\n")
        fid.write(" 3 12\n\n")

        fid.write("Atoms #bond\n\n")

        # Top cut flake: atom type 1
        for i in range(nT):
            x, y, z = gra_t_cut[i]
            fid.write(f"{i+1:5d} 1 {1:1d} {x:7.3f}  {y:7.3f}  {z:7.3f}\n")

        # Bottom graphene: atom type 2
        for i in range(nB):
            x, y, z = gra_b[i]
            fid.write(f"{i+1+nT:5d} 1 {2:1d} {x:7.3f}  {y:7.3f}  {z:7.3f}\n")

        # Stage copy of top flake: atom type 3
        for i in range(nT):
            x, y, z = gra_t_stage[i]
            fid.write(f"{i+1+nT+nB:5d} 1 {3:1d} {x:7.3f}  {y:7.3f}  {z:7.3f}\n")

        fid.write("\n Bonds\n\n")

        # Bond between each top flake atom and its staged counterpart
        for i in range(nBonds):
            atom1 = i + 1
            atom2 = i + 1 + nB + nT
            fid.write(f"{i+1:5d} {1:1d} {atom1:d}  {atom2:d}\n")

    print(f"LAMMPS data file written to: {filename}")
    print(f"lambda = {lambda_val:.10f}")


def rotate_z(theta: float) -> np.ndarray:
    """
    Construct a 2D rotation matrix corresponding to rotation around z-axis.

    Parameters
    ----------
    theta : float
        Rotation angle in radians.

    Returns
    -------
    np.ndarray
        2x2 rotation matrix.
    """
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])


def points_in_polygon(points: np.ndarray, polygon_vertices: np.ndarray) -> np.ndarray:
    """
    Determine whether points are inside or on the boundary of a polygon.

    MATLAB's inpolygon returns two arrays:
    - inside
    - on boundary
    Here we approximate the same behavior by:
    - checking interior points with Path.contains_points
    - checking boundary points using a small radius tolerance

    Parameters
    ----------
    points : np.ndarray
        Array of shape (N, 2), point coordinates.
    polygon_vertices : np.ndarray
        Array of shape (M, 2), polygon vertices.

    Returns
    -------
    np.ndarray
        Boolean mask of shape (N,), True for points inside or on boundary.
    """
    path = Path(polygon_vertices)

    # Strictly inside
    inside = path.contains_points(points)

    # Approximate "on boundary" by using a tiny positive and negative radius
    # This is a practical substitute for MATLAB's exact inpolygon boundary output.
    on_boundary = path.contains_points(points, radius=1e-9) ^ path.contains_points(points, radius=-1e-9)

    return inside | on_boundary



if __name__ == "__main__":
    main()