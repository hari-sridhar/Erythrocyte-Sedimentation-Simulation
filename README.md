
# Modeling and Simulation of the Dynamics of Erythrocyte Sedimentation and Plasma Counter-Flow

**Author:** Harish Sridhar

**Institute:** Department of Applied Mechanics, IIT Madras

---

## 📌 Overview

This MATLAB-based simulation project models the early dynamics of erythrocyte (red blood cell) sedimentation in blood plasma. This process is key to understanding inflammation and other pathological conditions in the body, often measured by ESR (Erythrocyte Sedimentation Rate). The simulation treats each erythrocyte as a discrete agent subjected to physical forces and implements graph-based current flow to represent plasma behavior around and between the particles.

---

## 🎯 Objectives

- Simulate sedimentation of erythrocytes under physiological forces.
- Model inter-particle and wall interactions using force-balance and lubrication theory.
- Use a graph-theoretic approach to simulate plasma flow via an electrical analogy.
- Visualize trajectories, velocities, and currents in 2D.
- Enable analysis from simple 2-particle to complex N-body systems under realistic boundary conditions.

---

## ⚙️ Methodology

### 1. Force Modeling

Each particle experiences a combination of:
- **Gravitational Force**: Downward pull based on mass.
- **Buoyant Force**: Upward push due to plasma displacement.
- **Drag Force**: Resists motion through viscous medium; modeled using Stokes' law.
- **Lubrication Force**: Short-range interaction force that prevents overlap and adds repulsion/attraction based on distance.

The system is assumed to be in **steady state**, i.e., the net force on each particle is zero:  
\[
ΣF = 0 {terminal velocity is quickly achieved}
\]

---

### 2. Particle and Graph Representation

- Each erythrocyte is modelled with **4 nodes** placed at its cardinal points (top, bottom, left, right).
- **Intra-particle edges** connect these nodes.
- **Inter-particle edges** are determined through Delaunay triangulation of mirrored and original nodes.
- Nodes are extended with **mirrored versions** to handle wall interactions.

---

### 3. Electrical Network Analogy

- The network formed by nodes and edges is interpreted as a **resistor network**.
- A **Laplacian matrix** is built from the resistances.
- Voltage sources/sinks are assigned based on particle configuration.
- The resulting voltage vector is solved using pseudo-inverse:  
  \[
  L*v = S
  \]
- Currents are then derived from voltages via Ohm's law.

---

### 4. Periodic Boundary Conditions

To simulate a continuous unbounded domain:
- Particles near the left/right walls are **mirrored** to interact across the boundary.
- This avoids edge artifacts and maintains realism in the system.

---

## 💻 Simulation Flow

1. **Initialization**
   - Define domain, number of particles, and physical parameters.
   - Assign initial positions and radii to particles.
   - Generate their 4 nodes and edges.

2. **Time Loop**
   - For each time step:
     - Recompute pairwise distances and interaction vectors.
     - Calculate force matrices: P (xx), Q (xy), S (yy).
     - Solve force balance equation for velocities.
     - Recompute resistance matrix and edge currents.
     - Update positions.

3. **Post-Processing**
   - Save data to `data.mat` file.
   - Plot velocity vs. time and position vs. time.
   - Animate streamline current flows over time.

---

## 📁 File Descriptions

| File | Purpose |
|------|---------|
| `og_main.m` | Primary driver script |
| `generateParticleNodes.m` | Generate particle nodes |
| `generateParticleEdges.m` | Create 4-node internal connections |
| `mirrorParticleNodes.m` | Mirror nodes across vertical boundaries |
| `compute_wall_edge_resistance.m` | Compute triangle-based resistance matrix |
| `generateSource.m` | Assign source and sink nodes per particle |
| `getEdgesAscending.m` | Sort edges by Y-coordinate |
| `computeNodeVoltages.m` | Compute voltages from Laplacian network |
| `computeEdgeCurrents.m` | Compute currents using Ohm's law |
| `wall_lubrication_correction.m` | Adjust drag near wall |
| `animateParticleStreamline.m` | Animate current arrows and positions |
| `plotParticleDynamics.m` | Plot trajectories and velocities |


---

## 📊 Outputs

- **Trajectory plots**: Shows particle movement over time.
- **Velocity plots**: X and Y components plotted over simulation duration.
- **Streamline animation**: Visualizes edge currents as arrows with color and magnitude.
- **Edge resistance matrix**: Printed to console and stored per timestep.

---

## 🔬 Scientific Insights

This simulation framework allows:
- Studying how erythrocytes aggregate (rouleaux formation).
- Analyzing how plasma flow counteracts or facilitates sedimentation.
- Observing effects of viscosity, particle spacing, and wall proximity.
- Extension to pathological cases (e.g., anisocytosis, abnormal plasma).

By combining agent-based modeling and electrical analogies, this simulation provides a scalable and realistic approach to biomedical sedimentation modeling.

---

## 🚀 How to Run

```bash
1. Open MATLAB.
2. Run Main.m
3. When prompted, enter the number of particles (example: 3 or 5).
4. View outputs: animation, velocity plots, trajectory graphs.
```

---

## 🧠 Future Work

- Introduce deformability of erythrocytes.
- Add plasma flow modeling via Navier-Stokes (hybrid approach).
- Integrate with clinical data for validation.
- Include aggregation dynamics based on fibrinogen or immunoglobulin interaction.

---

## 📎 License and Acknowledgements

Developed at the Department of Applied Mechanics, IIT Madras, under the guidance of Prof. Danny Raj M and Prof. Ramakrishnan S .  
Academic use only.

---

