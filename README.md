# Data-Driven Material Project

This repository contains a finite element method (FEM) simulation framework for analyzing in-plane loaded plates using quadrilateral elements under 2D plane stress conditions. The project explores the effects of varying material stiffness distributions on structural responses such as displacement, stress, strain, and compliance.

For a detailed overview of the problem being solved, please refer to the project statement:

📘 [Multi-Phase Material Design Project (PDF)](https://github.com/Waris-K/data_driven_mat_proj/blob/main/Multi-Phase-Material-Design-Project.pdf)

---

## 📌 Project Status

**Active**

---

## 🎯 Objective

The primary goal of this project is to develop a FEM-based simulation tool that can:

- Model plates with heterogeneous material properties
- Analyze the impact of different stiffness distributions (e.g., binary layouts)
- Evaluate key metrics such as average displacement, stress, strain, and compliance
- Visualize force vs. displacement and stress vs. strain relationships

---

## 🧰 Methods Used

- Finite Element Analysis (FEA)
- Linear Elasticity Theory
- Numerical Integration (Gaussian Quadrature)
- Matrix Assembly and Solution Techniques
- Data Visualization

---

## 🛠️ Technologies

- MATLAB
- Custom FEM Functions
- Visualization with MATLAB plots

---

## 📁 Project Structure

```bash
data_driven_mat_proj/
├── MeshRectanglularPlate.m
├── getGaussQuadrature.m
├── shapeFunctions.m
├── Jacobian.m
├── strainDisplacementMatrix.m
├── materialMatrix.m
├── plotBoundaryConditionsAndForces.m
├── plotDeformedShape.m
├── main.m
└── README.md
```

---

## 🚀 Getting Started

### Prerequisites

- MATLAB R2018b or later

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/Waris-K/data_driven_mat_proj.git
   ```

2. Navigate to the project directory:

   ```bash
   cd data_driven_mat_proj
   ```

3. Open MATLAB and add the project directory to the path.

4. Run the main simulation:

   ```matlab
   main
   ```

---

## 📊 Features

- Customizable material stiffness layouts
- Global stiffness matrix and force vector assembly
- Application of Dirichlet boundary conditions and external loads
- Nodal displacement and elemental stress/strain computation
- Compliance and average displacement calculations
- Force vs. displacement and stress vs. strain plots

---

## 📈 Results

Key outputs include:

- Nodal displacements
- Elemental stress and strain
- Compliance values
- Visual plots for structural response

---

## 🤝 Contributing

Contributions are welcome! Please fork the repository and open a pull request.

---

## 📄 License

This project is licensed under the MIT License.

---

## 📬 Contact

For questions or feedback, feel free to reach out via [GitHub Issues](https://github.com/Waris-K/data_driven_mat_proj/issues).


---



