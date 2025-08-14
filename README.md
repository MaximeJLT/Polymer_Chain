# Freely Jointed Chain (FJC) Polymer Model

## Objectives

This work aims to simulate the Freely Jointed Chain (FJC) polymer model and verify its validity by comparing characteristic quantities with theoretical predictions.  
In particular, we focus on the **end-to-end vector**:

$$
\vec{Q}
$$

and compare simulation results with analytical expectations.

---

## Introduction

A polymer is a long chain of atoms or molecules called **monomers**.  
In the FJC model, the polymer is represented as a chain of $N$ rigid segments of length $b$. Each segment can point in a **random direction**, independently of the others. This means there is:

- **No chain stiffness**
- **No steric hindrance**
- **No correlation** between neighboring segments

Because each monomer can adopt any spatial orientation, the polymer can take an extremely large number of possible conformations.

The position of monomer $n$ is given by:

$$
\vec{R_n} = \vec{R_1} + \sum_{i = 1}^{n-1} \vec{b_i}
$$

In our simulation, we set:

$$
\vec{R_1} = \vec{0}
$$

The **end-to-end vector** is then:

$$
\vec{Q} = \vec{R}_{N + 1} - \vec{R_1} = \sum_{i=1}^N \vec{b_i}
$$

Because the bond orientations are random, the first moment of $\vec{Q}$ is:

$$
\langle \vec{Q} \rangle = \vec{0}
$$

The useful information lies in the **mean squared end-to-end distance**:

$$
\begin{aligned}
\langle Q^2 \rangle &= \langle \vec{Q} \cdot \vec{Q} \rangle \\
&= \left\langle \left( \sum_{i=1}^N \vec{b_i} \right) \cdot \left( \sum_{j=1}^N \vec{b_j} \right) \right\rangle \\
&= \sum_{i=1}^N \sum_{j=1}^N \langle \vec{b_i} \cdot \vec{b_j} \rangle \\
&= \sum_{i=1}^N \sum_{j=1}^N b^2 \langle \cos(\theta_{ij}) \rangle \\
&= \sum_{i=1}^N \sum_{j=1}^N b^2 \delta_{ij} \\
&= N b^2
\end{aligned}
$$

---

## Simulation description

The simulation proceeds as follows:

1. **Parameter definition**
   - Number of bonds $N$
   - Bond length $b$
   - Number of conformations $T$

2. **Random orientation generation**  
   Bond vectors are generated using rejection sampling to ensure **uniform distribution** over the sphere:
   - Polar angle $\phi$ in $[0, 2\pi]$
   - Azimuthal angle $\theta$ in $[-\pi, \pi]$

3. **Conformation generation**  
   From these random orientations, the positions of all monomers are computed for each conformation.

---

## Results

### Figure 1: Polar angle distribution
![Polar angle distribution](TD7/Fig_2.pdf)  

Distribution of the polar angle $\phi$ obtained with rejection sampling.  
It matches the expected **theoretical distribution** in the range $[0, 2\pi]$.

---

### Figure 2: Azimuthal angle distribution
![Azimuthal angle distribution](TD7/Fig_3.pdf)  

Distribution of the azimuthal angle $\theta$ obtained with rejection sampling.  
It matches the expected **theoretical distribution** in the range $[-\pi, \pi]$.

---

### Figure 3: Example conformations (VMD)
![Conformation 1](TD7/1.pdf)  
![Conformation 2](TD7/2.pdf)  
![Conformation 3](TD7/3.pdf)  
![Conformation 4](TD7/4.pdf)  
![Conformation 5](TD7/5.pdf)  
![Conformation 6](TD7/6.pdf)  

Different polymer conformations visualized with **VMD**.  
The **dark blue** monomer is the first in the chain, and the **dark red** monomer is the last.

---

## Mean squared end-to-end distance

Two conditions should hold:

$$
\langle \vec{Q} \rangle = \vec{0}
$$

$$
\langle Q^2 \rangle = N b^2
$$

From the simulation, we obtain:

$$
\langle \vec{Q} \rangle =
\begin{bmatrix}
-0.38338 \\
-0.1347 \\
0.15766
\end{bmatrix}
\quad\text{and}\quad
\langle Q^2 \rangle = 891.73
$$

These values are **very close** to the theoretical expectations for $N = 100$ and $b = 3$.

---

## Scaling with $N$

To further validate the model, we vary $N$ from 100 to 1000 in steps of 100, and check the linear relationship:

$$
\langle Q^2 \rangle = N b^2
$$

### Figure 4: Linear regression of $\langle Q^2 \rangle$ vs $N$
![Linear regression](TD7/reg.pdf)  

From the slope, we estimate:

$$
b = 2.9537
$$

with a coefficient of determination:

$$
R^2 = 0.9975
$$

This confirms the correctness of the simulation.

---

## End-to-end distance distribution

We also analyze the distribution of $Q$ for $N = 100$.

### Figure 5: Theoretical distribution of $Q$
![Theoretical Q distribution](TD7/Fig_4.pdf)  

The expected shape is **Gaussian**.

### Figure 6: Computed distribution of $Q$
![Computed Q distribution](TD7/Fig_1.pdf)  

The simulated distribution matches the Gaussian profile predicted by theory.

---

## Conclusion

We have simulated the **Freely Jointed Chain (FJC)** model and validated it by:

1. Confirming that $\langle \vec{Q} \rangle \approx \vec{0}$.
2. Verifying that $\langle Q^2 \rangle \approx N b^2$.
3. Matching the simulated angular distributions to theoretical expectations.
4. Obtaining correct scaling of $\langle Q^2 \rangle$ with $N$.
5. Reproducing the Gaussian distribution of $Q$.

This model is widely used in polymer physics, including for the study of **DNA**.
