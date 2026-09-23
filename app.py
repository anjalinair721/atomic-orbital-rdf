import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from math import factorial
from scipy.special import genlaguerre

st.set_page_config(
    page_title="Atomic Orbital RDF Explorer",
    page_icon="⚛️",
    layout="centered"
)

st.title("⚛️ Atomic Orbital RDF Explorer")
st.write("Calculate and visualize the Radial Distribution Function (RDF) for a hydrogen-like atomic orbital.")

# Quantum-number controls
st.sidebar.header("Orbital Parameters")

n = st.sidebar.number_input(
    "Principal quantum number (n)",
    min_value=1,
    max_value=10,
    value=4,
    step=1
)

orbital_letters = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}

valid_l = list(range(n))
l = st.sidebar.selectbox(
    "Azimuthal quantum number (l)",
    options=valid_l,
    format_func=lambda x: f"{x} ({orbital_letters.get(x, f'l={x}')})"
)

orbital_name = f"{n}{orbital_letters.get(l, f'l={l}')}"

# Orbital analysis
angular_nodes = l
radial_nodes = n - l - 1
total_nodes = n - 1

# Average distance and plot boundary
r_avg = 0.5 * (3 * (n ** 2) - l * (l + 1))
r_max = 2.5 * r_avg

# Radial coordinate
r = np.linspace(0, r_max, 1000)

# Generalized radial wavefunction
rho = (2.0 * r) / n
norm = np.sqrt(
    (2.0 / n) ** 3
    * factorial(n - l - 1)
    / (2.0 * n * factorial(n + l))
)
laguerre = genlaguerre(n - l - 1, 2 * l + 1)

R_nl = norm * (rho ** l) * np.exp(-rho / 2.0) * laguerre(rho)

# Radial Distribution Function
rdf = (r ** 2) * (R_nl ** 2)

st.subheader(f"{orbital_name} Orbital Analysis")

c1, c2, c3 = st.columns(3)
c1.metric("Angular nodes", angular_nodes)
c2.metric("Radial nodes", radial_nodes)
c3.metric("Total nodes", total_nodes)

st.write(f"**Average electron distance ⟨r⟩:** {r_avg:.3f} a₀")
st.write(f"**Plot boundary (rₘₐₓ):** {r_max:.3f} a₀")

# Plot
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(r, rdf, label=f"{orbital_name} orbital", linewidth=2)
ax.set_title(f"Radial Distribution Function for {orbital_name}")
ax.set_xlabel("Distance from nucleus r (a₀)")
ax.set_ylabel("RDF")
ax.grid(True, linestyle="--", alpha=0.6)
ax.legend()
ax.set_xlim(0, 5 * (n ** 2))
ax.set_ylim(bottom=0)

st.pyplot(fig)
plt.close(fig)

st.caption(
    "Based on the equations and calculations in the supplied plotting_orbitals notebook."
)
