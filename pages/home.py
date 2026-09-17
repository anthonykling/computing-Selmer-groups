import streamlit as st
#from GaussianInt import GaussInt, gaussian_primes
#from Laplacian import Selmer_group

st.set_page_config(
    page_title="Gaussian Selmer Groups",
    page_icon="𝑖",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("Selmer Group Calculator over Q(i)")

st.write(r"""
This application implements the algorithm presented in https://arxiv.org/abs/2410.22714 
to compute the Selmer group of $E_b : y^2 = x^3 + bx$ over $\mathbb{Q}(i)$
where $b\in\mathbb{Z}[i]$.
""")

st.info(
    "Use the sidebar to choose between computing the Selmer group of a single curve or a family of curves"
)