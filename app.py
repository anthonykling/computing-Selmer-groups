import streamlit as st
from GaussianInt import GaussInt, gaussian_primes
from Laplacian import Selmer_group

import streamlit as st

st.set_page_config(
    page_title="Selmer Group",
    page_icon="𝑖",
    layout="wide"
)

st.title("Selmer Group Calculator over Q(i)")

st.write("""
This application provides tools for computing Selmer groups associated
to Gaussian integers and for exploring how Selmer group sizes vary
over families of Gaussian primes.
""")

st.info(
    "Use the sidebar to choose between the Selmer Group Calculator "
    "and the Gaussian Prime Explorer."
)