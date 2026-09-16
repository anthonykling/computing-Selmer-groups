import streamlit as st
import pandas as pd
from utils import parse_gaussian_integer

from GaussianInt import GaussInt, gaussian_primes
from Laplacian import Selmer_group



st.title("Selmer Group Calculator")

st.write(
    r"Compute the Selmer group of the elliptic curve $E_b/ \mathbb{Q}(i)$"
    r" associated to a Gaussian integer $b \in \mathbb{Z}[i]$"
    r" where $E_b : y^2 = x^3 + bx$."
)

input_type = st.radio(
    "How would you like to specify b?",
    ["Enter b directly", "Enter the factorization"]
)

show_factorization = st.toggle("Show factorization of $b$", value=False)


# ---------------------------------------------------------
# Direct input
# ---------------------------------------------------------

if input_type == "Enter b directly":

    st.write('*b* must be in the form *x+iy* or *x*')
    b_input = st.text_input(
        "Gaussian integer b",
        value="3 + 2i"
    )

    if st.button("Compute Selmer Group", type="primary"):

        try:
            b = parse_gaussian_integer(b_input)

            st.write(f"**b =** `{b}`")
            if show_factorization:
                with st.spinner("Factoring b..."):
                    f = b.factor()
                    result = {"Factor": f.keys(), "Exponent": f.values()}
                    st.subheader("Primary Decomposition")
                    st.dataframe(pd.DataFrame(result))
                

            with st.spinner("Computing Selmer group..."):
                S = Selmer_group(b)

            st.success("Computation complete!")

            st.metric("Size of Selmer group", len(S))

            st.write("Selmer group:")
            st.write(S)

        except Exception as e:
            st.error(f"Error: {e}")


# ---------------------------------------------------------
# Factorization input
# ---------------------------------------------------------

else:

    st.subheader("Specify the factorization of b")
    col1, col2 = st.columns(2)

    with col1:
        s_b = st.number_input(
            r"Exponent $s_b$ of $i$",
            min_value=0,
            max_value=3,
            value=0,
            step=1
        )

    with col2:
        t_b = st.number_input(
            r"Exponent $t_b$ of $(1+i)$",
            min_value=0,
            max_value=3,
            value=0,
            step=1
        )

    st.write(
        r"Enter the factors of $b$ "
        r"and their exponents."
    )

    num_factors = st.number_input(
        "Number of factors",
        min_value=1,
        max_value=20,
        value=1,
        step=1
    )

    prime_factors = []

    for j in range(num_factors):

        col1, col2 = st.columns(2)

        with col1:
            p_input = st.text_input(
                f"Factor {j+1}",
                key=f"prime_{j}",
                placeholder="e.g. 5 + 2i"
            )

        with col2:
            exponent = st.number_input(
                f"Exponent {j+1}",
                min_value=1,
                max_value=10,
                value=1,
                step=1,
                key=f"exp_{j}"
            )

        if p_input.strip():
            prime_factors.append((p_input, exponent))

    if st.button("Compute Selmer Group", type="primary"):

        try:
            b = GaussInt(0, 1) ** s_b
            b *= GaussInt(1, 1) ** t_b

            for p_input, exponent in prime_factors:
                p = parse_gaussian_integer(p_input)
                b *= p ** exponent

            st.write(f"**Constructed b =** `{b}`")
            
            if show_factorization:
                with st.spinner("Factoring b..."):
                    f = b.factor()
                    result = {"Factor": f.keys(), "Exponent": f.values()}
                    st.subheader("Primary Decomposition")
                    st.dataframe(pd.DataFrame(result))    

            with st.spinner("Computing Selmer group..."):
                S = Selmer_group(b)

            st.success("Computation complete!")

            st.metric("Length of Selmer group", len(S))

            st.write("Selmer group:")
            st.write(S)

        except Exception as e:
            st.error(f"Error: {e}")