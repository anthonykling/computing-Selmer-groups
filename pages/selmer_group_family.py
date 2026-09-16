import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from utils import parse_gaussian_integer

from GaussianInt import GaussInt, gaussian_primes
from Laplacian import Selmer_group


st.title("Gaussian Prime Explorer")

st.write(
    r"""
    Generate Gaussian primes satisfying
    \(N(p) < \text{norm\_size}\) and
    \(p \equiv a \pmod m\), then compute the size of the
    corresponding Selmer groups.
    """
)


col1, col2 = st.columns(2)

with col1:
    norm_size = st.number_input(
        "Norm bound",
        min_value=2,
        value=1000,
        step=100
    )

with col2:
    a_input = st.text_input(
        r"$a$",
        value="1"
    )

m_input = st.text_input(
    r"$m$",
    value="1 + i"
)


if st.button("Generate and Compute", type="primary"):

    try:
        a = parse_gaussian_integer(a_input)
        m = parse_gaussian_integer(m_input)

        with st.spinner("Generating Gaussian primes..."):
            primes = gaussian_primes(norm_size, a, m)

        st.write(f"Found **{len(primes)} Gaussian primes**.")

        if len(primes) == 0:
            st.warning("No Gaussian primes were found.")
            st.stop()

        results = []

        progress = st.progress(0)

        for j, p in enumerate(primes):

            S = Selmer_group(p)

            results.append({
                "prime": str(p),
                "norm": p.norm(),
                "Selmer group length": len(S)
            })

            progress.progress((j + 1) / len(primes))

        df = pd.DataFrame(results)

        st.subheader("Results")

        st.dataframe(df)

        # -------------------------------------------------
        # Bar graph
        # -------------------------------------------------

        st.subheader("Selmer Group Sizes")

        fig, ax = plt.subplots()

        ax.bar(
            df["prime"],
            df["Selmer group length"]
        )

        ax.set_xlabel("Gaussian prime")
        ax.set_ylabel("Length of Selmer group")
        ax.set_title("Selmer Group Lengths")

        plt.xticks(rotation=90)

        st.pyplot(fig)

    except Exception as e:
        st.error(f"Error: {e}")