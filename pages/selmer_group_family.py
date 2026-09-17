import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from utils import parse_gaussian_integer

from GaussianInt import GaussInt, gaussian_primes
from Laplacian import Selmer_group


st.title("Gaussian Prime Explorer")

st.write(
    r"""
    This will generate all Gaussian primes $p$ of norm less than the specified bound and
    $p \equiv a (\text{mod}\;m)$, then computes the size of the
    corresponding Selmer group of $E_p$.  The data is then presented as a frequency plot.
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
                "Selmer group size": len(S)
            })

            progress.progress((j + 1) / len(primes))

        df = pd.DataFrame(results)

        st.subheader("Results")

        st.dataframe(df)

        # -------------------------------------------------
        # Bar graph
        # -------------------------------------------------

        st.subheader("Frequency of Selmer Group Sizes")

        frequency = (
            df["Selmer group size"]
            .value_counts()
            .sort_index()
        )

        fig, ax = plt.subplots(figsize=(8, 5))

        ax.bar(
            frequency.index.astype(str),
            frequency.values,
            width=0.7
        )

        ax.set_xlabel("Selmer group size")
        ax.set_ylabel("Frequency")
        ax.set_title("Distribution of Selmer Group Sizes")

        # Add frequency labels above each bar
        for i, value in enumerate(frequency.values):
            ax.text(
                i,
                value,
                str(value),
                ha="center",
                va="bottom"
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        st.pyplot(fig)

    except Exception as e:
        st.error(f"Error: {e}")