import streamlit as st
from GaussianInt import GaussInt
from Laplacian import Selmer_group

st.title("Selmer Group Calculator")

st.write("Enter a Gaussian integer b = x + yi.")

x = st.number_input("Real part", value=1, step=1)
y = st.number_input("Imaginary part", value=0, step=1)

if st.button("Calculate"):
    st.write('Calculating Selmer group for the Gaussian integer...')
    b = GaussInt(int(x), int(y))

    result = Selmer_group(b)

    st.write("Gaussian integer:")
    st.latex(f"b = {x} + {y}i")

    st.write("Selmer group:")
    st.write(result)