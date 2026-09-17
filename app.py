import streamlit as st

pg = st.navigation([
    st.Page("pages/home.py", title="Home"),
    st.Page("pages/selmer_group.py", title="Selmer Group"),
    st.Page("pages/selmer_group_family.py", title="Selmer Group of Families"),
    ])

pg.run()