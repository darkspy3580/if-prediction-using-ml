import streamlit as st

# This must be the first Streamlit call in your script
st.set_page_config(
    page_title="IF Prediction App",  # Your custom page title
    page_icon="assets/Logov.png",  # Path to your favicon image
)

# Display a simple message
st.title("Welcome to IF Prediction App")
st.write("This is your Streamlit app with a custom title and favicon.")
