import streamlit as st
import base64

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
    return img_base64

# Load the local image and convert it to base64
image_path = 'assets/home_bg.png'  # Replace with the path to your image
try:
    image_base64 = image_to_base64(image_path)

    # Embed the image in CSS
    page_bg_img = f"""
    <style>
    .stApp {{
      background-color: #1E1E1E; 
      background-image: url("data:image/png;base64,{image_base64}");
      background-position: left;  /* Adjusts the position */
      background-size: 100%;  /* Zooms out the image */
      background-repeat: no-repeat;
      padding: 3rem;
      border-radius: 10px;
      text-align: center;  
    }}
    </style>
    """
    st.markdown(page_bg_img, unsafe_allow_html=True)
except FileNotFoundError:
    st.error(f"Image not found: {image_path}")
