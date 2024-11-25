import streamlit as st
import base64

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
    return img_base64

# Add custom HTML for page title and favicon
favicon_path = 'assets/favicon.ico'  # Replace with the path to your favicon
title = "Your Custom Title"  # Replace with your desired title

# Add custom HTML for the page meta (title and favicon)
try:
    favicon_base64 = image_to_base64(favicon_path)

    # Embed the title and favicon in the page's head
    page_meta = f"""
    <head>
        <title>{title}</title>
        <link rel="icon" href="data:image/x-icon;base64,{favicon_base64}" type="image/x-icon">
    </head>
    """
    st.markdown(page_meta, unsafe_allow_html=True)

except FileNotFoundError:
    st.error(f"Favicon not found: {favicon_path}")

# Load the local image and convert it to base64 for the background
image_path = 'assets/home_bg.png'  # Replace with the path to your image
try:
    image_base64 = image_to_base64(image_path)

    # Embed the image in CSS for background styling
    page_bg_img = f"""
    <style>
    .stApp {{
      background-color: #1E1E1E; 
      background-image: url("data:image/png;base64,{image_base64}");
      background-position: right;  /* Adjusts the position */
      background-size: 85%;  /* Zooms out the image */
      background-repeat: no-repeat;
      padding: 3rem;
      border-radius: 10px;
      text-align: center;  
      padding-left: 80px;  /* Add space for sidebar */
    }}

    /* Adjust the sidebar width (optional) */
    .css-1g1z59i {{
      width: 200px;  /* Increase the width of the sidebar */
    }}

    /* Adjust the main content area */
    .css-1y4v5i9 {{
      margin-left: 220px;  /* Shift main content to the right */
    }}
    </style>
    """
    st.markdown(page_bg_img, unsafe_allow_html=True)

except FileNotFoundError:
    st.error(f"Image not found: {image_path}")
