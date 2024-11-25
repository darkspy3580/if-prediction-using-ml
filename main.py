import streamlit as st
import base64

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
    return img_base64

# Function to embed title and favicon in the head section
def embed_title_and_favicon():
    title = "Your Custom Title"  # Set your custom title here
    favicon_path = '/mnt/data/45A20AEB-9AEA-43E2-810B-34C86E308FA5.png'  # Path to your favicon file
    try:
        # Convert favicon to base64
        with open(favicon_path, "rb") as favicon_file:
            favicon_base64 = base64.b64encode(favicon_file.read()).decode('utf-8')

        # Embed the title and favicon using HTML
        page_meta = f"""
        <head>
            <title>{title}</title>
            <link rel="icon" href="data:image/png;base64,{favicon_base64}" type="image/png">
        </head>
        """
        st.markdown(page_meta, unsafe_allow_html=True)
    except FileNotFoundError:
        st.error(f"Favicon not found: {favicon_path}")

# Add favicon and title to the page
embed_title_and_favicon()

# Load the background image
image_path = 'assets/home_bg.png'  # Replace with the path to your image

try:
    image_base64 = image_to_base64(image_path)

    # Embed the background image and custom styles
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
