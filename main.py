import streamlit as st
import base64

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
    return img_base64

# Load the local image and convert it to base64
image_path = 'assets/home_bg.png'  # Replace with your local PNG image path
image_base64 = image_to_base64(image_path)


# Add custom background and text styling using CSS
st.markdown(
    f"""
    <style>
    .main {{
        background-color: #1E1E1E;  /* Dark background color */
        background-image: url('data:image/png;base64,{image_base64}');  /* Local image as base64 */
        background-size: cover;
        background-position: center;
        color: #FFFFFF;  /* White text color */
        padding: 3rem;
        border-radius: 15px;
        text-align: center;
    }}
    
    h1, h2, h3 {{
        color: #FFD700;  /* Gold color for headings */
        font-family: 'Arial', sans-serif;
        font-weight: bold;
    }}
    
    .subheader {{
        font-size: 24px;
        color: #00BFFF;  /* Bright blue for subheader */
    }}
    </style>
    """, unsafe_allow_html=True
)

