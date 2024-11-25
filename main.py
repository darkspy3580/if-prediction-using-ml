import streamlit as st
import base64
import os

# Function to convert image to base64
def image_to_base64(image_path):
    try:
        with open(image_path, "rb") as img_file:
            img_base64 = base64.b64encode(img_file.read()).decode("utf-8")
        return img_base64
    except FileNotFoundError:
        st.error(f"Error: File not found at {image_path}")
        return None

# Set the page configuration
st.set_page_config(page_title="Custom Background App", layout="wide")

# Define the path to the background image
image_path = "assets/home_bg.png"  # Replace with your local PNG image path

# Ensure the file path works both locally and on Streamlit Cloud
current_dir = os.path.dirname(os.path.abspath(__file__))
absolute_image_path = os.path.join(current_dir, image_path)

# Check if the file exists and convert to base64
if os.path.exists(absolute_image_path):
    image_base64 = image_to_base64(absolute_image_path)
    if image_base64:
        st.success(f"Image successfully loaded from: {absolute_image_path}")
else:
    st.error(f"Image not found at: {absolute_image_path}")
    st.stop()  # Stop the app if the image isn't found

# Add custom background and text styling using CSS
st.markdown(
    f"""
    <style>
    body {{
        background-color: #1E1E1E; /* Dark background color */
        background-image: url('data:image/png;base64,{image_base64}');
        background-size: cover;
        background-position: center;
        color: #FFFFFF; /* White text color */
    }}
    
    h1, h2, h3 {{
        color: #FFD700; /* Gold color for headings */
        font-family: 'Arial', sans-serif;
        font-weight: bold;
    }}
    
    .subheader {{
        font-size: 24px;
        color: #00BFFF; /* Bright blue for subheader */
    }}
    </style>
    """,
    unsafe_allow_html=True,
)

# Add sample content to test the layout
st.title("Welcome to My Custom Background App")
st.subheader("This is a subheader with custom styling.")
st.write("Your background image should now be visible.")
